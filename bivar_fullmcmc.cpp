//// mcmc with full cond for regr parameters, sampled after spline
// intercept with regr parameters 
//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// MCMC FCNS

//simulate multivariate normal rvs
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols; 
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


//density of mv normal rv

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    //std::cout << x.row(i) << "\n" << mean;
    
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

double dmvnrm_arma(arma::vec x,  
                   arma::rowvec mean,  
                   arma::mat sigma, 
                   bool logd = false) { 
  //int n = x.n_rows;
  int xdim = x.size();
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  //for (int i=0; i < n; i++) {
  //std::cout << x.row(i) << "\n" << mean;
  
  arma::vec z = rooti * arma::trans( x.t() - mean) ;    
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  //}  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}


/*-------------------------------------------------------
# Generate Draws from an Inverse Wishart Distribution
# via the Bartlett Decomposition
#--------------------------------------------------------
# NOTE: output is identical to riwish from MCMCpack
#       provided the same random seed is used
#--------------------------------------------------------
#   n     number of samples
#   S     scale matrix 
#   v     degrees of freedom
#-------------------------------------------------------*/
// [[Rcpp::export]]
arma::cube rinvwish(int n, int v, arma::mat S){
  RNGScope scope;
  int p = S.n_rows;
  arma::mat L = chol(inv_sympd(S), "lower");
  arma::cube sims(p, p, n, arma::fill::zeros);
  for(int j = 0; j < n; j++){
    arma::mat A(p,p, arma::fill::zeros);
    for(int i = 0; i < p; i++){
      int df = v - (i + 1) + 1; //zero-indexing
      A(i,i) = sqrt(R::rchisq(df)); 
    }
    for(int row = 1; row < p; row++){
      for(int col = 0; col < row; col++){
        A(row, col) = R::rnorm(0,1);
      }
    }
    arma::mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
    sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}


//weighted sample function, equivalent of R sample()
// // [[Rcpp::export]]
int cpp_sample(int h, arma::vec probs){
  //int n = probs.size();
  int out = h;
  arma::vec cumprobs(h);
  double r = R::runif(0,1);
  cumprobs[0] = probs[0];
  if(r < cumprobs[0]){
    out = 1;
    return(out);
  }
  for(int i=1;i<h;i++){
    cumprobs[i] = cumprobs[i-1] + probs[i];
    if((r < cumprobs[i]) && (r > cumprobs[i-1])){
      out = i + 1;
      return(out);
    }
  } 
  return(out);
}

arma::vec cpp_sample(NumericVector options, NumericVector probs){
  int n = options.length();
  arma::vec out(1);
  NumericVector cumprobs(n);
  double r = R::runif(0,1);
  cumprobs[0] = probs[0];
  if(r < cumprobs[0]){
    out[0] = options[0];
    return(out);
  }
  for(int i=1;i<n;i++){
    cumprobs[i] = cumprobs[i-1] + probs[i];
    if((r < cumprobs[i]) && (r > cumprobs[i-1])){
      out[0] = options[i];
      return(out);
    }
  } 
  return(out);
}


//equivalent of R prod()
// // [[Rcpp::export]]
double cpp_prod(arma::vec v){
  int n = v.size(); 
  double out = 1.0;
  for(int i=0;i<n;i++){
    out *= (1.0 - v[i]);
  }
  return(out);
}

//samples zeta VECTOR, no outside for loop needed
// // [[Rcpp::export]]
arma::ivec sample_zeta(arma::mat x, arma::vec pi, arma::mat mu, arma::cube Sigma){
  int n = x.n_rows;
  int h = pi.size();
  arma::mat den(n,h);
  arma::mat probs(n,h);
  arma::rowvec sprobs(h);
  arma::ivec out(n);
  for(int k=0;k<h;k++){
    den.col(k) = dmvnrm_arma(x,mu.row(k),Sigma.slice(k),true);
    
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<h;j++){
      probs(i,j) = log(pi[j])+den(i,j);
    }

    sprobs = exp(probs.row(i)-log(sum(exp(probs.row(i)))));
    out[i] = cpp_sample(h, sprobs.t());
  }
  
  return(out);
}

double sample_v(int K, int h, int nh, double alpha, int nhp){
  if(K < h){
    return(R::rbeta(1+nh,alpha+nhp));
  }
  else{
    return(1);
  }
}

//calculates pi based on sampled V
//double calc_pi(int K, NumericVector V){
//  if(K > 1){
//    return(V[K]*cpp_prod(V[seq(0,K-1)]));
//  }
//  else{
//    return(V[K]);
//  }
//}


//samples mu
//NumericVector sample_mu()

//sample alpha
double sample_alpha(double a, double b, int h, arma::vec V){
  return(R::rgamma(a+h-1,1.0/(b-accu(log(1-V.subvec(0,h-2))))));
}

/*calculates the number of obs in each group h for a given iteration */
// [[Rcpp::export]]
arma::ivec calc_allnh(arma::ivec zeta, int h, int n){
  arma::ivec out = arma::zeros<arma::ivec>(h);
  for(int i=0;i<n;i++){
    out[zeta[i]-1]++;
  }
  return(out);
}

/*calculates value needed to update V in gibbs step */

// [[Rcpp::export]]
int calc_nhp(int k, int h, arma::ivec all){
  int out;
  if(k < h-1){
    out = int(arma::accu(all.subvec(k+1,h-1)));
  }
  else{
    out = 0;
  }
  return(out);
}

/*used to find which observations are in cluster h for a given iteration,
 * just subsets full data -- fulldata[zeta==k]
 */
arma::mat calc_groupdata(int n, int k, arma::ivec zeta, arma::vec xee, arma::vec xes){
  arma::mat out(n,2);
  int k3 = 0;
  for(int k2=0;k2<n;k2++){
    if(zeta[k2]==k+1){
      out(k3,0) = xee[k2];
      out(k3,1) = xes[k2];
      k3++;
    }
  }
  
  return(out);
}

/*bivariate normal likelihood used to calc acceptance probability 
 * for (Xee, Xes) proposal
 */
double cpp_log_bqx(double meanyee, double varyee, double meanyes, double varyes,
                   double varwee, double varwes, arma::rowvec meanx, arma::mat varx,
                   arma::vec x, arma::rowvec yee, arma::rowvec yes, arma::rowvec wee, arma::rowvec wes, int nr){

  //int n = yee.size();
  //int h = pivector.size();
  double ywt = 0;
  double xt = 0;
  for(int i=0;i<nr;i++){
  ywt += R::dnorm(yee[i],meanyee,sqrt(varyee),true) + //R::dnorm(yee[1],meanyee,sqrt(varyee),true) +
    R::dnorm(yes[i],meanyes,sqrt(varyes),true) + //R::dnorm(yes[1],meanyes,sqrt(varyes),true) +
    R::dnorm(wee[i],x[0],sqrt(varwee),true)+ //R::dnorm(wee[1],x[0],sqrt(varwee),true) +
    R::dnorm(wes[i],x[1],sqrt(varwes),true);//+ R::dnorm(wes[1],x[1],sqrt(varwes),true);
  }
  //for(int j=0;j<h;j++){
  xt = dmvnrm_arma(x,meanx,varx,true);
  //std::cout << dmvnrm_arma(x,meanmat.row(j),varmat.slice(j)) << "\n";
  //}
  return(ywt+xt);
}

// double cpp_log_bqx(arma::vec meanyee, arma::vec varyee, arma::vec meanyes, arma::vec varyes,
//                    arma::vec varwee, arma::vec varwes,
//                    arma::mat pivector, arma::mat meanx, arma::cube varx,
//                    arma::mat x, arma::mat yee, arma::mat yes, arma::mat wee, arma::mat wes){
//   
//   int n = yee.n_rows;
//   //int h = pivector.size();
//   double ywt = 0;
//   double xt = 0;
//   for(int i=0;i<n;i++){
//     ywt += R::dnorm(yee(i,0),meanyee[i],sqrt(varyee[i]),true) + R::dnorm(yee(i,1),meanyee[i],sqrt(varyee[i]),true) +
//       R::dnorm(yes(i,0),meanyes[i],sqrt(varyes[i]),true) + R::dnorm(yes(i,1),meanyes[i],sqrt(varyes[i]),true) +
//       R::dnorm(wee(i,0),x(i,0),sqrt(varwee[i]),true)+ R::dnorm(wee(i,1),x(i,0),sqrt(varwee[i]),true) +
//       R::dnorm(wes(i,0),x(i,1),sqrt(varwes[i]),true)+ R::dnorm(wes(i,1),x(i,1),sqrt(varwes[i]),true);
//   
//     xt += sum(dmvnrm_arma(x,meanx.row(i),varx.slice(i),true));
//     
//   }
// 
//   return(ywt+xt);
// }


// /*used to calculate acceptance probability for spline updates, 
//  * normal likelihood since assuming y~N(f(x),sigma2)
//  */
double cpp_log_q(arma::vec pmean, double sigma2, arma::vec y){
  int n = y.size();
  return(-n*log(sqrt(sigma2)) - 1/(2.0*sigma2)*(accu(pow((y-pmean),2))));
}

/*calculates sample covariance matrix for vectors x and y, output is 
 * thus a 2x2 cov matrix
 */
// arma::mat cpp_cov(arma::vec x, arma::vec y, int n){
//   arma::mat out(2,2);
//   double xbar = mean(x);
//   double ybar = mean(y);
//   double covar = 0.0;
//   out(0,0) = var(x);
//   out(1,1) = var(y);
//   for(int i=0;i<n;i++){
//     covar += (x[i]-xbar)*(y[i]-ybar);
//   }
//   out(0,1) = out(1,0) = covar/(double(n)-1.0);
//   return(out);
// }

/* calculates a block diagnoal covariance matrix, must have even number of cols!,
 * block diagnonals are 2x2 for EE and ES
 */
arma::mat cpp_cov(arma::mat x){
  int n = x.n_rows;
  int m = x.n_cols;
  arma::mat out = arma::zeros(m,m);
  arma::rowvec xbar = mean(x,0);
  double covar = 0.0;
  for(int i=0;i<m;i++){
    out(i,i) = var(x.col(i));
    if(out(i,i) == 0){
       out(i,i) = 100;
       std::cout << "Variance equaled 0, col = " << i << "\n";
    }
    if(i % 2 == 0){
      covar = 0.0;
      for(int j=0;j<n;j++){
        covar += (x(j,i)-xbar[i])*(x(j,i+1)-xbar[i+1]);
      }
      out(i,i+1) = out(i+1,i) = covar/(double(n)-1.0);
    }
  }
  return(out);
}

/*indicates which indicies of the data are not allowed to be selected to
 * become a new not for one of 2 reasons:
 * 1. it already is a knot
 * 2. it is within 3 positions of a current knot
 */
// IntegerVector cpp_pick_indicies(NumericVector fulldata, NumericVector sel,int nsamp){
//   IntegerVector matched;
//   
//   matched = match(sel,fulldata);
//   //int n = fulldata.size();
//   int k = matched.size();
//   IntegerMatrix intmatrix(7,k);
//   IntegerVector out(k*7+6);
//   
//   
//   for(int i=0;i<3;i++){
//     intmatrix(i,_) = i+1+matched; 
//     intmatrix(i+3,_) = -i-1+matched;
//   }
//   intmatrix(6,_) = matched;
//   
//   for(int i=0;i<7;i++){
//     for(int j=0;j<k;j++){
//       if(intmatrix(i,j)>0 && intmatrix(i,j)<nsamp){
//         out[i*k+j]=intmatrix(i,j);
//       }
//       else{
//         out[i*k+j]=0;
//       }
//     }
//   }
//   out[k*7] = 0;
//   out[k*7+1] = nsamp-1;
//     out[k*7+2] = 1;
//   out[k*7+3] = nsamp-2;
//       out[k*7+4] = 2;
//   out[k*7+5] = nsamp-3;
//   return(unique(out));
// }

IntegerVector cpp_pick_indicies(NumericVector fulldata, NumericVector sel,int nsamp){
  int n = fulldata.size();
  int m = sel.size();
  //int nm = n-m;
  int count1 = 0;
  //arma::ivec out;
  IntegerVector res(m+16,0);
  //full = arma::sort(full);
  //sel = arma::sort(sel);
  
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      if(fulldata[i]==sel[j]){
        res[count1] = i;
        count1++;
        break;
      }
    }
  }
  res[m] = 0;
  res[m+1] = n-1;
  res[m+2] = 1;
  res[m+3] = n-2;
  res[m+4] = 2;
  res[m+5] = n-3;
  res[m+6] = 3;
  res[m+7] = n-4;
  res[m+8] = 4;
  res[m+9] = n-5;
  res[m+10] = 5;
  res[m+11] = n-6;
  res[m+12] = 6;
  res[m+13] = n-7;
  res[m+14] = 7;
  res[m+15] = n-8;
  //std::cout << unique(res) << "\n";
  return(unique(res));
}

/* calls my_bs() from R which is a slightly modified version of 
 * bs() from R splines library */
arma::mat call_my_bs(Function f, NumericVector x, NumericVector knots) {
  NumericMatrix res = f(x,knots=knots);
  return as<arma::mat>(res);
}

/*computes linear regression coefficients and predicted values using
those coefficients on the given x values, ie predict(m1) */
List my_lm(const arma::vec & y, const arma::mat & X) {
  
  int n = X.n_rows; //, k = X.n_cols;
  arma::vec intercept = arma::ones(n);
  arma::mat modelmat = join_rows(intercept,X);
    
  arma::colvec coef = arma::solve(modelmat, y); 
  //arma::colvec resid = y - X*coef; 
  arma::colvec preds = modelmat*coef;
  
  return List::create(Named("coefficients") = coef,
                      Named("preds")       = preds);
}

List my_lm_2(const arma::vec & y, const arma::mat & X) {
  
  //int n = X.n_rows; //, k = X.n_cols;
  //arma::vec intercept = arma::ones(n);
  //arma::mat modelmat = join_rows(intercept,X);
  
  arma::colvec coef = arma::solve(X, y); 
  //arma::colvec resid = y - X*coef; 
  arma::colvec preds = X*coef;
  
  return List::create(Named("coefficients") = coef,
                      Named("preds")       = preds);
}

List my_lm_qp(const arma::vec & y, const arma::mat & X, Function f){
  int k = X.n_cols;
  int n = X.n_rows;
  arma::vec coef(k);
  arma::vec preds(n);
  arma::rowvec dvec(k);
  arma::mat amat = arma::zeros(k-1,k);
  for(int i=0;i<k-1;i++){
    amat(i,i) = -1;
    amat(i,i+1) = 1;
  }
  
  dvec = y.t()*X;
  arma::vec bvec = arma::zeros(k-1);
  //std::cout << "chol decomp " << chol(trans(X)*X) << "\n"
  //          << "dvec = " << dvec << "\n" << "amat = " << amat << "\n";
  
  coef = as<arma::vec>(f(arma::inv(arma::chol(trans(X)*X)), dvec, amat.t(), bvec,0,true));
  preds = X*coef;
  
  return List::create(Named("coefficients") = coef,
                      Named("preds")       = preds);
}


/*removes the one knot chosen to be removed from the vector of knot 
 * locations */
arma::vec cpp_remove_one(arma::vec fullvec, double remove){
  int n = fullvec.size();
  int k = 0;
  arma::vec out(n-1);
  for(int i=0;i<n;i++){
    if(fullvec[i] != remove){
      out[k] = fullvec[i];
      k++;
    }
  }
  return(out);
}

//equivalent to x[-y] in R
arma::vec rmv_by_ind(arma::vec x, arma::vec y) {
  int nx = x.size();
  int ny= y.size();
  int count1 = ny - 1;
  int count2 = nx - ny - 1;
  int n = nx-ny;
  arma::vec res(n);
  x = arma::sort(x);
  y= arma::sort(y);
  for(int i=nx-1;i>=0;i--){
    if(count2 < 0 && i!=y[count1]){
      stop("Count 2 is < 0 which is not possible");
      break;
    }
    if(count1<0){
      res[count2]=x[i];
      count2--;
    }
    else{
      if(i==y[count1]){
        count1--;
      }
      else{
        res[count2]=x[i];
        count2--;
      }
    }
  }
  return(res);
}

/*---------------------------------------------------
 *-------------------------------------------
 *----------------------------------------------
 *----------------------------------------------------
 *---------------------------------------------------
 *--------------------------------------------------
 *---------------------------------------------*/


// [[Rcpp::export]]
List mcmc_full(
    const arma::mat yee,
    const arma::mat yes,
    const arma::mat wee,
    const arma::mat wes,
    const arma::mat Z,
    List initial_values,
    const List prior,
    const int nreps,
    const int burn,
    const int h,
    const int maxknot,
    Function spline,
    Function my_qp
){
  int n = yee.n_rows;
  const int l = 3; //number of cts derivatives, set at 3
  const int np = Z.n_cols;
  const int nr = yee.n_cols;

  arma::vec yeeb = mean(yee,1);
  arma::vec yesb = mean(yes,1);
  arma::vec yeebdiff(n);
  arma::vec yesbdiff(n);

  
  // Set initial values
  int currentkee                = as<int>(initial_values["currentkee"]);
  int currentkes                = as<int>(initial_values["currentkes"]);
  double ck                     = as<double>(initial_values["ck"]);
  arma::vec knotsee             = (as<arma::vec>(initial_values["knotsee"]));
  arma::vec knotses             = (as<arma::vec>(initial_values["knotses"]));
  arma::vec currentxee          = (as<arma::vec>(initial_values["currentxee"]));
  arma::vec currentxes          = (as<arma::vec>(initial_values["currentxes"]));
  arma::mat currentx(n,2);
  currentx.col(0)               = currentxee;
  currentx.col(1)               = currentxes;
  arma::mat currentxt           = currentx.t();

  arma::vec currentv            = (as<arma::vec>(initial_values["currentv"]));
  arma::vec currentpi           = (as<arma::vec>(initial_values["currentpi"]));
  double currentalpha           = (as<double>(initial_values["currentalpha"]));
  arma::ivec currentzeta         = (as<arma::ivec>(initial_values["currentzeta"]));
  arma::vec currentpredee       = (as<arma::vec>(initial_values["currentpredee"]));
  arma::vec currentpredes       = (as<arma::vec>(initial_values["currentpredes"]));
  arma::cube currentSigma(2,2,h); // = clone(arma::cube(initial_values["currentSigma"]));
  arma::cube tune(2,2,n);
  
  arma::vec currentmuee         = as<arma::vec>(initial_values["currentmuee"]);
  arma::vec currentmues         = as<arma::vec>(initial_values["currentmues"]);
  arma::mat currentmu(h,2);
  currentmu.col(0)              = currentmuee;
  currentmu.col(1)              = currentmues;
  arma::mat currentbetaee       = as<arma::mat>(initial_values["currentbetaee"]);
  arma::mat currentbetaes       = as<arma::mat>(initial_values["currentbetaes"]);
  

  double currentsigma2ee        = as<double>(initial_values["currentsigma2ee"]);
  double currentsigma2es        = as<double >(initial_values["currentsigma2es"]);
  double currentsigma2ve        = as<double >(initial_values["currentsigma2ve"]);
  double currentsigma2vs        = as<double >(initial_values["currentsigma2vs"]);

  arma::vec currentsigma2x      = as<arma::vec>(initial_values["currentsigma2x"]);
  arma::vec tunevar             = as<arma::vec>(initial_values["tunevar"]);
  double tunecor             = as<double>(initial_values["tunecor"]);
  

  for(int i=0;i<h;i++){
      currentSigma.slice(i).diag() = currentsigma2x;
      currentSigma(0,1,i) = currentSigma(1,0,i) = 0;
  }

  tune = arma::zeros(2,2,n);
  for(int i=0;i<n;i++){
    tune.slice(i).diag() = tunevar;
    tune(0,1,i)=tune(1,0,i)=sqrt(tunevar[0]*tunevar[1])*tunecor;
  }


  
  //allocate storage
  arma::mat latentxee(nreps-burn,n);
  arma::mat latentxes(nreps-burn,n);
  arma::mat muee(nreps-burn,h);
  arma::mat mues(nreps-burn,h);
  
  arma::mat sigma2xee(nreps-burn,h);
  arma::mat sigma2xes(nreps-burn,h);
  arma::mat corrx(nreps-burn,h);
  
  arma::mat pi(nreps-burn,h);
  arma::imat zeta(nreps-burn,n);
  arma::mat meanfcnee(nreps-burn,n);
  arma::mat meanfcnes(nreps-burn,n);
  arma::vec alpha(nreps-burn);

  //  predsee <- matrix(0,nrow=nreps-burn,ncol=3)
  //  predses <- matrix(0,nrow=nreps-burn,ncol=3)
  arma::vec sigma2eee(nreps-burn);
  arma::vec sigma2ees(nreps-burn);
  arma::vec sigma2vee(nreps-burn);
  arma::vec sigma2ves(nreps-burn);
  arma::vec kee(nreps-burn);
  arma::vec kes(nreps-burn);
  
  arma::mat betaee(nreps-burn,np);
  arma::mat betaes(nreps-burn,np);
  

  // Set prior values
  double lambda     = as<double>(prior["lambda"]   ); //prior for number of kntos
  double ae         = as<double>(prior["ae"]   ); //priors for sigmae
  double be         = as<double>(prior["be"]   );
  double av         = as<double>(prior["av"]   ); //priors for sigmav
  double bv         = as<double>(prior["bv"]   ); 
  //double c1    = as<double>(prior["c1"]   ); // tau ~ Ca+(0,c)
  //double d1    = as<double>(prior["d1"]   ); // tau ~ Ca+(0,c)
  double a_alp      = as<double>(prior["a_alp"]);
  double b_alp      = as<double>(prior["b_alp"]); // 
  double d          = as<double>(prior["d"]);   // 
  double Mb         = as<double>(prior["Mb"]);   // 
  double Vb         = as<double>(prior["Vb"]);   // 

  
  arma::vec m       = as<arma::vec>(prior["m"]   ); //prior mean for mus
  arma::mat v2      = as<arma::mat>(prior["v2"]   ); //prior var for mus
  arma::mat psi     = as<arma::mat>(prior["psi"]   ); //prior var for mus


  arma::ivec allnh(h); allnh.zeros();
  arma::vec xhee(n);
  arma::vec xhes(n);
  arma::mat xh(n,2);
  arma::vec xhb(2);
  arma::mat S(n,n);
  arma::vec mup(2);
  arma::mat v2p(2,2);
  arma::mat temp(1,2);
  arma::vec centeree(n);
  arma::vec centeres(n);
  arma::mat propx(n,2);
  int nhp;
  double nu;
  
  double acceptprob;
  double acceptprob1ee;
  double acceptprob2ee;
  double acceptprob3ee;
  double acceptprob1es;
  double acceptprob2es;
  double acceptprob3es;
  

  double bkee;
  double dkee;
  double bkes;
  double dkes;
  double u;
  arma::vec newee(1);
  arma::vec newes(1);
  arma::vec rmvee(1);
  arma::vec rmves(1);
  arma::vec chgee(1);
  arma::vec chges(1);
  arma::vec reduced;
  arma::vec compare(2);
  arma::vec knotsandxee;
  arma::vec knotsandxes(n+maxknot+3);
  arma::vec knotspropee(n+maxknot+3);
  arma::vec knotspropes(n+maxknot+3);
  arma::vec tempee(n);
  arma::vec tempes(n);
  IntegerVector unavailableee = rep(0,n);
  IntegerVector unavailablees = rep(0,n);
  arma::vec availableee;
  arma::vec availablees;
  arma::mat bmatrixee(n,maxknot);
  arma::mat bmatrixes(n,maxknot);
  List modelee;
  List modeles;
  List modelee_prop;
  List modeles_prop;
  List modelee_alt;
  List modeles_alt;
  arma::mat Vee(np,np);
  arma::mat Ves(np,np);
  arma::vec Mee(np);
  arma::vec Mes(np);
  arma::vec tempxee(n);
  arma::vec tempxes(n);
  arma::vec proppredee(n);
  arma::vec proppredes(n);
  
  double zbee;
  double zbes;
  
  double e2ee;
  double e2es;
  double v2ee;
  double v2es;
    
  arma::mat latentxee2(burn,n);
  arma::mat latentxes2(burn,n);
  arma::imat zeta2(burn,n);
  arma::mat acceptvec(nreps,n);
  arma::mat muee2(burn,h);
  arma::mat mues2(burn,h);
  arma::mat sigma2xee2(burn,h);
  arma::mat sigma2xes2(burn,h);
  
  //to check acceptance rates of spline step
  arma::vec accept_rate(2);
  accept_rate[0] = 0;
  accept_rate[1] = 0;
  
  //std::cout << "start\n";
  
  modelee = my_lm_qp(yeeb,call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotsee.begin(),knotsee.end())),my_qp);
  modeles = my_lm_qp(yesb,call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotses.begin(),knotses.end())),my_qp);
  
  //std::cout << "0\n";
  //RNGScope rngScope;
  //run mcmc 
  for(int i=0;i<nreps;i++){
    //sample zeta
    //std::cout << "1\n" << currentSigma << "\n";
    currentzeta = sample_zeta(currentx,currentpi,currentmu,currentSigma);
    allnh = calc_allnh(currentzeta,h,n);

    //loop to sample pi, Sigma, mus, alpha
    for(int k=0;k<h;k++){
      nhp = calc_nhp(k,h,allnh);
 
      //calc mean of group k
      if(allnh[k]==0){
        xhb = arma::zeros<arma::vec>(2);
        S = arma::zeros<arma::mat>(2,2) + psi;
      }
      else{
        //for every zeta k, make subgroup of data that is in group k
        xh = calc_groupdata(n, k, currentzeta, currentx.col(0), currentx.col(1));
        
        xhee = xh.submat(0,0,allnh[k]-1,0);//xh.col(0).rows(0,int(allnh[k])-1);
        xhes = xh.submat(0,1,allnh[k]-1,1);//xh.col(1).rows(0,int(allnh[k])-1);
        //only want first allnh[k] values for mean since vector is length n
        xhb[0] = mean(xhee.subvec(0,allnh[k]-1));
        xhb[1] = mean(xhes.subvec(0,allnh[k]-1));
        S = arma::trans(arma::join_rows(xhee.subvec(0,allnh[k]-1)-currentmu(k,0),xhes.subvec(0,allnh[k]-1)-currentmu(k,1)))*arma::join_rows(xhee.subvec(0,allnh[k]-1)-currentmu(k,0),xhes.subvec(0,allnh[k]-1)-currentmu(k,1)) + psi ;
      }
      
      //sample v
      if(k<h-1){
        currentv[k] = R::rbeta(1+allnh[k],currentalpha+nhp);
      }
      else{
        currentv[k] = 1;
      }
      //calc pi
      if(k>0){
        currentpi[k] = currentv[k]*cpp_prod(currentv.subvec(0,k-1));
      }
      else{
        currentpi[k] = currentv[k];
      }
      
      
      //std::cout << "2\n";
      nu = d + double(allnh[k]);
      //S = arma::trans(arma::join_rows(xhee-currentmu(k,0),xhes-currentmu(k,1)))*arma::join_rows(xhee-currentmu(k,0),xhes-currentmu(k,1)) + psi;
      currentSigma.slice(k) = rinvwish(1,nu,S);
      
      
      v2p = arma::inv(v2.i() + allnh[k]*arma::inv(currentSigma.slice(k)));
      mup = v2p*(v2.i()*m + arma::inv(currentSigma.slice(k))*xhb*allnh[k]);
      currentmu.row(k) = mvrnormArma(1,mup,v2p);
      
    }
    //currentalpha = sample_alpha(a_alp,b_alp,h,currentv);
    
    e2ee = 0;
    e2es = 0;
    v2ee = 0;
    v2es = 0;
    for(int ii=0;ii<nr;ii++){
      e2ee += accu(pow(yee.col(ii)-currentpredee-Z*currentbetaee.t(),2))/2.0;
      e2es += accu(pow(yes.col(ii)-currentpredes-Z*currentbetaes.t(),2))/2.0;
      v2ee += accu(pow(wee.col(ii)-currentx.col(0),2))/2.0;
      v2es += accu(pow(wes.col(ii)-currentx.col(1),2))/2.0;
    }
    currentsigma2ee = 1/R::rgamma(ae+nr*n/2.0,1.0/(be+e2ee));
    currentsigma2es = 1/R::rgamma(ae+nr*n/2.0,1.0/(be+e2es));
    currentsigma2ve = 1/R::rgamma(av+nr*n/2.0,1.0/(bv+v2ee));
    currentsigma2vs = 1/R::rgamma(av+nr*n/2.0,1.0/(bv+v2es));

    //std::cout << "3\n";
    //sample x
    for(int g=0;g<n;g++){
      //std::cout << "3a, iteration = " << i << "g = " << g << "\n";
      //std::cout << "i= " << i << "\n tuning = " << tune.slice(g) << "\n";
      propx.row(g) = mvrnormArma(1,currentx.row(g).t(),2.88*tune.slice(g));
      //std::cout << "3b\n";
      
      zbee = 0.0;
      zbes = 0.0;
      for(int gg=0;gg<np;gg++){
        zbee += Z(g,gg)*currentbetaee(0,gg);
        zbes += Z(g,gg)*currentbetaes(0,gg);
      }
      
      tempxee = currentx.col(0);
      tempxes = currentx.col(1);
      tempxee[g] = propx(g,0);
      tempxes[g] = propx(g,1);
       
      proppredee = call_my_bs(spline,NumericVector(tempxee.begin(),tempxee.end()),NumericVector(knotsee.begin(),knotsee.end()))*as<arma::vec>(modelee["coefficients"]); 
      proppredes = call_my_bs(spline,NumericVector(tempxes.begin(),tempxes.end()),NumericVector(knotses.begin(),knotses.end()))*as<arma::vec>(modeles["coefficients"]); 


      acceptprob = cpp_log_bqx(proppredee[g]+zbee,currentsigma2ee,proppredes[g]+zbes,currentsigma2es,currentsigma2ve,currentsigma2vs,currentmu.row(currentzeta[g]-1),currentSigma.slice(currentzeta[g]-1),trans(propx.row(g)),yee.row(g),yes.row(g),wee.row(g),wes.row(g),nr)-
        cpp_log_bqx(currentpredee[g]+zbee,currentsigma2ee,currentpredes[g]+zbes,currentsigma2es,currentsigma2ve,currentsigma2vs,currentmu.row(currentzeta[g]-1),currentSigma.slice(currentzeta[g]-1),trans(currentx.row(g)),yee.row(g),yes.row(g),wee.row(g),wes.row(g),nr);
      
      acceptvec(i,g) = acceptprob;
      
      if(log(R::runif(0,1)) < acceptprob){
        currentx.row(g) = propx.row(g);
      }

      if((i < burn) && (i > 300) && (i%20==0)){
        //tune.slice(g) = cpp_cov(latentxee.rows(0,i).col(g),latentxes.rows(0,i).col(g),i+1);
        tune.slice(g) = cpp_cov(arma::join_rows(latentxee.rows(0,i-1).col(g),latentxes.rows(0,i-1).col(g)));
      }
      //
      // if(i==burn+1){
      //   std::cout << sqrt(tune.slice(g)[0]) << " ";
      //   std::cout << sqrt(tune.slice(g)[3]) << " " << tune.slice(g)[1]/(sqrt(tune.slice(g)[0])*sqrt(tune.slice(g)[3])) << "\n";
      // 
      // }
    }
    
    currentxee = currentx.col(0);
    currentxes = currentx.col(1);
    
    //std::cout << "4\n";
    //spline step
    //EE
    compare = arma::zeros(2);
    compare[0] = 1;
    compare[1] = R::dpois(currentkee+1,lambda,false)/R::dpois(currentkee,lambda,false);
    bkee = ck*compare.min();
    compare[1] = R::dpois(currentkee,lambda,false)/R::dpois(currentkee+1,lambda,false);
    dkee = ck*compare.min();
    
    
     
    if((knotsee.min() < currentxee.min())){
      knotsee[which_min(NumericVector(knotsee.begin(),knotsee.end()))] = currentxee.min()+10;
      std::cout << "knot too small ee\n";
    } 
    if((knotsee.max() > currentxee.max())){
      knotsee[which_max(NumericVector(knotsee.begin(),knotsee.end()))] = currentxee.max()-10;
      std::cout << "knot too large ee\n";
    } 

    if((knotsee.min() < currentxee.min())){
      knotsee[which_min(NumericVector(knotsee.begin(),knotsee.end()))] = currentxee.min()+30;
      std::cout << "knot too small ee\n";
    } 
    if((knotsee.max() > currentxee.max())){
      knotsee[which_max(NumericVector(knotsee.begin(),knotsee.end()))] = currentxee.max()-30;
      std::cout << "knot too large ee\n";
    } 
    
    if((knotsee.min() < currentxee.min())){
      knotsee[which_min(NumericVector(knotsee.begin(),knotsee.end()))] = currentxee.min()+50;
      std::cout << "knot too small ee\n";
    } 
    if((knotsee.max() > currentxee.max())){
      knotsee[which_max(NumericVector(knotsee.begin(),knotsee.end()))] = currentxee.max()-50;
      std::cout << "knot too large ee\n";
    } 
    
    knotsandxee = arma::join_cols(currentxee,knotsee);
    knotsandxee = arma::sort(knotsandxee);
    
    //unavailableee = cpp_pick_indicies(as<NumericVector>(wrap(knotsandxee)),NumericVector(knotsee.begin(),knotsee.end()),n);    
    unavailableee = cpp_pick_indicies(as<NumericVector>(wrap(arma::sort(currentxee))),NumericVector(knotsee.begin(),knotsee.end()),n);    
    availableee = rmv_by_ind(arma::sort(currentxee),arma::sort(as<arma::vec>(unavailableee)));
        
    u = R::runif(0,1);

        
    //birth step
    if(u <= bkee || knotsee.size()==0){
      newee = cpp_sample(NumericVector(availableee.begin(),availableee.end()),NumericVector(rep(1.0/double(n-unavailableee.size()),n-unavailableee.size())));
      knotspropee = arma::sort(arma::join_cols(knotsee,newee));
      bmatrixee = call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotspropee.begin(),knotspropee.end()));                
      //modelee = my_lm_qp(yeeb,bmatrixee);

      modelee_prop = my_lm_qp(yeeb,bmatrixee,my_qp);

      acceptprob1ee = cpp_log_q(as<arma::vec>(modelee_prop["preds"]),currentsigma2ee/2.0,yeeb)-
        cpp_log_q(currentpredee,currentsigma2ee/2.0,yeeb) -0.5*log(n) +
        log(double(n) - (2.0*(double(l)+1.0) + (currentkee)*(2.0*double(l)+1.0))) - log(double(n));
                      //std::cout << "18a" << "\n";

      if(log(u) < acceptprob1ee){

        modelee = modelee_prop; //
        currentpredee = as<arma::vec>(modelee["preds"]);
        knotsee = knotspropee;
        
        accept_rate[0]++;
      }
      else{
        //modelee_alt = my_lm_qp(yeeb,call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotsee.begin(),knotsee.end())));
        //currentpredee = as<arma::vec>(modelee_alt["preds"]);
        currentpredee = call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotsee.begin(),knotsee.end()))*as<arma::vec>(modelee["coefficients"]);
        
      }
    }
    //death step
    else if(u <= bkee+dkee && knotsee.size()>1){
      rmvee = cpp_sample(NumericVector(knotsee.begin(),knotsee.end()),NumericVector(rep(1.0/double(knotsee.size()),knotsee.size())));
          //std::cout << "14b"<< "\n";

      knotspropee = arma::sort(cpp_remove_one(knotsee,rmvee[0]));
          //std::cout << "15b"<< "\n";

      bmatrixee = call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotspropee.begin(),knotspropee.end()));
                //std::cout << "16b"<< "\n";
               
                
      //modelee = my_lm_qp(yeeb,bmatrixee);
      modelee_prop = my_lm_qp(yeeb,bmatrixee,my_qp);


      acceptprob2ee = cpp_log_q(as<arma::vec>(modelee_prop["preds"]),currentsigma2ee/2.0,yeeb)-
        cpp_log_q(currentpredee,currentsigma2ee/2.0,yeeb) -0.5*log(n) -
        log(double(n) - (2.0*(double(l)+1.0) + (currentkee)*(2.0*double(l)+1.0))) + log(double(n));

          //std::cout << "18b"<< "\n";

      if(log(u) < acceptprob2ee){

        modelee = modelee_prop; //
        currentpredee = as<arma::vec>(modelee["preds"]);
        knotsee = knotspropee;
        
        accept_rate[0]++;
      }
      else{
        //modelee_alt = my_lm_qp(yeeb,call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotsee.begin(),knotsee.end())));
        //currentpredee = as<arma::vec>(modelee_alt["preds"]);
        currentpredee = call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotsee.begin(),knotsee.end()))*as<arma::vec>(modelee["coefficients"]);
        
      }

    }
    //move step
    else{
      chgee = cpp_sample(NumericVector(knotsee.begin(),knotsee.end()),NumericVector(rep(1.0/double(knotsee.size()),knotsee.size())));
      reduced = cpp_remove_one(knotsee,chgee[0]);
      knotspropee = arma::join_cols(reduced,cpp_sample(NumericVector(availableee.begin(),availableee.end()),rep(1.0/double(availableee.size()),availableee.size())));//sample(sort(currentxee)[-unavailableee],1)
      knotspropee = arma::sort(knotspropee); 
      bmatrixee = call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotspropee.begin(),knotspropee.end()));           
      //modelee = my_lm_qp(yeeb,bmatrixee);
      modelee_prop = my_lm_qp(yeeb,bmatrixee,my_qp);

      acceptprob3ee = cpp_log_q(as<arma::vec>(modelee_prop["preds"]),currentsigma2ee/2.0,yeeb)-
        cpp_log_q(currentpredee,currentsigma2ee/2.0,yeeb) -0.5*log(n);

      if(log(u) < acceptprob3ee){

        modelee = modelee_prop; //
        currentpredee = as<arma::vec>(modelee["preds"]);
        knotsee = knotspropee;
        
        accept_rate[0]++;
      }
      else{
        //modelee_alt = my_lm_qp(yeeb,call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotsee.begin(),knotsee.end())));
        //currentpredee = as<arma::vec>(modelee_alt["preds"]);
        currentpredee = call_my_bs(spline,NumericVector(currentxee.begin(),currentxee.end()),NumericVector(knotsee.begin(),knotsee.end()))*as<arma::vec>(modelee["coefficients"]);
        
      }
    }
    
    //std::cout << "5\n";
    //ES spline step
      compare = arma::zeros(2);
    compare[0] = 1;
    compare[1] = R::dpois(currentkes+1,lambda,false)/R::dpois(currentkes,lambda,false);
    bkes = ck*compare.min();
    compare[1] = R::dpois(currentkes,lambda,false)/R::dpois(currentkes+1,lambda,false);
    dkes = ck*compare.min();
    
    if((knotses.min() < currentxes.min())){
      knotses[which_min(NumericVector(knotses.begin(),knotses.end()))] = currentxes.min()+5;
      std::cout << "knot too small es\n";
    } 
    if((knotses.max() > currentxes.max())){
      knotses[which_max(NumericVector(knotses.begin(),knotses.end()))] = currentxes.max()-5;
      std::cout << "knot too large es\n";
    } 
    
    if((knotses.min() < currentxes.min())){
      knotses[which_min(NumericVector(knotses.begin(),knotses.end()))] = currentxes.min()+15;
      std::cout << "knot too small es\n";
    } 
    if((knotses.max() > currentxes.max())){
      knotses[which_max(NumericVector(knotses.begin(),knotses.end()))] = currentxes.max()-15;
      std::cout << "knot too large es\n";
    } 
    
    if((knotses.min() < currentxes.min())){
      knotses[which_min(NumericVector(knotses.begin(),knotses.end()))] = currentxes.min()+25;
      std::cout << "knot too small es\n";
    } 
    if((knotses.max() > currentxes.max())){
      knotses[which_max(NumericVector(knotses.begin(),knotses.end()))] = currentxes.max()-25;
      std::cout << "knot too large es\n";
    } 
    
    knotsandxes = arma::join_cols(currentxes,knotses);
    knotsandxes = arma::sort(knotsandxes);
    
    //unavailablees = cpp_pick_indicies(as<NumericVector>(wrap(knotsandxes)),NumericVector(knotses.begin(),knotses.end()),n);
    unavailablees = cpp_pick_indicies(as<NumericVector>(wrap(arma::sort(currentxes))),NumericVector(knotses.begin(),knotses.end()),n);    
    availablees = rmv_by_ind(arma::sort(currentxes),arma::sort(as<arma::vec>(unavailablees)));
    
    u = R::runif(0,1);
    //birth step
    if(u <= bkes || knotses.size()==0){
      newes = cpp_sample(NumericVector(availablees.begin(),availablees.end()),NumericVector(rep(1.0/double(n-unavailablees.size()),n-unavailablees.size())));
      knotspropes = arma::sort(arma::join_cols(knotses,newes));
      bmatrixes = call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotspropes.begin(),knotspropes.end()));
      //modeles = my_lm_qp(yesb,bmatrixes);
      modeles_prop = my_lm_qp(yesb,bmatrixes,my_qp);
      
 
      acceptprob1es = cpp_log_q(as<arma::vec>(modeles_prop["preds"]),currentsigma2es/2.0,yesb)-
        cpp_log_q(currentpredes,currentsigma2es/2.0,yesb) -0.5*log(n) +
        log(double(n) - (2.0*(double(l)+1.0) + (currentkes)*(2.0*double(l)+1.0))) - log(double(n));

      if(log(u) < acceptprob1es){
        modeles = modeles_prop; //
        currentpredes = as<arma::vec>(modeles["preds"]);
        knotses = knotspropes;
        
        accept_rate[1]++;
      }
      else{
        //modeles_alt = my_lm_qp(yesb,call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotses.begin(),knotses.end())));
        //currentpredes = as<arma::vec>(modeles_alt["preds"]);
        currentpredes = call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotses.begin(),knotses.end()))*as<arma::vec>(modeles["coefficients"]);
        
      }
    }
    //death step
    else if(u <= bkes+dkes && knotses.size()>1){
      rmves = cpp_sample(NumericVector(knotses.begin(),knotses.end()),NumericVector(rep(1.0/double(knotses.size()),knotses.size())));
      knotspropes = arma::sort(cpp_remove_one(knotses,rmves[0]));
      bmatrixes = call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotspropes.begin(),knotspropes.end()));              
      //modeles = my_lm_qp(yesb,bmatrixes);
      modeles_prop = my_lm_qp(yesb,bmatrixes,my_qp);
      

      acceptprob2es = cpp_log_q(as<arma::vec>(modeles_prop["preds"]),currentsigma2es/2.0,yesb)-
        cpp_log_q(currentpredes,currentsigma2es/2.0,yesb) -0.5*log(n) -
        log(double(n) - (2.0*(double(l)+1.0) + (currentkes)*(2.0*double(l)+1.0))) + log(double(n));

      if(log(u) < acceptprob2es){
        modeles = modeles_prop; //
        currentpredes = as<arma::vec>(modeles["preds"]);
        knotses = knotspropes;
        
        accept_rate[1]++;
      }
      else{
        //modeles_alt = my_lm_qp(yesb,call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotses.begin(),knotses.end())));
        //currentpredes = as<arma::vec>(modeles_alt["preds"]);
        currentpredes = call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotses.begin(),knotses.end()))*as<arma::vec>(modeles["coefficients"]);
        
      }

    }
    //move step
    else{
      chges = cpp_sample(NumericVector(knotses.begin(),knotses.end()),NumericVector(rep(1.0/double(knotses.size()),knotses.size())));
      reduced = cpp_remove_one(knotses,chges[0]);
      knotspropes = arma::join_cols(reduced,cpp_sample(NumericVector(availablees.begin(),availablees.end()),rep(1.0/double(availablees.size()),availablees.size())));//sample(sort(currentxee)[-unavailableee],1)
      knotspropes = arma::sort(knotspropes); 
      bmatrixes = call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotspropes.begin(),knotspropes.end()));           
      //modeles = my_lm_qp(yesb,bmatrixes);
      modeles_prop = my_lm_qp(yesb,bmatrixes,my_qp);
      

      acceptprob3es = cpp_log_q(as<arma::vec>(modeles_prop["preds"]),currentsigma2es/2.0,yesb)-
        cpp_log_q(currentpredes,currentsigma2es/2.0,yesb) -0.5*log(n);

      if(log(u) < acceptprob3es){
        modeles = modeles_prop; //
        currentpredes = as<arma::vec>(modeles["preds"]);
        knotses = knotspropes;
        
        accept_rate[1]++;
      }
      else{
        //modeles_alt = my_lm_qp(yesb,call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotses.begin(),knotses.end())));
        //currentpredes = as<arma::vec>(modeles_alt["preds"]);
        currentpredes = call_my_bs(spline,NumericVector(currentxes.begin(),currentxes.end()),NumericVector(knotses.begin(),knotses.end()))*as<arma::vec>(modeles["coefficients"]);
        
      }
    }
    
    
    yeebdiff = yeeb - currentpredee;
    yesbdiff = yesb - currentpredes; 
            //sample beta ee
    Vee = arma::inv((1/currentsigma2ee)*Z.t()*Z + (1/Vb)*arma::eye(np,np));
    Mee = Vee*((1/currentsigma2ee)*Z.t()*(yeebdiff)+Mb/Vb);
    currentbetaee.row(0) = mvrnormArma(1,Mee,Vee);

    //sample beta es
    Ves = arma::inv((1/currentsigma2es)*Z.t()*Z + (1/Vb)*arma::eye(np,np));
    Mes = Ves*((1/currentsigma2es)*Z.t()*(yesbdiff)+Mb/Vb);
    currentbetaes.row(0) = mvrnormArma(1,Mes,Ves);


    currentkee = knotsee.size();
    currentkes = knotses.size();
    
    if(i >= burn){
          //store samples
      kee[i-burn] = knotsee.size();  
      kes[i-burn] = knotses.size();
      meanfcnee.row(i-burn) = currentpredee.t();
      meanfcnes.row(i-burn) = currentpredes.t();
      
      sigma2eee[i-burn] = currentsigma2ee;
      sigma2ees[i-burn] = currentsigma2es;
      sigma2vee[i-burn] = currentsigma2ve;
      sigma2ves[i-burn] = currentsigma2vs;
      
      //predsee.rows(i,i) = currentpredee[1:3] ;
      //predses.rows(i,i) = currentpredes[1:3]; 
      //currentkee = kee[i-burn]; 
      //currentkes = kes[i-burn];
    
      zeta.row(i-burn) = currentzeta.t();
      pi.row(i-burn) = currentpi.t(); 
      muee.row(i-burn) = trans(currentmu.col(0));
      mues.row(i-burn) = trans(currentmu.col(1));
      
      betaee.row(i-burn) = currentbetaee.row(0);
      betaes.row(i-burn) = currentbetaes.row(0);
    
    
      for(int u=0;u<h;u++){
        sigma2xee(i-burn,u) = currentSigma.subcube(0,0,u,0,0,u)[0];
        sigma2xes(i-burn,u) = currentSigma.subcube(1,1,u,1,1,u)[0];
        corrx(i-burn,u) = currentSigma.subcube(0,1,u,0,1,u)[0]/(sqrt(sigma2xee(i-burn,u))*sqrt(sigma2xes(i-burn,u)));
      }

      latentxee.row(i-burn) = trans(currentx.col(0));
      latentxes.row(i-burn) = trans(currentx.col(1));
      
      alpha[i-burn] = currentalpha;

    }
    else{ //need to update using burnin to tune random walk, overide
    //once i>=burnin
      latentxee.row(i) = trans(currentx.col(0));
      latentxes.row(i) = trans(currentx.col(1));
      
      // latentxee2.row(i) = trans(currentx.col(0));
      // latentxes2.row(i) = trans(currentx.col(1));
      // zeta2.row(i) = currentzeta.t();
      // muee2.row(i) = trans(currentmu.col(0));
      // mues2.row(i) = trans(currentmu.col(1));
      // for(int u=0;u<h;u++){
      //   sigma2xee2(i,u) = currentSigma.subcube(0,0,u,0,0,u)[0];
      //   sigma2xes2(i,u) = currentSigma.subcube(1,1,u,1,1,u)[0];
      // }
      
    }
    
    // if(i % 1000==0){
    //     std::cout << "i= " << i << "\n";
    // }
    //std::cout << i << "\n";
  }

  
  return List::create(
    Named("latentxee")  = latentxee,
    Named("latentxes")  = latentxes,
    Named("muee")       = muee,
    Named("mues")       = mues,
    //Named("muee2")       = muee2,
    //Named("mues2")       = mues2,
    Named("sigma2xee")  = sigma2xee,
    Named("sigma2xes")  = sigma2xes,
    //Named("sigma2xee2")  = sigma2xee2,
    //Named("sigma2xes2")  = sigma2xes2,
    Named("corrx")      = corrx,
    Named("pi")         = pi,
    //Named("acceptprob")   = acceptvec,
    
    Named("zeta")       = zeta,
    //Named("zeta2")       = zeta2,
    Named("meanfcnee")  = meanfcnee,
    Named("meanfcnes")  = meanfcnes,
    Named("alpha")      = alpha,
    //Named("predsee")  = predsee,
    //Named("predses")  = predses,
    Named("sigma2eee")  = sigma2eee,
    Named("sigma2ees")  = sigma2ees,
    Named("sigma2vee")  = sigma2vee,
    Named("sigma2ves")  = sigma2ves,
    Named("kee")        = kee,
    Named("kes")        = kes,
    //Named("latentxee2")  = latentxee2,
    //Named("latentxes2")  = latentxes2,
    Named("betaee")     = betaee,
    Named("betaes")     = betaes);
    //Named("acceptance_rate") = accept_rate/double(nreps));
 

}

