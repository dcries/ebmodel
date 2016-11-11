//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <string>
//#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}



NumericVector quantileCpp(NumericVector x, NumericVector probs) {
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  int npr = probs.size();
  NumericVector ans(npr);
  for(int i=0; i<npr; i++){
    ans[i] = as<double>(quantile(x, probs[i]));
  }
return ans;
}


// [[Rcpp::export]]
List pp_check(List sample, arma::vec truth,
              arma::vec wsum, arma::vec ysum, arma::vec xee,
              arma::vec xes, arma::vec yee, arma::vec yes){
                
                //std::cout << "1a" <<"\n";
                arma::mat muee = as<arma::mat>(sample["muee"]);
                arma::mat mues = as<arma::mat>(sample["mues"]);
                arma::mat zeta = as<arma::mat>(sample["zeta"]);
                arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
                arma::mat meanfcnes = trans(as<arma::mat>(sample["meanfcnes"]));
                arma::mat latentxee = as<arma::mat>(sample["latentxee"]);
                arma::mat latentxes = as<arma::mat>(sample["latentxes"]);
                arma::mat sigma2xee = as<arma::mat>(sample["sigma2xee"]);
                arma::mat sigma2xes = as<arma::mat>(sample["sigma2xes"]);
                arma::vec sigma2eee = as<arma::vec>(sample["sigma2eee"]);
                arma::vec sigma2ees = as<arma::vec>(sample["sigma2ees"]);
                arma::vec sigma2vee = as<arma::vec>(sample["sigma2vee"]);
                arma::vec sigma2ves = as<arma::vec>(sample["sigma2ves"]);
                arma::mat corrx = as<arma::mat>(sample["corrx"]);
                
                //arma::mat mu(nreps,2);
                //mu.col(0) = muee;
                //mu.col(1) = mues;

                int n =  meanfcnee.n_rows;
                int nreps = sigma2eee.size();
                //std::cout << "1c" <<"\n";

                NumericVector quants = NumericVector::create(0.05,0.1,0.25,0.75,0.9,0.95);
                                //std::cout << "1d" <<"\n";

                //allocate storage
                arma::mat pmse(nreps,2);
                arma::mat checkx(nreps,truth.size());
                arma::mat checkw(nreps,wsum.size());
                arma::mat checky(nreps,wsum.size());
                arma::mat data(n,2);
                                //std::cout << "1e" <<"\n";

                arma::mat covmat(2,2);
                arma::vec meanvec = arma::zeros(2);
                arma::vec datawee(n);
                arma::vec datawes(n);
                arma::vec datayee(n);
                arma::vec datayes(n);
                arma::vec datasum(n);
                arma::vec a(n);
                int indx;
                
                //std::cout << 1 <<"\n";
                
                for(int i=0;i<nreps;i++){
                  for(int k=0;k<n;k++){
                    indx = zeta(i,k)-1;
                    meanvec[0] = muee(i,indx);
                    meanvec[1] = mues(i,indx);
                    covmat(0,0) = sigma2xee(i,indx);
                    covmat(1,1) = sigma2xes(i,indx);
                    covmat(1,0) = covmat(0,1) = sqrt(sigma2xee(i,indx)*sigma2xes(i,indx))*corrx(i,indx);
                    data.row(k) = mvrnormArma(1,meanvec,covmat);

                    datawee[k] = R::rnorm(latentxee(i,k),sqrt(sigma2vee[i]));
                    datawes[k] = R::rnorm(latentxes(i,k),sqrt(sigma2ves[i]));
                    datayee[k] = R::rnorm(meanfcnee(k,i),sqrt(sigma2eee[i]));
                    datayes[k] = R::rnorm(meanfcnes(k,i),sqrt(sigma2ees[i]));

                  }
                                  //std::cout << 3 <<"\n";
                  datasum = data.col(0);
                  a = data.col(1);

                  
                  checkx(i,0) = min(data.col(0));
                  checkx.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datasum.begin(),datasum.end()),quants));
                  checkx(i,7) = max(data.col(0));
                  checkx(i,8) = min(data.col(1));
                  checkx.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(a.begin(),a.end()),quants));
                  checkx(i,15) = max(data.col(1));
                                  //std::cout << 4 <<"\n";
                  checkw(i,0) = min(datawee);
                  checkw.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datawee.begin(),datawee.end()),quants));
                  checkw(i,7) = max(datawee);
                  checkw(i,8) = min(datawes);
                  checkw.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datawes.begin(),datawes.end()),quants));
                  checkw(i,15) = max(datawes);
                    //std::cout << 5 <<"\n";
                  checky(i,0) = min(datayee);
                  checky.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datayee.begin(),datayee.end()),quants));
                  checky(i,7) = max(datayee);
                  checky(i,8) = min(datayes);
                  checky.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datayes.begin(),datayes.end()),quants));
                  checky(i,15) = max(datayes);
                    //std::cout << 6 <<"\n";
                    
                  pmse(i,0) = mean(pow(meanfcnee.col(i)-yee,2));
                  pmse(i,1) = mean(pow(meanfcnes.col(i)-yes,2));


                  if(i % 1000==0){
                    std::cout << i << "\n";
                  }
                }

      
                return(List::create(
                  Named("pmse")     = pmse,
                  Named("checkx")   = checkx,
                  Named("checkw")   = checkw,
                  Named("checky")   = checky));

              }
              
// [[Rcpp::export]]
List pp_jags_dp(List sample, arma::vec truth, 
              arma::vec wsum, arma::vec ysum, arma::vec xee,
              arma::vec xes, arma::vec yee, arma::vec yes, 
              arma::vec zg, arma::vec zb, arma::vec za){
                
                arma::mat muee = as<arma::mat>(sample["muee"]);
                arma::mat mues = as<arma::mat>(sample["mues"]);
                arma::mat zeta = as<arma::mat>(sample["zeta"]);
                //arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
                //arma::mat meanfcnes = trans(as<arma::mat>(sample["meanfcnes"]));
                arma::mat latentxee = as<arma::mat>(sample["latentxee"]);
                arma::mat latentxes = as<arma::mat>(sample["latentxes"]);
                arma::mat sigmaxee = as<arma::mat>(sample["sigmaxee"]);
                arma::mat sigmaxes = as<arma::mat>(sample["sigmaxes"]);
                arma::vec sigmaeee = as<arma::vec>(sample["sigmaeee"]);
                arma::vec sigmaees = as<arma::vec>(sample["sigmaees"]);
                arma::vec sigmavee = as<arma::vec>(sample["sigmavee"]);
                arma::vec sigmaves = as<arma::vec>(sample["sigmaves"]);
                arma::mat corrx = as<arma::mat>(sample["corrx"]);

                arma::vec be0 = as<arma::vec>(sample["be0"]);
                arma::vec be1 = as<arma::vec>(sample["be1"]);
                arma::vec bs0 = as<arma::vec>(sample["bs0"]);
                arma::vec bs1 = as<arma::vec>(sample["bs1"]);
                arma::vec geg = as<arma::vec>(sample["geg"]);
                arma::vec geb = as<arma::vec>(sample["geb"]);
                arma::vec gea = as<arma::vec>(sample["gea"]);
                arma::vec gig = as<arma::vec>(sample["gig"]);
                arma::vec gib = as<arma::vec>(sample["gib"]);
                arma::vec gia = as<arma::vec>(sample["gia"]);


                //arma::mat mu(nreps,2);
                //mu.col(0) = muee;
                //mu.col(1) = mues;

                int n =  zeta.n_cols;
                int nreps = sigmaeee.size();
                //std::cout << "1c" <<"\n";

                NumericVector quants = NumericVector::create(0.05,0.1,0.25,0.75,0.9,0.95);
                                //std::cout << "1d" <<"\n";

                //allocate storage
                arma::mat pmse(nreps,2);
                arma::mat checkx(nreps,truth.size());
                arma::mat checkw(nreps,wsum.size());
                arma::mat checky(nreps,wsum.size());
                arma::mat data(n,2);

                arma::mat covmat(2,2);
                arma::vec meanvec = arma::zeros(2);
                arma::vec datawee(n);
                arma::vec datawes(n);
                arma::vec datayee(n);
                arma::vec datayes(n);
                arma::vec datasum(n);
                arma::vec a(n);
                int indx;

                arma::vec meanfcnee(n);
                arma::vec meanfcnes(n);
                
                //std::cout << 1 <<"\n";
                
                for(int i=0;i<nreps;i++){
                  //std::cout << "1g" <<"\n";
                  
                  meanfcnee = be0[i]+be1[i]*trans(latentxee.row(i)) + geg[i]*zg + geb[i]*zb + gea[i]*za;
                  meanfcnes = bs0[i]+bs1[i]*trans(latentxes.row(i)) + gig[i]*zg + gib[i]*zb + gia[i]*za;
                  for(int k=0;k<n;k++){
                    //std::cout << "1h" <<"\n";
                    
                    indx = zeta(i,k)-1;
                    meanvec[0] = muee(i,indx);
                    meanvec[1] = mues(i,indx);
                    covmat(0,0) = pow(sigmaxee(i,indx),2);
                    covmat(1,1) = pow(sigmaxes(i,indx),2);
                    covmat(1,0) = covmat(0,1) = sigmaxee(i,indx)*sigmaxes(i,indx)*corrx(i,indx);
                    data.row(k) = mvrnormArma(1,meanvec,covmat);
                    //std::cout << "1a" <<"\n";
                    
                    datawee[k] = R::rnorm(latentxee(i,k),sigmavee[i]);
                    datawes[k] = R::rnorm(latentxes(i,k),sigmaves[i]);
                    datayee[k] = R::rnorm(meanfcnee(k),sigmaeee[i]);
                    datayes[k] = R::rnorm(meanfcnes(k),sigmaees[i]);
                    //std::cout << "1b" <<"\n";
                    
                  }
                                  //std::cout << 3 <<"\n";
                                  datasum = data.col(0);
                  a = data.col(1);
                  
                  
                  checkx(i,0) = min(data.col(0));
                  checkx.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datasum.begin(),datasum.end()),quants));
                  checkx(i,7) = max(data.col(0));
                  checkx(i,8) = min(data.col(1));
                  checkx.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(a.begin(),a.end()),quants));
                  checkx(i,15) = max(data.col(1));
                  //std::cout << 4 <<"\n";
                  checkw(i,0) = min(datawee);
                  checkw.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datawee.begin(),datawee.end()),quants));
                  checkw(i,7) = max(datawee);
                  checkw(i,8) = min(datawes);
                  checkw.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datawes.begin(),datawes.end()),quants));
                  checkw(i,15) = max(datawes);
                    //std::cout << 5 <<"\n";
                  checky(i,0) = min(datayee);
                  checky.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datayee.begin(),datayee.end()),quants));
                  checky(i,7) = max(datayee);
                  checky(i,8) = min(datayes);
                  checky.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datayes.begin(),datayes.end()),quants));
                  checky(i,15) = max(datayes);
                    //std::cout << 6 <<"\n";
                    
                  pmse(i,0) = mean(pow(meanfcnee-yee,2));
                  pmse(i,1) = mean(pow(meanfcnes-yes,2));
                  //std::cout << 7 <<"\n";
                  

                  if(i % 1000==0){
                    std::cout << i << "\n";
                  }
                }

      
                return(List::create(
                  Named("pmse")     = pmse,
                  Named("checkx")   = checkx,
                  Named("checkw")   = checkw,
                  Named("checky")   = checky));

              }              
              
// [[Rcpp::export]]
List pp_jags_lm(List sample, arma::vec truth,
              arma::vec wsum, arma::vec ysum, arma::vec xee,
              arma::vec xes, arma::vec yee, arma::vec yes, 
              arma::vec zg, arma::vec zb, arma::vec za){
                
                //std::cout << "1a" <<"\n";
                arma::vec muee = as<arma::vec>(sample["muee"]);
                arma::vec mues = as<arma::vec>(sample["mues"]);
                //arma::mat zeta = as<arma::mat>(sample["zeta"]);
                //arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
                //arma::mat meanfcnes = trans(as<arma::mat>(sample["meanfcnes"]));
                arma::mat latentxee = as<arma::mat>(sample["latentxee"]);
                arma::mat latentxes = as<arma::mat>(sample["latentxes"]);
                arma::vec sigmaxee = as<arma::vec>(sample["sigmaxee"]);
                arma::vec sigmaxes = as<arma::vec>(sample["sigmaxes"]);
                arma::vec sigmaeee = as<arma::vec>(sample["sigmaeee"]);
                arma::vec sigmaees = as<arma::vec>(sample["sigmaees"]);
                arma::vec sigmavee = as<arma::vec>(sample["sigmavee"]);
                arma::vec sigmaves = as<arma::vec>(sample["sigmaves"]);
                arma::vec corrx = as<arma::vec>(sample["corrx"]);
                
                arma::vec be0 = as<arma::vec>(sample["be0"]);
                arma::vec be1 = as<arma::vec>(sample["be1"]);
                arma::vec bs0 = as<arma::vec>(sample["bs0"]);
                arma::vec bs1 = as<arma::vec>(sample["bs1"]);
                arma::vec geg = as<arma::vec>(sample["geg"]);
                arma::vec geb = as<arma::vec>(sample["geb"]);
                arma::vec gea = as<arma::vec>(sample["gea"]);
                arma::vec gig = as<arma::vec>(sample["gig"]);
                arma::vec gib = as<arma::vec>(sample["gib"]);
                arma::vec gia = as<arma::vec>(sample["gia"]);


                //arma::mat mu(nreps,2);
                //mu.col(0) = muee;
                //mu.col(1) = mues;

                int n =  latentxee.n_cols;
                int nreps = sigmaeee.size();
                //std::cout << "1c" <<"\n";

                NumericVector quants = NumericVector::create(0.05,0.1,0.25,0.75,0.9,0.95);
                                //std::cout << "1d" <<"\n";

                //allocate storage
                arma::mat pmse(nreps,2);
                arma::mat checkx(nreps,truth.size());
                arma::mat checkw(nreps,wsum.size());
                arma::mat checky(nreps,wsum.size());
                arma::mat data(n,2);
                                //std::cout << "1e" <<"\n";

                arma::mat covmat(2,2);
                arma::vec meanvec = arma::zeros(2);
                arma::vec datawee(n);
                arma::vec datawes(n);
                arma::vec datayee(n);
                arma::vec datayes(n);
                arma::vec datasum(n);
                arma::vec a(n);
                
                arma::vec meanfcnee(n);
                arma::vec meanfcnes(n);
                
                //std::cout << 1 <<"\n";
                
                for(int i=0;i<nreps;i++){
                  meanfcnee = be0[i]+be1[i]*trans(latentxee.row(i)) + geg[i]*zg + geb[i]*zb + gea[i]*za;
                  meanfcnes = bs0[i]+bs1[i]*trans(latentxes.row(i)) + gig[i]*zg + gib[i]*zb + gia[i]*za;
                  meanvec[0] = muee[i];
                  meanvec[1] = mues[i];
                  covmat(0,0) = pow(sigmaxee[i],2);
                  covmat(1,1) = pow(sigmaxes[i],2);
                  covmat(1,0) = covmat(0,1) = sigmaxee[i]*sigmaxes[i]*corrx[i];
                  
                  for(int k=0;k<n;k++){
                    data.row(k) = mvrnormArma(1,meanvec,covmat);

                    datawee[k] = R::rnorm(latentxee(i,k),sigmavee[i]);
                    datawes[k] = R::rnorm(latentxes(i,k),sigmaves[i]);
                    datayee[k] = R::rnorm(meanfcnee[k],sigmaeee[i]);
                    datayes[k] = R::rnorm(meanfcnes[k],sigmaees[i]);

                  }

                  
                                  //std::cout << 3 <<"\n";
                                  datasum = data.col(0);
                  a = data.col(1);
                  
                  
                  checkx(i,0) = min(data.col(0));
                  checkx.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datasum.begin(),datasum.end()),quants));
                  checkx(i,7) = max(data.col(0));
                  checkx(i,8) = min(data.col(1));
                  checkx.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(a.begin(),a.end()),quants));
                  checkx(i,15) = max(data.col(1));
                  checkw(i,0) = min(datawee);
                  checkw.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datawee.begin(),datawee.end()),quants));
                  checkw(i,7) = max(datawee);
                  checkw(i,8) = min(datawes);
                  checkw.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datawes.begin(),datawes.end()),quants));
                  checkw(i,15) = max(datawes);
                    //std::cout << 5 <<"\n";
                  checky(i,0) = min(datayee);
                  checky.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datayee.begin(),datayee.end()),quants));
                  checky(i,7) = max(datayee);
                  checky(i,8) = min(datayes);
                  checky.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datayes.begin(),datayes.end()),quants));
                  checky(i,15) = max(datayes);
                    //std::cout << 6 <<"\n";
                    
                  pmse(i,0) = mean(pow(meanfcnee-yee,2));
                  pmse(i,1) = mean(pow(meanfcnes-yes,2));


                  if(i % 1000==0){
                    std::cout << i << "\n";
                  }
                }

      
                return(List::create(
                  Named("pmse")     = pmse,
                  Named("checkx")   = checkx,
                  Named("checkw")   = checkw,
                  Named("checky")   = checky));

              }              
              
              
              // [[Rcpp::export]]
List pp_jags_lm2(List sample, arma::vec truth, 
              arma::vec wsum, arma::vec ysum, arma::vec xee,
              arma::vec xes, arma::vec yee, arma::vec yes){
                
                arma::vec muee = as<arma::vec>(sample["muee"]);
                arma::vec mues = as<arma::vec>(sample["mues"]);
                //arma::mat zeta = as<arma::mat>(sample["zeta"]);
                //arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
                //arma::mat meanfcnes = trans(as<arma::mat>(sample["meanfcnes"]));
                arma::mat latentxee = as<arma::mat>(sample["latentxee"]);
                arma::mat latentxes = as<arma::mat>(sample["latentxes"]);
                arma::vec sigmaxee = as<arma::vec>(sample["sigmaxee"]);
                arma::vec sigmaxes = as<arma::vec>(sample["sigmaxes"]);
                arma::vec sigmaeee = as<arma::vec>(sample["sigmaeee"]);
                arma::vec sigmaees = as<arma::vec>(sample["sigmaees"]);
                arma::vec sigmavee = as<arma::vec>(sample["sigmavee"]);
                arma::vec sigmaves = as<arma::vec>(sample["sigmaves"]);
                arma::vec corrx = as<arma::vec>(sample["corrx"]);

                arma::vec be0 = as<arma::vec>(sample["be0"]);
                arma::vec be1 = as<arma::vec>(sample["be1"]);
                arma::vec bs0 = as<arma::vec>(sample["bs0"]);
                arma::vec bs1 = as<arma::vec>(sample["bs1"]);


                //arma::mat mu(nreps,2);
                //mu.col(0) = muee;
                //mu.col(1) = mues;

                int n =  latentxee.n_cols;
                int nreps = sigmaeee.size();
                //std::cout << "1c" <<"\n";

                NumericVector quants = NumericVector::create(0.05,0.1,0.25,0.75,0.9,0.95);
                                //std::cout << "1d" <<"\n";

                //allocate storage
                arma::mat pmse(nreps,2);
                arma::mat checkx(nreps,truth.size());
                arma::mat checkw(nreps,wsum.size());
                arma::mat checky(nreps,wsum.size());
                arma::mat data(n,2);
                                //std::cout << "1e" <<"\n";

                arma::mat covmat(2,2);
                arma::vec meanvec = arma::zeros(2);
                arma::vec datawee(n);
                arma::vec datawes(n);
                arma::vec datayee(n);
                arma::vec datayes(n);
                arma::vec datasum(n);
                arma::vec a(n);

                arma::vec meanfcnee(n);
                arma::vec meanfcnes(n);
                

                for(int i=0;i<nreps;i++){
                  meanfcnee = be0[i]+be1[i]*trans(latentxee.row(i)) ;
                  meanfcnes = bs0[i]+bs1[i]*trans(latentxes.row(i));
                  meanvec[0] = muee[i];
                  meanvec[1] = mues[i];
                  covmat(0,0) = pow(sigmaxee[i],2);
                  covmat(1,1) = pow(sigmaxes[i],2);
                  covmat(1,0) = covmat(0,1) = sigmaxee[i]*sigmaxes[i]*corrx[i];

                  for(int k=0;k<n;k++){
                    data.row(k) = mvrnormArma(1,meanvec,covmat);

                    datawee[k] = R::rnorm(latentxee(i,k),sigmavee[i]);
                    datawes[k] = R::rnorm(latentxes(i,k),sigmaves[i]);
                    datayee[k] = R::rnorm(meanfcnee[k],sigmaeee[i]);
                    datayes[k] = R::rnorm(meanfcnes[k],sigmaees[i]);

                  }

                  //std::cout << "1c" <<"\n";
                  
                                  //std::cout << 3 <<"\n";
                                  datasum = data.col(0);
                  a = data.col(1);
                  
                  
                  checkx(i,0) = min(data.col(0));
                  checkx.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datasum.begin(),datasum.end()),quants));
                  checkx(i,7) = max(data.col(0));
                  checkx(i,8) = min(data.col(1));
                  checkx.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(a.begin(),a.end()),quants));
                  checkx(i,15) = max(data.col(1));
                  //std::cout << 4 <<"\n";
                  checkw(i,0) = min(datawee);
                  checkw.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datawee.begin(),datawee.end()),quants));
                  checkw(i,7) = max(datawee);
                  checkw(i,8) = min(datawes);
                  checkw.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datawes.begin(),datawes.end()),quants));
                  checkw(i,15) = max(datawes);
                    //std::cout << 5 <<"\n";
                  checky(i,0) = min(datayee);
                  checky.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datayee.begin(),datayee.end()),quants));
                  checky(i,7) = max(datayee);
                  checky(i,8) = min(datayes);
                  checky.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datayes.begin(),datayes.end()),quants));
                  checky(i,15) = max(datayes);
                    //std::cout << 6 <<"\n";
                    
                  pmse(i,0) = mean(pow(meanfcnee-yee,2));
                  pmse(i,1) = mean(pow(meanfcnes-yes,2));


                  if(i % 1000==0){
                    std::cout << i << "\n";
                  }
                }

      
                return(List::create(
                  Named("pmse")     = pmse,
                  Named("checkx")   = checkx,
                  Named("checkw")   = checkw,
                  Named("checky")   = checky));

              }              

// [[Rcpp::export]]
List pp_full(List sample, arma::vec truth,
              arma::vec wsum, arma::vec ysum, arma::vec xee,
              arma::vec xes, arma::vec yee, arma::vec yes,
              arma::mat Z){
  
  //std::cout << "1a" <<"\n";
  arma::mat muee = as<arma::mat>(sample["muee"]);
  arma::mat mues = as<arma::mat>(sample["mues"]);
  arma::mat zeta = as<arma::mat>(sample["zeta"]);
  arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
  arma::mat meanfcnes = trans(as<arma::mat>(sample["meanfcnes"]));
  arma::mat latentxee = as<arma::mat>(sample["latentxee"]);
  arma::mat latentxes = as<arma::mat>(sample["latentxes"]);
  arma::mat sigma2xee = as<arma::mat>(sample["sigma2xee"]);
  arma::mat sigma2xes = as<arma::mat>(sample["sigma2xes"]);
  arma::vec sigma2eee = as<arma::vec>(sample["sigma2eee"]);
  arma::vec sigma2ees = as<arma::vec>(sample["sigma2ees"]);
  arma::vec sigma2vee = as<arma::vec>(sample["sigma2vee"]);
  arma::vec sigma2ves = as<arma::vec>(sample["sigma2ves"]);
  arma::mat corrx     = as<arma::mat>(sample["corrx"]);
  arma::mat betaee       = as<arma::mat>(sample["betaee"]);
  arma::mat betaes       = as<arma::mat>(sample["betaes"]);

  
  
  //arma::mat mu(nreps,2);
  //mu.col(0) = muee;
  //mu.col(1) = mues;
  
  int n =  meanfcnee.n_rows;
  int nreps = sigma2eee.size();
  int p = Z.n_cols;
  //std::cout << "1c" <<"\n";
  
  NumericVector quants = NumericVector::create(0.05,0.1,0.25,0.75,0.9,0.95);
  //std::cout << "1d" <<"\n";
  
  //allocate storage
  arma::mat pmse(nreps,2);
  arma::mat checkx(nreps,truth.size());
  arma::mat checkw(nreps,wsum.size());
  arma::mat checky(nreps,wsum.size());
  arma::mat data(n,2);
  //std::cout << "1e" <<"\n";
  
  arma::mat covmat(2,2);
  arma::vec meanvec = arma::zeros(2);
  arma::vec datawee(n);
  arma::vec datawes(n);
  arma::vec datayee(n);
  arma::vec datayes(n);
  arma::vec datasum(n);
  arma::vec a(n);
  int indx;
  double regmeanee;
  double regmeanes;
  
  //std::cout << 1 <<"\n";
  
  for(int i=0;i<nreps;i++){
    for(int k=0;k<n;k++){
      regmeanee = 0;
      regmeanes = 0;
      indx = zeta(i,k)-1;
      meanvec[0] = muee(i,indx);
      meanvec[1] = mues(i,indx);
      covmat(0,0) = sigma2xee(i,indx);
      covmat(1,1) = sigma2xes(i,indx);
      covmat(1,0) = covmat(0,1) = sqrt(sigma2xee(i,indx)*sigma2xes(i,indx))*corrx(i,indx);
      data.row(k) = mvrnormArma(1,meanvec,covmat);
      
      for(int l=0;l<p;l++){
        regmeanee += Z(k,l)*betaee(i,l);
        regmeanes += Z(k,l)*betaes(i,l);
      }
      
      datawee[k] = R::rnorm(latentxee(i,k),sqrt(sigma2vee[i]));
      datawes[k] = R::rnorm(latentxes(i,k),sqrt(sigma2ves[i]));
      datayee[k] = R::rnorm(meanfcnee(k,i)+regmeanee,sqrt(sigma2eee[i]));
      datayes[k] = R::rnorm(meanfcnes(k,i)+regmeanes,sqrt(sigma2ees[i]));
      
    }
    //std::cout << 3 <<"\n";
    datasum = data.col(0);
    a = data.col(1);
    
    
    checkx(i,0) = min(data.col(0));
    checkx.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datasum.begin(),datasum.end()),quants));
    checkx(i,7) = max(data.col(0));
    checkx(i,8) = min(data.col(1));
    checkx.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(a.begin(),a.end()),quants));
    checkx(i,15) = max(data.col(1));
    //std::cout << 4 <<"\n";
    checkw(i,0) = min(datawee);
    checkw.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datawee.begin(),datawee.end()),quants));
    checkw(i,7) = max(datawee);
    checkw(i,8) = min(datawes);
    checkw.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datawes.begin(),datawes.end()),quants));
    checkw(i,15) = max(datawes);
    //std::cout << 5 <<"\n";
    checky(i,0) = min(datayee);
    checky.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datayee.begin(),datayee.end()),quants));
    checky(i,7) = max(datayee);
    checky(i,8) = min(datayes);
    checky.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datayes.begin(),datayes.end()),quants));
    checky(i,15) = max(datayes);
    //std::cout << 6 <<"\n";
    
    pmse(i,0) = mean(pow(meanfcnee.col(i)+Z*trans(betaee.row(i))-yee,2));
    pmse(i,1) = mean(pow(meanfcnes.col(i)+Z*trans(betaes.row(i))-yes,2));
    
    
    if(i % 1000==0){
      std::cout << i << "\n";
    }
  }
  
  
  return(List::create(
      Named("pmse")     = pmse,
      Named("checkx")   = checkx,
      Named("checkw")   = checkw,
      Named("checky")   = checky));
  
}


// [[Rcpp::export]]
List pp_bvn(List sample, arma::vec truth,
             arma::vec wsum, arma::vec ysum, arma::vec xee,
             arma::vec xes, arma::vec yee, arma::vec yes,
             arma::mat Z){
  
  //std::cout << "1a" <<"\n";
  arma::vec muee = as<arma::vec>(sample["muee"]);
  arma::vec mues = as<arma::vec>(sample["mues"]);
  arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
  arma::mat meanfcnes = trans(as<arma::mat>(sample["meanfcnes"]));
  arma::mat latentxee = as<arma::mat>(sample["latentxee"]);
  arma::mat latentxes = as<arma::mat>(sample["latentxes"]);
  arma::vec sigma2xee = as<arma::vec>(sample["sigma2xee"]);
  arma::vec sigma2xes = as<arma::vec>(sample["sigma2xes"]);
  arma::vec sigma2eee = as<arma::vec>(sample["sigma2eee"]);
  arma::vec sigma2ees = as<arma::vec>(sample["sigma2ees"]);
  arma::vec sigma2vee = as<arma::vec>(sample["sigma2vee"]);
  arma::vec sigma2ves = as<arma::vec>(sample["sigma2ves"]);
  arma::vec corrx     = as<arma::vec>(sample["corrx"]);
  arma::mat betaee       = as<arma::mat>(sample["betaee"]);
  arma::mat betaes       = as<arma::mat>(sample["betaes"]);
  
  
  
  //arma::mat mu(nreps,2);
  //mu.col(0) = muee;
  //mu.col(1) = mues;
  
  int n =  meanfcnee.n_rows;
  int nreps = sigma2eee.size();
  int p = Z.n_cols;
  //std::cout << "1c" <<"\n";
  
  NumericVector quants = NumericVector::create(0.05,0.1,0.25,0.75,0.9,0.95);
  //std::cout << "1d" <<"\n";
  
  //allocate storage
  arma::mat pmse(nreps,2);
  arma::mat checkx(nreps,truth.size());
  arma::mat checkw(nreps,wsum.size());
  arma::mat checky(nreps,wsum.size());
  arma::mat data(n,2);
  //std::cout << "1e" <<"\n";
  
  arma::mat covmat(2,2);
  arma::vec meanvec = arma::zeros(2);
  arma::vec datawee(n);
  arma::vec datawes(n);
  arma::vec datayee(n);
  arma::vec datayes(n);
  arma::vec datasum(n);
  arma::vec a(n);
  double regmeanee;
  double regmeanes;
  
  //std::cout << 1 <<"\n";
  
  for(int i=0;i<nreps;i++){
    for(int k=0;k<n;k++){
      regmeanee = 0;
      regmeanes = 0;
      meanvec[0] = muee[i];
      meanvec[1] = mues[i];
      covmat(0,0) = sigma2xee[i];
      covmat(1,1) = sigma2xes[i];
      covmat(1,0) = covmat(0,1) = sqrt(sigma2xee[i]*sigma2xes[i])*corrx[i];
      data.row(k) = mvrnormArma(1,meanvec,covmat);
      
      for(int l=0;l<p;l++){
        regmeanee += Z(k,l)*betaee(i,l);
        regmeanes += Z(k,l)*betaes(i,l);
      }
      
      datawee[k] = R::rnorm(latentxee(i,k),sqrt(sigma2vee[i]));
      datawes[k] = R::rnorm(latentxes(i,k),sqrt(sigma2ves[i]));
      datayee[k] = R::rnorm(meanfcnee(k,i)+regmeanee,sqrt(sigma2eee[i]));
      datayes[k] = R::rnorm(meanfcnes(k,i)+regmeanes,sqrt(sigma2ees[i]));
      
    }
    //std::cout << 3 <<"\n";
    datasum = data.col(0);
    a = data.col(1);
    
    
    checkx(i,0) = min(data.col(0));
    checkx.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datasum.begin(),datasum.end()),quants));
    checkx(i,7) = max(data.col(0));
    checkx(i,8) = min(data.col(1));
    checkx.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(a.begin(),a.end()),quants));
    checkx(i,15) = max(data.col(1));
    //std::cout << 4 <<"\n";
    checkw(i,0) = min(datawee);
    checkw.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datawee.begin(),datawee.end()),quants));
    checkw(i,7) = max(datawee);
    checkw(i,8) = min(datawes);
    checkw.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datawes.begin(),datawes.end()),quants));
    checkw(i,15) = max(datawes);
    //std::cout << 5 <<"\n";
    checky(i,0) = min(datayee);
    checky.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datayee.begin(),datayee.end()),quants));
    checky(i,7) = max(datayee);
    checky(i,8) = min(datayes);
    checky.row(i).cols(9,14) = as<arma::rowvec>(quantileCpp(NumericVector(datayes.begin(),datayes.end()),quants));
    checky(i,15) = max(datayes);
    //std::cout << 6 <<"\n";
    
    pmse(i,0) = mean(pow(meanfcnee.col(i)+Z*trans(betaee.row(i))-yee,2));
    pmse(i,1) = mean(pow(meanfcnes.col(i)+Z*trans(betaes.row(i))-yes,2));
    
    
    if(i % 1000==0){
      std::cout << i << "\n";
    }
  }
  
  
  return(List::create(
      Named("pmse")     = pmse,
      Named("checkx")   = checkx,
      Named("checkw")   = checkw,
      Named("checky")   = checky));
  
}


// [[Rcpp::export]]
List pp_uni(List sample, arma::vec truth,
             arma::vec wsum, arma::vec ysum, arma::vec xee,
             arma::vec yee,
             arma::mat Z){
  
  //std::cout << "1a" <<"\n";
  arma::mat muee = as<arma::mat>(sample["muee"]);
  arma::mat zeta = as<arma::mat>(sample["zeta"]);
  arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
  arma::mat latentxee = as<arma::mat>(sample["latentxee"]);
  arma::mat sigma2xee = as<arma::mat>(sample["sigma2xee"]);
  arma::vec sigma2eee = as<arma::vec>(sample["sigma2eee"]);
  arma::vec sigma2vee = as<arma::vec>(sample["sigma2vee"]);
  arma::mat betaee       = as<arma::mat>(sample["betaee"]);

  
  
  //arma::mat mu(nreps,2);
  //mu.col(0) = muee;
  //mu.col(1) = mues;
  
  int n =  meanfcnee.n_rows;
  int nreps = sigma2eee.size();
  int p = Z.n_cols;
  //std::cout << "1c" <<"\n";
  
  NumericVector quants = NumericVector::create(0.05,0.1,0.25,0.75,0.9,0.95);
  //std::cout << "1d" <<"\n";
  
  //allocate storage
  arma::vec pmse(nreps);
  arma::mat checkx(nreps,truth.size());
  arma::mat checkw(nreps,wsum.size());
  arma::mat checky(nreps,wsum.size());
  arma::vec data(n);
  //std::cout << "1e" <<"\n";
  
  arma::vec datawee(n);
  arma::vec datayee(n);
  arma::vec datasum(n);
  arma::vec a(n);
  int indx;
  double regmeanee;

  //std::cout << 1 <<"\n";
  
  for(int i=0;i<nreps;i++){
    for(int k=0;k<n;k++){
      regmeanee = 0;
      indx = zeta(i,k)-1;
      data[k] = R::rnorm(muee(i,indx),sqrt(sigma2xee(i,indx)));
      
      for(int l=0;l<p;l++){
        regmeanee += Z(k,l)*betaee(i,l);
      }
      
      datawee[k] = R::rnorm(latentxee(i,k),sqrt(sigma2vee[i]));
      datayee[k] = R::rnorm(meanfcnee(k,i)+regmeanee,sqrt(sigma2eee[i]));

    }
    //std::cout << 3 <<"\n";

    
    checkx(i,0) = min(data);
    checkx.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(data.begin(),data.end()),quants));
    checkx(i,7) = max(data);
    //std::cout << 4 <<"\n";
    checkw(i,0) = min(datawee);
    checkw.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datawee.begin(),datawee.end()),quants));
    checkw(i,7) = max(datawee);
    //std::cout << 5 <<"\n";
    checky(i,0) = min(datayee);
    checky.row(i).cols(1,6) = as<arma::rowvec>(quantileCpp(NumericVector(datayee.begin(),datayee.end()),quants));
    checky(i,7) = max(datayee);
    //std::cout << 6 <<"\n";
    
    pmse[i] = mean(pow(meanfcnee.col(i)+Z*trans(betaee.row(i))-yee,2));

    
    if(i % 1000==0){
      std::cout << i << "\n";
    }
  }
  
  
  return(List::create(
      Named("pmse")     = pmse,
      Named("checkx")   = checkx,
      Named("checkw")   = checkw,
      Named("checky")   = checky));
  
}

// [[Rcpp::export]]
arma::mat pp_pmse(List sample,
            arma::vec yee, arma::vec yes,
            arma::mat Z){
  
  arma::mat meanfcnee = trans(as<arma::mat>(sample["meanfcnee"]));
  arma::mat meanfcnes = trans(as<arma::mat>(sample["meanfcnes"]));
  arma::mat betaee       = as<arma::mat>(sample["betaee"]);
  arma::mat betaes       = as<arma::mat>(sample["betaes"]);
  
  int nreps = betaee.n_rows;
  //std::cout << "1c" <<"\n";
  
  //allocate storage
  arma::mat pmse(nreps,2);
  //std::cout << "1e" <<"\n";

  
  for(int i=0;i<nreps;i++){
    
    pmse(i,0) = mean(pow(meanfcnee.col(i)+Z*trans(betaee.row(i))-yee,2));
    pmse(i,1) = mean(pow(meanfcnes.col(i)+Z*trans(betaes.row(i))-yes,2));
    
    
    if(i % 1000==0){
      std::cout << i << "\n";
    }
  }
  return(pmse);
}
