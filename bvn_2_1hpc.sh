#!/bin/bash

module load R
R CMD BATCH ebmodel/bvn_sim_2_1.R

# Instructions for new hpc-class users:
#  To use this script:
#   1) Save this script as a file named myjob on hpc-class
#   2) On hpc-class, Issue                   
#       qsub myjob    to submit the job 
#        Use qstat -a to see job status, 
#         Use qdel jobname to delete one of your jobs
#         jobnames are of the form 1234.hpc-class 

# This script has cd  as the first command. 
# qsub command was executed. This is what most users want. Change
# that command if you want something else.
###########################################
# Output goes to file BATCH_OUTPUT.
# Error output goes to file BATCH_ERRORS.
# If you want the output to go to another file, change BATCH_OUTPUT 
# or BATCH_ERRORS in the following lines to the full path of that file. 

#PBS  -o bvn_2_1  
#PBS  -e bvn_2_1_errors 

#PBS -lnodes=1:ppn=1:compute,walltime=48:00:00

# Change to directory from which qsub was executed 
   cd $PBS_O_WORKDIR

  