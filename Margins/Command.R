#########################################################
##                                                     ##
## This code launches the jobs in SubJobs_eCDF_excmu.R ##
##                                                     ##
#########################################################
module load R/3.4.2/intel-2017

############
### Test ###
############
Rscript SubJobs_eCDF_excmu.R "Replic <- 1" "ncpus <- 1" "jobtime <- '00:45:00'" # the first location
Rscript SubJobs_eCDF_excthresh.R "Replic <- 1" "ncpus <- 1"

###########
### Run ###
###########
Rscript SubJobs_eCDF_excmu.R "Replic <- 16703" "ncpus <- 16"
Rscript SubJobs_eCDF_excthresh.R "Replic <- 16703" "ncpus <- 20"
