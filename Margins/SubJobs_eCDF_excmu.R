########################################################################################
##                                                                                    ##
## This code submits jobs for all locations using MarginalTransformation_eCDF_excmu.R ##
##                                                                                    ##
########################################################################################

#!/usr/bin/Rscript
args <- commandArgs(TRUE)
for (arg in args) eval(parse(text = arg))
rm(arg, args)

##################
# Default values #
##################
MAXCPU = 16
omp = if(exists("omp") && omp>1) omp else 1
inicio = if(exists("inicio") && omp>0) inicio else 1
partition = if(exists("partition")) partition else "batch"
defaulttime='1-0' 
jobtime = if( exists("jobtime") ) jobtime else defaulttime

#####################################
# Folder to save .sub and time file #
#####################################
options(warn = -1)
rutasub = '~/EVAChallenge2019/Margins/OutputsByLocMu/Sub/'
dir.create(rutasub, recursive = TRUE)
options(warn = 0)

###########################
# Number of cpus and jobs #
###########################
# Replic = 16703
ncpus = if(exists("ncpus")) ncpus else MAXCPU
njobs = if(exists("njobs")) njobs else ceiling(Replic/ncpus) # This is the total number of jobs to be submitted


# Check
if (inicio > njobs ) stop ( paste("[ERROR] Inicio[",inicio,"] cannot be greater than njobs[",njobs,"]",sep="") )

for (i in inicio:njobs) {
  resto = Replic%%ncpus
  tmp = if ( i == njobs && resto > 0 ) resto else ncpus
  Rs <- (ncpus)*(i-1) + c(1:(tmp))
  ncpus = tmp
  
  nombre <- "excmu_"
  sbatch_line <- paste0("#!/bin/bash \n",
                  paste0("#SBATCH --job-name=",nombre, i,"\n"),
		  "#SBATCH -p ",partition,"\n",
		  # "#SBATCH -A stsda\n",
                  "#SBATCH -o Err_Outmu/%j.out\n",
                  "#SBATCH -e Err_Outmu/%j.err\n",
                  "#SBATCH --nodes=1 \n",
                  "#SBATCH -c ",ncpus,"\n",
		              "#SBATCH --mem=256000 \n",
                  "#SBATCH --time=",jobtime,"\n\n")
  r_line <- "module load R/3.4.2/intel-2017\n\n"
  time_line <- paste0("echo \"",i," $SLURM_JOBID `date`\" | tee -a ", '~/EVAChallenge2019/Margins/', "Runtime_eCDF_excmu.txt\n\n")
  code_line <- paste0("OMP_NUM_THREADS=",omp," Rscript /home/castroda/EVAChallenge2019/Margins/MarginalTransformation_eCDF_excmu.R \"Rs<-c(",
                   paste(Rs,collapse=","),
                   ")\" \"mc.cores<-",ncpus,
                   "\" \n")
  # .sub filename
  filename <- paste0(rutasub,"Application_",i,".sub")
  sub <- paste0(sbatch_line, r_line, time_line, code_line, time_line)
  
  if (file.exists(filename)) {
    unlink(filename)
  }
  cat(sub, file = filename)
  slurmjobid <- system( paste0("sbatch ",filename," | grep -oE '[^ ]+$' ") , intern=TRUE )
  message( paste0( "[",i,"]\tJOBID:", slurmjobid) )
  Sys.sleep(1)
  # system( paste("scontrol update JobId=",slurmjobid," Partition=defaultq Features=\"\"",sep="") )
}


