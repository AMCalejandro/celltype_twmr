#!/usr/bin/env Rscript


# Load libraries and utils
library(tidyverse)
library(purrr)
library(data.table)
library(here)
library(parallel) #To send tasks to the processing cores in parallel
source(here::here("../R", "matrixcreator.R"))
source(here::here("../R", "cojogenefiles.R"))
source(here::here("../R", "gwasqc.R"))
source(here::here("../R", "cojoget_expand.R"))
source(here::here("../R", "workdirops.R"))
source(here::here("../R", "compute_alpha.R"))

# ARGS
cmd_args=commandArgs(trailingOnly = TRUE)
workdirpath = paste0(cmd_args[2], "/")
gwaseqtl_sign = fread(cmd_args[3])
eqtlname = cmd_args[4] # Only eqtlgen and metabrain for now
ngwas = cmd_args[5]
neqtl = cmd_args[6]




# Run Cojo inference parallel job
cat("Starting cojo job \n")
start=Sys.time()
cojores = cojocompute(workdirpath = workdirpath,
                         gwaseqtl_sign = gwaseqtl_sign,
                         eqtlname = eqtlname)
end=Sys.time()
cat("Time running cojo job: ", end - start, "\n")




#Run matrix expansion job
cat("Starting expansion job \n")
start=Sys.time()
expandcojores = matrix_expansion(workdirpath = workdirpath,
                                gwaseqtl_sign = gwaseqtl_sign,
                                eqtlname = eqtlname)
end=Sys.time()
cat("Time running matrix expansion job: ", end - start, "\n")



# Calculate alpha for all genes under the working directory
alphaget(wd = workdirpath,
         ngwas = ngwas,
         neqtl = neqtl)
mergeres = mergealphares(wd = workdirpath)
cat("Number of unique genes: ", length(unique(mergeres$gene)))
cat("\n")
cat("Total entries in alpha: ", dim(mergeres)[1])


# Clean up the whole working directory
#cleanworkdir(wd = workdirpath, remove_all = TRUE)
