#!/usr/bin/env Rscript


# Load libraries and utils
library(tidyverse)
library(purrr)
library(data.table)
library(here)
library(AnnotationDbi) #BiocManager::install("AnnotationDbi")
library(org.Hs.eg.db) #BiocManager::install("org.Hs.eg.db")
library(parallel) #To send tasks to the processing cores in parallel
source(here::here("../R", "matrixcreator.R"))
source(here::here("../R", "cojogenefiles.R"))
source(here::here("../R", "gwasqc.R"))


# ARGS
cmd_args=commandArgs(trailingOnly = TRUE)
workdirpath = cmd_args[2]
gwaseqtl_sign = fread(cmd_args[3])
eqtlname = cmd_args[4] # Only eqtlgen and metabrain for now



# Using 5x 10 cores




# Getting ready to launch the nested loops to run cojo
allchrs = list.files(workdirpath)

cojodata_prep =
  mclapply(allchrs, mc.cores = 2, FUN = function(chrdir) {
    chrnumber = as.numeric(gsub(".*?([0-9]+)", "\\1", chrdir))
    allgenes = list.files(paste0(workdirpath, chrdir))
    
    # Get the eqtldata for the current chromosome
    d_filtered = readeqtlchunk(eqtlname, chr = chrnumber)

    path_amppd = paste0("/home/rstudio/DATA/AMPPD/QC_DATA/AMPPD-SPLIT/AMPPD_RSID_CHR",chrnumber)
    newdirpath = paste0(workdirpath, chrdir, "/")
    
    savesignsnps = 
      mclapply(1:length(allgenes), mc.cores = 2, FUN = function(index) {
        gene = allgenes[index]
        gwaseqtl_sign_snps = gwaseqtl_sign %>%
          dplyr::filter(ENSEMBL == .env[["gene"]]) %>%
          pull(SNP)
        
        if (length(gwaseqtl_sign_snps) == 1) {
          singleinst_gene = oneinst_process(gwaseqtl = gwaseqtl_sign,
                                            workdir = newdirpath,
                                            gene = gene,
                                            snp = gwaseqtl_sign_snps)
          # TODO
          #I think in the case that a single SNP has multiple exposures,
          # I would be missing an extra step to try to capture any extra SNP for third GENES
        } else {
          path_snpsextract = paste0(newdirpath, gene, "/", gene, ".cojosnplist")
          data.table::fwrite(data.frame(SNPS = gwaseqtl_sign_snps), path_snpsextract)
          path_cojodf = getcojogene(eqtldf = d_filtered,
                                    gene = gene,
                                    pathwrite = newdirpath)
          
          cojoutput = runcojo(path_ref = path_amppd,
                              cojo = path_cojodf,
                              snps = path_snpsextract,
                              chr = chrnumber,
                              out = path_cojodf)
          
          captureerr = paste(cojoutput, collapse = " ") %>%
            str_extract_all(".*error.*", simplify = TRUE)
          
          if (!purrr::is_empty(captureerr)) {
            newpath_out = paste0(tools::file_path_sans_ext(path_cojodf), "_1000G.txt")
            path_g1000 = paste0("~/celltyping_wd/ref_panel/g1000_eur/g1000_split/g1000_chr", chrnumber)
            cojoutput = runcojo(path_ref = path_g1000,
                                cojo = path_cojodf,
                                snps = path_snpsextract,
                                chr = chrnumber,
                                out = newpath_out)
          }
        }
        gwaseqtl_sign_snps
      })
    savesignsnps
  })