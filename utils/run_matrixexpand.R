#!/usr/bin/env Rscript

# Load libraries and utilr
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


allchrs = list.files(workdirpath)

cojodata_prep =
  parallel::mclapply(allchrs, mc.cores = 2, FUN = function(chrdir) {
    
    chrnumber = as.numeric(gsub(".*?([0-9]+)", "\\1", chrdir))
    allgenes = list.files(paste0(workdirpath, chrdir))
    
    d_filtered = readeqtlchunk(eqtlname, chr = chrnumber)
    
    path_amppd = paste0("/home/rstudio/DATA/AMPPD/QC_DATA/AMPPD-SPLIT/AMPPD_RSID_CHR",chrnumber) 
    newdirpath = paste0(workdirpath, chrdir, "/")
    copy_d_filtered = d_filtered
    

    getmatrices = 
      parallel::mclapply(1:length(allgenes), mc.cores = 10, FUN = function(index) {
        
        d_filtered = copy_d_filtered
        gene = allgenes[index]
        
        
        getref_focal = inferRefPanel(paste(newdirpath,gene, sep = "/"))
        geneinst = jmaprocessor(gene = gene,
                                chrpath = newdirpath,
                                refpanel = getref_focal,
                                focalgene = TRUE)
        
        twmr_matrix = NULL
        if (!is.null(geneinst)) {
          if (!file.exists(paste0(newdirpath, gene, "/", gene, ".matrix"))) {
            
            
            newgenes = exp_expander(geneinst$SNP, eqtldf = d_filtered)
            newgenes = base::setdiff(newgenes, gene)
            newgenes = base::intersect(newgenes, allgenes)
            
            # extrageneinst = map_df(newgenes, function(mygene) {
            #   getref = inferRefPanel(paste(newdirpath,mygene, sep = "/"))
            #   extraexp_jmaprocessor(gene = mygene,
            #                         focalsnps = geneinst$SNP,
            #                         chrpath = newdirpath,
            #                         refpanel = getref)
            # })
            
            extrageneinst = map_df(newgenes, function(mygene) {
              getref = inferRefPanel(paste(newdirpath,mygene, sep = "/"))
              jmaprocessor(gene = mygene,
                           chrpath = newdirpath,
                           refpanel = getref,
                           focalgene = FALSE,
                           focalsnps = geneinst$SNP)
            })
            
            if (!nrow(extrageneinst)==0) {
              overlapsnps = base::intersect(geneinst$SNP, extrageneinst$SNP)
              allgeneinst = rbind(geneinst, extrageneinst)
              twmr_matrix = matrixCreate(eqtl = allgeneinst, gwas = gwaseqtl_sign,
                                         overlapsnps = overlapsnps)
            } else {
              twmr_matrix = matrixCreate(eqtl = geneinst, gwas = gwaseqtl_sign)
              overlapsnps = NULL
            }
            
            
            
            if (!nrow(twmr_matrix)==0) {
              
              data.table::fwrite(data.frame(SNPS = twmr_matrix$GENE),
                                 paste0(newdirpath, gene, "/", gene, ".snps"))
              
              # Get LD and prune snps
              pathgene = paste0(newdirpath, gene, "/", gene)
              ldmatrix = computeld(genepath = pathgene, chr = chrnumber, refpanel = getref_focal)  
              
              if (!is.null(overlapsnps)) {
                snpsarranged = c(overlapsnps, setdiff(twmr_matrix$GENE, overlapsnps))
              } else {
                snpsarranged = twmr_matrix$GENE
              }
              
              snpspruned = ldpruner(ldmatrix = ldmatrix,
                                    snps = snpsarranged)
              
              get_matrices = matrices_proc(genematrix = twmr_matrix,
                                           ldmatrix = ldmatrix,
                                           snps = snpspruned)
              
              
              data.table::fwrite(get_matrices[[1]],
                                 paste0(newdirpath, gene, "/", gene, ".matrix"))
              write.table(get_matrices[[2]],
                          paste0(newdirpath, gene, "/", gene, ".ld"), 
                          col.names = FALSE, row.names = FALSE)
              twmr_matrix = get_matrices
              
            }
          } else {
            message("Matrix already exists... Skipping gene ", gene)
            twmr_matrix = NULL
          }
        } else {
          message("Gene focal matrix is empty. Skipping gene: ", gene)
          twmr_matrix = NULL
        }
        rm(d_filtered)
        gc()
        return(twmr_matrix)
      })
    return(getmatrices)
  })


# Save results
cat("Saving results")
saveRDS(cojodata_prep, "ALLRESULTS_MATRIXEXPAND.rds")

