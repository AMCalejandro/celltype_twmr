#' construct_fullmatrix
#'
#' This function takes .jma file as input, and removes variants in LD (LD_r > 0.1)
#' @param gene gene we want to load the .jma file from
#' @param workdir The path to the chromosome chunk working dir
#' @param eqtldf The full eqtl data frame filtered by chromosome chunk
#' @param allgenes List of all genes in chromosome 1
#' @param gwas 
#' @param pthresh The threshold to consider significant instruments
#' @return The pruned jma file for the input gene
#' 
#' 
# Get focal gene SNPs and expand the matrix
construct_fullmatrix = function(gene, geneinst, eqtldf, gwas, allgenes, workdir) {
  
  # getref_focal = inferRefPanel(paste(workdir,gene, sep = "/"))
  # geneinst = jmaprocessor(gene = gene,
  #                         chrpath = workdir,
  #                         refpanel = getref_focal,
  #                         focalgene = TRUE)
  
  newgenes = exp_expander(geneinst$SNP, eqtldf = eqtldf)
  newgenes = base::setdiff(newgenes, gene)
  newgenes = base::intersect(newgenes, allgenes)
  
  
  extrageneinst = map_df(newgenes, function(mygene) {
    getref = inferRefPanel(paste(workdir,mygene, sep = "/"))
    jmaprocessor(gene = mygene,
                 chrpath = workdir,
                 refpanel = getref,
                 focalgene = FALSE,
                 focalsnps = geneinst$SNP)
  })
  
  if (!nrow(extrageneinst)==0) {
    overlapsnps = base::intersect(geneinst$SNP, extrageneinst$SNP)
    allgeneinst = rbind(geneinst, extrageneinst)
    twmr_matrix = matrixCreate(eqtl = allgeneinst, gwas = gwas,
                               overlapsnps = overlapsnps)
  } else {
    twmr_matrix = matrixCreate(eqtl = geneinst, gwas = gwas)
    overlapsnps = NULL
    
  }
  
  return(list(overlapsnps, 
              twmr_matrix))
  
}


#' jmaprocessor
#'
#' This function takes .jma file as input, and removes variants in LD (LD_r > 0.1)
#' @param gene gene we want to load the .jma file from
#' @param chrpath The path to the chromosome foldr containing target genes
#' @param pthresh The threshold to consider significant instruments
#' @param refpanel The name of the refpanel to use
#' @param focalgene If we are doing to the COJO results processing for the focal gene
#' @param focalsnps A vector of the focal SNPs
#' @return The pruned jma file for the input gene
#' 
jmaprocessor = function(gene, chrpath,
                        pthresh =  0.05 / 16000,
                        refpanel = NULL, 
                        focalgene = TRUE,
                        focalsnps = c()) {
  
  if (is.null(refpanel)) {
    if (file.exists(paste0(chrpath, gene, "/", gene, ".tmpmatrix"))) {
      onegenematrix = fread(paste0(chrpath, gene, "/",gene, ".tmpmatrix"))
      res = gather(onegenematrix, key = GENE, value = b, 
                   -c(SNP, BETA_GWAS)) %>%
        dplyr::select(SNP, b, GENE) %>%
        dplyr::filter(GENE == .env[["gene"]])
      
    } else {
      #stop(".matrix expected but not found")
      warning("Skipping gene as no cojo output found")
      res = NULL
    }
  } else {
    if (refpanel == "AMPPD") {
      if (file.exists(paste0(chrpath, gene, "/", gene, ".txt.jma.cojo"))) {
        jmagene = fread(paste0(chrpath,
                               gene, "/",
                               gene,
                               ".txt.jma.cojo"))
      } else {
        #stop("1000G as ref panel, but jma not found")
        jmagene = NULL
      }
    } else if (refpanel == "1000G") {
      if (file.exists(paste0(chrpath, gene, "/", gene, "_1000G.txt.jma.cojo"))) {
        jmagene = fread(paste0(chrpath,
                               gene, "/",
                               gene,
                               "_1000G.txt.jma.cojo"))
      } else {
        #stop("AMMPD as ref panel, but jma not found")
        jmagene = NULL
      }
    }
    
    if (!is.null(jmagene)) {
      jmagene = subset(jmagene, p < pthresh)
      jmagene$GENE = gene
      
      res = jmagene %>%
        dplyr::select(SNP, b, GENE)
    } else {
      res = NULL
    }
  }
  if ( (!focalgene) && (!is.null(res)) ) { # If it is not the focal gene, make sure it has some overlapping snps with the focalgene
    if (is.null(focalsnps)) {
      stop("Processing xtra exposure ", gene, " but focalsnps missing")
    }
    res_tmp = subset(res, SNP %in% focalsnps)
    if (nrow(res_tmp) == 0) {
      res = NULL # Extra instrument does not share independent signifcantly associated SNPs
    }
  }
  return(res)
}


#' ldpruner
#'
#' This function prunes an input LD matrix
#' @param ldmatrix The SNPs LD squared matrix
#' @param snps P-value ordered snp vector to iterate on
#' @return a vector of genes
#' 
ldpruner = function(ldmatrix, snps) {
  n = length(snps)
  to_keep = character()
  to_remove = character()
  
  ldmatrix = as.data.frame(ldmatrix)
  
  for (snpind in 1:length(snps)) {
    
    if (sum(length(to_keep), length(to_remove)) == n) {
      break
    }
    
    snp = snps[snpind]
    #print(snp)
    
    if (snpind == 1) {
      to_keep = snp
    } else {
      if (snp %in% to_remove) {
        next
      } else {
        to_keep = append(to_keep, snp)
      }
    }
    
    snp_ld = subset(ldmatrix, SNP == snp)[!colnames(ldmatrix) %in% 'SNP']
    
    # Get snps to remove
    snpsremove = colnames(snp_ld)[which((snp_ld > 0.1) & (snp_ld < 1))]
    snpsremove = setdiff(snpsremove, c(to_keep, to_remove))
    to_remove = append(to_remove, snpsremove)
  }
  return(to_keep)
}


#' oneinst_process
#'
#' This function generates the ld and gene matrices when a gene only has 1 eqtl
#' @param gwaseqtl data frame of significant instrument and exposures
#' @param workdir The path to the gene being processed
#' @param gene The gene to process
#' @param snp the single eqtl for the gene being processed
#' @return a vector of genes
#' 
oneinst_process = function(gwaseqtl, workdir, gene, snp) {
  genematrix = gwaseqtl %>%
    dplyr::filter(SNP == snp) %>%
    dplyr::select(SNP, ENSEMBL, BETA_EQTL, BETA_GWAS) %>%
    spread(ENSEMBL, BETA_EQTL) %>%
    dplyr::relocate(BETA_GWAS, .after = last_col())
  
   
  data.table::fwrite(genematrix, paste0(workdir, gene, "/", gene, ".tmpmatrix"))
  ldmatrix = as.data.frame(matrix(nrow = nrow(genematrix),ncol = nrow(genematrix)))
  ldmatrix = base::replace(ldmatrix, is.na(ldmatrix), 1)
  data.table::fwrite(ldmatrix, paste0(workdir, gene, "/", gene, ".ld"),
                     col.names = F, row.names = F )
}



#' matrixCreate
#'
#' Generate the matrix combining eqtl and GWAS betas
#' @param eqtl All the eqtl data generated in long format from COJO
#' @param gwas GWAS df
#' @param overlapsnps Get the instruments shared between the focal gene and the extra eposures
#' @return a matrix in wide format ready to use it in TWMR
#' 
matrixCreate = function(eqtl, gwas, overlapsnps = c()) {
  # Process matrix
  twmr_matrix_tmp = eqtl %>%
    dplyr::select(GENE, SNP, b) %>%
    tidyr::spread(GENE, b) %>%
    dplyr::rename(GENE = SNP) 
  
  # Filter qced gwas by the genes present in twmr
  gwas_proc = gwas %>% dplyr::distinct(SNP, .keep_all = TRUE)
  
  
  # Merge the data
  twmr_matrix = twmr_matrix_tmp %>%
    dplyr::inner_join(
      gwas_proc %>% dplyr::select(SNP, BETA_GWAS), 
      by = c("GENE"="SNP")) %>%
    base::replace(., is.na(.), 0)
  
  # This is some code that removes SNPs that are significantly associated with the exposure, but the 
  # SNP estimate in relation to the outcome is 0.
  # I think we shuold still keep this SNPs
  
  # Checking there are no unexpected SNPs with no effect on the phenotype
  # if (!length(which(twmr_matrix$BETA_GWAS == 0)) == 0) {
  #   dropindex = which(twmr_matrix$BETA_GWAS == 0)
  #   snps = twmr_matrix[dropindex, GENE]
  #   
  #   if (snps %in% overlapsnps) {
  #     stop(paste0("Some focal gene extra exposures overlapping genes \
  #                 do not have an effect on the phenotype. \n
  #                 Check: ", snps))
  #   #dropindex = which(twmr_matrix$BETA_GWAS == 0)
  #   #twmr_matrix = twmr_matrix[-dropindex, ]
  #   } else {
  #     twmr_matrix = twmr_matrix[-dropindex, ]
  #   }
  # }
  return(twmr_matrix)
}

  


#' inferRefPanel
#'
#' This function determines which ref panel was used to run cojo
#' @param pathfolder This is the path to the gene folder
#' @return a boolean indicating whether 1000G was used as reference panel or not
#' 
inferRefPanel = function(pathfolder) {
  whichref = NULL
  ref <- list.files(pathfolder, pattern=".*_1000G.*")
  
  if (length(ref) > 0 ) {
    whichref = "1000G" 
  } else {
    ref <- list.files(pathfolder, pattern=".txt.jma.cojo")
    if (length(ref) > 0 ) {
      whichref = "AMPPD"
    }  
  }
  return(whichref)
}


#' computeld
#'
#' This function generates an squared ld matrix, using the smae reference panel than for cojo
#' @param  genepath Path to gene folder
#' @param  chrnumber THe chromosome number
#' @param  refpanel Which ref panel was used to compute
#' @return a boolean indicating whether 1000G was used as reference panel or not
#' 
computeld = function(genepath, chr, refpanel) {
  
  if (is.null(refpanel)) {
    path_ref = paste0("/home/rstudio/DATA/AMPPD/QC_DATA/AMPPD-CONTROLS-SPLIT/AMPPD_RSID_CHR",chr)
  } else {
    if (refpanel == "AMPPD") {
      path_ref = paste0("/home/rstudio/DATA/AMPPD/QC_DATA/AMPPD-CONTROLS-SPLIT/AMPPD_RSID_CHR",chr)
    } else if (refpanel == "1000G") {
      path_ref = paste0("/home/rstudio/DATA/REF_PANEL/g1000_eur/g1000_split/g1000_chr", chr)
    }
  }
  
  ldplinkcmd =
    paste("~/software/plink",
          #"--bfile /home/rstudio/celltyping_wd/ref_panel/g1000_eur/g1000_eur",
        paste0("--bfile ", path_ref),
        "--r2 square",
        paste0("--extract ", paste0(genepath, ".snps")),
        paste0("--out ", genepath),
        sep = " ")
  ldout = system(ldplinkcmd, intern = TRUE)
  ldpath = ldout[which(str_detect(ldout, "\\.ld$"))]
  snps = fread(paste0(genepath, ".snps"))$SNPS
  ldmatrix = fread(ldpath, col.names = snps)
  ldmatrix$SNP = snps
  ldmatrix = as.data.frame(ldmatrix)
  
  return(ldmatrix)
}


#' matrices_proc
#'
#' This function generates an squared ld matrix, using the smae reference panel than for cojo
#' @param genematrix The full gene matrix with all significant instrument for the focal gene and extra exposures
#' @param  ldmatrix The LD squared matrix of the SNPs that make up the gene matrix
#' @param  snps Arranged SNPs to do the pruning
#' @return a list containing the pruned gene matrix and the corresponding LD suqared matrix
#' 
matrices_proc = function(genematrix, ldmatrix, snps) {
  
  matrix_pruned = subset(genematrix, GENE %in% snps)
  
  ldmatrix_pruned = subset(ldmatrix, SNP %in% snps) %>%
    as.data.frame()
  
  ldmatrix_pruned = ldmatrix_pruned[, (names(ldmatrix_pruned) %in% snps)]
  colnames(ldmatrix_pruned) <- NULL

  return(list(matrix_pruned, ldmatrix_pruned))
}









#' jmaprocessor
#'
#' This function takes .jma file as input, and removes variants in LD (LD_r > 0.1)
#' @param gene gene we want to load the .jma file from
#' @param chrpath The path to the chromosome foldr containing target genes
#' @param pthresh The threshold to consider significant instruments
#' @return The pruned jma file for the input gene
#' 
# jmaprocessor = function(gene, chrpath,
#                         pthresh =  0.05 / 16000,
#                         refpanel = NULL) {
#   if (is.null(refpanel)) {
#     if (file.exists(paste0(chrpath, gene, "/", gene, ".matrix"))) {
#       onegenematrix = fread(paste0(chrpath, gene, "/",gene, ".matrix"))
#       res = gather(onegenematrix, key = GENE, value = b, 
#                    -c(SNP, BETA_GWAS)) %>%
#         dplyr::select(SNP, b, GENE) %>%
#         dplyr::filter(GENE == .env[["gene"]])
#     } else {
#       stop(".matrix expected but not found")
#     }
#   } else {
#     if (refpanel == "AMPPD") {
#       if (file.exists(paste0(chrpath, gene, "/", gene, ".txt.jma.cojo"))) {
#         jmagene = fread(paste0(chrpath,
#                                gene, "/",
#                                gene,
#                                ".txt.jma.cojo"))
#         ldgene = fread(paste0(chrpath,
#                               gene, "/",
#                               gene,
#                               ".txt.ldr.cojo"), header = T, sep = "\t")
#         ldgene = ldgene[,1:(ncol(ldgene) - 1)]
#         # Try with getting the abs LD scores
#         ldgene[,2:ncol(ldgene)] = abs(ldgene[,2:ncol(ldgene)])
#       } else {
#         stop("1000G as ref panel, but jma not found")
#       }
#     } else if (refpanel == "1000G") {
#       if (file.exists(paste0(chrpath, gene, "/", gene, "_1000G.txt.jma.cojo"))) {
#         jmagene = fread(paste0(chrpath,
#                                gene, "/",
#                                gene,
#                                "_1000G.txt.jma.cojo"))
#         ldgene = fread(paste0(chrpath,
#                               gene, "/",
#                               gene,
#                               "_1000G.txt.jma.cojo"), header = T, sep = "\t")
#       } else {
#         stop("AMMPD as ref panel, but jma not found")
#       }
#     }
#     
#     # Use ldpruner to process ldgene and get SNP in LE (r2<0.1) only
#     snps = jmagene %>% arrange(p) %>% pull(SNP)
#     getsnps = ldpruner(ldmatrix = ldgene, snps = snps)
#     jmagene = subset(jmagene, SNP %in% getsnps)
#     jmagene = subset(jmagene, p < pthresh) # This should no longer be necessary after adding the extract step when running COJO
#     jmagene$GENE = gene
#   
#     res = jmagene %>%
#       dplyr::select(SNP, b, GENE)
#   }
#   return(res)
# }

# jmaprocessor = function(gene, chrpath,
#                         pthresh =  0.05 / 16000,
#                         refpanel = NULL, focalgene = FALSE) {
#   if (is.null(refpanel)) {
#     if (file.exists(paste0(chrpath, gene, "/", gene, ".matrix"))) {
#       onegenematrix = fread(paste0(chrpath, gene, "/",gene, ".matrix"))
#       res = gather(onegenematrix, key = GENE, value = b, 
#                    -c(SNP, BETA_GWAS)) %>%
#         dplyr::select(SNP, b, GENE) %>%
#         dplyr::filter(GENE == .env[["gene"]])
#       
#     } else {
#       stop(".matrix expected but not found")
#     }
#   } else {
#     if (refpanel == "AMPPD") {
#       if (file.exists(paste0(chrpath, gene, "/", gene, ".txt.jma.cojo"))) {
#         jmagene = fread(paste0(chrpath,
#                                gene, "/",
#                                gene,
#                                ".txt.jma.cojo"))
#       } else {
#         #stop("1000G as ref panel, but jma not found")
#         jmagene = NULL
#       }
#     } else if (refpanel == "1000G") {
#       if (file.exists(paste0(chrpath, gene, "/", gene, "_1000G.txt.jma.cojo"))) {
#         jmagene = fread(paste0(chrpath,
#                                gene, "/",
#                                gene,
#                                "_1000G.txt.jma.cojo"))
#       } else {
#         #stop("AMMPD as ref panel, but jma not found")
#         jmagene = NULL
#       }
#     }
#     
#     if (!is.null(jmagene)) {
#       jmagene = subset(jmagene, p < pthresh)
#       jmagene$GENE = gene
#       
#       res = jmagene %>%
#         dplyr::select(SNP, b, GENE)
#     } else {
#       res = NULL
#     }
#   }
#   return(res)
# }


#' jmaprocessor
#'
#' This function takes the jma matrix and processses it according to the input gene
#' @param gne The input gene
#' @param chrpath The path to the input gene working directory
#' @param pthresh The threshold to consider an instrument independently and 
#                 significantly asssociated with the epxosure
#' @param refpanel The reference panel that was used to get the jma matrix
#' @param focalgene A boolean. is the input gene the focal gene or a candidate
#'                    extra exposure
#' @param focalsnps The snp vector of the focal gene
#' @return a vector of genes
#'



# extraexp_jmaprocessor = function(gene, focalsnps, chrpath,
#                         pthresh =  0.05 / 16000,
#                         refpanel = NULL) {
#   if (is.null(refpanel)) {
#     if (file.exists(paste0(chrpath, gene, "/", gene, ".matrix"))) {
#       onegenematrix = fread(paste0(chrpath, gene, "/",gene, ".matrix"))
#       res = gather(onegenematrix, key = GENE, value = b, 
#                    -c(SNP, BETA_GWAS)) %>%
#         dplyr::select(SNP, b, GENE) %>%
#         dplyr::filter(GENE == .env[["gene"]])
#       
#     } else {
#       stop(".matrix expected but not found")
#     }
#   } else {
#     if (refpanel == "AMPPD") {
#       if (file.exists(paste0(chrpath, gene, "/", gene, ".txt.jma.cojo"))) {
#         jmagene = fread(paste0(chrpath,
#                                gene, "/",
#                                gene,
#                                ".txt.jma.cojo"))
#       } else {
#         #stop("AMPPD as ref panel, but jma not found")
#         jmagene = NULL
#       }
#     } else if (refpanel == "1000G") {
#       if (file.exists(paste0(chrpath, gene, "/", gene, "_1000G.txt.jma.cojo"))) {
#         jmagene = fread(paste0(chrpath,
#                                gene, "/",
#                                gene,
#                                "_1000G.txt.jma.cojo"))
#       } else {
#         #stop("1000G as ref panel, but jma not found")
#         jmagene = NULL
#       }
#     }
#     
#     if (!is.null(jmagene)) {
#       jmagene = subset(jmagene, p < pthresh)
#       jmagene$GENE = gene
#       res = jmagene %>%
#         dplyr::select(SNP, b, GENE)
#     } else {
#       res = NULL
#     }
#   }
#   res_tmp = subset(res, SNP %in% focalsnps) # Checking the extra exposure share instruments with the focal gene
#   if (nrow(res_tmp) == 0) { # Extra instrument does not share independent signifcantly associated SNPs
#     res = NULL
#   }
#   return(res)
# }



#' exp_expander
#'
#' It takes the independent snps associated with a gene expression and it finds other regulated genes
#' @param snps the independent snps for a certain gene
#' @param eqtldf the filtered eqtl data frame
#' @return a vector of genes
#' 
exp_expander = function(snps, eqtldf) {
  exp = subset(eqtldf, SNP %in% snps) %>%
    pull(ENSEMBL) %>% unique()
}






