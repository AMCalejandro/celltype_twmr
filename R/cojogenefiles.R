#' getcojogene
#'
#' This function generates the input data to run COJO for a given gene
#' 
getcojogene = function(eqtldf, gene, pathwrite) {
  eqtldf_gene = eqtldf %>%
  dplyr::filter(ENSEMBL == .env[["gene"]]) %>%
    dplyr::select(SNP, A1, A2, freq = MAF,
                  b = Beta, se = SE, p = P, N)

  path_cojodf = paste0(pathwrite,
                       gene, "/",
                       gene, ".txt")
  
  data.table::fwrite(eqtldf_gene, path_cojodf, 
                     sep = "\t", quote = F,
                     col.names = T, row.names = F)
  
  return(path_cojodf)
}


#' runcojo
#'
#' This function runs cojo-slct analysis using system R function
#' @param path_ref Path to reference panel
#' @param cojo Path to cojo file we are going to use for indpependent associations analysis
#' @param snps Path to the candidate snps which are significant instruments
#' @param chr The chromosome number
#' @return The log of the cojo run
#' 
runcojo = function(path_ref = NULL,
                   cojo = path_cojodf,
                   snps = path_snpsextract,
                   chr = chrnumber,
                   out = path_cojodf) {

    # gctacmd = paste("~/software/gcta-1.94.1",
    #                 paste0("--bfile ", path_ref),
    #                 paste0("--cojo-file ", cojo),
    #                 paste0("--extract ", snps),
    #                 paste0("--chr ", chr),
    #                 "--cojo-slct", 
    #                 "--cojo-p 1e-3",
    #                 paste0("--out ", out),
    #                 sep = " ")
    gctacmd = paste("~/software/gcta-1.94.1",
                    paste0("--bfile ", path_ref),
                    paste0("--cojo-file ", cojo),
                    paste0("--extract ", snps),
                    paste0("--chr ", chr),
                    "--cojo-slct",
                    "--cojo-p 1e-3",
                    paste0("--out ", out),
                    sep = " ")
    
    
    getlog = system(gctacmd, intern = T)
    
    return(getlog)
}

# gctacmd = paste("~/software/gcta-1.94.1",
#                  paste0("--bfile ", path_amppd),
#                  paste0("--cojo-file ", path_cojodf),
#                  paste0("--extract ", path_snpsextract),
#                  paste0("--chr ", chrnumber),
#                  "--cojo-slct", 
#                  "--cojo-p 1e-3",
#                  paste0("--out ", path_cojodf),
#                  sep = " ")
# getlog = system(gctacmd, intern=TRUE)
# 
# captureerr = paste(getlog, collapse = " ") %>% 
#   str_extract_all(".*error.*", simplify = TRUE)
# 
# if (!purrr::is_empty(captureerr)) { # If I detect an error, I run using 1000G as the ref panel
#   path_ref = "~/celltyping_wd/ref_panel/g1000_eur/g1000_eur"
#   gctacmd = paste("~/software/gcta-1.94.1",
#                    paste0("--bfile ", path_ref),
#                    paste0("--cojo-file ", path_cojodf),
#                    paste0("--extract ", path_snpsextract),
#                    paste0("--chr ", chrnumber),
#                    "--cojo-slct", 
#                    "--cojo-p 1e-3",
#                    paste0("--out ", path_cojodf),
#                    sep = " ")              
#  system(gctacmd)

