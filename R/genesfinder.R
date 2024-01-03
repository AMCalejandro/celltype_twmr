#' Genes finder
#'
#' Given a pvalue threshold and a gwas, it maps snps to genes based on distance alone
#' @param gwas gwas data.
#' @param generef the ncbi gene reference panel.
#' @param pvalue threshold to filter gwas snps.
#' @param parallel Whether to run the function on a parallel backend.
#' @param nthreads Number of threads to use if running on a parallel backend.
#' @return A list of vectors.
#' 
genesfinder = function(gwas, generef, pvalue = 5e-8,
                       parallel = FALSE, nthreads = 5) {
  gwasfilt = gwas %>%
    dplyr::filter(P < .env[["pvalue"]]) %>%
    dplyr::select(SNP, CHR, BP)
  
  generef = generef %>%
    dplyr::mutate(START = ifelse(START - 1500000 < 0, 0, START - 1500000),
                  END = END + 1500000)
  
  if (parallel) {
    getgenes = mclapply(c(1:22), 
                        mc.cores = nthreads,
                        FUN = function(chr) {
      gwasfiltchunk = gwasfilt %>%
        filter(CHR == chr)
      
      res = c()
      if (!(nrow(gwasfiltchunk)==0)) {
        generefchunk = 
          subset(generef, CHR == chr)
        
        for (mybp in gwasfiltchunk$BP) {
          genesnpmapping = generefchunk %>%
            filter((START < mybp) & (END > mybp)) %>%
            pull(ENSEMBL)
          res <- c(res, genesnpmapping)
        }
        res <- unique(res)
      }
      res
    })
  } else {
    getgenes = lapply(c(1:22), function(chr) {
      gwasfiltchunk = gwasfilt %>%
        filter(CHR == chr)
      
      res = c()
      if (!(nrow(gwasfiltchunk)==0)) {
        generefchunk = 
          subset(generef, CHR == chr)
        
        for (mybp in gwasfiltchunk$BP) {
          genesnpmapping = generefchunk %>%
            filter((START < mybp) & (END > mybp)) %>%
            pull(ENSEMBL)
          res <- c(res, genesnpmapping)
        }
        res <- unique(res)
      }
      res
    })
  }
  getgenes = getgenes[sapply(getgenes, Negate(is.null))] %>% unlist()
  # TODO
  # Approach to filter genes whose START BP is lower than X and 
  # whose END BP is higher than X
  # THE ONE ON TOP WORKS BUT IT COULD BE IMPROVED
  genechr = generef %>%
    dplyr::filter(ENSEMBL %in% getgenes) %>%
    dplyr::select(ENSEMBL, CHR)
  
  return(genechr)
}

#' createWorkdir
#'
#' This function does some minor QC on an input gene reference panel
#' @param generef gene reference panel df containing at least ENTREZ, chromosome , and start and end position
#' @param ensembl_id lgl whether the ensembl_id is present or not
#' @return gene reference panel on an harmonised format that will be used on downstream operations
#' 
generefqc = function(generef, ensembl_id = F) {
  
  generef$ENTREZ <- as.character(generef$ENTREZ)
  
  if (!ensembl_id) {
    entrez_ensembl <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
    entrez_ensembl_entrez <- entrez_ensembl %>% count(gene_id) %>% filter(n==1)
    entrez_ensembl_ens <- entrez_ensembl %>% count(ensembl_id) %>% filter(n==1)
    
    entrez_ensembl <-
      filter(entrez_ensembl,
             (gene_id %in% entrez_ensembl_entrez$gene_id) & 
               (ensembl_id %in% entrez_ensembl_ens$ensembl_id)) %>%
      dplyr::rename(ENTREZ = gene_id, ENSEMBL = ensembl_id)
    
    generef = 
      inner_join(generef, entrez_ensembl) 
  }
  
  generef = generef %>%
    dplyr::select(ENSEMBL, ENTREZ, HGNC, CHR, START, END)
  
  return(generef)
  
}


#' Genes finder
#'
#' Given a pvalue threshold and a gwas, it maps snps to genes based on distance alone
#' @param gwas gwas data.
#' @param generef the ncbi gene reference panel.
#' @param pvalue threshold to filter gwas snps.
#' @return A list of vectors.
#' 
# genesfinder = function(gwas, generef, pvalue = 5e-8) {
#   
#   gwasfilt = gwas %>%
#     dplyr::filter(P < .env[["pvalue"]]) %>%
#     dplyr::select(SNP, CHR, BP)
#   
#   generef = generef %>%
#     dplyr::mutate(START = ifelse(START - 1000000 < 0, 0, START - 1000000),
#                   END = END + 1000000)
# 
#   
#   
#   getgenes = lapply(c(1:22), function(chr) {
#     gwasfiltchunk = gwasfilt %>%
#       filter(CHR == chr)
#     
#     res = c()
#     if (!(nrow(gwasfiltchunk)==0)) {
#       generefchunk = 
#         subset(generef, CHR == chr)
#       
#       for (mybp in gwasfiltchunk$BP) {
#         genesnpmapping = generefchunk %>%
#           filter((START < mybp) & (END > mybp)) %>%
#           pull(ENSEMBL)
#         res <- c(res, genesnpmapping)
#       }
#       res <- unique(res)
#     }
#     res
#     })
#   getgenes = getgenes[sapply(getgenes, Negate(is.null))] %>% unlist()
#   # TODO
#   # Approach to filter genes whose START BP is lower than X and 
#   # whose END BP is higher than X
#   # THE ONE ON TOP WORKS BUT IT COULD BE IMPROVED
#   
#   genechr = generef %>%
#     dplyr::filter(ENSEMBL %in% getgenes) %>%
#     dplyr::select(ENSEMBL, CHR)
#   
#   return(genechr)
# }






