#' eqtlprocessor
#'
#' TGiven a chr chunk of a eqtl study dataset, it does some processing 
#' according to the target analysis
#' @param eqtlstudy a character with the eqtl study name
#' @param chrnumber The chrnumber of the input chunk
#' @param outputformat The output format the eqtl is being processed for
#' @return An harmonised eqtl df for the input chr chunk
#' 
eqtlprocessor = function(eqtlstudy, chrnumber, outputformat) {
  
  if (eqtlstudy == "eqtlgen") {
    allcolnames = c("p.value", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                    "OtherAllele", "Zscore", "Gene", "GeneSymbol", "GeneChr", 
                    "GenePos", "NrCohorts", "NrSamples", "FDR", "BonferroniP",
                    "maf" , "eQTL", "BETA", "SE")
    d = 
      fread(
        input = paste0("~/celltyping_wd/celltype_twmr/test_epqtl/tmp/heightgwas_eqtlgen_twmr/eqtl_filtered/eqtlfiltered_chr",
                       chrnumber,
                       ".txt"),
        col.names = allcolnames)
    
    if (outputformat == "cojo") {
      d_filtered = d %>%
        dplyr::select(SNP, CHR = SNPChr, BP = SNPPos, ENSEMBL = Gene,
                      HGNC = GeneSymbol, A1 = AssessedAllele, A2 =OtherAllele, 
                      MAF = maf, P = `p.value`, Beta = BETA, SE, N = NrSamples, eQTL)
    }
    
  } else if (eqtlstudy == "amppd") {
    #TODO
  }
  
  
  return(d_filtered)  
}
