#' readeqtlchunk
#'
#' It takes a GWAS QCed by mungeSumStats and returns the data on COJO format
#' @param eqtlname The eqtl study we are using on the nalysis
#' @param chrnumber The chromoome number we are going to get a chunk for the eqtl data
#' @return d_filtered A data frame in the desire format ready for further processing
#'
readeqtlchunk = function(eqtlname = "eqtlgen", chr = chrnumber) {
  if (eqtlname == "eqtlgen") {
    d_filtered =
      fread(
        input = paste0("/home/rstudio/DATA/eQTLdata/eQTLGen/eQTLGen_split_tmp/eqtlgen_chr",
                       chr,
                       ".txt")) %>%
      dplyr::select(SNP, CHR = SNPChr, BP = SNPPos, ENSEMBL = Gene,
                    HGNC = GeneSymbol, A1 = AssessedAllele, A2 =OtherAllele,
                    MAF = maf, P = `p.value`, Beta = BETA, SE, N = NrSamples, eQTL)
    } else if (eqtlname == "metabrain") {
  
        d_filtered =
          fread(
            input = paste0("/home/rstudio/DATA/eQTLdata/METABRAIN_QC_V2/CORTEX/QCed_2020-05-26-Cortex-EUR-",
                           chr,
                           "-biogenformat.txt.gz")) %>%
          
          dplyr::mutate(A2 = ifelse(AlleleAssessed == AltAllele, 
                                    RefAllele, 
                                    AltAllele)) %>%
          dplyr::select(SNP, CHR = SNPChr, BP = SNPChrPos, 
                        ENSEMBL = ProbeName, HGNC = HGNCName, 
                        A1 = AlleleAssessed, A2, MAF, P = PValue, Beta = `Meta-Beta`,
                        SE = `Meta-SE`, N = NrSamples, eQTL)
    } else {
      stop("Unknown eqtl study")
    }
  return(d_filtered)
}




#' df_harmonise
#'
#' It takes a GWAS QCed by mungeSumStats or a eqtl study and returns the data on COJO format
#' @param df data.table a gwas df coming from MungeSumStats
#' @param eqtl bool if df is a eqtl study or not
#' @param eqtlname chr Name od the the eqtl study to process |eQTLGen or Metabrain|
#' @param get_sign bool if TRUE, df is filtered to extract significant eQTLs
#' @return data.table a gwas df COJO formatted
#' 
df_harmonise = function(df, eqtl = FALSE, eqtlname = NULL, 
                        get_sign = FALSE) {
  # if (!eqtl) {
  #   
  #   if ("N_CAS" %in% colnames(df)) {
  #     df$N =  df$N_CAS +  df$N_CON
  #   } else if ("NSTUDY" %in% colnames(df)) {
  #     df = df %>% dplyr::rename(N = NSTUDY)
  #   } else {
  #     df[ , 'N'] <- NA
  #   }
  #   
  #   if (!"FRQ" %in% colnames(df)) {
  #     df[ , 'FRQ'] <- NA
  #   }
  if (!eqtl) {
    
    if ("N_CAS" %in% colnames(df)) {
      df$N =  df$N_CAS +  df$N_CON
    } 
    if ("NSTUDY" %in% colnames(df)) {
      df = df %>% dplyr::rename(N = NSTUDY)
    }
    if (!"N" %in% colnames(df)) {
      df[ , 'N'] <- NA
    }
    if (!"FRQ" %in% colnames(df)) {
      df[ , 'FRQ'] <- NA
    }
    
    df = df %>%
      dplyr::select(SNP, A1, A2, freq = FRQ,
                    b = BETA, se = SE,
                    p = P, N)
  } else {
    if (eqtlname == "eqtlgen") {
      df = df %>%
        dplyr::mutate(eQTL = "eQTLGen") %>%
        dplyr::select(SNP, CHR = SNPChr, 
                      BP = SNPPos,
                      ENSEMBL = Gene, HGNC = GeneSymbol,
                      A1 = AssessedAllele, 
                      A2 = OtherAllele, MAF = maf,
                      P = Pvalue, Beta = BETA, SE,
                      N = NrSamples, eQTL)
      
    } else if (eqtlname == "metabrain"){
      df = df %>%
        dplyr::mutate(A2 = ifelse(AlleleAssessed == AltAllele,
                                  RefAllele,
                                  AltAllele)) %>%
        dplyr::select(SNP, CHR = SNPChr, BP = SNPChrPos, 
                      ENSEMBL = ProbeName, HGNC = HGNCName, 
                      A1 = AltAllele, A2, MAF, P = PValue, 
                      Beta = `Meta-Beta`, SE = `Meta-SE`, 
                      N = NrSamples, eQTL)
    } else {
      stop("Unrecognised eqtl dataset", eqtlname)
    }
    
    if (get_sign) {
      df = df %>%
        dplyr::filter(P < (0.05 / 16000))
    }
  }
  return(df)
}


#' cojoformat_qc
#'
#' This function finds the overlap between significant eQTLs and GWAS SNPs
#' @param gwas a gwas data frame
#' @param eqtl The eqtl data frame
#' @param gwasname The name of the gwas being procesed
#' @param eqtlname The name of the eqtl data being processed 
#' @param writepath Path to write the data
#' @return a gwas data.table COJO formatted
#' 
inst_sign_overlap = function(gwas, eqtl, 
                             gwasname, eqtlname,
                             writepath = NULL) {
  
  
  gwas = gwas %>%
   dplyr::select(SNP, BETA_GWAS = b, A1_GWAS = A1)
  eqtl = eqtl %>%
   dplyr::select(SNP, CHR, ENSEMBL, HGNC,
                 BETA_EQTL = Beta, A1_eqtl = A1)

  # Inner join significant data
  gwaseqtl_sign = gwas %>%
    inner_join(eqtl) %>% 
    dplyr::relocate(CHR, .after = SNP)

  
  if (is.null(writepath)) {
    get_tmp = tempfile(pattern = "file", tmpdir = tempdir())
    fwrite(gwaseqtl_sign, get_tmp, quote = F, sep = "\t",
           col.names = T, row.names = F)
    
  } else {
    fwrite(gwaseqtl_sign, writepath, quote = F, sep = "\t",
           col.names = T, row.names = F)
  }
  return(gwaseqtl_sign)
}



# FUNCTION BORROWED FROM  https://github.com/RHReynolds/colochelpR
#' Convert rs ids to CHR:BP locations.
#'
#' @param df dataframe. Dataframe containing SNPs as rs ids.
#' @param SNP_column chr. Name of column (in quotation marks) containing rs ids
#'   in dataframe.
#' @param dbSNP BS genome reference snps (choose appropriate dbSNP build
#'   dependent on genome build).
#'
#' @return Dataframe with rs ids and CHR:BP locations.
#' @export
#'

convert_rs_to_loc <- function(df, SNP_column, dbSNP){
  
  rs <- BSgenome::snpsById(dbSNP, df[[SNP_column]], ifnotfound = "drop") %>%
    as.data.frame() %>%
    tidyr::unite(col = "loc", seqnames, pos, sep = ":", remove = T) %>%
    dplyr::rename(rs = RefSNP_id) %>%
    dplyr::select(rs, loc)
  
  filter_vector <- c("rs")
  names(filter_vector) <- SNP_column
  
  df <- df %>%
    dplyr::inner_join(rs, by = filter_vector)
  
  rs_split = stringr::str_split_fixed(string = df$loc, pattern = ":", n = 2) %>%
    as.data.frame() %>%
    mutate(across(everything(), as.numeric))
  colnames(rs_split) = c("CHR", "BP")
  
  df = cbind(df, rs_split) %>% 
    dplyr::select(-loc)
  
  return(df)
  
}


# FUNCTION BORROWED FROM  https://github.com/RHReynolds/colochelpR
#' Convert CHR:BP locations to rs ids.
#'
#' @description Function will convert genomic co-ordinates to rs ids.
#'
#' @section Warning:
#' \itemize{
#'   \item Some CHR:BP locations have more than one
#'   associated rs id, thus some filtering for duplicates may have to occur
#'   after conversion. We leave this to the user to decide how they wish to
#'   filter.
#'   \item Some CHR:BP locations may not have an associated rs id --
#'   these will be represented by NA in the SNP column after conversion.
#'   }
#'
#' @param df dataframe. Dataframe containing SNPs as genomic locations. Must
#'   contain 2 columns labelled \code{CHR} and \code{BP}, with chromosome and
#'   base pair positions, respectively. If the dataframe contains additional
#'   columns included these will be preserved.
#' @param dbSNP BS genome reference snps (choose appropriate dbSNP build
#'   dependent on genome build).
#'
#' @return Dataframe with rs ids and CHR:BP locations.
#' @export
#'

convert_loc_to_rs <- function(df, dbSNP){
  
  # If df CHR column has "chr" in name, remove
  if(stringr::str_detect(df$CHR[1], "chr")){
    
    df <-
      df %>%
      dplyr::mutate(CHR = stringr::str_replace(CHR, "chr", ""))
    
  }
  
  # If columns CHR are not correct format, this can cause problems with later join
  df <-
    df %>%
    dplyr::mutate(CHR = as.factor(CHR),
                  BP = as.integer(BP))
  
  # Convert df to GRanges object
  df_gr <-
    GenomicRanges::makeGRangesFromDataFrame(df,
                                            keep.extra.columns = FALSE,
                                            ignore.strand = TRUE,
                                            seqinfo = NULL,
                                            seqnames.field = "CHR",
                                            start.field = "BP",
                                            end.field = "BP",
                                            starts.in.df.are.0based = FALSE)
  
  # Genomic position object as dataframe with SNP locations converted to RS id.
  df_gr <-
    BSgenome::snpsByOverlaps(dbSNP, df_gr, minoverlap = 1L) %>%
    # Note that the default value for minoverlap is 0 which means that, by default, in addition to the SNPs that are
    # located within the genomic regions specified thru the ranges argument, snpsByOverlaps also returns SNPs that are
    # adjacent to these regions. Use minoverlap=1L to omit these SNPs.
    as.data.frame()
  
  combined <-
    df_gr %>%
    dplyr::rename(
      SNP = RefSNP_id,
      CHR = seqnames,
      BP = pos,
      # Rename strand column just in case this is a column in the inputted df
      gr_strand = strand) %>%
    dplyr::right_join(df, by = c("CHR", "BP")) %>%
    dplyr::select(-gr_strand, -alleles_as_ambig)
  
  return(combined)
  
}