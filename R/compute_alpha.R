#' alphaget
#'
#' This function takes the wd path, and it calculates alpha assoc of exposure 
#' outcome for each gene across chromosomes
#' @param wd working direcotry to look for run alpha in
#' @return a data.frame with all alpha results merged
#'
alphaget = function(wd, ngwas, neqtl) { 
  runalpha = system(paste0("bash ",
              " /home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/utils/run_twmr.sh ",
              wd, " ", ngwas, " ", neqtl,  " > /home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/res/",
              basename(wd), "_alphacalculate.log"))
  return(" Alpha successfully inferred ")
}


#' mergealphares
#'
#' It takes the working directory, look for alpha files, and store them under the results wd
#' #' @param wd working direcotry to look for run alpha in
#' @return a data.frame with all alpha results merged
mergealphares = function(wd) {

  gwaseqtlname = basename(wd)
  files = purrr::map(list.files(wd, full.names = T), function(chr) {
    genefiles = list.files(chr, full.names = T)
    res = list.files(genefiles, pattern = "*\\.alpha", full.names = T)
    res
  }) %>% purrr::compact()
  
  alldata = mclapply(files, mc.cores = 11, FUN = function(files_chr) {
     filecontent =  map_df(files_chr, function(x) {
       
       data = fread(x) 
       if (is_empty(data)) {
         return(NULL)
       }
       
       data = data %>%
         mutate(study = gwaseqtlname,
                chr = as.numeric(str_match(x, "chr(.*?)/")[,2]),
                alpha = as.numeric(alpha),
                SE = as.numeric(SE),
                P = as.numeric(P),
                Nsnps = as.integer(Nsnps),
                Ngene = as.integer(Ngene))
       return(data)
     })
     return(filecontent)
  }) %>% bind_rows() %>% as.data.frame()
  
  # Storing results
  storealpha = writeres(wd = wd, alpha = alldata)
  
  return(alldata)
}