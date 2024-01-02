#' createworkdir
#'
#' This function generates the working directory to generate tmp files and run TWMR
#' @param workdir workind directory
#' @param gwasname name of your gwas study.
#' @param qtlname name of your qtl study.
#' @param genename the name of a gene or a vector of genes
#' @return newdir The path to working directory
createworkdir = function(workdir = tempdir(), 
                         gwasname = NULL, 
                         qtlname = NULL,
                         genename = NULL) {
  
  if (is.null(genename)) {
    
    newdir = paste0(workdir, paste0(gwasname, "_", qtlname, "/"))
    
    if(!file.exists(newdir)) { 
      dir.create(newdir)
    } else {
      warning("Folder already exists")
    }
  } else {
    genename = unique(genename)
    lapply(genename, function(gene) {
     newdir = paste0(workdir, gene)
     dir.create(newdir)
    })
    newdir = workdir
  }
  return(newdir)
}



#' cleanworkdir
#'
#' This function scans the workding directory and it recursively removes all temporary folders
#' @param workdir workind directory
#' 
cleanworkdir = function(wd, remove_all = FALSE) {
  if (remove_all) {
    cleanupcmd =
      paste0("find ", 
             wd,
             " -type f ", 
             "\\( -iname \\* \\)",
             " -delete")
  } else {
    cleanupcmd = 
      paste0("find ", 
             wd,
             " -type f ", "\\(",
             " -iname \\*.cojosnplist",
             " -o -iname \\*.\\*snps", 
             " -o -iname \\*.cojo",
             " -o -iname \\*.log",
             " -o -iname \\*.txt \\)",
             " -delete")
  }
  system(cleanupcmd)
}


# removeworkdir = function(workdir) {
#   removecmd = 
#     paste0("rm -rf ", workdir, "/*")
#   system(removecmd)
# }


writeres = function(wd, alpha) {
  
  gwaseqtlname = basename(wd)
  print("Writing results to res/ directory")
  filename = paste0("/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/res/", gwaseqtlname, ".res")
  
  if (file.exists(filename)) {
    file.remove(filename)
  }
  
  write.table(alpha, 
              paste0("/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/res/", gwaseqtlname, ".res"), 
              col.names = T, row.names = F, sep = "\t", quote = F)
  
}

writeresv2 = function(wd, df, alpha = T) {
  
  
  gwaseqtlname = basename(wd)
  print("Writing results to res/ directory")
  if (alpha) {  
    filename = paste0("/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/res/", gwaseqtlname, ".res")
    if (file.exists(filename)) {
      file.remove(filename)
    }
  } else {
    filename = paste0("/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/res_enrichment/", gwaseqtlname, "_enrich.res")
    if (file.exists(filename)) {
      file.remove(filename)
    }
    
  }
  write.table(df, filename, col.names = T, row.names = F, sep = "\t", quote = F)
}

