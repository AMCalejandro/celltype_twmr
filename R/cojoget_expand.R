matrix_expansion = function(workdirpath, 
                            gwaseqtl_sign, 
                            eqtlname, 
                            ncorechr = 5, 
                            ncoregenes = 15) {
  
  allchrs = list.files(workdirpath)
  cojodata_prep =
    parallel::mclapply(allchrs, mc.cores = ncorechr, FUN = function(chrdir) {
      
      chrnumber = as.numeric(gsub(".*?([0-9]+)", "\\1", chrdir))
      allgenes = list.files(paste0(workdirpath, chrdir))
      d_filtered = readeqtlchunk(eqtlname, chr = chrnumber)
      path_amppd = paste0("/home/rstudio/DATA/AMPPD/QC_DATA/AMPPD-CONTROLS-SPLIT/AMPPD_RSID_CHR",chrnumber)
      newdirpath = paste0(workdirpath, chrdir, "/")
      copy_d_filtered = d_filtered
      
      
      getmatrices = 
        parallel::mclapply(1:length(allgenes), mc.cores = ncoregenes, FUN = function(index) {
          
          d_filtered = copy_d_filtered
          gene = allgenes[index]
          
          getref_focal = inferRefPanel(paste(newdirpath, gene, sep = "/"))
          geneinst = jmaprocessor(gene = gene,
                                  chrpath = newdirpath,
                                  refpanel = getref_focal,
                                  focalgene = TRUE)

          twmr_matrix = NULL
          if (!is.null(geneinst)) {
            
            #if (!file.exists(paste0(newdirpath, gene, "/", gene, ".matrix"))) {
              
              # Get focal gene SNPs and expand the matrix
              matrix_metadata = construct_fullmatrix(gene=gene,
                                                     geneinst = geneinst,
                                                     workdir = newdirpath, 
                                                     eqtldf = d_filtered, 
                                                     allgenes = allgenes, 
                                                     gwas = gwaseqtl_sign)
              
              twmr_matrix = matrix_metadata[[2]]
              if (!nrow(twmr_matrix)==0) {
                
                data.table::fwrite(data.frame(SNPS = twmr_matrix$GENE),
                                   paste0(newdirpath, gene, "/", gene, ".snps"))
                
                # Get LD and prune snps
                pathgene = paste0(newdirpath, gene, "/", gene)
                ldmatrix = computeld(genepath = pathgene, chr = chrnumber, refpanel = getref_focal)  
    
                overlapsnps = matrix_metadata[[1]]
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
                
              } else {
                message("matrix for gene: ", gene, "is null")
                twmr_matrix = NULL
              }
            
            #} else {
            #  message("Matrix already exists... Skipping gene ", gene)
            #  twmr_matrix = NULL
            #}
          
          } else {
            message("Gene focal matrix is empty. Skipping gene: ", gene)
            twmr_matrix = NULL
          }
          rm(d_filtered)
          gc()
          return(twmr_matrix)
        })
      
      # Run the cleanwd function to get rid of intermediary files
      cleanworkdir(wd = newdirpath, remove_all = FALSE)
      return(getmatrices)
    })
  return(cojodata_prep)
}





cojocompute = function(workdirpath, 
                       gwaseqtl_sign, 
                       eqtlname, 
                       ncorechr = 5, 
                       ncoregenes = 15) {
  
  allchrs = list.files(workdirpath)
  matrix_expand =
    mclapply(allchrs, mc.cores = ncorechr, FUN = function(chrdir) {
      
      chrnumber = as.numeric(gsub(".*?([0-9]+)", "\\1", chrdir))
      allgenes = list.files(paste0(workdirpath, chrdir))
      d_filtered = readeqtlchunk(eqtlname, chr = chrnumber)
      
      path_amppd = paste0("/home/rstudio/DATA/AMPPD/QC_DATA/AMPPD-CONTROLS-SPLIT/AMPPD_RSID_CHR",chrnumber)
      path_g1000 = paste0("/home/rstudio/DATA/REF_PANEL/g1000_eur/g1000_split/g1000_chr", chrnumber)
      
      newdirpath = paste0(workdirpath, chrdir, "/")
      
      # TEST using all overlapping snps to run COJO
      # After that, we construct the matrices
      # gwaseqtl_sign_chr = gwaseqtl_sign %>% filter(CHR == 1) 
      # path_snpsextract = paste0(newdirpath, "ALL.cojosnplist")
      # data.table::fwrite(data.frame(SNPS = gwaseqtl_sign_chr$SNP), path_snpsextract)
      # 
      # 
      # eqtldf_gene = d_filtered %>%
      #   dplyr::filter(CHR == 1) %>%
      #   dplyr::select(SNP, A1, A2, freq = MAF,
      #                 b = Beta, se = SE, p = P, N)
      # 
      # path_cojodf = paste0(newdirpath, "ALLEQTL.txt")
      # 
      # data.table::fwrite(eqtldf_gene, path_cojodf, 
      #                    sep = "\t", quote = F,
      #                    col.names = T, row.names = F)
      # 
      # path_cojodf = getcojogene(eqtldf = d_filtered,
      #                           gene = gene,
      #                           pathwrite = newdirpath)
      
      # 
      # 
      # newpath_out = paste0(tools::file_path_sans_ext(path_cojodf), "_1000G.txt")
      # 
      # cojoutput = runcojo(path_ref = path_g1000,
      #                     cojo = path_cojodf,
      #                     snps = path_snpsextract,
      #                     chr = chrnumber,
      #                     out = newpath_out)
      
      savesignsnps = 
        mclapply(1:length(allgenes), mc.cores = ncoregenes, FUN = function(index) {
          
          gene = allgenes[index]
          gwaseqtl_sign_snps = gwaseqtl_sign %>%
            dplyr::filter(ENSEMBL == .env[["gene"]]) %>%
            pull(SNP) %>% unique()
          
          if (length(gwaseqtl_sign_snps) == 1) {
            singleinst_gene = oneinst_process(gwaseqtl = gwaseqtl_sign,
                                              workdir = newdirpath,
                                              gene = gene,
                                              snp = gwaseqtl_sign_snps)
            # TODO
            #I think in the case that a single SNP has multiple exposures,
            # I would be missing an extra step to try to capture any extra SNP for third GENES
          } else {
            
            # path_snpsextract = paste0(newdirpath, gene, "/", gene, ".cojosnplist")
            # data.table::fwrite(data.frame(SNPS = gwaseqtl_sign_snps), path_snpsextract)
            # path_cojodf = getcojogene(eqtldf = d_filtered,
            #                           gene = gene,
            #                           pathwrite = newdirpath)
            # newpath_out = paste0(tools::file_path_sans_ext(path_cojodf), "_1000G.txt")
            # 
            # cojoutput = runcojo(path_ref = path_g1000,
            #                     cojo = path_cojodf,
            #                     snps = path_snpsextract,
            #                     chr = chrnumber,
            #                     out = newpath_out)
            # 
            # captureerr = paste(cojoutput, collapse = " ") %>%
            #   str_extract_all(".*error.*", simplify = TRUE)
            # 
            # 
            # if (!purrr::is_empty(captureerr)) {
            #   cojoutput = runcojo(path_ref = path_amppd,
            #                       cojo = path_cojodf,
            #                       snps = path_snpsextract,
            #                       chr = chrnumber,
            #                       out = path_cojodf)
            # }
            
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
              cojoutput = runcojo(path_ref = path_g1000,
                                  cojo = path_cojodf,
                                  snps = path_snpsextract,
                                  chr = chrnumber,
                                  out = newpath_out)
             }
          }
          return(gwaseqtl_sign_snps)
        })
      return(savesignsnps)
    })
  return(matrix_expand)
}