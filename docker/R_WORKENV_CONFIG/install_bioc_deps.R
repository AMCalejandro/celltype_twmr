#!/usr/bin/env r

# Installing bioconductor deps
BiocManager::install(c("rtracklayer", 
                       "BSgenome.Hsapiens.1000genomes.hs37d5", 
                       "SNPlocs.Hsapiens.dbSNP155.GRCh37",
                       "AnnotationDbi",
                       "org.Hs.eg.db",
                       "MungeSumstats"), ask = F)
