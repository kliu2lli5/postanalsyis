#GX_PUBLIC directory
dir.gx <- "/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC"

#gene table
genetable <- read.table(gzfile("/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/UCSC-genome-tables/hg19.GENCODE19.gz"),
                      header = TRUE, comment.char = "", sep = "\t", as.is = TRUE, check.names = F)
save(genetable, file = "/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/hg19.GENCODE19.RData")


#recombination map
rm(list=ls())
recomb<- list()
for(chr in c(1:22, "X", "X_par1", "X_par2")) {
  infile <- paste("/GWD/appbase/projects/statgen/RD-MDD-GX_PUBLIC/HapMap/recombination/2011-01_phaseII_B37", 
                  paste("genetic_map_GRCh37_chr", chr, ".txt.gz", sep =""), sep = "/")
  recomb[[paste("chr", chr, sep ="")]] <- read.table(gzfile(infile), header = TRUE, check.names = F, 
                                                                    comment.char = "", sep = "\t", as.is = TRUE)
}
save(recomb, file = "/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/recombHapMapII_GRCh37.RData")


