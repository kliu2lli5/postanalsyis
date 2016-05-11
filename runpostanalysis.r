#cmd: /GWD/appbase/projects/statgen/GXapp/R3.0.0/R-3.0.0/bin/R --vanilla --args /GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/config_postanalysis.txt </GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/runpostanalysis.r >test.log

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >=1) {
  configFile <- args[1]
} else {
  stop("Please input [config file with full path]")
}
# configFile<-"/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/config_postanalysis.txt"

source("/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/pipeline.postanalysis.r")
postanalysis.pipeline(configFile)
