#cmd: /GWD/appbase/projects/statgen/GXapp/R3.0.0/R-3.0.0/bin/R --vanilla --args [analysis directory] [modelFile] </GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/GCandContrast.r >test.log


args <- commandArgs(trailingOnly = TRUE)
if (length(args) >=2) {
  adir0 <- args[1]
  modelFile <- args[2]
} else {
  stop("Please input [model file with column analysis and contrasts]")
}


adir0<- "/GWD/appbase/projects/statgen3/PGx7656_sirukumab/Analysis/ARA3003"
setwd(adir0)
modelFile<- "/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/GDCgtx/package_update/lili/model.contrast.txt"

#required columns: analysis, contrasts
#significant threshold?
models<- read.table(modelFile, na.strings=c('NA',''),colClasses="character", sep = "\t", 
                      as.is = T, header = T, strip.white = TRUE, quote = "", fill = TRUE)
if(!"analysis" %in% names(models))
  stop("Missing column analysis in ", modelFile)
if(!"contrasts" %in% names(models) & !"groups" %in% names(models))
  stop("Please specify column groups and/or contrasts  in ", modelFile)


gtxloc = "/GWD/appbase/projects/statgen/GXapp/G-P_assoc_pipeline/R-packages/x86_64-unknown-linux-gnu-library/3.0"
.libPaths(gtxloc)
library(gtx)
library(data.table) 

for( idx in 1:nrow(models)){
  m<- models[idx, "analysis"]
  agroups <- NULL
  acontrast <- NULL
  if("groups" %in% names(models)) {
    agroups <- tokenise.whitespace(models[idx, "groups"])
  }
  if("contrasts" %in% names(models)){
    acontrasts <- tokenise.whitespace(models[idx, "contrasts"])
    acontrasts.bad <- sapply(acontrasts, function(contrast1) return(length(unlist(strsplit(contrast1, "/"))) != 2))
    if (any(acontrasts.bad)) {
      stop('Model "', m, '" has contrasts(s) ',
         paste(acontrasts[acontrasts.bad], collapse = ", "),' not in format group1/group2')
    }
  }
  
  #obtain meta data for each group:modelname, groupname, modelcall, phenotype, association results
  agroups <- unique(c(agroups, unlist(strsplit(acontrasts, "/"))))
  
  if(length(agroups)<1) {
    message("Analysis ", m,  " is ignored due to missing group and/or contrast info")
    next
  }

  resg <- lapply(agroups, function(agroup1) {
    adir <- file.path(adir0,"analyses", m, agroup1)
    file.out <- file.path(adir, "ALL.out.txt.gz")
    file.out.gc<- file.path(adir, "ALL.out.gc.txt.gz")
    if(!file.exists(file.out) & !file.exists(file.out.gc))
      stop("No results for model ", m, " in group ", agroup1)
    message(Sys.time(), ": Collating results for model ", m, " in group ", agroup1)
    if(file.exists(file.out.gc)){
       res1 <- read.table(gzfile(file.out.gc), quote = "", comment.char = "", header = TRUE, stringsAsFactors = FALSE)
       names(res1) <- gsub("X.", "", names(res1), fixed = T)
       res1 <- data.table(res1[!is.na(res1$pvalue), ])
       setkey(res1, SNP)
    }
    else {
      res1 <- read.table(gzfile(file.out), quote = "", comment.char = "", header = TRUE, stringsAsFactors = FALSE)
      ## Confirm expected columns present
      stopifnot(all(c("SNP", "pvalue", "beta", "SE") %in% names(res1)))
      names(res1) <- gsub("X.", "", names(res1), fixed = T)
      ## Convert to data table to improve efficiency retaining only needed columns and rows with non-missing pvalue
      res1 <- data.table(res1[!is.na(res1$pvalue),])
      ## Sort by SNP
      setkey(res1, SNP)
      ## Makes sense to apply GC and not store redundant (non-GCed) results here
      ## pvalue from 'best' method, LRT or F test
      lambda <- res1[ , gclambda(pvalue)]
      setattr(res1, "lambda", lambda)
      invisible(res1[ , pvalue.GC := pvalue])
      invisible(if (lambda > 1.) {
        res1[ , pvalue.GC := pchisq(qchisq(pvalue, df = 1, lower.tail = FALSE)/lambda, df = 1, lower.tail = FALSE)]
      })
  
      ## Apply GC to SEs using pvalues from Wald test
      lambdaWald <- res1[ , gclambda(pchisq((beta/SE)^2, df = 1, lower.tail = FALSE))]
      setattr(res1, "lambdaWald", lambdaWald)
      invisible(res1[ , SE.GC := SE])
      invisible(if (lambdaWald > 1.) {
        res1[ , SE.GC := SE*sqrt(lambdaWald)] # overwrite, save memory
      })
      ## Set SE.GC to NA when mis-calibrated,
      ## working definition being Wald chisquare statistic > 2. times the LRT chisquare statistic
      res1[(beta/SE.GC)^2 / qchisq(pvalue.GC, df = 1, lower.tail = FALSE) > 2., SE.GC := NA]   
   
      for(c in c("SE", "SE.GC", "beta", "pvalue","pvalue.GC", "analysed.Freq1", "analysed.Rsq" ))
        setnames(res1, c, paste(c, agroup1, sep = "."))
      message(Sys.time(), ": Writing gc result for model ", m, " in group ", agroup1, " to ", file.out.gc)
      write.table(res1, file=gzfile(file.out.gc), sep = "\t", row.names = F, quote = F)
    }
    if("P_threshold_group" %in% names(models) & !is.na(models[idx, "P_threshold_group"])) {
      file.out.gc.sig<- file.path(adir, "ALL.out.gc.sig.txt")
      message(Sys.time(), ": Writing significant gc result for model ", m, " in group ", agroup1, " to ", file.out.gc.sig)
      gcP <- paste("pvalue.GC", agroup1, sep = ".")
      write.table(res1[!is.na(res1[[gcP]])&res1[[gcP]] <= models[idx, "P_threshold_group"], ], 
                  file=file.out.gc.sig, sep = "\t", row.names = F, quote = F)
    }
    return(res1)
  })
  names(resg) <- agroups

  if(is.null(acontrasts)) next
  
  resc <- lapply(acontrasts, function(contrast1) {
    message(Sys.time(), ": Collating results for model ", m, " for contrast ", contrast1)
    groups <- unlist(strsplit(contrast1, "/"))
    group1 <- groups[1]
    group2 <- groups[2]
  
    ## More efficient to work with chi2 statistics, monotonic in P and using bespoke GC calculation
    ## ce = contrast expression
    ce <- parse(text = paste("list(chi2 = (beta.", group1, " - beta.", group2, ")^2/",
                           "(SE.GC.", group1, "^2 + SE.GC.", group2, "^2))", sep = ""))
    res1 <- resg[[group1]][resg[[group2]], eval(ce)]
    lambda <- res1[ , median(chi2, na.rm = TRUE)/qchisq(0.5, df = 1, lower.tail = FALSE)]
    setattr(res1, "lambda", lambda) # This is a Wald-like test but it's the primary one for the contrast
    invisible(if (lambda > 1.) {
      res1[ , pvalue.GC := pchisq(chi2/lambda, df = 1, lower.tail = FALSE)]
    } else {
      res1[ , pvalue.GC := pchisq(chi2, df = 1, lower.tail = FALSE)]
    })
    res1[ , chi2 := NULL]
    contrastP <- paste("pvalue.GC_", group1,"_vs_",  group2, sep ="")
    setnames(sig, "pvalue.GC", contrastP)
             
    file.contrast<- file.path(adir0,"analyses", m, paste("ALL.", group1, "_vs_", group2, ".out.gc.txt.gz", sep= ""))
    message(Sys.time(), ": Writing contrast result for model ", m, " for contrast ", contrast1, " to ", file.contrast)
    write.table(res1, file=gzfile(file.contrast),  sep = "\t", row.names = F, quote = F)
    
    if("P_threshold_contrast" %in% names(models) & !is.na(models[idx, "P_threshold_contrast"])) {
      file.out.gc.sig<- file.path(adir0,"analyses", m, paste("ALL.", group1, "_vs_", group2, ".out.gc.sig.txt", sep= ""))
      sig<- res1[!is.na(res1[["contrastP"]])&res1[["contrastP"]]<= models[idx, "P_threshold_contrast"], ]
      lapply(c(group1, group2), function(g){sig[resg[[g]][sig$SNP, ]]})
      message(Sys.time(), ": Writing significant contrast result for model ", m, " for contrast ", contrast1, " to ", file.out.gc.sig)
      write.table(sig, 
                  file=file.out.gc.sig, sep = "\t", row.names = F, quote = F)
    }
    return(res1)
  })
  names(resc) <- contrasts1

}

