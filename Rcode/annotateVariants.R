library(tidyverse)


annotate_variants <- function(sampleName, inputFolder, variantList){

  
  data <- read_tsv(file.path(inputFolder,sampleName, paste0(sampleName,'.txt')), col_names = FALSE)
  colnames(data)[1] <- 'CHROM'
  colnames(data)[2] <- 'POS'
  colnames(data)[3] <- 'REF'
  colnames(data)[4] <- 'ALT'
  
  
  
  # Read VCF file to extract column names
  file <- file.path(inputFolder, 'filtered', 'vcf_files_annotated', paste0(sampleName, '.ann.vcf'))
  lines <- readLines(file, warn = FALSE)
  vcf_names <- strsplit(lines[grep("^#CHROM", lines)], "\t")[[1]]
  
  # Read VCF file into a data frame
  vcf <- read.table(file.path(inputFolder, 'filtered', 'vcf_files_annotated', paste0(sampleName, '.ann.vcf')),
                    comment.char = '#', sep = "\t", header = FALSE, col.names = vcf_names)
  colnames(vcf)[1] <- '#CHROM'
  # Extract functional annotations
  include <- rep('NONE',nrow(vcf))
  for (i in seq_along(vcf$INFO)) {
    functionalAnnotation <- unlist(strsplit(strsplit(vcf$INFO[i], ';')[[1]][2], ','))
    if(any(sapply(strsplit(functionalAnnotation, '\\|'), "[[", 3) == 'MODERATE')){
      impact <- 'MODERATE'
      include[i] <- impact
    }
    if(any(sapply(strsplit(functionalAnnotation, '\\|'), "[[", 3) == 'HIGH')){
      impact <- 'HIGH'
      include[i] <- impact
    }
  }
  
  
  
  # Check and filter based on functional annotation
  includeFunctionalAnnotation <- logical(nrow(data))
  for (i in seq_len(nrow(data))) {
    subset_rows <- vcf[vcf$'#CHROM' == data$'CHROM'[i] & vcf$POS == data$POS[i], ]
    
    if (nrow(subset_rows) != 1) {
      print(subset_rows)
      stop('More than one hit in the annotation file. ERROR')
    } else {
      includeFunctionalAnnotation[i] <- include[which(vcf$'#CHROM' == data$'CHROM'[i] & vcf$POS == data$POS[i])]
    }
  }
  
  data$relevant <- includeFunctionalAnnotation
  sum(data$relevant == 'FALSE')
  
  data <- data %>% mutate(variantName = paste(CHROM, POS, sep = '_'))
  
  return(data[, c('variantName', 'REF', 'ALT', 'relevant')])
}
