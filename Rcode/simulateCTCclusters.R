source('functions.R')
library(viridis)
library(VGAM)
library(pscl)
library(MASS)

############
#Config
############
inputFolder <- "../../input_folder"
treeName <- "Br7"


input <- load_data(inputFolder, treeName)



############
#Exploratory data analysis
############


#input <- load_data(inputFolder, treeName)
#totalReadCounts <- input$totalReadCounts
#sampleDescription <- input$sample_description



#totalReadCountVector <- totalReadCounts %>% unlist()



#fit1 <- glm(totalReadCountVector ~ 1, family = poisson(link = 'log'))
#fit2 <- glm.nb(totalReadCountVector ~ 1)
#fit3 <- zeroinfl(totalReadCountVector ~1, dist = 'negbin')
#fit4 <- zeroinfl(totalReadCountVector ~1, dist = 'poisson')

#summary(fit1)
#summary(fit2)
#exp(coef(fit))
#coef(fit2)


#exp(summary(fit3)$coefficients$zero[1])


#parameterOfNegBinom <- exp(summary(fit3)$coefficients$count[,1])


#exp(coef(fit3))
#summary(fit4)


#simNew <- ifelse(rbinom(length(totalReadCountVector), size = 1, prob = exp(coef(fit3))[2]) > 0,
#                 0, rnegbin(length(totalReadCountVector), exp(coef(fit3))[1], theta = exp(-0.76961)))

#sim <- data.frame(sim = vector(), run = vector())

#for(i in 1:100){
  
#  simNew <- ifelse(rbinom(length(totalReadCountVector), size = 1, prob = exp(coef(fit3))[2]) > 0,
#                   0, rnegbin(length(totalReadCountVector), exp(coef(fit3))[1], theta = exp(-0.76961)))
  
#  sim <- rbind(sim, data.frame(sim = simNew, run = i))
#}


#sim <- rbind(sim, data.frame(sim = totalReadCountVector, run = 0))

#sim %>%
 # ggplot(aes(x = sim, group = run)) + 
  #geom_histogram(data = sim[sim$run == 0,], alpha = 0.4, color = 'darkseagreen', fill = 'darkseagreen') +
#  geom_freqpoly(data = sim[sim$run != 0,], aes(x = sim), color = 'red', position = 'identity', alpha = 0.4)





#' Fits a zero inflated negative binomial distribution to the total read count data.
#'
#' @param input The loaded dataset
#'
#' @return The parameters of the distribution.
#' @export
#'
#' @examples
fitReadCountDistribution <- function(input){
  totalReadCounts <- input$totalReadCounts
  sampleDescription <- input$sample_description
  totalReadCountVector <- totalReadCounts %>% unlist()
  fit <- zeroinfl(totalReadCountVector ~1, dist = 'negbin')
  return(list(zeroProb = exp(summary(fit)$coefficients$zero[1]), theta = exp(summary(fit)$coefficients$count[2,1]), expValue = exp(summary(fit)$coefficients$count[1,1]) ))
}








#comparing different models it looks like a zero-inflated beta binomial model is
#appropriate to simulate read counts The coefficients are determined in fit3

#' From a number of wildtyoe and mutated genotypes, read counts are simulated as
#' follows:
#' 1. Each of the alleles drops out at constant rate "dropoutRate".
#' 2. The read count distribution is estimated from the data using a zero-inflated
#' negative-binomial model. This distribution is then used to estimate the total number
#' of read counts.
#' 3. The multiple-displacement amplification is modelled using a beta-binomial
#' model, given the total read counts sampled in step 2.
#' 4. Each allele mayflip its genotype at rate "errorRate".
#'
#' @param nWildtypeAlleles 
#' @param nMutatedAlleles 
#' @param dropoutRate 
#' @param errorRate 
#' @param mu 
#' @param theta 
#'
#' @return
#' @export
#'
#' @examples
simulateReads <- function(nWildtypeAlleles,nMutatedAlleles, dropoutRate, errorRate, readCountFit){
  #draw from a binomial model to simulate dropouts
  nWildtypeAlleles <- rbinom(1, size = nWildtypeAlleles,prob = (1-dropoutRate))
  nMutatedAlleles <- rbinom(1, size = nMutatedAlleles, prob = (1-dropoutRate))
  
  

  
  #draw from a negative-binomial to simulate the total read count
  isZero <- rbinom(1, size = 1, p = readCountFit$zeroProb)
  if(isZero == 0){
    nReads <- 0
  }
  else{
    nReads <- rnegbin(1, mu = readCountFit$expValue, theta = readCountFit$theta)
  }
  
  
  #draw from a beta-binomial to simulate overdispersion through multiple-
  #displacement amplification
  nMutatedReads <- rbetabinom.ab(n = 1, size = nReads, shape1 = nWildtypeAlleles, shape2 = nMutatedAlleles)
  
  #randomly flip the genotypes of reads with a certain error rate
  falsePositives <- rbinom(1, size = nReads-nMutatedReads, prob = errorRate)
  falseNegatives <- rbinom(1, size = nMutatedReads, prob = errorRate)
  
  nMutatedReads <- nMutatedReads + falsePositives - falseNegatives
  
  return(list(read_counts = c(nReads- nMutatedReads, nMutatedReads)))
}








#' Calls genotypes of single cells based on the CTC-SCITE algorithm 
#'
#' @param nTreeSamplingEvents number of sampled trees. Appricimated postserior gets better the higher this number is.
#' @param input The loaded data.
#'
#' @return returns a data frame in long format that gives the genotype and the
#' posterior genotype probability for each cell and sample. 
#' @export
#'
#' @examples
getGenotypeMatrix <- function(nTreeSamplingEvents = 100, input){
  
  postSampling <- input$postSampling
  nCells <- input$nCells
  nMutations <- input$nMutations
  nClusters <- input$nClusters
  alleleCount <- input$alleleCount
  ClusterID <- input$clusterID
  mutatedReadCounts <- input$mutatedReadCounts
  totalReadCounts <- input$totalReadCounts
  
  
  desired_values <- sample(1:length(postSampling), size = nTreeSamplingEvents, replace = FALSE) %>% sort()
  postSampling <- postSampling[desired_values]  
  postSamplingTrees <- lapply(postSampling, FUN = function(entry){return(entry$Tree)})

  
  logGenotypes <- getProbabilityOfBeingMutated(postSampling, nCells, nMutations, nClusters,
                                               alleleCount, ClusterID, mutatedReadCounts, totalReadCounts,
                                               rep(0,nCells))
  
  genotypes_wide <- lapply(logGenotypes, FUN = exp)
  genotypes_wide <- data.frame(do.call(cbind,genotypes_wide))
  
  genotypes <-
    genotypes_wide %>%
    as_tibble() %>%
    rownames_to_column("Mutation") %>%
    pivot_longer(-Mutation,names_to = "Sample", values_to = "Posterior")
  

  ggplot(genotypes, aes(Mutation, Sample)) +
    geom_tile(aes(fill = Posterior)) +
    scale_fill_viridis()
  
  genotypes$WBC <- input$sample_description$WBC[(genotypes$Sample %>% substr(start = 2, stop = nchar(.)) %>% as.numeric())]
  
  genotypes %>% mutate(WBC = as.factor(WBC)) %>%
    filter(genotypes$WBC ==1) %>%
    ggplot(mapping = aes(x = Posterior, alpha = 0.6)) +
    geom_histogram(position = "identity", binwidth = 0.005)
  
  
  
  genotypes <- genotypes %>% mutate(Mutation = as.numeric(Mutation), Genotype = as.integer(Posterior > 0.5))
  
  genotypes %>% filter(Sample %in% 
                         paste0("X",which(input$sample_description$single_cell == TRUE & input$sample_description$WBC == FALSE))) %>%
    ggplot(aes(Mutation, Sample)) +
    geom_tile(aes(fill = Genotype)) +
    scale_fill_viridis()
  
  return(genotypes)
}



#' Creates the input dataset for CTC-SCITE run with simulated CTC-clusters.
#'
#' @param samplingSize number of trees to determine the genotype of individual cells.
#' To be passed to getGenotypeMatrix.
#' @param clusterSizeVector A numeric vector that indicates how many clusters of which size
#' are to be simulated. number i at position j means that i many clusters with j
#' cells are to be simulated.
#' @param input the loaded dataset
#' @param output_directory Directory to write the simulated input files for
#' the CTC-SCITE run to.
#' @param dropoutRate The dropout rate to assume for the simulation
#' @param errorRate The error rate to assume for the simulation
#' @param seed Set a seed for reproducibility
#'
#' @return No return, but a "samples_nodeDescription.tsv and .txt file are written to
#' disk.
#' @export
#'
#' @examples
simulateCTCclusters <- function(samplingSize, clusterSizeVector, input, output_directory, dropoutRate = 0.3, errorRate = 0.001, seed = 123){
  set.seed(seed)
  color_palette <- list("orchid", "orchid1", "orchid2", "orchid3", "orchid4", "darkorchid1","darkorchid2", "darkorchid3", "darkorchid4", "purple", "purple1", "purple2", "purple3", "purple4")
  
  
  fit <- fitReadCountDistribution(input)
  
  print("Calling genotypes")
  genotypes <- getGenotypeMatrix(nTreeSamplingEvents = samplingSize, input = input)
  
  if(sum(clusterSizeVector) >= length(unique(genotypes$Sample))){
    stop('You want to sample more genotypes than can be provided')
  }
  
  
  genotypesOutputFormat <- data.frame(matrix(0, nrow =input$nMutations, ncol = 0))
  sampleDescriptionOutputFormat <- data.frame(matrix(0, nrow = 0, ncol = 5))
  colnames(sampleDescriptionOutputFormat) <- c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description") 
  
  cellIDs <- paste0("X",1:nrow(input$sample_description))
  cells <- sample(size = sum(clusterSizeVector), x = cellIDs, replace = FALSE)
  iterator <- 0
  
  
  for (clusterSize in 1:length(clusterSizeVector)){
    clustersBySize <- 1
    while(clustersBySize <= clusterSizeVector[clusterSize]){
      print(paste("Simulating CTC cluster", iterator))
      print(paste("Number of cells:", clusterSize))
      genotype <- genotypes %>%
        filter(Sample == cells[clustersBySize]) %>%
        arrange(Mutation)
      
      genotype <- pull(genotype, Genotype)
      
      nMutatedAlleles <- clusterSize * genotype
      nAllelesTotal <- clusterSize * rep(2, length(genotype))
      nWildtypeAlleles <- nAllelesTotal - nMutatedAlleles
      data <- data.frame(nWildtypeAlleles = nWildtypeAlleles, nMutatedAlleles = nMutatedAlleles)
      print("Starting simulation of read counts")
      reads <- apply(data, FUN = function(x){return(simulateReads(x[1], x[2],  dropoutRate, errorRate, fit)$read_counts)}, MARGIN = 1) %>% t()
      genotypesOutputFormat <- cbind(genotypesOutputFormat, reads)
      print("Done")
      
      newSample <- data.frame(sample_name = paste0(input$sampleName, '_sim', iterator),
                        total_number_cells = clusterSize, tumor_cells = clusterSize,
                        WBCs = 0,
                        description = 
                          paste0('[color=', color_palette[[iterator+1]], ',label=', input$sampleName, '_sim', iterator , ',fillcolor=', color_palette[[iterator+1]], ',image="../CTC-cluster-icons/cluster_', clusterSize,'-0.png"]') )
      sampleDescriptionOutputFormat <- rbind(sampleDescriptionOutputFormat, newSample)
      
      iterator <- iterator + 1
      clustersBySize <- clustersBySize + 1
    }
  }
  print("Writing output files")
  description_data <- read_delim(file.path(input$directory, input$sampleName, paste0(input$sampleName, '_samples_nodeDescription.tsv')), delim = '\t', col_names = FALSE, quote = "none")
  colnames(description_data) <- c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description") 
  description_data <- rbind(description_data, sampleDescriptionOutputFormat)
  write_delim(x = description_data, file = file.path(output_directory,input$sampleName, paste0(input$sampleName, '_samples_nodeDescription.tsv')), delim = '\t', col_names = FALSE, quote = "none", escape = "none")
  
  read_data <- read_delim(file.path(input$directory, input$sampleName, paste0(input$sampleName, '.txt')), delim = '\t', col_names = FALSE, escape_backslash	= TRUE)
  read_data <- cbind(read_data, genotypesOutputFormat)
  write_delim(x = read_data, file = file.path(output_directory,input$sampleName, paste0(input$sampleName, '.txt')), delim = '\t', col_names = FALSE, quote = "none", escape = "none")
}






simulateCTCclusters(samplingSize = 100, clusterSizeVector = c(0,2,1), input, output_directory = "../../input_folder/test", dropoutRate =  0.35, errorRate = 0.0015, seed = 123)
  
  
