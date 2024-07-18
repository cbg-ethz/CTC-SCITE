source("functions.R")
library(viridis)
library(VGAM)
library(pscl)
library(MASS)
library(boot)
library("optparse")

############
# Config
############

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input-file"),
  type = "character",
  default = "~/Documents/projects/CTC_backup/input_folder", help = "Path to the folder containing all input files"
)
parser <- add_option(parser, c("-n", "--name-of-tree"),
  type = "character",
  default = "Br23", help = "Name of the tree for which to simulate CTC-clusters"
)
parser <- add_option(parser, c("-s", "--simulation-cluster-size"),
  type = "character",
  default = "2", help = "Number of cells in the simulated clusters"
)
parser <- add_option(parser, c("-o", "--output-folder"),
  type = "character",
  default = "~/Documents/projects/CTC_backup/simulations/simulation2", help = ""
)
args <- parse_args(parser, args = c("--input-folder", "--name-of-tree", "--simulation-cluster-size", "--output-folder"))


inputFolder <- dirname(args$input - file)
treeName <- args$name_of_tree
clusterSize <- args$simulation - cluster - size
outputFolder <- args$output - folder

# inputFolder <- "~/Documents/projects/CTC_backup/input_folder"
# treeName <- "Br23"


input <- load_data(inputFolder, treeName)


#
# ############
# # Exploratory data analysis
# ############
#
#
# # input <- load_data(inputFolder, treeName)
# # totalReadCounts <- input$totalReadCounts
# # sampleDescription <- input$sample_description
#
#
#
# # totalReadCountVector <- totalReadCounts %>% unlist()
#
#
# sum(totalReadCountVector == 0) / length(totalReadCountVector)
#
# # fit1 <- glm(totalReadCountVector ~ 1, family = poisson(link = 'log'))
# fit2 <- glm.nb(totalReadCountVector ~ 1)
# fit3 <- zeroinfl(totalReadCountVector ~ 1, dist = "negbin")
# # fit4 <- zeroinfl(totalReadCountVector ~1, dist = 'poisson')
#
# # summary(fit1)
# summary(fit2)
# summary(fit3)
# # summary(fit4)
# # exp(coef(fit))
# # coef(fit2)
#
#
# # exp(summary(fit3)$coefficients$zero[1])
#
#
# # parameterOfNegBinom <- exp(summary(fit3)$coefficients$count[,1])
#
#
# # exp(coef(fit3))
# # summary(fit4)
#
#
# # sim() <-
#
#
# simNew <- ifelse(rbinom(length(totalReadCountVector), size = 1, prob = exp(coef(fit3))[2]) > 0,
#   0, rnegbin(length(totalReadCountVector), exp(coef(fit3))[1], theta = exp(-0.76961))
# )
#
# sim <- data.frame(sim = vector(), run = vector())
#
# for (i in 1:100) {
#   simNew <- ifelse(rbinom(length(totalReadCountVector), size = 1, prob = exp(coef(fit3))[2]) > 0,
#     0, rnegbin(length(totalReadCountVector), exp(coef(fit3))[1], theta = exp(-0.76961))
#   )
#
#   # simNew <- rnegbin(length(totalReadCountVector), exp(coef(fit2)), theta = 0.9222)
#   sim <- rbind(sim, data.frame(sim = simNew, run = i))
# }
#
#
# sim <- rbind(sim, data.frame(sim = totalReadCountVector, run = 0))
#
# sim %>%
#   ggplot(aes(x = sim, group = run)) +
#   geom_histogram(data = sim[sim$run == 0, ], alpha = 0.4, color = "darkseagreen", fill = "darkseagreen") +
#   geom_freqpoly(data = sim[sim$run != 0, ], aes(x = sim), color = "red", position = "identity", alpha = 0.4)
#




#' Fits a zero inflated negative binomial distribution to the total read count data.
#'
#' @param input The loaded dataset
#' @param zeroInfl If this boolean value is FALSE, then a negative binomial is fit to the data
#'
#'
#' @return The parameters of the distribution. If zeroInfl is false, then the zero probability
#' is set to 0.
#' @export
#'
#' @examples
fitReadCountDistribution <- function(input, zeroInfl = TRUE) {
  totalReadCounts <- input$totalReadCounts
  sampleDescription <- input$sample_description
  totalReadCountVector <- totalReadCounts %>% unlist()

  if (zeroInfl == TRUE) {
    fit <- zeroinfl(totalReadCountVector ~ 1, dist = "negbin")
    return(list(zeroProb = inv.logit(summary(fit)$coefficients$zero[1]), theta = exp(summary(fit)$coefficients$count[2, 1]), expValue = exp(summary(fit)$coefficients$count[1, 1])))
  } else {
    fit <- glm.nb(totalReadCountVector ~ 1)
    return(list(zeroProb = 0, theta = summary(fit)$theta, expValue = exp(coef(fit))))
  }
}






# comparing different models it looks like a zero-inflated beta binomial model is
# appropriate to simulate read counts The coefficients are determined in fit3

#' From a number of wildtyoe and mutated genotypes, read counts are simulated as
#' follows:
#' 1. Each of the alleles drops out at constant rate "dropoutRate".
#' 2. The read count distribution is estimated from the data using a zero-inflated
#' negative-binomial model. This distribution is then used to estimate the total number
#' of read counts.
#' 3. The multiple-displacement amplification is modelled using a beta-binomial
#' model, given the total read counts sampled in step 2.
#' 4. Each allele may flip its genotype at rate "errorRate".
#'
#' @param nWildtypeAlleles
#' @param nMutatedAlleles
#' @param dropoutRate
#' @param errorRate
#' @param mu
#' @param theta
#'
#' @return A pair of read counts; the first one being the total number of reads
#' and the second one being the number of mutated reads.
#' @export
#'
#' @examples
simulateReads <- function(nWildtypeAlleles, nMutatedAlleles, dropoutRate, errorRate, readCountFit) {
  # draw from a binomial model to simulate dropouts
  nWildtypeAlleles <- rbinom(1, size = nWildtypeAlleles, prob = (1 - dropoutRate))
  nMutatedAlleles <- rbinom(1, size = nMutatedAlleles, prob = (1 - dropoutRate))




  # draw from a negative-binomial to simulate the total read count
  isZero <- rbinom(1, size = 1, p = readCountFit$zeroProb)
  if (isZero == TRUE) {
    nReads <- 0
  } else {
    nReads <- rnegbin(1, mu = readCountFit$expValue, theta = readCountFit$theta)
  }


  # draw from a beta-binomial to simulate overdispersion through multiple-
  # displacement amplification
  nWildtypeReads <- rbetabinom.ab(n = 1, size = nReads, shape1 = nWildtypeAlleles, shape2 = nMutatedAlleles)

  nMutatedReads <- nReads - nWildtypeReads

  # randomly flip the genotypes of reads with a certain error rate
  falsePositives <- rbinom(1, size = nReads - nMutatedReads, prob = errorRate)
  falseNegatives <- rbinom(1, size = nMutatedReads, prob = errorRate)

  nMutatedReads <- nMutatedReads + falsePositives - falseNegatives

  return(list(read_counts = c(nReads, nMutatedReads)))
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
getGenotypeMatrix <- function(nTreeSamplingEvents = 1000, input) {
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
  postSamplingTrees <- lapply(postSampling, FUN = function(entry) {
    return(entry$Tree)
  })


  logGenotypes <- getProbabilityOfBeingMutated(
    postSampling, nCells, nMutations, nClusters,
    alleleCount, ClusterID, mutatedReadCounts, totalReadCounts,
    rep(0, nCells)
  )

  genotypes_wide <- lapply(logGenotypes, FUN = exp)
  genotypes_wide <- data.frame(do.call(cbind, genotypes_wide))

  genotypes <-
    genotypes_wide %>%
    as_tibble() %>%
    rownames_to_column("Mutation") %>%
    pivot_longer(-Mutation, names_to = "Sample", values_to = "Posterior")


  ggplot(genotypes, aes(Mutation, Sample)) +
    geom_tile(aes(fill = Posterior)) +
    scale_fill_viridis()

  genotypes$WBC <- input$sample_description$WBC[(genotypes$Sample %>% substr(start = 2, stop = nchar(.)) %>% as.numeric())]

  genotypes %>%
    mutate(WBC = as.factor(WBC)) %>%
    filter(genotypes$WBC == 1) %>%
    ggplot(mapping = aes(x = Posterior, alpha = 0.6)) +
    geom_histogram(position = "identity", binwidth = 0.005)



  genotypes <- genotypes %>% mutate(Mutation = as.numeric(Mutation), Genotype = as.integer(Posterior > 0.5))

  genotypes %>%
    filter(Sample %in%
      paste0("X", which(input$sample_description$single_cell == TRUE & input$sample_description$WBC == FALSE))) %>%
    ggplot(aes(Mutation, Sample)) +
    geom_tile(aes(fill = Genotype)) +
    scale_fill_viridis()

  return(genotypes)
}



#' Creates the input dataset for CTC-SCITE run with simulated CTC-clusters.
#'
#' For the simulation, the follwing steps were performed:
#' 1. A zero-inflated negative binomial distribution is fit to the total read counts of a sample.
#' 2. For a new cell cluster, total read counts for each genomic position are sampled from the distribution fit in (1).
#' 3. For each mutation size:
#'    a) the total number of alleles is set to 2*(the number of cells in simulated cluster). The number of
#'       of mutated alleles is set to 0 (unmutated) or number of cells in the clustser (one mutated allele in each cell)
#'    b) Each of the alleles is removed at the dropout rate.
#'    c) The number of mutated reads are drawn from a beta-binomial distribution with alpha=mutated alleles, beta=non-mutated alleles, and n = total read count
#'    d) Each of the mutated reads is changed to non-mutated and vice versa at the error rate.
#'
#'
#'
#' @param samplingSize number of trees to determine the genotype of individual cells.
#' To be passed to getGenotypeMatrix.
#' @param clusterSizeVector A number that indicates the cluster complexity to by simulated
#' (i.e. the number of cells in the cluster)
#' @param input the loaded dataset
#' @param output_directory Directory to write the simulated input files for
#' the CTC-SCITE run to.
#' @param dropoutRate The dropout rate to assume for the simulation
#' @param errorRate The error rate to assume for the simulation
#' @param seed Set a seed for reproducibility
#' @param zeroInflated If this boolean vector is false, then the total read count will be
#' sampled from a negative binomial and not a zero-inflated negative binomial
#'
#' @return No return, but a "samples_nodeDescription.tsv and .txt file are written to
#' disk.
#' @export
#'
#' @examples
simulateCTCclusters <- function(samplingSize, clusterSizeVector, input, output_directory, output_label, dropoutRate = 0.3, errorRate = 0.001, seed = 123, zeroInflated = TRUE) {
  set.seed(seed)
  color_palette <- list("orchid", "orchid1", "orchid2", "orchid3", "orchid4", "darkorchid", "darkorchid1", "darkorchid2", "darkorchid3", "darkorchid4", "purple", "purple1", "purple2", "purple3", "purple4")


  fit <- fitReadCountDistribution(input, zeroInfl = zeroInflated)

  print("Calling genotypes")
  # Output data frame in long format. This is essentially a cell x mutation genotype matrix.
  # This represents to pool of genotypes from which I can now sample for the simulation.

  genotypes <- getGenotypeMatrix(nTreeSamplingEvents = samplingSize, input = input)


  genotypesOutputFormat <- data.frame(matrix(0, nrow = input$nMutations, ncol = 0))
  sampleDescriptionOutputFormat <- data.frame(matrix(0, nrow = 0, ncol = 5))
  colnames(sampleDescriptionOutputFormat) <- c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description")

  ## Sample as many genotypes as there should be simulated clusters.
  cellIDs <- paste0("X", 1:nrow(input$sample_description))



  if (sum(clusterSizeVector) > length(unique(genotypes$Sample))) {
    cells <- sample(size = length(unique(genotypes$Sample)), x = cellIDs, replace = FALSE)
    stop("You want to sample more genotypes than can be provided")
  } else {
    cells <- sample(size = sum(clusterSizeVector), x = cellIDs, replace = FALSE)
  }


  iterator <- 0

  # iterating over the size of the clusters to be simulated
  for (clusterSize in 1:length(clusterSizeVector)) {
    # iterating over the number of clusters of the same size to be simulated. Here not a for loop, to avoid backwards counting in R.
    clustersBySize <- 1
    while (clustersBySize <= clusterSizeVector[clusterSize]) {
      print(paste("Simulating CTC cluster ", iterator))
      print(paste("Number of cells: ", clusterSize))
      genotype <- genotypes %>%
        filter(Sample == cells[clustersBySize]) %>%
        arrange(Mutation)

      genotype <- pull(genotype, Genotype)

      nMutatedAlleles <- clusterSize * genotype
      nAllelesTotal <- clusterSize * rep(2, length(genotype))
      nWildtypeAlleles <- nAllelesTotal - nMutatedAlleles
      data <- data.frame(nWildtypeAlleles = nWildtypeAlleles, nMutatedAlleles = nMutatedAlleles)
      print("Starting simulation of read counts")
      reads <- apply(data, FUN = function(x) {
        return(simulateReads(x[1], x[2], dropoutRate, errorRate, fit)$read_counts)
      }, MARGIN = 1) %>% t()
      genotypesOutputFormat <- cbind(genotypesOutputFormat, reads)
      print("Done")

      newSample <- data.frame(
        sample_name = paste0(input$sampleName, "_sim", iterator),
        total_number_cells = clusterSize, tumor_cells = clusterSize,
        WBCs = 0,
        description =
          paste0("[color=", color_palette[[iterator + 1]], ',label="', input$sampleName, "_sim", iterator, '",fillcolor=', color_palette[[iterator + 1]], ',image="../CTC-cluster-icons/cluster_', clusterSize, '-0.png"]')
      )
      sampleDescriptionOutputFormat <- rbind(sampleDescriptionOutputFormat, newSample)

      iterator <- iterator + 1
      clustersBySize <- clustersBySize + 1
    }
  }
  print("Writing output files")


  dir.create(file.path(output_directory, paste(input$sampleName, output_label, sep = "_")), recursive = TRUE)
  description_data <- read_delim(file.path(input$directory, input$sampleName, paste0(input$sampleName, "_samples_nodeDescription.tsv")), delim = "\t", col_names = FALSE, quote = "none")
  colnames(description_data) <- c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description")
  description_data <- rbind(description_data, sampleDescriptionOutputFormat)
  write_delim(x = description_data, file = file.path(output_directory, paste(input$sampleName, output_label, sep = "_"), paste0(input$sampleName, "_", output_label, "_samples_nodeDescription.tsv")), delim = "\t", col_names = FALSE, quote = "none", escape = "none")

  read_data <- read_delim(file.path(input$directory, input$sampleName, paste0(input$sampleName, ".txt")), delim = "\t", col_names = FALSE, escape_backslash = TRUE)
  read_data <- cbind(read_data, genotypesOutputFormat)
  write_delim(x = read_data, file = file.path(output_directory, paste(input$sampleName, output_label, sep = "_"), paste0(input$sampleName, "_", output_label, ".txt")), delim = "\t", col_names = FALSE, quote = "none", escape = "none")
}





# for (tree in c("Br11", "Br16_AC_max2", "Br16_AC_max3", "Br16_AC_max4", "Br16_B_max1", "Br16_B_max2", "Br16_B_max3", "Br16_B_max4", "Br16_C_max1", "Br16_C_max2", "Br16_C_max3", "Br23", "Br26", "Br30", "Br37", "Br38", "Br39", "Br44", "Br45", "Br46", "Br53", "Br57", "Brx50", "Lu2", "Lu7", "Ov8", "Pr6", "Pr9")) {}
clusterSizeVector <- c(0, 4, 3, 2, 2, 2, 2, 2, 2)

keep <- rep(0, length(clusterSizeVector))
keep[clusterSize] <- 1
clusterSizeVector[keep == 0] <- 0



print(paste("Running simulation for", tree))

simulateCTCclusters(
  samplingSize = 100, clusterSizeVector = clusterSizeVector, input,
  output_directory = outputFolder, output_label = output_label,
  dropoutRate = 0.35, errorRate = 0.0015, seed = 124,
  zeroInflated = TRUE
)
