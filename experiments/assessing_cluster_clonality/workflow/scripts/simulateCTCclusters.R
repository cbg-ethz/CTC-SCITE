library(viridis)
library(VGAM)
library(pscl)
library(MASS)
library(boot)
source("~/work/CTC-SCITE/experiments/assessing_cluster_clonality/workflow/resources/functions.R")
library("optparse")


############
# Config
############


color_palette <-
  list(
    "orchid", "orchid1", "orchid2", "orchid3", "orchid4", "darkorchid",
    "darkorchid1", "darkorchid2", "darkorchid3", "darkorchid4", "purple",
    "purple1", "purple2", "purple3", "purple4"
  )



parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input-folder"),
  type = "character",
  default = "~/Documents/projects/CTC_backup/input_folder", help = "Path to the folder containing all input files"
)
parser <- add_option(parser, c("-n", "--name-of-tree"),
  type = "character",
  default = "Br23", help = "Name of the tree for which to simulate CTC-clusters"
)
parser <- add_option(parser, c("-s", "--simulation-cluster-size"),
  type = "numeric",
  default = "2", help = "Number of cells in the simulated clusters"
)
parser <- add_option(parser, c("-o", "--output-folder"),
  type = "character",
  default = "~/Documents/projects/CTC_backup/simulations/simulation3", help = ""
)
parser <- add_option(parser, c("-m", "--monoclonal"),
  type = "logical",
  default = TRUE, help = ""
)


#' Fits a zero inflated negative binomial distribution
#' to the total read count data.
#'
#' @param input The loaded dataset
#' @param zeroInfl If this boolean value is FALSE,
#' then a negative binomial is fit to the data
#'
#'
#' @return The parameters of the distribution. If zeroInfl is false,
#' then the zero probability is set to 0.
#' @export
#'
#' @examples
fitReadCountDistribution <- function(input, zeroInfl = TRUE) {
  totalReadCounts <- input$totalReadCounts
  sampleDescription <- input$sample_description
  totalReadCountVector <- totalReadCounts %>% unlist()

  if (zeroInfl == TRUE) {
    fit <- zeroinfl(totalReadCountVector ~ 1, dist = "negbin")
    return(
      list(
        zeroProb = inv.logit(summary(fit)$coefficients$zero[1]),
        theta = exp(summary(fit)$coefficients$count[2, 1]),
        expValue = exp(summary(fit)$coefficients$count[1, 1])
      )
    )
  } else {
    fit <- glm.nb(totalReadCountVector ~ 1)
    return(
      list(zeroProb = 0, theta = summary(fit)$theta, expValue = exp(coef(fit)))
    )
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
simulateReads <-
  function(
      nWildtypeAlleles, nMutatedAlleles, dropoutRate, errorRate, readCountFit) {
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
    nWildtypeReads <-
      rbetabinom.ab(
        n = 1, size = nReads, shape1 = nWildtypeAlleles, shape2 = nMutatedAlleles
      )

    nMutatedReads <- nReads - nWildtypeReads

    # randomly flip the genotypes of reads with a certain error rate
    falsePositives <- rbinom(1, size = nReads - nMutatedReads, prob = errorRate)
    falseNegatives <- rbinom(1, size = nMutatedReads, prob = errorRate)

    nMutatedReads <- nMutatedReads + falsePositives - falseNegatives

    return(list(read_counts = c(nReads, nMutatedReads)))
  }








#' Calls genotypes of single cells based on the CTC-SCITE algorithm
#'
#' @param nTreeSamplingEvents number of sampled trees. Approximated posterior
#' gets better the higher this number is.
#' @param input The loaded data.
#'
#' @return returns a data frame in long format that gives the genotype and the
#' posterior genotype probability for each cell and sample.
#' @export
#'
#' @examples
call_genotypes <- function(n_tree_sampling_events = 1000, input) {
  postSampling <- input$post_sampling
  nCells <- input$n_cells
  nMutations <- input$n_mutations
  nClusters <- input$n_clusters
  alleleCount <- input$allele_count
  ClusterID <- input$cluster_id
  mutatedReadCounts <- input$mutated_read_counts
  totalReadCounts <- input$total_read_counts


  desired_values <-
    sample(
      1:length(postSampling),
      size = n_tree_sampling_events, replace = FALSE
    ) %>%
    sort()
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

  genotypes$WBC <-
    input$sample_description$WBC[
      (genotypes$Sample %>% substr(start = 2, stop = nchar(.)) %>% as.numeric())
    ]

  genotypes %>%
    mutate(WBC = as.factor(WBC)) %>%
    filter(genotypes$WBC == 1) %>%
    ggplot(mapping = aes(x = Posterior, alpha = 0.6)) +
    geom_histogram(position = "identity", binwidth = 0.005)



  genotypes <-
    genotypes %>%
    mutate(
      Mutation = as.numeric(Mutation), Genotype = as.integer(Posterior > 0.5)
    )

  genotypes %>%
    filter(
      Sample %in%
        paste0(
          "X",
          which(
            input$sample_description$single_cell == TRUE &
              input$sample_description$WBC == FALSE
          )
        )
    ) %>%
    ggplot(aes(Mutation, Sample)) +
    geom_tile(aes(fill = Genotype)) +
    scale_fill_viridis()

  return(genotypes)
}


#' Samples a a specified number of genotypes
#'
#' @param input the generic input data from the CTC-SCITE tree sampling
#' @param genotypes called genotypes for each cell in long format
#' @param sampling_size a vector that indicates how many
#'
#' @return
#' @export
#'
#' @examples
sample_genotypes <- function(input, genotypes, sampling_size) {
  cellIDs <- paste0("X", 1:nrow(input$sample_description))

  if (sampling_size > length(unique(genotypes$Sample))) {
    cells <-
      sample(
        size = length(unique(genotypes$Sample)), x = cellIDs, replace = FALSE
      )
    stop("You want to sample more genotypes than can be provided")
  } else {
    cells <- sample(size = sampling_size, x = cellIDs, replace = FALSE)
  }
  return(cells)
}



#' Appends the simulated data to the original data and writes to new files
#'
#' @param output_directory
#' @param input The generic tree sampling data
#' @param output_label a name for the simulated dataset files, e.g. the number
#' of clusters
#' @param simulated_sample_description The lines for the sample description file
#' adding the simulated CTC clusters
#' @param genotypes_output_format The read counts for the simulated data
#'
#' @return No return, but writes the files to disk
#' @export
#'
#' @examples
create_simulated_output <-
  function(output_directory, input, output_label, simulated_sample_description,
           genotypes_output_format) {
    print("Writing output files")

    dir.create(
      file.path(
        output_directory, paste(input$sampleName, output_label, sep = "_")
      ),
      recursive = TRUE
    )
    description_data <-
      read_delim(
        file.path(
          input$directory,
          input$sampleName,
          paste0(input$sampleName, "_samples_nodeDescription.tsv")
        ),
        delim = "\t",
        col_names = FALSE,
        quote = "none"
      )
    colnames(description_data) <-
      c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description")

    description_data <-
      rbind(description_data, simulated_sample_description)
    write_delim(
      x = description_data,
      file = file.path(
        output_directory,
        paste(input$sampleName, output_label, sep = "_"),
        paste0(
          input$sampleName, "_",
          output_label,
          "_samples_nodeDescription.tsv"
        )
      ),
      delim = "\t",
      col_names = FALSE,
      quote = "none",
      escape = "none"
    )

    read_data <-
      read_delim(
        file.path(
          input$directory,
          input$sampleName,
          paste0(input$sampleName, ".txt")
        ),
        delim = "\t",
        col_names = FALSE,
        escape_backslash = TRUE
      )

    read_data <- cbind(read_data, genotypes_output_format)

    write_delim(
      x = read_data,
      file = file.path(
        output_directory,
        paste(input$sampleName, output_label, sep = "_"),
        paste0(input$sampleName, "_", output_label, ".txt")
      ),
      delim = "\t",
      col_names = FALSE,
      quote = "none",
      escape = "none"
    )
  }


#' Create input files for the CTC SCITE algorithm with
#' one simulated oligoclonal cluster
#'
#' @param input The generic posterior sampling from CTC-SCITE
#' @param number_of_cells The size of the output CTC-cluster
#' @param output_directory The directory to write the output to
#'
#' @return
#' @export
#'
#' @examples
simulate_oligoclonals <- function(input, output_directory, number_of_cells, sampling_size = 100) {
  read_data <-
    read_delim(
      file.path(
        input$directory,
        input$sampleName,
        paste0(input$sampleName, ".txt")
      ),
      delim = "\t",
      col_names = FALSE,
      escape_backslash = TRUE
    )

  description_data <-
    read_delim(
      file.path(
        input$directory,
        input$sampleName,
        paste0(input$sampleName, "_samples_nodeDescription.tsv")
      ),
      delim = "\t",
      col_names = FALSE,
      quote = "none"
    )
  colnames(description_data) <-
    c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description")




  cluster_identity_of_cells <- c()
  idx <- 1
  for (cell_count in description_data$total_number_cells) {
    for (i in 1:cell_count) {
      cluster_identity_of_cells <- c(cluster_identity_of_cells, idx)
    }
    idx <- idx + 1
  }



  genotypes <- call_genotypes(nTreeSamplingEvents = sampling_size, input = input)
  genotypes_wide <-
    genotypes %>%
    dplyr::select(c(Sample, Genotype, Mutation)) %>%
    pivot_wider(names_from = "Sample", values_from = Genotype)

  rownames(genotypes_wide) <- genotypes_wide$Mutation
  genotypes_wide <- genotypes_wide[, 2:ncol(genotypes_wide)]


  single_cell_indices <- which(description_data$tumor_cells == 1 & description_data$total_number_cells == 1)

  hamming_distance <- function(x, y) {
    return(sum(x != y))
  }

  it <- 0
  while (TRUE) { # I sample until I get a cluster where at least two cells are distinct
    # Sample cells to merge
    cells_to_merge <-
      sample(single_cell_indices, number_of_cells, replace = FALSE)
    cell_identity <- which(cluster_identity_of_cells %in% cells_to_merge)

    ## Check if all cells have the same genotype:
    ## If not, the cluster is oligoclonal and we are fine.

    sum_of_distances <- 0
    for (j in 1:length(cell_identity)) {
      for (k in 1:length(cell_identity)) {
        sum_of_distances <-
          sum_of_distances +
          hamming_distance(
            genotypes_wide[, cell_identity[j]],
            genotypes_wide[, cell_identity[k]]
          )
      }
    }
    if (sum_of_distances > 0) {
      break
    }
    it <- 1 + it
    if (it > 100) {
      break
    }
  }


  aggregated_data_ref <- rep(0, dim(read_data)[1])
  aggregated_data_alt <- rep(0, dim(read_data)[1])

  for (cell in cells_to_merge) {
    aggregated_data_ref <- aggregated_data_ref + read_data[, 4 + 2 * cell - 1]
    aggregated_data_alt <- aggregated_data_alt + read_data[, 4 + 2 * cell]
  }

  columns_to_remove <- c(4 + 2 * cells_to_merge - 1, 4 + 2 * cells_to_merge)

  read_data <- cbind(read_data, aggregated_data_ref, aggregated_data_alt)
  read_data <- read_data[, -columns_to_remove]

  output_label <- paste(cells_to_merge, collapse = "_")

  newSample <- data.frame(
    sample_name = paste(input$sampleName, "sim", output_label, sep = "_"),
    total_number_cells = number_of_cells, tumor_cells = number_of_cells,
    WBCs = 0,
    description =
      paste0(
        "[color=", color_palette[[1]],
        ',label="', input$sampleName, "_sim",
        '",fillcolor=',
        color_palette[[1]],
        ',image="../CTC-cluster-icons/cluster_',
        number_of_cells,
        '-0.png"]'
      )
  )

  description_data_output_format <-
    rbind(description_data, newSample)
  description_data_output_format <-
    description_data_output_format[-cells_to_merge, ]



  dir.create(
    file.path(
      output_directory, paste(input$sampleName, output_label, sep = "_")
    ),
    recursive = TRUE
  )

  write_delim(
    x = read_data,
    file = file.path(
      output_directory,
      paste(input$sampleName, output_label, sep = "_"),
      paste0(input$sampleName, "_", output_label, ".txt")
    ),
    delim = "\t",
    col_names = FALSE,
    quote = "none",
    escape = "none"
  )

  write_delim(
    x = description_data_output_format,
    file = file.path(
      output_directory,
      paste(input$sampleName, output_label, sep = "_"),
      paste0(
        input$sampleName, "_",
        output_label,
        "_samples_nodeDescription.tsv"
      )
    ),
    delim = "\t",
    col_names = FALSE,
    quote = "none",
    escape = "none"
  )
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
#' To be passed to call_genotypes
#' @param cluster_size_vector A number that indicates the cluster complexity to be simulated
#' (i.e. the number of cells in the cluster)
#' @param input the loaded dataset
#' @param output_directory Directory to write the simulated input files for
#' the CTC-SCITE run to.
#' @param output_label The number of cells in the simulated cluster
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
simulateCTCclusters <- function(
    samplingSize,
    cluster_size_vector,
    input,
    output_directory,
    output_label,
    dropoutRate = 0.3,
    errorRate = 0.001,
    seed = 123,
    zeroInflated = TRUE) {
  set.seed(seed)
  #  color_palette <-
  #    list(
  #      "orchid", "orchid1", "orchid2", "orchid3", "orchid4", "darkorchid",
  #      "darkorchid1", "darkorchid2", "darkorchid3", "darkorchid4", "purple",
  #      "purple1", "purple2", "purple3", "purple4"
  #      )


  print("Calling genotypes")
  # Output data frame in long format.
  # This is essentially a cell x mutation genotype matrix.
  # This represents to pool of genotypes from which
  # I can now sample for the simulation.

  genotypes <-
    call_genotypes(nTreeSamplingEvents = samplingSize, input = input)

  fit <- fitReadCountDistribution(input, zeroInfl = zeroInflated)



  cells <-
    sample_genotypes(
      input = input,
      genotypes = genotypes,
      sampling_size = sum(cluster_size_vector)
    )

  genotypes_output_format <-
    data.frame(matrix(0, nrow = input$nMutations, ncol = 0))
  sample_description_output_format <- data.frame(matrix(0, nrow = 0, ncol = 5))
  colnames(sample_description_output_format) <-
    c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description")


  iterator <- 0
  # iterating over the size of the clusters to be simulated
  for (size_of_cluster in 1:length(cluster_size_vector)) {
    # iterating over the number of clusters of the same size to be simulated.
    # Here not a for loop, to avoid backwards counting in R.
    clustersBySize <- 1
    while (clustersBySize <= cluster_size_vector[size_of_cluster]) {
      print(paste("Simulating CTC cluster ", iterator))
      print(paste("Number of cells: ", size_of_cluster))
      genotype <- genotypes %>%
        filter(Sample == cells[clustersBySize]) %>%
        arrange(Mutation)

      genotype <- pull(genotype, Genotype)

      nMutatedAlleles <- size_of_cluster * genotype
      nAllelesTotal <- size_of_cluster * rep(2, length(genotype))
      nWildtypeAlleles <- nAllelesTotal - nMutatedAlleles
      data <- data.frame(
        nWildtypeAlleles = nWildtypeAlleles, nMutatedAlleles = nMutatedAlleles
      )
      print("Starting simulation of read counts")
      reads <- apply(data, FUN = function(x) {
        return(
          simulateReads(x[1], x[2], dropoutRate, errorRate, fit)$read_counts
        )
      }, MARGIN = 1) %>% t()
      genotypes_output_format <- cbind(genotypes_output_format, reads)
      print("Done")

      newSample <- data.frame(
        sample_name = paste0(input$sampleName, "_sim", iterator),
        total_number_cells = size_of_cluster, tumor_cells = size_of_cluster,
        WBCs = 0,
        description =
          paste0(
            "[color=", color_palette[[iterator + 1]],
            ',label="', input$sampleName, "_sim",
            iterator,
            '",fillcolor=',
            color_palette[[iterator + 1]],
            ',image="../CTC-cluster-icons/cluster_',
            size_of_cluster,
            '-0.png"]'
          )
      )
      sample_description_output_format <-
        rbind(sample_description_output_format, newSample)

      iterator <- iterator + 1
      clustersBySize <- clustersBySize + 1
    }

    if (cluster_size_vector[size_of_cluster] > 0) {
      create_simulated_output(
        output_directory,
        input,
        size_of_cluster,
        sample_description_output_format,
        genotypes_output_format
      )
    }
  }
}






main <- function(){
  
  args <- parse_args(parser)
  
  
  input_folder <- args$"input-folder"
  tree_name <- args$"name-of-tree"
  cluster_size <- args$"simulation-cluster-size"
  output_folder <- args$"output-folder"
  monoclonal <- args$monoclonal
  
   input_folder <- "~/Documents/projects/CTC_backup/input_folder"
   tree_name <- "Br23"
  
  input <- load_data(input_folder, tree_name)
  print("Input data successfully loaded.")
  
  
  
  
  
  cluster_size_vector <- c(0, 3, 3, 3, 3, 3, 3, 3, 3)
  
  
  print(paste("Running simulation for", tree_name))
  
  all_cluster_sizes <- input$sample_description %>%
    filter(WBC == 0 & color != "gray93") %>%
    group_by(color) %>%
    filter(n() > 1) %>%
    summarize(cluster_size = n()) %>%
    dplyr::select("cluster_size") %>%
    unique()
  
  
  if (monoclonal == TRUE) {
    keep <- rep(0, length(cluster_size_vector))
    keep[cluster_size] <- 1
    cluster_size_vector[keep == 0] <- 0
    print("Simulating monoclonal clusters.")
    simulateCTCclusters(
      samplingSize = 100, cluster_size_vector = cluster_size_vector, input,
      output_directory = output_folder, output_label = output_label,
      dropoutRate = 0.35, errorRate = 0.0015, seed = 124,
      zeroInflated = TRUE
    )
  } else {
    for (idx in 1:nrow(all_cluster_sizes)) {
      cluster_size <- all_cluster_sizes$cluster_size[idx]
      print("Simulating oligoclonal clusters.")
      for (idx2 in 1:cluster_size_vector[cluster_size]) {
        simulate_oligoclonals(input, output_folder, cluster_size, sampling_size = 100)
      }
      print("Success.")
    }
  }
}


