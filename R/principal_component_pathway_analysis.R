library(fields)
library(Matrix)
library(dplyr)
library(rlist)
library(stringr)
library(fgsea)
library(pbapply)
#' @description
#' @param mass_matrix is a matrix like object with each row corresponding to a pixel, and each column corresponding to a m/z value (Colnames = m/z values), support sparse/dense matrix,m each cell is the intensity of at the pixel; Note that the mass matrix should be half-processed, which means not filtered against any background signals, but undergoes normalization and baseline correction
#' @param width is the width of each imaging dimension, in pixels
#' @param height is the height of each imaging dimension, in pixels
#' @param ppm_error is the parts-per-million error tolerance of matching m/z value with potential metabolites
#' @param ion_mode is only needed when ppm_error is not specified, goes to function estimate_mz_resolution_error will be used to access the ppm_error
#' @param matrix_molecule is only needed when ppm_error is not specified, goes to function estimate_mz_resolution_error will be used to access the ppm_error
#' @param tof_resolution is the tof resolution of the instrument used for MALDI run, calculated by ion (ion mass,m/z)/(Full width at half height)
#' @param input_mz is used when mass_matrix doesn't have the column names to be the m/z value, it a list of m/z values with one-to-one correspondence with each column of the mass_matrix
#' @param num_retained_component is an integer value to indicated preferred number of PCs to retain
#' @param resampling_factor is a numerical value >0, indicate how you want to resample the size of roginal matrix
#' @param p_val_threshold is the p value threshold for pathways to be significant
#' @param biplot is a boolean to indcate whether biplot is included

#' @return pca is the pca information of original data
#' @return pathway_enrichment_pc is the pathway enrichment results for each PC




principal_component_pathway_analysis = function(mass_matrix,
                                                width,
                                                height,
                                                ppm_error = NULL,
                                                matrix_molecule = NULL,
                                                ion_mode = NULL,
                                                tof_resolution = 30000,
                                                input_mz = NULL,
                                                num_retained_component = NULL,
                                                variance_explained_threshold = 0.9,
                                                resampling_factor = 2,
                                                p_val_threshold = 0.05,
                                                biplot = F) {
  rescale <- function(x, newrange = range(x)) {
    xrange <- range(x)
    mfac <- (newrange[2] - newrange[1]) / (xrange[2] - xrange[1])
    newrange[1] + (x - xrange[1]) * mfac
  }

  ResizeMat <- function(mat, ndim = dim(mat)) {
    if (!require(fields))
      stop("`fields` required.")

    # input object
    odim <- dim(mat)
    obj <- list(x = 1:odim[1],
                y = 1:odim[2],
                z = mat)

    # output object
    ans <- matrix(NA, nrow = ndim[1], ncol = ndim[2])
    ndim <- dim(ans)

    # rescaling
    ncord <-
      as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
    loc <- ncord
    loc[, 1] = rescale(ncord[, 1], c(1, odim[1]))
    loc[, 2] = rescale(ncord[, 2], c(1, odim[2]))

    # interpolation
    ans[ncord] <- interp.surface(obj, loc)
    ans
  }
  # PCA analysis
  print("Scaling original matrix")
  mass_matrix = Matrix(as.matrix(mass_matrix), sparse = T)
  if (!is.null(resampling_factor)) {
    print("Running matrix resampling")
    pb = txtProgressBar(
      min = 0,
      max = ncol(mass_matrix),
      initial = 0,
      style = 3
    )
    if (!is.numeric(resampling_factor)) {
      stop("Please enter correct resampling_factor")
    }
    new.width = as.integer(width / resampling_factor)
    new.height = as.integer(height / resampling_factor)

    resampled_mat = matrix(nrow =  new.height * new.width)
    for (i in 1:ncol(mass_matrix)) {
      temp_mz_matrix = matrix(mass_matrix[, i],
                              ncol = width,
                              nrow = height,
                              byrow = F)
      resampled_temp = ResizeMat(temp_mz_matrix, c(new.width,
                                                   new.height))
      resampled_mat = cbind(resampled_mat, as.vector(resampled_temp))
      setTxtProgressBar(pb, i)
    }
    close(pb)
    resampled_mat = resampled_mat[, -1]
    colnames(resampled_mat) = colnames(mass_matrix)
    print("Resampling finished!")
    gc()
  }
  print("Running the principal component analysis (can take some time)")
  # Runing PCA

  resampled_mat_standardised =   as.matrix(t(t(resampled_mat) - Matrix::colSums(resampled_mat)/nrow(resampled_mat)))
  print("Computing the covariance")
  cov_mat <- cov(resampled_mat_standardised)
  print("Computing the eigenvalue/eigenvectors")
  eigen_result <- eigen(cov_mat)
  gc()
  # Extract eigenvectors and eigenvalues
  eigenvectors <- eigen_result$vectors
  eigenvalues <- eigen_result$values

    print("Computing PCA")
    pc = pbapply::pblapply(1:ncol(resampled_mat_standardised),function(i){
      temp = resampled_mat_standardised[,1]*eigenvectors[1,i]
      for(j in 2:ncol(resampled_mat_standardised)){
        temp =temp+resampled_mat_standardised[,j]*eigenvectors[j,i]
      }
      return(temp)
    })
    pc = do.call(cbind,pc)

  # make pca object
  pca = list(sdev = sqrt(eigenvalues),
             rotation = eigenvectors,
             center = Matrix::colSums(resampled_mat)/nrow(resampled_mat),
             scale = FALSE,
             x = pc)
  print("PCA finished")

  gc()
  rm(mass_matrix)
  eigenvalues = pca$sdev ^ 2
  # Step 5: Compute Principal Components
  # Choose number of principal components, k
  # if not input, use scree test to help find retained components
  if (is.null(num_retained_component)) {
    if (!is.null(variance_explained_threshold)) {
      tryCatch({
        cumulative_variance = cumsum(eigenvalues) / sum(eigenvalues)

        # Plot cumulative proportion of variance explained
        plot(
          cumulative_variance,
          type = 'b',
          main = "Cumulative Variance Explained",
          xlab = "Number of Principal Components",
          ylab = "Cumulative Proportion of Variance Explained"
        )

        # Add a horizontal line at the desired threshold
        threshold = variance_explained_threshold  # Example threshold
        abline(h = threshold,
               col = "red",
               lty = 2)

        # Find the number of principal components to retain based on the threshold
        retained =  which(cumulative_variance >= threshold)[1] - 1
      },
      error = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      },
      warning = function(cond) {
        stop(
          "Check if correct variance threshold for principle components are inputted, should be numeric value between 0 and 1"
        )
      })

    } else{
      # if threshold not inputted, use Kaiser's criterion
      print(
        "Both variance_explained_threshold and num_retained_component not inputted, use Kaiser's criterion for determination"
      )
      plot(
        eigenvalues,
        type = 'b',
        main = "Scree Plot",
        xlab = "Principal Component",
        ylab = "Eigenvalue"
      )

      # Add a horizontal line at 1 (Kaiser's criterion)
      abline(h = 1,
             col = "red",
             lty = 2)

      # Add a vertical line at the elbow point
      elbow_point <- which(diff(eigenvalues) < 0)[1]
      abline(v = elbow_point,
             col = "blue",
             lty = 2)
      retained = length(which(eigenvalues >= 1))
    }
  } else{
    retained = as.integer(num_retained_component)
    if (is.na(retained)) {
      stop("Please enter correct number of principle components to retain")
    }
  }

  if (!is.null(input_mz)) {
    if ((length(input_mz) != ncol(mass_matrix)) |
        (!is.numeric(input_mz))) {
      stop(
        "Please ensure input_mz has one-to-one correspondence with each column of the mass_matrix"
      )
    } else{
      input_mz = input_mz
    }
  } else{
    tryCatch({
      input_mz = data.frame(cbind(
        row_id = 1:length(colnames(mass_matrix)),
        mz = as.numeric(colnames(mass_matrix))
      ))
    },
    error = function(cond) {
      stop(
        "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
      )
    },
    warning = function(cond) {
      stop(
        "Check whether column names of the input matrix is correctly labelled as the m/z ratio"
      )
    })
  }

  #
  source(paste0(dirname(system.file(package = "SpaMTP")),"/R/fct_db_adduct_filter.R"))
  source(paste0(dirname(system.file(package = "SpaMTP")),"/R/fct_formula_filter.R"))
  source(paste0(dirname(system.file(package = "SpaMTP")),"/R/fct_proc_db.R"))
  #### Load the Cleaned and summarized DB ####
  Chebi_db     = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/inst/db_files/Chebi_1_names.rds"))
  HMDB_db      = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/inst/db_files/HMDB_1_names.rds"))

  # Set the db that you want to search against
  db = rbind(HMDB_db, Chebi_db)
  # set which adducts you want to search for
  adduct_file = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/inst/adduct_file.rds"))
  if (ion_mode == "positive") {
    test_add = sub(" ","",adduct_file$adduct_name[which(adduct_file$charge >= 0)])
  } else if (ion_mode == "negative") {
    test_add = sub(" ","",adduct_file$adduct_name[which(adduct_file$charge <= 0)])
  } else{
    stop("Please enter correct polarity")
  }
  # Using Chris' pipeline for annotation
  # 1) Filter DB by adduct.
  db_1 = db_adduct_filter(db, test_add, polarity = ifelse(ion_mode == "positive",
                                                          "pos", "neg"))

  # 2) only select natural elements
  db_2 = formula_filter(db_1)

  # 3) search db against mz df return results
  # Need to specify ppm error
  # If ppm_error not specified, use function to estimate

  if (is.null(ppm_error)) {
    if(is.null(matrix_molecule)){
      stop("Please enter correct matrix molecule, either in molecular mass (unit au) or entries")
    }else{
      ppm_error_df = estimate_mz_resolution_error(
        mass_matrix = resampled_mat,
        ion_mode = ion_mode,
        matrix_molecule = matrix_molecule,
        mass_resolution = tof_resolution,
        plot = F
      ) 
    }
    if (nrow(ppm_error_df) == 0) {
      stop(
        "No close metabolites find around given matrix mass, please specify a ppm_error, e.g 20"
      )
    } else{
      print("Successfully found peaks by given matrix molecule, ppm_error estimated")
      ppm_error_mean = mean(ppm_error_df$error_ppm)
    }
  }
  # Set error tolerance
  ppm_error = 1e6 / tof_resolution / sqrt(2 * log(2))
  db_3 = proc_db(input_mz, db_2, ppm_error) %>% mutate(entry = str_split(Isomers,
                                                                         pattern = "; "))
  print("Query necessary data and establish pathway database")
  input_id = lapply(db_3$entry, function(x) {
    x = unlist(x)
    index_hmdb = which(grepl(x, pattern = "HMDB"))
    x[index_hmdb] = paste0("hmdb:", x[index_hmdb])
    index_chebi = which(grepl(x, pattern = "CHEBI"))
    x[index_chebi] = tolower(x[index_chebi])
    return(x)
  })
  chem_props =readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/chem_props.rds"))
  db_3 = db_3 %>% mutate(inputid = input_id)
  rampid = c()
  chem_source_id = unique(chem_props$chem_source_id)
  pb = txtProgressBar(
    min = 0,
    max = nrow(db_3),
    initial = 0,
    style = 3
  )
  for (i in 1:nrow(db_3)) {
    rampid[i] = (chem_props$ramp_id[which(chem_source_id %in% db_3$inputid[i])])[1]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  db_3 = cbind(db_3, rampid)
  print("Query finished")
  ####################################################################################################
  get_analytes_db = function(input_id,analytehaspathway,
                             chem_props,pathway) {

    rampid = chem_props$ramp_id[which(chem_props$chem_source_id %in% unique(input_id))]
    #
    pathway_ids = analytehaspathway$pathwayRampId[which(rampid %in% analytehaspathway$rampId)]

    analytes_db = lapply(pathway_ids, function(x) {
      content = analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == x)]
      content = content[which(grepl(content, pattern = "RAMP_C"))]
      return(content)
    })
    analytes_db_name = unlist(lapply(pathway_ids, function(x) {
      name = pathway$pathwayName[which(pathway$pathwayRampId == x)]
      return(name)
    }))
    names(analytes_db) = analytes_db_name
    return(analytes_db)
  }
  # get rank pathway database
  print("Getting reference pathways")
  analytehaspathway = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/analytehaspathway.rds"))
  pathway = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/pathway.rds"))
  source = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/source.rds"))
  pathway_db = get_analytes_db(input_id,analytehaspathway,
                               chem_props,pathway)
  pathway_db = pathway_db[which(!duplicated(names(pathway_db)))]
  # get names for the ranks
  name_rank = lapply(input_mz$mz, function(x) {
    return(unique(na.omit(db_3$rampid[which(db_3$observed_mz == x)])))
  })

  #Set progress bar
  pb_new = txtProgressBar(
    min = 0,
    max = retained,
    initial = 0,
    style = 3
  )
  print("Runing set enrichment analysis")
  pca_sea_list = list()

  for (i in 1:retained) {
    # get the absolute value and sign of the loading
    loading = data.frame(cbind(
      abs_loading = abs(pca[["rotation"]][, i]),
      sign_loading = sign(pca[["rotation"]][, i])
    )) %>% arrange(desc(abs_loading))
    # run set enrichment analysis

    ranks = unlist(lapply(1:length(pca[["rotation"]][, i]), function(x) {
      pc_new = rep(pca[["rotation"]][x, i], times = length(name_rank[[x]]))
      names(pc_new) = name_rank[[x]]
      return(pc_new)
    }))
    ranks = ranks[which(!duplicated(names(ranks)))]
    suppressWarnings({
      gsea_result = fgsea(
        pathways =  pathway_db,
        stats = ranks,
        minSize = 5,
        maxSize = 500
      ) %>% filter(pval <= p_val_threshold) %>% mutate(principle_component = paste0("PC", i)) %>%
        mutate(leadingEdge_metabolites = lapply(leadingEdge, function(x) {
          temp = unlist(x)
          metabolites_name = unique(tolower(chem_props$common_name[which(chem_props$ramp_id %in% temp)]))
          return(metabolites_name)
        }))
    })
    setTxtProgressBar(pb_new, i)
    # Make sure sign of loading is positive to make it positively correlate with the PC
    pca_sea_list = list.append(pca_sea_list,
                               gsea_result)
  }
  close(pb_new)
  names(pca_sea_list) = paste0("PC", 1:retained)

  par(mfrow = c(2, 2))
  par(mar = c(2, 2, 1, 1))
  if (biplot == T) {
    biplot(pca, choices = c(1, 2), cex = c(0.05, 0.8))
    biplot(pca, choices = c(2, 3), cex = c(0.05, 0.8))
    biplot(pca, choices = c(1, 3), cex = c(0.05, 0.8))
  }
  return(list(pca = pca,
                     pathway_enrichment_pc = pca_sea_list,
                     new.width = as.integer(width/as.numeric(resampling_factor)),
                     new.height = as.integer(height/as.numeric(resampling_factor)))
}
