#' @estimate_mz_resolution_error 
# This is an application estimate the actual mass accuracy of your maldi data
#' @param mass_matrix is a matrix like object with each row corresponding to a pixel, and each column corresponding to a m/z value (Colnames = m/z values), support sparse/dense matrix,m each cell is the intensity of at the pixel; Note that the mass matrix should be half-processed, which means not filtered against any background signals, but undergoes normalization and baseline correction
#' @param matrix_molecule is either a single numerical value, representing the molecular mass of the matrix, or the chemical ID (chebi,HMDB,CAS,pubchem,Kegg, chemspider). If mutiple matrix is used as mixed, please use this function for individual ones
#' @param ion_mode is a char variable selected from either "positive" or "negative"
#' @param sparse is a boolean variable indicates whether sparse matrix is used to save memory
#' @param use_colnames_as_mz is a Boolean variable indicates whether use the colnames of mass_matrix as list of m/z ratioes
#' @param mz_vector is only needed when use_colnames_as_mz is set as FALSE, indicates the m/z values of the column vectors in mass_matrix
#' @param mass_resolution is a parameter of the instrument used for MALDI run, calculated by ion (ion mass,m/z)/(Full width at half height), set to be 30000 as default
#' @param Quantile is a numeric value between 0 and 1, is the typical quantile of intensity of matrix molecule, as in generally MALDI, matrix molecule is extensively applied, this value is close to 1

#' @return A dataframe contains the potential matches with the matrix molecule and corresponding mass errors.

# mass_matrix_dense = read.csv("/stornext/Bioinf/data/lab_brain_cancer/projects/tme_spatial/venture_multi_omics/venture_pt1/Ven1A_information_matrix_unfiltered.csv")
# mass_matrix_dense = mass_matrix_dense[,-1]
# colnames(mass_matrix_dense) = sub("X","",colnames(mass_matrix_dense))

estimate_mz_resolution_error = function(mass_matrix,
                                        matrix_molecule,
                                        ion_mode,
                                        sparse = T,
                                        use_colnames_as_mz = T,
                                        mass_resolution = 30000,
                                        mz_vector = NULL,
                                        quantile = 0.9){
  #(1) Converting matrix to sparse matrix when necessary
  require(Matrix)
  if(sparse == T){
    print("Converting dense matrix to sparse matrix")
    mass_matrix = Matrix(as.matrix(mass_matrix), sparse = TRUE)
    print("Conversion finished")
  }
  #(2) Get lsit of m/z values
  if(use_colnames_as_mz == T){
      tryCatch({
        mz_list = as.numeric(colnames(mass_matrix))
      },
      error = function(cond){
        stop("Error occurred: ", conditionMessage(cond),
             message("Please check whether the input matrix has correctnumeric  m/z ratio as colnames,
              or try to set 'use_colnames_as_mz = T', and put in a vector of m/z ratios 'mz_vector' in the same order as the columns of the input matrix"))
      },
      warning = function(cond){
        stop("Error occurred: ", conditionMessage(cond),
             message("Please check whether the input matrix has correctnumeric  m/z ratio as colnames,
              or try to set 'use_colnames_as_mz = T', and put in a vector of m/z ratios 'mz_vector' in the same order as the columns of the input matrix"))
      }
      )
  }else if(!is.null(mz_vector)){
    mz_list = as.numeric(mz_vector)
  }else{
    stop("Cannot find list of screened m/z ratioes.
         Check if m/z ratio vector or colnames of input matrix is correctly input!")
  }
  #(3) query for all possible m/z for matrix molecule
  source("./R/fct_db_adduct_filter.R")
  source("./R/fct_formula_filter.R")
  source("./R/fct_proc_db.R")
  #### Load the Cleaned and summarized DB ####
  # require krish's DB
  print("Loading query databases...")
  adduct_file = readRDS(".QIMR_project/inst/adduct_file.rds")
  load(file = "./data/chem_props.rda")
  print("Database loading finished.")
  
  #(4) Determine the polarity for adducts
  
  if(is.numeric(matrix_molecule)){
    matrix_molecule = matrix_molecule
  }else{
    # query for the give id
    matrix_molecule = as.numeric(chem_props$monoisotop_mass[which(grepl(chem_props$chem_source_id,
                       pattern = matrix_molecule,
                       ignore.case = T))][1])
    if(length(matrix_molecule) == 0 ){
      stop("Check if correct matrix molecule id is input")
    }
  }
    if(ion_mode == "positive"){
      positive_adducts = adduct_file[which(adduct_file$pol == "pos"),]
      formulas = gsub("([0-9])([A-Za-z])", "\\1*\\2",  positive_adducts$ion.mass)
      M = matrix_molecule
      potential_mz = c(sapply(formulas, function(x) eval(parse(text = x))))
      names(potential_mz) = positive_adducts$adduct_name
    }else if(ion_mode == "negative"){
      negative_adducts = adduct_file[which(adduct_file$pol == "neg"),]
      formulas = gsub("([0-9])([A-Za-z])", "\\1*\\2",  negative_adducts$ion.mass)
      M = matrix_molecule
      potential_mz = c(sapply(formulas, function(x) eval(parse(text = x))))
      names(potential_mz) = negative_adduct$adduct_name
    }else{
      print("Please enter the correct ion mode from following: 'positive' or 'negative'")
    }
  closest_peaks = lapply(potential_mz, function(x){
    id = which.min(abs(x-mz_list))
    closest_mz = mz_list[id]
    error_mz = (x-mz_list)[id]
    error_ppm = abs(error_mz)/matrix_molecule*1e6
    closest_peak_intensity = mean(mass_matrix[,id])
    closest_return = c(matched_mz_id= id,
                      closest_mz = closest_mz,
                      error_mz = error_mz,
                      error_ppm = error_ppm,
                      closest_peak_intensity =  closest_peak_intensity)
    return(closest_return)
  })
  # Calculate full_width at half height
  fwhm = matrix_molecule/mass_resolution
  intensities = colSums(mass_matrix)
  potential_matrix_peak_indices <- which(intensities >= quantile(intensities, probs =quantile))
  
  filtered_closest_peaks = as.data.frame(do.call(rbind,lapply(closest_peaks, function(x){
    if((x[1] %in% potential_matrix_peak_indices) & abs(x[3])<=3*fwhm){
      return(x)
    }else{
      return(NULL)
    }
  })))
  
  if(nrow(filtered_closest_peaks) == 0){
    stop("Warnning: no peaks found with given conditions consider increase mass accuracy or decrease quantile.
         Make sure you enter the correct matrix as well")
  }else{
    # print result
    require(ggplot2)
    require(fitdistrplus)
    # Plot histogram with fitted distribution
    print("Plotting matched adducted molecules with their corresponding mass error")
    plotdf  = cbind(filtered_closest_peaks,
                    adduct_status = rownames(filtered_closest_peaks))
    plot = ggplot(plotdf, aes(x = error_mz, y = closest_peak_intensity)) +
      geom_histogram(stat = "identity", fill = "skyblue", color = "black", bins = 10) +
      labs(x = "Error (m/z, unit: u)", y = "Peak Intensity") +
      geom_text(aes(y = closest_peak_intensity, label = adduct_status),position = "stack", vjust = -0.5, size = 3.5, color = "black")+
      geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
      theme_minimal()
    suppressWarnings(print(plot))
    print("Finished")
    return(plotdf)
  }
}





