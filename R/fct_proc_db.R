#' @proc_db
#' This is the main function. It does a few things. Probably best if you just read the comments below.
#' @observed_df The dataframe of experimental data. To be searched against, i.e., row_id and m/z
#' @reference_df A dataframe of a database to be aligned. There is a very specific format to this database.
#' @ppm_threshold The ppm threshold that is to be used to search through the database.
#'

#### Search an mz list against a data base ####
proc_db <- function(observed_df,
                    reference_df,
                    ppm_threshold = 10) {
  # create an empty list to store matched results.
  result_list <- list()

  # Check if observed_df has only one row
  if (nrow(observed_df) == 1) {
    # Create a dummy entry with all values set to 0
    dummy_row <-
      data.frame(row_id = 0, mz = 0)  # Modify this line based on your column names

    # Combine the dummy entry with the original observed_df
    observed_df <- rbind(observed_df, dummy_row)
  } else if (nrow(observed_df) == 0) {
    # Handle the case when observed_df is empty
    message(
      "Warning: No entries detected in input mz list.\n Please check format of input list,\n it must contain row_id and mz as the headers"
    )
    return(NULL)
  }

  # Calculate mz range of observed_df
  calculate_bounds <- function(input_df) {
    lower_bound <- min(input_df$mz, na.rm = TRUE)
    upper_bound <- max(input_df$mz, na.rm = TRUE)
    bounds <-
      list(lower_bound = lower_bound, upper_bound = upper_bound)
    return(bounds)
  }
  # extract out bounds
  lower_bound <- calculate_bounds(observed_df)$lower_bound
  upper_bound <- calculate_bounds(observed_df)$upper_bound

  # Function to calculate ppm error as a valve
  ppm_error <- function(observed_mz, reference_mz, ppm) {
    abs_diff_ppm <-
      abs(observed_mz - reference_mz) / abs(reference_mz) * 1e6
    if (abs_diff_ppm <= ppm) {
      return(abs_diff_ppm)
    }
    else{
      return("Out")
    }
  }

  # Function to calculate ppm range and check if mz values are within the range
  # Returns TRUE if match is found and false if no match.
  ppm_range_match <- function(observed_mz, reference_mz, ppm) {
    abs_diff_ppm <-
      abs(observed_mz - reference_mz) / abs(reference_mz) * 1e6
    abs_diff_ppm <= ppm
  }

  #  For loop to go through each column of the reference_df that is provided.
  #  probably a good idea to filter reference_df to only adducts that you want before putting it into this function.
  # the -c(1:4) essentially makes it loop over the numeric portions of the DBs.

  for (col_name in names(reference_df)[-c(1:5)]) {
    result_col <- list()
    for (i in seq_len(nrow(observed_df))) {
      # Check if the reference mz is within the range
      within_range <-
        reference_df[[col_name]] >= lower_bound &
        reference_df[[col_name]] <= upper_bound

      # Condition 1: Only proceed if there are values within the range
      if (!any(within_range)) {
        next
      }

      # check for matches
      matches <- ppm_range_match(
        observed_mz = observed_df$mz[i],
        reference_mz = reference_df[[col_name]],
        ppm = ppm_threshold
      )

      # Condition 2: skip to the next iteration if no match
      if (!any(matches)) {
        next
      }

      # Condition 3: index the matches in the reference df
      # which keeps things that are TRUE.
      # This allows referencing to be quicker below.
      matching_indices <- which(matches)
      # print(matching_indices)

      # Calculate ppm error for each match.
      error <- mapply(
        ppm_error,
        observed_mz = observed_df$mz[i],
        reference_mz = reference_df[[col_name]][matching_indices],
        ppm = ppm_threshold
      )

      # Extract out relevant info for each match.
      filtered_matches <- data.frame(
        ID = observed_df$row_id[i],
        Match = matches[matching_indices],
        observed_mz = observed_df$mz[i],
        Reference_mz = reference_df[[col_name]][matching_indices],
        Error = error,
        Adduct = col_name,
        # DB_ID = reference_df$compound_id[matching_indices],
        Formula = reference_df$formula[matching_indices],
        Exactmass = reference_df$exactmass[matching_indices],
        Isomers = reference_df$isomers[matching_indices],
        InchiKeys = reference_df$isomers_inchikey[matching_indices],
        IsomerNames = reference_df$isomers_names[matching_indices]
      )
      # store each results df for each mz value in observed_df
      result_col[[i]] <- filtered_matches
    }
    # store results df in a list for each mz value for each adduct
    result_list[[col_name]] <- result_col
  }

  # Combine the individual matches per adduct df from the list into a dataframe
  combined_results <-
    do.call(rbind, lapply(result_list, function(result_col)
      do.call(rbind, result_col)))
}
