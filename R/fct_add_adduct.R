#' @description
#' This function pretty much gets the db ready for use in the data pipeline.
# It calculates the adducts from the exact mass value (molecular weight). Then appends them to the cleaned database dataframe. The adduct_charge can be set between 1-3, the adduct_df is the file that exist in the inst directory as adduct_file.rds. This file contains 47 common adducts accross positive and negative polarity.
#'
#' @param db standard format of a db table that is generated by the CompoundDb package. See inst/db_files/nativeDBfiles for examples. It's probably a good idea to clean up the database i.e., remove 0 values and NA's before proceeding with the calculation below. Use the clean_data function to do this.
#' @param adduct_df this is a file of possible adducts, see inst/adduct_file.rds for file. For original reference;
#' "Kind, T., Fiehn, O. Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry. BMC Bioinformatics 8, 105 (2007). https://doi.org/10.1186/1471-2105-8-105" for original data source.
#' @param adduct_charge This is a number betwee 1-3. representing the lowest to highest charge states of adducts that can be calculated from the below function.

#load libraries
library(roxygen2)
library(dplyr)
library(tidyr)
library(data.table)

# function
Add_Adduct <- function(db, adduct_df, adduct_charge) {
  # Make sure the adduct_charge is numeric
  if (!is.numeric(adduct_charge)) {
    stop("Input must be a numeric value.")
  }
  # Set the max value adduct_charge can be
  if (adduct_charge > 3) {
    stop("Input exceeds the maximum value, must be an absolute integer between 1-3.")
  }

  db_strip <- function(data_base) {
    df <- data_base %>%
      select(compound_id, formula, exactmass, inchikey, name)
    return (df)
  }

  # strip db to only take compound_id, formula, exactmass, inchikey and name.

  stripped_db <- db_strip(db)

  # print(head(stripped_db))
  # Filter the adduct to charge value
  filtered_df <- adduct_df %>%
    filter(abs(charge) <= adduct_charge)

  # Get every possible combination of data base exactmass value to adduct
  cross_df <- stripped_db %>%
    crossing(filtered_df)

  # print(head(cross_df))

  # Calculate the adduct masses from the exactmass values
  adductmass_df <- cross_df %>%
    mutate(adduct_mz = (exactmass * mult) + add_mass)

  # print(head(adductmass_df))
  # append the adduct values to the exactmass values
  # retain db ID and formula exactmass values of original entry.
  adduct_df <- adductmass_df %>%
    pivot_wider(
      id_cols = c(compound_id, formula, exactmass, inchikey, name),
      names_from = adduct_name,
      values_from = adduct_mz
    )

  # print(head(adduct_df))

  iso_grp_fn_isomers <- function(adduct_df) {
    grouped_data <-  adduct_df %>%
      select(-inchikey,-name) %>%
      group_by(across(-compound_id)) %>%
      summarise(isomers = paste(compound_id, collapse = "; ")) %>%
      ungroup()
    return(grouped_data)
  }

  iso_grp_fn_inchikey <- function(adduct_df) {
    library(dplyr)

    grouped_data <-  adduct_df %>%
      select(-compound_id,-name) %>%
      group_by(across(-inchikey)) %>%
      summarise(isomers_inchikey = paste(inchikey, collapse = "; ")) %>%
      ungroup()
    return(grouped_data)
  }

  iso_grp_fn_names <- function(adduct_df) {
    library(dplyr)

    grouped_data <-  adduct_df %>%
      select(-compound_id,-inchikey) %>%
      group_by(across(-name)) %>%
      summarise(isomers_names = paste(name, collapse = "; ")) %>%
      ungroup()
    return(grouped_data)
  }

  iso_df <- iso_grp_fn_isomers(adduct_df)
  inchi_df <- iso_grp_fn_inchikey(adduct_df)
  names_df <- iso_grp_fn_names(adduct_df)

  fin_df <- merge(iso_df, inchi_df)
  fin_df <- merge(fin_df, names_df)

  fin_df %>% select()

  df_out <- fin_df %>%
    select(formula,
           exactmass,
           isomers,
           isomers_inchikey,
           isomers_names,
           everything())
}
