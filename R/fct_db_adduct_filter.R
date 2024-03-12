#' @db_adduct_filter
#' filters a db and only includes adducts set by the user
#' @param db this takes in a db file as defined by CompoundDB
#' @param adduct takes a vector of adducts of the form c("M+H", "M+NH4", etc)
#' @param polarity can be set to "neg" or "pos" depending on acquisition mode.



#### Filter DB by polarity and adducts ####
db_adduct_filter <- function(db, adduct, polarity = "neg") {
  # only include adducts from either neg or pos polarity
  if (polarity == "neg") {
    neg_adducts_1 <-
      c(
        "M-H2O-H",
        "M-H",
        "M+Na-2H",
        "M+Cl",
        "M+K-2",
        "M+FA-H",
        "M+Hac-H",
        "M+Br",
        "M+TFA-H",
        "2M-H",
        "2M+FA-H",
        "2M+Hac-H",
        "3M-H"
      )
    pol <- neg_adducts_1
  }
  else if (polarity == "pos") {
    pos_adducts_1 <-
      c(
        "M+H",
        "M+NH4",
        "M+Na",
        "M+CH3OH+H",
        "M+K",
        "M+ACN+H",
        "M+2Na-H",
        "M+IsoProp+H",
        "M+ACN+Na",
        "M+2K+H",
        "M+DMSO+H",
        "M+2ACN+H",
        "M+IsoProp+Na+H",
        "2M+H",
        "2M+NH4",
        "2M+Na",
        "2M+K",
        "2M+ACN+H",
        "2M+ACN+Na"
      )
    pol <- pos_adducts_1
  } else{
    stop("Invalid polarity. Choose 'neg' or 'pos'.")
  }

  # get rid of spaces in the adduct names
  # in col names
  db <- db %>%
    rename_all( ~ gsub(" ", "", .))

  # Filter the db by polarity
  db <- db %>%
    dplyr::select(formula, exactmass, isomers, isomers_inchikey, isomers_names, pol)

  # in adduct entry
  adduct <- gsub(" ", "", adduct)

  check_and_truncate_adduct_vector <- function(adduct, db) {
    element_exists <- adduct %in% colnames(db)
    missing_elements <- adduct[!element_exists]
    if (length(missing_elements) > 0) {
      for (missing_element in missing_elements) {
        cat(
          "Adduct",
          missing_element,
          "is not in the DB, it has been removed from the search.\n"
        )
      }
      truncated_adduct <- adduct[element_exists]
      return(truncated_adduct)
    } else {
      return(adduct)
    }
  }

  adduct <- check_and_truncate_adduct_vector(adduct, db)

  db_filtered <- db %>%
    dplyr::select(formula,
           exactmass,
           isomers,
           isomers_inchikey,
           isomers_names,
           any_of(adduct)) %>%
    as.data.frame()
  return(db_filtered)
}
