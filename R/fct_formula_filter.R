#' @formula_filter
#' The Formula Filter function
#' #' filters out elements that are deemed not natural, or only allows elements set by the user.
#'  see the publication  "Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry"
#'  Citation & url: Kind, T., Fiehn, O. Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry. BMC Bioinformatics 8, 105 (2007). https://doi.org/10.1186/1471-2105-8-105. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-105
#' @param df database in the form of a dataframe to be searched against and filtered.
#' @param elements a vector of elements to include in the filtering. By default will exclude anything not in this list.

#### Filter by formula ####
formula_filter <- function(df, elements = NULL) {
  if (is.null(elements)) {
    elements <- c("H", "C", "N", "O", "S", "Cl", "Br",
                  "F", "Na", "P", "I")
  }


  # Elements to allow
  allowed_elements <- elements

  # Function to check if a formula contains only allowed elements
  is_formula_valid <- function(formula) {
    elements <-
      stringr::str_extract_all(formula, "[A-Z][a-z]*")[[1]] # defines an element as a Uppercase immediately followed by none or more lowercase letters.
    all(elements %in% allowed_elements) # Then checks if they are in the vector of allowed elements.
  }

  # Filter rows based on allowed elements
  filtered_df <- df %>%
    filter(sapply(formula, is_formula_valid))

  return(filtered_df)

}
