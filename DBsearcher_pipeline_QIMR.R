#' @DataBaseSearchingPipeline
#' Authored by: Christopher CJ Fitzgerald, July 2023.
#' For use of this pipeline in publications please ask Christopher Fitzgerald, Bioinformatician @MetabolomicsAustralia. github - https://github.com/ChemCharles

#' The below code is a test pipeline used to develop the functions in the db explorer application. This pipeline allows the user to search a list of m/z values derived from HRAM experiments against a database.
#'
#' The pipeline is made up of three functions. These are;
#'
#' The db_adduct_filter function and formula_filter function are aimed at filtering the databases. This makes the matching more efficient as it reduces the size of the databases by removing things that the user doesn't want to search. The final function proc_db does the matching between the databases and returns a list of matches, if any, to your input list.
#'
#' The structure of the databases are normalised by grouping all isomers (defined by identical formula and exactmass values), this is done by the associated Clean_db_pipeline and originally CompoundDB. The auddcts for any given exact mass are then calculated. This database searching pipeline is not designed to give a 1-1 match to any given mz value, however, may give a 1-many match . Put in terms of the levels of confidence in metabolite identification as originally defined by the International Metabolomics Standards Initiative, This pipeline would give a level 4 (very low confidence) match.

# load libraries and functions.
library(dplyr)
library(data.table)
library(tidyr)

source("./R/fct_db_adduct_filter.R")
source("./R/fct_formula_filter.R")
source("./R/fct_proc_db.R")

#### Load the Cleaned and summarized DB ####
Chebi_db     <- readRDS("./inst/db_files/Chebi_1_names.rds")
HMDB_db      <- readRDS("./inst/db_files/HMDB_1_names.rds")
LIPIDMAPS_db <- readRDS("./inst/db_files/Lipidmaps_1_names.rds")

##### DB searching Pipeline ####

# This first pipeline lets you query one database at a time.

# File of m/z values to query against

test_mz_df <- read.csv("./inst/Test_input/data_peaks.csv")

# Set the db that you want to search against
db <- HMDB_db

# set which adducts you want to search for
test_add_pos <- c("M+H")
# Note; "M+NH4","M+Na","M+CH3OH+H","M+K" etc. Look at the formula filter func to get the rest of the possible adducts.

# set ppm error/threshold for searching
ppm_error <- 5

# Three main steps relates to the three main functions

# Steps 1) & 2) are aimed at condensing the databases by applying 1) a filter to only consider the adducts that the user specifies. 2) Filtering the molecular formulas to contain only elements that the user specifies. # Step 3) This last function then does the database matching and searching.

# 1) Filter DB by adduct.
db_1 <- db_adduct_filter(db, test_add_pos, polarity = "pos")

# 2) only select natural elements
db_2 <- formula_filter(db_1)

# 3) search db against mz df return results
db_3 <- proc_db(test_mz_df, db_2, ppm_error)

# export the results
# write.csv(db_3, "test.csv")





