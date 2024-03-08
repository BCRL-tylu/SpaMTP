#' @CleanDB_pipeline
#' @description
#'
#'

#load libraries
library(CompoundDb)
library(dplyr)
library(tidyr)
library(data.table)
library(gtools)

# load helper functions
source("./R/fct_add_adduct.R")
source("./R/fct_clean_data.R")

#Step 1) ####
#use the compound_tbl_sdf function from CompoundDb to extract databases into a uniform format. Please see the documentation for CompoundDb. WARNING THIS FUNCTION TAKES ALOT OF COMPUTE TIME! so do it over lunch for larger databases!.
# @https://bioconductor.org/packages/release/bioc/html/CompoundDb.html
# use BiocManager::install("CompoundDb") to install.

# Note you need to download the SDF files from database websites. A good resource for this is the fiehn group's repository https://mona.fiehnlab.ucdavis.edu/downloads or going directly to the db website.

# save the db table as a rds file.
# saveRDS(chebi_file, file = "chebi_db.rds")
loaded_data <- readRDS("./inst/db_files/nativeDBfiles/HMDB_db.rds")

# Step 2) ####
# Clean the databases

adduct_df <- readRDS("./inst/adduct_file.rds")

db <- loaded_data

#### Clean up data frames with the below func. ####
# Gets rid of any NA entries and where exactmass is equal to zero.
# use the clean data function
db_clean <- clean_data(db)

#### This function pretty much gets the db ready for use in the data pipeline.
# It calculates the adducts from the exact mass value. Then appends them to the cleaned database dataframe. The adduct_charge can be set between 1-3, the adduct_df is the file that exist in the Adduct_file. This file contains 47 common adducts accross positive and negative polarity.
# Use add_adduct function to add adduct columns to db
# Below I have entered the cleaned up db from above, the adduct data file and wanted to calculate all adducts of charge 1.
db_w_adducts <- Add_Adduct(db_clean,adduct_df, 1)

#Save the cleaned database file
# This data file are the ones used in the DBsearcher_pipeline script and in the shiny application. The naming structure I have given for clarity is the Database_Charge_ImportantMetaData.

# Step 3)
# Export cleaned DB to use in pipeline####
saveRDS(db_w_adducts, file = "GNPS_1_names.rds")


