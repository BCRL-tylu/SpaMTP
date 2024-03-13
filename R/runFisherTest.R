#' @runFisherTest
#' @description
#' This the function used to compute the exact fisher test for over-represntation based pathway analysis
#'

#' @param Analytes IDs in following format: "db source: ID", for example ""ensembl:ENSG00000135679", "hmdb:HMDB0000064""
#' @param analyte_type = "metabolites" or "genes" or "both"
#' @param max_path_size The max number of metabolites in a specific pathway
#' @param min_path_size The min number of metabolites in a specific pathway
#' @param alternative The hypothesis of the fisher exact test

#' @return a dataframe with the relevant pathway information
runFisherTest = function (analytes, analyte_type = "metabolites",
           alternative = "greater", min_path_size = 5, max_path_size = 150)
{
  require(dplyr)
  now <- proc.time()
  print("Fisher Testing ......")
  pathwayRampId <- rampId <- c()
  # Get the RaMP ids for metabolites/genes

  print("Loading files ......")
  analytehaspathway = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/analytehaspathway.rds"))
  pathway = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/pathway.rds"))
  source = readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/source.rds"))
  chem_props =readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/chem_props.rds"))
  analyte =readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/analyte.rds"))
  pathway =readRDS(paste0(dirname(system.file(package = "SpaMTP")),"/data/pathway.rds"))

  print("Loading files finished!")

  if(analyte_type == "metabolites"){
    source = source[which(grepl(source$rampId, pattern = "RAMP_C") == T),]
    analytehaspathway = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_C") == T),]
    analyte = analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),]
  }else if(analyte_type == "genes"){
    source = source[which(grepl(source$rampId, pattern = "RAMP_G") == T),]
    analytehaspathway = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_G") == T),]
    analyte = analyte[which(grepl(analyte$rampId, pattern = "RAMP_G") == T),]
  } else {
    stop("analyte_type was not specified correctly.  Please specify one of the following options: metabolites, genes")
  }

  # Pathway enrichment
  if(analyte_type == "metabolites"){
    #####################################
    #####Metabolic pathway analysis######
    #####################################
  print("Begin metabolic pathway analysis ......")
  analytes_rampids = unique(source$rampId[which(source$sourceId %in% analytes)])
  # (1) Get candidate pathways
  # Get all analytes and number of analytes within a specific pathway

  pathway_rampids = analytehaspathway[which(analytehaspathway$rampId %in% analytes_rampids),] %>% rowwise() %>%
    mutate(ananlytes_id_list = list(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathwayRampId)])) %>%
    count(pathwayRampId,
          ananlytes_id_list,
          sort = T,
          name = "analytes_in_pathways")

  pathway_rampids = pathway_rampids %>%  mutate(screened_analytes = list(ananlytes_id_list[which(ananlytes_id_list %in% analytes_rampids)]))
  # (2) Get total analytes in each pathway
  total_in_pathways = analytehaspathway[which(analytehaspathway$pathwayRampId %in%  pathway_rampids$pathwayRampId),] %>%
    count(pathwayRampId,
          name = "total_in_pathways")
  # (3) Creare a df that store the enrichment square for each pathways
  enrichment_df = merge(total_in_pathways,
                        pathway_rampids,
                        by = "pathwayRampId") %>%
    mutate(total_analytes = length(unique(analytes_rampids)))
  # (4) Conduct pathway enrichment
  total_in_selected_pathways = length(unique(analytehaspathway$rampId))
  print("Calculating p value......")
  enrichment_df = enrichment_df %>% rowwise() %>%mutate(p_val = fisher.test(matrix(c(analytes_in_pathways,# Detected metabolites in pathway
                                                                        total_analytes-analytes_in_pathways,# Detected metabolites not in pathway
                                                                        total_in_pathways-analytes_in_pathways, # Pathway elements not detected
                                                                        total_in_selected_pathways-total_in_pathways-total_analytes+analytes_in_pathways), 2,2),
                                                                        alternative = alternative)$p.value)
  enrichment_df = cbind(enrichment_df,
                  fdr= p.adjust(enrichment_df$p_val, method="fdr"))%>% mutate(background_analytes_number = total_in_selected_pathways)
  print("P value obtained")
  # (5) Append pathway information to the original df
  enrichment_df_with_info = merge(enrichment_df,
        pathway[which(pathway$pathwayRampId %in% enrichment_df$pathwayRampId),],
        by = "pathwayRampId") %>% filter(!duplicated(pathwayName))

  # (6) Append metabolites information to the original df
  # Paste back the original Ids
  print("Formatting outputs ......")
  metabolites_information  = do.call(rbind,lapply(enrichment_df_with_info$ananlytes_id_list, function(x){
    temp_list = unique(chem_props[which(chem_props$ramp_id %in% x),]$chem_source_id)
    return(paste(temp_list[which(temp_list %in% analytes)],
                 collapse = ";"))
  }))
  # (7) Reduce the dataframe with respected to the User input pathway size
  enrichment_df_with_both_info = cbind(enrichment_df_with_info,
                                       metabolites_information) %>% filter(total_in_pathways<= max_path_size &
                                                                           total_in_pathways>= min_path_size)
  print("Done")
  return(enrichment_df_with_both_info)
  }else if(analyte_type == "genes"){
    #####################################
    ########Gene pathway analysis########
    #####################################
    print("Begin gene pathway analysis ......")
    analytes_rampids = unique(source$rampId[which(source$sourceId %in% analytes)])
    # (1) Get candidate pathways
    # Get all analytes and number of analytes within a specific pathway
    pathway_rampids = analytehaspathway[which(analytehaspathway$rampId %in% analytes_rampids),] %>% rowwise() %>%
      mutate(ananlytes_id_list = list(analytehaspathway$rampId[which(analytehaspathway$pathwayRampId == pathwayRampId)])) %>%
      count(pathwayRampId,
            ananlytes_id_list,
            sort = T,
            name = "analytes_in_pathways")
    pathway_rampids = pathway_rampids %>%  mutate(screened_analytes = list(ananlytes_id_list[which(ananlytes_id_list %in% analytes_rampids)]))
    # (2) Get total analytes in each pathway
    total_in_pathways = analytehaspathway[which(analytehaspathway$pathwayRampId %in%  pathway_rampids$pathwayRampId),] %>%
      count(pathwayRampId,
            name = "total_in_pathways")
    # (3) Creare a df that store the enrichment square for each pathways
    enrichment_df = merge(total_in_pathways,
                          pathway_rampids,
                          by = "pathwayRampId") %>%
      mutate(total_analytes = length(unique(analytes_rampids)))
    # (4) Conduct pathway enrichment
    total_in_selected_pathways = length(unique(analytehaspathway$rampId))
    print("Calculating p value......")
    enrichment_df = enrichment_df %>% rowwise() %>%mutate(p_val = fisher.test(matrix(c(analytes_in_pathways,# Detected metabolites in pathway
                                                                        total_analytes-analytes_in_pathways,# Detected metabolites not in pathway
                                                                        total_in_pathways-analytes_in_pathways, # Pathway elements not detected
                                                                        total_in_selected_pathways-total_in_pathways-total_analytes+analytes_in_pathways), 2,2),
                                                                        alternative = alternative)$p.value)
    enrichment_df = cbind(enrichment_df,
                          fdr= p.adjust(enrichment_df$p_val, method="fdr")) %>% mutate(background_analytes_number = total_in_selected_pathways)
    print("P value obtained")
    # (5) Append pathway information to the original df
    enrichment_df_with_info = merge(enrichment_df,
                                    pathway[which(pathway$pathwayRampId %in% enrichment_df$pathwayRampId),],
                                    by = "pathwayRampId") %>% filter(!duplicated(pathwayName))

    # (6) Append metabolites information to the original df
    print("Formatting outputs ......")
    genes_information  =  do.call(rbind,lapply(enrichment_df_with_info$ananlytes_id_list, function(x){
      temp_list = unique(source[which(source$rampId %in% x),]$sourceId)
      gene_symbol = sub("gene_symbol:","",temp_list[which(grepl(templist,
                                          pattern = "gene_symbol"))])
      return(paste(temp_list[which(temp_list %in% analytes)],
                   collapse = ";"))
    }))
    # (7) Reduce the dataframe with respected to the User input pathway size
    enrichment_df_with_both_info = cbind(enrichment_df_with_info,
                                         metabolites_information)%>% filter(total_in_pathways<= max_path_size &
                                                                           total_in_pathways>= min_path_size)
    print("Done")
    return(enrichment_df_with_both_info)
    }
}
