# This function converts the probes to ensemble id
# df = df with the column "Probe"
# 
convertProbesToEnsemble <- function(df, affy_attribute){
  if (!exists("ensembl")){
    ensembl <<- useEnsembl(biomart = "ensembl", 
                           dataset = 'hsapiens_gene_ensembl')
  }
  if (!"Probe" %in% colnames(df)){
    print("The column 'Probe' is missing from df")
    return()
  }
  return_list = list()
  mapping_table =  select(ensembl, 
                          keys= df$Probe, 
                          columns= c(affy_attribute,
                                     "ensembl_gene_id"),
                          keytype= affy_attribute)
  return_list[["mapping_table"]] = mapping_table
  # Ensure the each affy attribute maps to only 1 ensemble id 
  # (the smallest one) 
  # (Ensure that affy->ensemble is a one-to-one mapping)
  mapping_table = mapping_table[
    order(mapping_table[[affy_attribute]],
          mapping_table[["ensembl_gene_id"]]),]
  mapping_table = mapping_table[!duplicated(
    mapping_table[[affy_attribute]]),]
  
  # Remove all probes with no corresponding ensemble id
  df = df[df$Probe %in% mapping_table[[affy_attribute]],]
  # Replace rownames with the ensemble gene id
  rownames(mapping_table) = mapping_table[[affy_attribute]]
  df$ensembl_ID = mapping_table[df$Probe,"ensembl_gene_id"]
  return_list[["before_dedup"]] = df
  df$Probe = NULL
  # ensure order of variables stay the same for checking purposes
  current_order = unique(df$ensembl_ID)
  # Get the median signal of probes with the same ensemble id
  df = aggregate(df[,-ncol(df)], by = list(df[,"ensembl_ID"]) , median)
  rownames(df) = df$Group.1
  df$Group.1 = NULL
  df = df[current_order, ]
  return_list[["returned_df"]] = df
  return(return_list)
}
# Integrates file loading, replacing probes with ensemble id and caching of cleaned_data
getDataFile<- function(filename, affy_attribute, batch_name){
  cleaned_filename = file.path(cleaned_data_dir, basename(filename))
  if (!file.exists(cleaned_filename)){
    temp_df = read.csv(filename)
    temp_df$PROT = NULL # Remove PROT column
    # Convert probes to ensemble ids
    conversion_results <<- convertProbesToEnsemble(
      df = temp_df,
      affy_attribute = affy_attribute)
    temp_df = conversion_results$returned_df
    # Add H to the end so that we can differentiate it later
    colnames(temp_df) = paste(colnames(temp_df), batch_name, sep = "_")
    write.csv(temp_df, cleaned_filename)
    return(temp_df)
  } else {
    temp_df = read.csv(cleaned_filename, row.names = 1)
    return(temp_df)
  }
}
# Generates meta data table from an existing dataframe (Assumes each row is a different patient)
getMetaDataDataframe <-function(df, batch_name){
  meta = data.frame(row.names = colnames(df))
  meta$Patient_ID = rownames(meta)
  meta$Batch = batch_name
  meta$Class = substr(meta$Patient_ID, 1, 3)
  return(meta)
}

getGolubDataFile <- function(){
  clean_golub_path = file.path(cleaned_data_dir, "Leuk-GolubData.csv")
  golub_path = file.path(data_dir, "Leuk-GolubData.csv")
  if (!file.exists(clean_golub_path)){
    leuk_G = read.csv(golub_path)
    # Different chip from the rest, not found in biomart
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    if (!"hu6800.db" %in% installed.packages()){
      BiocManager::install("hu6800.db")
    }
    library(hu6800.db)
    mapping_table = select(hu6800.db, leuk_G$Probe, c("ENSEMBL"), keytype="PROBEID")
    mapping_table = data.frame(mapping_table)
    # Remove probes with no ensemble id
    mapping_table = mapping_table[!is.na(mapping_table$ENSEMBL),]
    # Ensure the each affy attribute maps to only 1 ensemble id 
    # (the smallest one) 
    # (Ensure that affy->ensemble is a one-to-one mapping)
    mapping_table = mapping_table[
      order(mapping_table[["PROBEID"]],
            mapping_table[["ENSEMBL"]]),]
    mapping_table = mapping_table[!duplicated(
      mapping_table[["PROBEID"]]),]
    # Remove variables that have no ensemble id
    leuk_G = leuk_G[leuk_G$Probe %in% mapping_table$PROBEID,]
    rownames(mapping_table) = mapping_table$PROBEID
    # Map Probes to ensembl id
    leuk_G$ensembl_ID = mapping_table[leuk_G$Probe,"ENSEMBL"]
    leuk_G$Probe = NULL
    # ensure order of variables stay the same for checking purposes
    current_order = unique(leuk_G$ensembl_ID)
    # Get the median signal of probes with the same ensemble id
    leuk_G = aggregate(leuk_G[,-ncol(leuk_G)], by = list(leuk_G[,"ensembl_ID"]) , median)
    rownames(leuk_G) = leuk_G$Group.1
    leuk_G$Group.1 = NULL
    leuk_G = leuk_G[current_order, ]
    # Rename cols so that we can differentiate the batches
    colnames(leuk_G) = paste(colnames(leuk_G), "G", sep="_")
    write.csv(leuk_G, clean_golub_path)
  } else {
    leuk_G = read.csv(clean_golub_path, row.names = 1)
  }
  return(leuk_G)
}