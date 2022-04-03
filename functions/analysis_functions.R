
dfSetDifference <- function(df1, df2){
  setdifference = c()
  for (i in 1:nrow(df1)){
    if (!any(df2$Sample1 == df1[i, "Sample1"] & 
             df2$Sample2 == df1[i, "Sample2"])){
      setdifference = c(setdifference, i)
    }
  }
  return(df1[setdifference,])
}
getDoppelnNonDoppelSamples <- function(doppel_result, metadata, batchname){
  samples = unique(doppel_result$PPCC_df[doppel_result$PPCC_df$DoppelgangerLabel=="Doppelganger", "Sample1"])
  samples = c(samples, unique(doppel_result$PPCC_df[doppel_result$PPCC_df$DoppelgangerLabel=="Doppelganger", "Sample2"]))
  all_samples_of_batch = rownames(metadata[metadata$Batch==batchname,])
  doppel_samples = intersect(all_samples_of_batch, samples)
  non_doppel_samples = setdiff(all_samples_of_batch, doppel_samples)
  # Output in df form
  return_list = list(doppel_samples, non_doppel_samples)
  ## Compute maximum length
  max_length = max(sapply(return_list, length))
  ## Add NA values to list elements
  return_list = lapply(return_list, function(v) { 
    c(v, rep(NA, max_length-length(v)))
  })
  # Column bind them
  return_df = do.call(cbind, return_list)
  colnames(return_df) = c("Doppelganger Samples", "Non-Doppelganger Samples")
  return(return_df)
}

# Oversample the raw and meta data automatically
oversample_batch <- function(raw_data, meta_data, seed=2021){
  return_list = list()
  # get batch sizes
  batch_sizes = table(meta_data$Batch)
  # Get number of samples to over sample
  add_num_samples = max(batch_sizes) - min(batch_sizes)
  # Get batch to over sample
  min_batch = names(which.min(batch_sizes))
  # Get min_batch sample names
  min_batch_samples = rownames(meta_data[meta_data$Batch==min_batch, ])
  # Over sample some samples
  set.seed(seed)
  if (add_num_samples < length(min_batch_samples)){
    new_sample_names = sample(min_batch_samples, add_num_samples, 
                              replace = FALSE) 
  }
  else {
    print("Error: This function does not support oversampling of over twice the smaller batch's size")
    return()
  }
  
  new_samples_meta = meta_data[new_sample_names, ]
  rownames(new_samples_meta) = paste(rownames(new_samples_meta), "dup", sep="_")
  return_list[["meta_data"]] = as.data.frame(rbind(meta_data,new_samples_meta))
  
  new_samples_raw = raw_data[, new_sample_names]
  colnames(new_samples_raw) = paste(colnames(new_samples_raw), "dup", sep="_")
  return_list[["raw_data"]] = as.data.frame(cbind(raw_data,new_samples_raw))
  
  return(return_list)
}

# Removes added duplicate samples and sample pairs
remove_all_dup <-function(doppel_result){
  doppel_return = doppel_result
  # Get column names with "_dup"
  duplicates = grep('_dup', colnames(doppel_return$Batch_corrected), value=TRUE)
  doppel_return$Batch_corrected = doppel_return$Batch_corrected[, !(colnames(doppel_return$Batch_corrected) %in% duplicates)]
  doppel_return$PPCC_matrix = doppel_return$PPCC_matrix[, !(colnames(doppel_return$PPCC_matrix) %in% duplicates)]
  doppel_return$PPCC_matrix = doppel_return$PPCC_matrix[!(rownames(doppel_return$PPCC_matrix) %in% duplicates),]
  doppel_return$PPCC_df = doppel_return$PPCC_df[!(doppel_return$PPCC_df$Sample1 %in% duplicates) & !(doppel_return$PPCC_df$Sample2 %in% duplicates), ]
  return(doppel_return)
}
