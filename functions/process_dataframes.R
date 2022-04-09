# Generates meta data data frame from raw data data frame
generate_meta_data <- function(df){
  sample_names = colnames(df)
  classes = sapply(sample_names, function(sample_name){
    split_name = strsplit(sample_name, "_")[[1]]
    class_vec = split_name[2:length(split_name)]
    return(paste(class_vec, collapse="_"))
  })
  meta_df = data.frame(Class  = classes,
                       Patient_ID = sample_names,
                       Batch = "ccle"
  )
  row.names(meta_df) = sample_names
  return(meta_df)
}

# Sets the first row of a data frame to header
set_first_row_as_header <- function(df) {
  names(df) = as.character(unlist(df[1,]))
  return(df[-1,])
}


# Converts a dataframe to a list of vectors (each row is a new vec)
convert_df_to_list <- function(df_to_be_converted){
  return_list = list()
  for (rowname in rownames(df_to_be_converted)){
    vec_row = as.character(as.vector(df_to_be_converted[rowname,]))
    vec_row = vec_row[!is.na(vec_row)]
    return_list[[rowname]] = vec_row
  } 
  return(return_list)
}



