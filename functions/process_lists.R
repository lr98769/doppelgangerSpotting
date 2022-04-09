# Converts a list of vectors to a data frame (each vec is a row)
convert_list_to_df <- function(list_to_be_converted){
  max_length_of_vec = max(sapply(list_to_be_converted, length))
  list_to_be_converted = lapply(list_to_be_converted, pad, 
                                n=max_length_of_vec)
  list_to_be_converted_df =  do.call(rbind,list_to_be_converted)
  return(data.frame(list_to_be_converted_df))
}

# Sorts a list by the length of it's elements
sort_list_by_length_of_el <- function(list_to_sort, descreasing=TRUE){
  length_list = sapply(list_to_sort, length)
  list_to_sort = list_to_sort[order(length_list, decreasing = TRUE)]
  return(list_to_sort)
}

# Converts experiment plan (list of list) to a dataframe
convert_experiment_plan_to_df <- function(experiment_plan){
  train_valid_list = experiment_plan[["train_valid_sets"]]
  
  # Change from nested list to non-nested-list
  new_list_of_sets = list()
  for (train_valid_set in names(train_valid_list)){
    train = train_valid_list[[train_valid_set]][["train"]]
    valid = train_valid_list[[train_valid_set]][["valid"]]
    
    train_name = paste(train_valid_set, "train", sep=".")
    valid_name = paste(train_valid_set, "valid", sep=".")
    new_list_of_sets[[train_name]] = train
    new_list_of_sets[[valid_name]] = valid
  }
  
  # Pad each vector with NA
  max_length = max(sapply(new_list_of_sets, length))
  new_list_of_sets = lapply(new_list_of_sets, pad, n=max_length)
  
  # Convert list to dataframe
  return_df = as.data.frame(new_list_of_sets)
  
  return(return_df)
}