# Added NA to the end of a vector 
pad <- function(vec, n){
  length_of_vec = length(vec)
  return(c(vec, rep(NA, n-length_of_vec)))
}

# Sort a vector by it's last n characters
sort_by_last_n_char <- function(vec_to_sort, n=1){
  vec_to_sort = vec_to_sort[order(
    substr(vec_to_sort, nchar(vec_to_sort)-n+1, nchar(vec_to_sort))
  )]
  return(vec_to_sort)
}

# If out if classes are balanced
is_balanced <- function(sample_vec, meta_data_df){
  class_vec = meta_data_df[sample_vec, "Class"]
  class_table = table(class_vec)
  max_length = max(as.numeric(class_table))
  return(all(class_table==max_length) & (length(class_table)==2))
}