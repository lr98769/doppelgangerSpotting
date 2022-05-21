# Returns a list of important dataframes for doppelganger verification
# experiment planning
# 1. Meta data dataframe
# 2. Doppel Map dataframe (a dataframe with a rowname of all
#    doppelganger samples and each subsequent column is its doppelganger
#    partner)
generate_planning_dataframes <- function(meta_data_df,
                                    doppel_results){
  return_list = list()
  # Sort meta data
  meta_data_df = meta_data_df[
    sort_by_last_n_char(row.names(meta_data_df)),]
  
  ppcc_df = doppel_results$PPCC_df
  
  # Data frame with only doppelgangers 
  doppelganger_df = ppcc_df[
    (ppcc_df$DoppelgangerLabel =="Doppelganger"),]
  
  # Rearrange data frame to list:
  # - Sample_1: [Sample_2, Sample_3]
  all_doppelganger_samples = unique(
    c(doppelganger_df$Sample1,
      doppelganger_df$Sample2))
  all_doppelganger_samples = sort_by_last_n_char(all_doppelganger_samples)
  rearranged_doppel = list()
  for (doppelganger_sample in all_doppelganger_samples){
    is_doppelgangers_with_sample2 = doppelganger_df[
      doppelganger_df$Sample1 == doppelganger_sample,
      "Sample2"]
    is_doppelgangers_with_sample1 = doppelganger_df[
      doppelganger_df$Sample2 == doppelganger_sample,
      "Sample1"]
    is_doppelgangers_with = unique(c(is_doppelgangers_with_sample2,
                                     is_doppelgangers_with_sample1))
    rearranged_doppel[[doppelganger_sample]] = is_doppelgangers_with
  }
  
  # Convert Sample_1: [Sample_2, Sample_3] to dataframe
  rearranged_doppel_df = convert_list_to_df(rearranged_doppel)
  colnames(rearranged_doppel_df) = rep("DoppelPartner", ncol(rearranged_doppel_df))
  
  # Reorder the rows according to number of doppel partners
  rearranged_doppel_df$num_doppel_partners = apply(rearranged_doppel_df, 1, function(row){
    return(sum(!is.na(row)))
  })
  rearranged_doppel_df = rearranged_doppel_df[order(rearranged_doppel_df$num_doppel_partners, 
                                                    decreasing = TRUE), ]
  rearranged_doppel_df$num_doppel_partners = NULL
  
  # Output relevant data frames
  return_list[["meta_data"]] = meta_data_df
  return_list[["doppel_mapping"]] = rearranged_doppel_df
  
  return(return_list)
}

# Returns list of doppelganger samples in validation:
# - These doppelganger samples maximize the number of doppelganger samples
#   in validation while also ensuring each of these samples are doppelgangers
#   with at least one other sample in training
# - The doppelganger samples not in this list would be added to training
find_largest_valid_doppel <- function(planning_list){
  # Convert doppel mapping df to a sorted list of doppel mappings
  doppel_mapping = planning_list$doppel_mapping
  doppel_mapping_list = convert_df_to_list(doppel_mapping)
  doppel_mapping_list = sort_list_by_length_of_el(doppel_mapping_list)
  
  # Keep track of doppel samples in training set
  doppel_sample_in_training = c()
  
  # Keep track of doppel samples in validation set
  doppel_sample_in_valid = c()
  
  # Starting from doppelganger samples with many matches
  for (doppel_sample in names(doppel_mapping_list)){
    # Just a check 
    # (this should not be triggered since names are unique)
    if (doppel_sample %in% doppel_sample_in_training){
      next
    }
    # If current sample is already in the validation set:
    # - We do not add it to training (that would cause duplicates) 
    # - We do not add it to validation (we have already added it)
    if (doppel_sample %in% doppel_sample_in_valid){
      next
    }
    # If we have not added this sample to training or validation set:
    # - We add it to the training set
    # - All its doppelganger sample partners will be added to the
    #   validation set
    doppel_sample_in_training = c(
      doppel_sample_in_training, 
      doppel_sample
    )
    doppel_sample_in_valid = c(
      doppel_sample_in_valid,
      doppel_mapping_list[[doppel_sample]]
    )
  }
  
  # Remove duplicates
  doppel_sample_in_valid = unique(doppel_sample_in_valid)
  
  # Sort by class
  doppel_sample_in_valid = sort_by_last_n_char(doppel_sample_in_valid)
  return(doppel_sample_in_valid)
}


# Separate a vector of samples into a lists of vectors
# where the name of each vector is the class
# Eg. {
# "class_1": c("sample1", "sample2"),
# "class_2": c("sample3", "sample4")
#}
create_class_list <- function(sample_vec, meta_data_df){
  meta_data_df = data.frame(meta_data_df)
  meta_data_df = meta_data_df[sample_vec,]
  
  class_list = list()
  unique_classes = unique(meta_data_df$Class)
  for (class in unique_classes){
    class_list[[class]] = rownames(
      meta_data_df[(meta_data_df$Class == class),]
    )
  }
  
  return(class_list)
}


# Find samples of the same class to replace doppel samples
find_valid_non_doppel <- function(all_non_doppel,
                                  valid_doppel, 
                                  meta_data_df){
  valid_non_doppel = c()
  
  # Get classes of valid doppels
  valid_doppel_class = meta_data_df[valid_doppel, "Class"]
  
  # Group non doppels that we can select from by class
  non_doppel_class_list = create_class_list(
    sample_vec = all_non_doppel,
    meta_data_df = meta_data_df
  )
  
  class_counters = list()
  
  #Iterate through all valid doppel classes
  for (class in valid_doppel_class){
    if (!(class %in% names(class_counters))){
      class_counters[[class]] = 1
    }
    if (class_counters[[class]] > length(non_doppel_class_list[[class]])){
      print("Exceeded number of samples")
      next
    }
    # Add the first non doppel sample of this class
    valid_non_doppel = append(
      valid_non_doppel,
      non_doppel_class_list[[class]][class_counters[[class]]]
    )
    # Increment counter
    class_counters[[class]] = class_counters[[class]] + 1
  }
  
  return (valid_non_doppel)
}

# Find non doppel samples (that are not in valid_non_doppel)
# to balance the classes in the validation set
find_valid_non_doppel_extra <- function(all_non_doppel,
                                        valid_doppel,
                                        valid_non_doppel,
                                        meta_data_df){
  
  
  # Check if it is balanced
  if (is_balanced(valid_doppel, meta_data_df)){
    # If it is already balanced,
    # no need extra samples to balance classes in validation
    valid_non_doppel_extra = c()
  }
  else {
    # Not balanced hence we need to get non doppel samples for balancing
    
    # Remaining non doppel samples we can use for balancing
    remaining_non_doppel = setdiff(all_non_doppel, valid_non_doppel)
    
    # Separate the remain non doppel samples by class
    non_doppel_class_list = create_class_list(
      sample_vec = remaining_non_doppel,
      meta_data_df = meta_data_df
    )
    
    # 3e) Get the number of validation samples in the major class
    valid_classes = table(meta_data_df[valid_doppel, "Class"])
    major_class_size = max(valid_classes)
    
    # 3f) Get minor class in validation
    if (length(valid_classes) == 2){
      minor_class_size = min(valid_classes)
      minor_class = names(valid_classes[valid_classes == minor_class_size])
    } else{
      minor_class_size = 0
      major_class = names(valid_classes)
      minor_class = setdiff(unique(meta_data_df$Class), major_class)
    }

    # 3g) Add in extra samples to balance classes in validation
    num_extra_samples = major_class_size - minor_class_size
    valid_non_doppel_extra = non_doppel_class_list[[minor_class]][
      1:num_extra_samples
    ]
  } 
  
  return(valid_non_doppel_extra)
}

create_experiment_plan <- function(experiment_plan,
                                   train_doppel,
                                   train_non_doppel,
                                   valid_doppel,
                                   valid_non_doppel,
                                   valid_non_doppel_extra,
                                   num_samples_added){
  
  experiment_plan[["train_valid_sets"]] = list()
  
  max_num_valid_doppel = length(valid_doppel)
  current_num_valid_doppel = 0
  
  while (TRUE){
    # Ensure we do not exceed maximum number of doppelgangers in valid
    if (current_num_valid_doppel > max_num_valid_doppel){
      current_num_valid_doppel = max_num_valid_doppel
    }
    
    # Get training and validation samples samples
    train = union(train_doppel, train_non_doppel)
    valid = valid_non_doppel_extra
    
    # For Doppel_0 case
    if (current_num_valid_doppel == 0){
      # Don't add any val doppel to the validation set
      current_valid_doppel = c() 
      # Add all val non doppel to the validation set
      current_valid_non_doppel = valid_non_doppel 
    } 
    # For all other doppel cases
    else {
      current_num_valid_non_doppel = (
        max_num_valid_doppel - current_num_valid_doppel)
      
      # Add the first current_valid_doppel-th to validation
      current_valid_doppel = head(valid_doppel,
                                  current_num_valid_doppel)
      # Add last num_valid_non_doppel samples to validation
      current_valid_non_doppel = tail(valid_non_doppel,
                                      current_num_valid_non_doppel)
    }
    
    # Add the remaining of the val doppel to the training set
    current_train_doppel = setdiff(valid_doppel,
                                   current_valid_doppel)
    # Add remaining of the val non doppel to the training set
    current_train_non_doppel = setdiff(valid_non_doppel,
                                       current_valid_non_doppel)
    
    # Add non doppel and doppel to valid and train sets with union
    train = union(current_train_doppel, train)
    train = union(current_train_non_doppel, train)
    
    valid = union(current_valid_non_doppel, valid)
    valid = union(current_valid_doppel, valid)
    
    # Add train-valid sets to experiment plan list
    train_valid_name = paste("Doppel", 
                             current_num_valid_doppel, sep="_")
    experiment_plan[["train_valid_sets"]][[train_valid_name]] = list()
    experiment_plan[["train_valid_sets"]][[train_valid_name]][["train"]] = train
    experiment_plan[["train_valid_sets"]][[train_valid_name]][["valid"]] = valid
    
    # Stop once we reach max
    if (current_num_valid_doppel == max_num_valid_doppel){
      break;
    }
    
    # Increase number of validation doppel samples
    current_num_valid_doppel = (
      current_num_valid_doppel + num_samples_added)
  }
  
  # Add in positive control
  
  # Get training and validation samples samples
  train = union(train_doppel, train_non_doppel)
  valid = valid_non_doppel_extra
  
  # Resembles the final doppel case but instead of doppelgangers,
  # we add duplicates
  
  # Add valid_non_doppel to both train and valid (this is the duplicate)
  train = union(valid_non_doppel, train)
  valid = union(valid_non_doppel, valid)
  
  # Add train-valid sets to experiment plan list
  train_valid_name = paste("Pos_Con", 
                           current_num_valid_doppel, sep="_")
  experiment_plan[["train_valid_sets"]][[train_valid_name]] = list()
  experiment_plan[["train_valid_sets"]][[train_valid_name]][["train"]] = train
  experiment_plan[["train_valid_sets"]][[train_valid_name]][["valid"]] = valid
  
  return(experiment_plan)
}

# Adds into the experiment plan a list containing the sample names
# in each of the 5 subsets:
# 1. train_doppel: Doppelganger samples in training
# 2. train_non_doppel: Non doppelgager samples in training
# 3. valid_doppel: Doppelganger samples in validation (in nth Doppel)
# 4. valid_non_doppel: Non Doppelganger samples in validation (in 0 Doppel)
# 5. valid_non_doppel_extra: Samples added to validation to balance classes
add_sample_types_to_exp_plan <- function(experiment_plan,
                                         train_doppel,
                                         train_non_doppel,
                                         valid_doppel,
                                         valid_non_doppel,
                                         valid_non_doppel_extra){
  experiment_plan[["types_of_samples"]] = list()
  experiment_plan[["types_of_samples"]][["train_doppel"]] = train_doppel
  experiment_plan[["types_of_samples"]][["train_non_doppel"]] = train_non_doppel
  experiment_plan[["types_of_samples"]][["valid_doppel"]] = valid_doppel
  experiment_plan[["types_of_samples"]][["valid_non_doppel"]] = valid_non_doppel
  experiment_plan[["types_of_samples"]][["valid_non_doppel_extra"]] = valid_non_doppel_extra
  return(experiment_plan)
}

get_class_distribution <- function(vec, meta_data_df, class_order){
  classes = meta_data_df[vec, "Class"]
  class_table = table(classes)
  return(as.numeric(class_table[class_order]))
}

get_num_doppel_pairs <- function(doppel_mapping_df,
                                 train,
                                 valid){
  total_num_doppel_pairs = 0
  doppel_mapping = convert_df_to_list(doppel_mapping_df)
  for (sample in valid){
    doppel_partners = doppel_mapping[[sample]]
    num_doppel_pairs = length(intersect(doppel_partners, train))
    total_num_doppel_pairs = total_num_doppel_pairs + num_doppel_pairs
  }
  return(total_num_doppel_pairs)
}

get_num_val_doppel_samples <- function(doppel_mapping_df,
                                       train,
                                       valid){
  num_val_doppel_samples = 0
  doppel_mapping = convert_df_to_list(doppel_mapping_df)
  for (sample in valid){
    doppel_partners = doppel_mapping[[sample]]
    num_doppel_pairs = length(intersect(doppel_partners, train))
    if (num_doppel_pairs > 0){
      num_val_doppel_samples = num_val_doppel_samples + 1
    }
  }
  return(num_val_doppel_samples)
}

# Adds in training and validation stats 
add_train_valid_stats <- function(experiment_plan, 
                                  meta_data_df,
                                  doppel_mapping_df){
  
  train_valid_sets = names(experiment_plan[["train_valid_sets"]])
  
  # Get unique classes
  classes = unique(meta_data_df[,"Class"])
  classes = classes[order(classes)]
  
  # Instantiate stats df
  stats_df = data.frame(matrix(ncol = 6, nrow = 0))
  stats_colnames = c()
  for (class in classes){
    stats_colnames = append(stats_colnames, paste("train", class, sep="_"))
    stats_colnames = append(stats_colnames, paste("valid", class, sep="_"))
  }
  stats_colnames = stats_colnames[order(stats_colnames)]
  stats_colnames = append(stats_colnames, "num_doppel_pairs")
  stats_colnames = append(stats_colnames, "num_val_doppel_samples")
  colnames(stats_df) = stats_colnames
  
  for (train_valid in train_valid_sets){
    return_row = c()
    
    train_dist = get_class_distribution(
      experiment_plan[["train_valid_sets"]][[train_valid]][["train"]],
      meta_data_df,
      classes)
    
    valid_dist = get_class_distribution(
      experiment_plan[["train_valid_sets"]][[train_valid]][["valid"]],
      meta_data_df,
      classes)
    
    num_doppel_pairs = get_num_doppel_pairs(
      doppel_mapping_df,
      experiment_plan[["train_valid_sets"]][[train_valid]][["train"]],
      experiment_plan[["train_valid_sets"]][[train_valid]][["valid"]]
    )
    
    num_val_doppel_sample = get_num_val_doppel_samples(
      doppel_mapping_df,
      experiment_plan[["train_valid_sets"]][[train_valid]][["train"]],
      experiment_plan[["train_valid_sets"]][[train_valid]][["valid"]]
    )
    
    return_row = c(return_row, train_dist, valid_dist, num_doppel_pairs, num_val_doppel_sample)
    stats_df[train_valid, ] = return_row
    
  }
  
  experiment_plan[["train_valid_stats"]] = stats_df
  
  return(experiment_plan)
}



# Generates an experiment plan with the planning list and val_doppel given by the user
generate_experiment_plan <- function(planning_list,
                                     valid_doppel,
                                     num_samples_added,
                                     shuffle_seed=2022){
  experiment_plan = list()
  
  # 1. Load dataframes
  print("1. Loading data frames")
  meta_data_df = planning_list$meta_data
  doppel_mapping_df = planning_list$doppel_mapping
  
  # 2. Classify samples into doppel and non doppel
  print("2. Generating sample subsets")
  # 2a) Get all samples
  all_samples = unique(row.names(meta_data_df))
  # 2b) Get all doppel samples
  all_doppel = unique(row.names(doppel_mapping_df))
  # 2c) Get all non doppel samples
  all_non_doppel = setdiff(all_samples, all_doppel)
  
  print("3. Splitting data set into 5 key groups")
  # 3. Split the doppel and non doppel samples into the following groups:
  #   i. train_non_doppel = Samples that remain in the train set for 
  #      all train-valid sets
  #   ii. train_doppel = Samples that have to remain in the train set 
  #       for samples in val_doppel to be considered doppelganger 
  #       samples
  #   iii. valid_doppel = Samples that are incrementally added into the
  #        validation set to increase number of doppelgangers
  #   iv. valid_non_doppel = Samples that are incrementally replaced
  #       in the validation set to increase number of doppelgangers
  #   v. valid_non_doppel_extra = Non doppel samples in the validation
  #      used to balance the classes in validation
  
  print("- Generating train_doppel")
  # 3a) Get training doppel samples
  train_doppel = setdiff(all_doppel, valid_doppel)
  
  print("- Generating valid_doppel")
  # 3b) Shuffle the validation doppel samples
  set.seed(shuffle_seed)
  valid_doppel = sample(valid_doppel)
  
  print("- Generating valid_non_doppel")
  # 3c) Get valid non doppel
  valid_non_doppel = find_valid_non_doppel(
    all_non_doppel = all_non_doppel,
    valid_doppel = valid_doppel, 
    meta_data_df = meta_data_df
  )
  
  print("- Generating valid_non_doppel_extra")
  # 3d) Get valid non doppel extra
  valid_non_doppel_extra = find_valid_non_doppel_extra(
    all_non_doppel = all_non_doppel,
    valid_doppel = valid_doppel,
    valid_non_doppel = valid_non_doppel,
    meta_data_df = meta_data_df
  )
  
  print("- Generating train_non_doppel")
  # 3e) Get non doppelganger samples in training
  train_non_doppel = setdiff(
    all_non_doppel, union(valid_non_doppel, valid_non_doppel_extra))
  
  print("4. Creating the experiment plan list")
  # 4. Make the experiment plan
  experiment_plan = create_experiment_plan(
    experiment_plan = experiment_plan,
    train_doppel = train_doppel,
    train_non_doppel = train_non_doppel,
    valid_doppel = valid_doppel,
    valid_non_doppel = valid_non_doppel,
    valid_non_doppel_extra = valid_non_doppel_extra,
    num_samples_added = num_samples_added
  )
  
  print("5. Adding sample type information")
  # 5. Add in all 5 types of samples
  experiment_plan = add_sample_types_to_exp_plan(
    experiment_plan = experiment_plan,
    train_doppel = train_doppel,
    train_non_doppel = train_non_doppel,
    valid_doppel = valid_doppel,
    valid_non_doppel = valid_non_doppel,
    valid_non_doppel_extra = valid_non_doppel_extra
  )
  
  experiment_plan = add_train_valid_stats(
    experiment_plan = experiment_plan,
    doppel_mapping_df = doppel_mapping_df,
    meta_data_df = meta_data_df
  )
  
  print("Experiment plan generated!")
  return(experiment_plan)
}
