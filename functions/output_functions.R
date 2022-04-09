if (!("openxlsx" %in% installed.packages())){
  instaall.packages("openxlsx")
}
library(openxlsx)
options("openxlsx.maxWidth" = 250) 

# Output a list of dataframes to a single excel file
output_list_of_df_to_excel <- function(list_to_save, filename){
  # Make excel file
  wb = createWorkbook(filename)
  i = 1
  # For each dataframe
  for (df_name in names(list_to_save)){
    # add sheet
    addWorksheet(wb, df_name)
    # add data
    writeData(wb, 
              sheet = df_name, 
              list_to_save[[df_name]], rowNames = TRUE)
    setColWidths(wb, 
                 df_name, 
                 cols = 1:ncol(list_to_save[[df_name]]), 
                 widths = "auto")
    i = i+1
  }
  # Save the excel file
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

# This function generates a highlighted planning xlsx file
# - Highlighted samples are the doppelganger samples we are using in 
#   the validation set
output_highlighted_planning_file <- function(highlighted_filename,
                                             planning_list,
                                             doppel_samples_in_valid){
  # Save planning list to excel
  output_list_of_df_to_excel(planning_list, highlighted_filename)
  
  # Load workbook
  wb = loadWorkbook(highlighted_filename)
  
  ## Style indicating doppelganger
  doppel_style = createStyle(fgFill = "#d0adff")
  
  # 1. Colour doppelgangers meta_data sheet
  
  # a) Load meta_data df
  meta_data_df = read.xlsx(wb, 
                           sheet = "meta_data",
                           colNames = TRUE,
                           rowNames = FALSE)
  
  # b) Loop through rownames in sheet and colour doppelgangers
  for (i in 1:nrow(meta_data_df)){
    if (meta_data_df[i, 1] %in% doppel_samples_in_valid){
      addStyle(wb, 
               sheet = "meta_data", 
               doppel_style, rows = i+1, cols = 1)
    }
  }
  # 2. Colour doppelgangers in doppel_mapping sheet
  
  # a) Load doppel mapping df
  doppel_mapping_df = read.xlsx(wb, 
                                sheet = "doppel_mapping",
                                colNames = TRUE,
                                rowNames = FALSE)
  
  # b) Loop through all cells in the sheet and colour doppelgangers
  for (i in 1:nrow(doppel_mapping_df)){
    for (j in 1:ncol(doppel_mapping_df)){
      value = doppel_mapping_df[i, j]
      if (is.na(value)){
        next
      }
      if (value %in% doppel_samples_in_valid){
        addStyle(wb, 
                 sheet = "doppel_mapping", 
                 doppel_style, rows = i+1, cols = j)
      }
    }
  }
  
  # Save the excel file
  saveWorkbook(wb, file = highlighted_filename, overwrite = TRUE)
}

# Highlights sheets in a experiment plan workbook according to highlight style list
highlight_sheets_in_ex_plan_wb <- function(wb,
                                           sample_types,
                                           highlight_style_list){
  
  for (sheet_name in names(wb)){
    # Load this sheet
    sheet_df = read.xlsx(wb, 
                         sheet = sheet_name,
                         colNames = TRUE,
                         rowNames = FALSE)
    
    # Go through cells in this sheet
    for (i in 1:nrow(sheet_df)){
      for (j in 1:ncol(sheet_df)){
        value = sheet_df[i, j]
        if (is.na(value)){
          next
        }
        for (sample_type in names(sample_types)){
          if (value %in% sample_types[[sample_type]]){
            style_of_sample = highlight_style_list[[sample_type]]
            addStyle(wb, 
                     sheet = sheet_name, 
                     style_of_sample, 
                     rows = i+1, 
                     cols = j)
            break; #since we already found the true colour of this cell
          }
        }
      }
    }
  }
  
  # Add in colour legend sheet
  colour_df = data.frame(sample_type = names(highlight_style_list),
                         colour = rep(NA, length(names(highlight_style_list))))
  addWorksheet(wb, "colour_legend")
  writeData(wb, 
            sheet = "colour_legend", 
            colour_df, 
            colNames = TRUE)
  setColWidths(wb, 
               sheet = "colour_legend", 
               cols = 1:ncol(colour_df), 
               widths = "auto")
  j = 1
  for (sample_type in names(highlight_style_list)){
    style_of_sample = highlight_style_list[[sample_type]]
    addStyle(wb, 
             sheet = "colour_legend", 
             style_of_sample, 
             rows = j+1, 
             cols = 2)
    j = j + 1
  }
  
  return(wb)
}

# Output an experiment plan (list of list)
output_experiment_plan <- function(experiment_plan,
                                   planning_list,
                                   ex_plan_csv_fp,
                                   ex_plan_xlsx_fp,
                                   highlight_style_list){
  
  # Get experiment plan in data frame
  experiment_plan_df = convert_experiment_plan_to_df(experiment_plan)
  
  # Get types of samples list
  sample_types = experiment_plan[["types_of_samples"]]
  stats_df = experiment_plan[["train_valid_stats"]]
  
  # 1. Output plan to csv 
  print("1. Outputing plan to csv...")
  write.table(experiment_plan_df, 
              ex_plan_csv_fp,
              na = "",
              row.names = FALSE,
              col.names = TRUE,
              append = FALSE,
              sep = ",")  
  
  # 2. Output plan to xlsx
  print("2. Outputing plan to xlsx")
  # a) Generating planning xlsx
  print("  a) Generating  planning xlsx")
  output_list_of_df_to_excel(
    list_to_save = planning_list,
    filename = ex_plan_xlsx_fp
  )
  
  # b) Load workbook
  print("  b) Loading workbook")
  wb = loadWorkbook(ex_plan_xlsx_fp)
  
  # c) Add experiment plan to a new sheet in the wb
  print("  c) Add in experiment plan")
  addWorksheet(wb, "experiment_plan")
  writeData(wb, 
            sheet = "experiment_plan", 
            experiment_plan_df, 
            colNames = TRUE)
  setColWidths(wb, 
               sheet = "experiment_plan", 
               cols = 1:ncol(experiment_plan_df), 
               widths = "auto")
  
  # c) Add stats to a new sheet in the wb
  print("  d) Add in stats")
  addWorksheet(wb, "train_valid_stats")
  writeData(wb, 
            sheet = "train_valid_stats", 
            stats_df, 
            rowNames = TRUE,
            colNames = TRUE)
  setColWidths(wb, 
               sheet = "train_valid_stats", 
               cols = 1:ncol(stats_df), 
               widths = "auto")
  
  # d) Highlight all sheets in the wb
  print("  e) Highlighting all sheets in the workbook")
  wb = highlight_sheets_in_ex_plan_wb(
    wb,
    highlight_style_list = highlight_style_list,
    sample_types = sample_types
  )
  
  # e) Save the excel file
  print("  f) Saving the excel file")
  saveWorkbook(wb, file = ex_plan_xlsx_fp, overwrite = TRUE)
}



