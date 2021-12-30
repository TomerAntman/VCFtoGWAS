Save_GWAS <- function(GWAS,
                      threshold_type,
                      nstrains,
                      rand = FALSE,
                      dir_results = getwd(),
                      results_name = NA){
  step_name = "Step2.2-GWAS"
  if (is.na(results_name)){
    results_name = name_by_time()
    }

  #* get directory where results will be saved
  results_directory <- create_directory(called_from = step_name,
                                          dir_results = dir_results,
                                          results_name = results_name)

  message("\nSaving GWAS RDS file to: ", results_directory,"\n")
  file_name <- paste0("GWAS_",threshold_type,"-",nstrains,"_strains")
  message(sprintf("File name to be saved:\n %s \n ",file_name))
  GWAS_as_list <- list(GWAS)
  names(GWAS_as_list) <- file_name
  Save_as_RDS(GWAS_as_list, directory = results_directory)

}
