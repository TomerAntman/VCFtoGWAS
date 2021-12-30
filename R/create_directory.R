
create_directory <- function(called_from,
                             dir_results = getwd(),
                             results_name = name_by_time()){
    #* The results_name parameter dictates the folder (within dir_results folder) that will be created
    #* to save the results.
    #* If the results name isn't changed from the default (name_by_time()),
    #* then a name based on time will be constructed.
    #* Nevertheless, the "step" of the process is always added to the name.
    #* Examples:
    #* "151221_10.58_Step1.1-Upload_VCF-results"
    #* "OnlyHomozyg_Step1.1-Upload_VCF-results"

    results_name<- paste0(results_name,"_", called_from) #deparse(substitute(called_from))
    results_directory <- paste0(dir_results,"/",results_name)
    dir.create(results_directory)

    #* check if directory was created and if not try something different:
    warn <- names(warnings())
    if (!dir.exists(results_directory) &&
        !is.null(warn) &&
        grepl('No such file or directory',warn)){
      message("\n",results_name, " folder couldn't be created because\n",dir_results,
              " directory doesn't exist.\nTrying to create ", results_name, " folder in\n",
              getwd(),"\n")
      results_directory2 <- paste0(getwd(),"/",results_name)
      dir.create(results_directory2)
      warn <- names(warnings())
      if (!dir.exists(results_directory2) &&
          !is.null(warn) &&
          grepl('No such file or directory',warn)){
        stop("\nCouldn't create a directory to save the results")
      }else{
        results_directory <- results_directory2
        message("\nResults will be saved into: \n",results_directory,"\n")
      }

    #* If the directory already exists, no new directory will be created
    }else if (!is.null(warn) && grepl('already exists',warn)){
        message("\nThe results directory already existed, the results will be saved in:\n",
                results_directory,"\n")

    #* If is wrong... just print a simple message
    }else{ message("\nResults directory created:\n",
                   results_directory,"\n")}

    #* return the directory
    return(results_directory)
}
