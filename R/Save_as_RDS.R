Save_as_RDS <- function(variables_list,
                        directory = getwd()){

  for (i in 1:length(variables_list)){
    var = variables_list[[i]]
    saveRDS(var, file = paste0(directory,"/",names(variables_list)[i],".RDS"))
    message(Sys.time(), " - ",names(variables_list)[i]," saved to: ", directory)
  }




}
