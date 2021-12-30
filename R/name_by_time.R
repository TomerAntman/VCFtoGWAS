name_by_time <- function(){
    results_name <- gsub(":",".",Sys.time())
    results_name <- gsub(" ","_",results_name)
    results_name <- substr(results_name,1,nchar(results_name)-3)
    year <- substr(strsplit(results_name,"-")[[1]][1],3,4)
    month <- strsplit(results_name,"-")[[1]][2]
    day <- strsplit(strsplit(results_name,"-")[[1]][3],"_")[[1]][1]
    time <- strsplit(strsplit(results_name,"-")[[1]][3],"_")[[1]][2]
    final_results_name <- paste0(day,month,year,"_",time)
    # return format: DDMMYY_hh.mm
    return(final_results_name)
}
