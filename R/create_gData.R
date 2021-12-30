create_gData <- function(GWAS_mat,
                         map,
                         phenotype,
                         trial=2,
                         features=3,
                         covariate = NULL,
                         do_save = TRUE,
                         dir_results = getwd(),
                         results_name = NA
                         ){
  if(do_save){
    step_name = "Step2.1-gData_creation"
    if (is.na(results_name)){
      results_name = name_by_time()
    }

    #* get directory where results will be saved
    results_directory <- create_directory(called_from = step_name,
                                          dir_results = dir_results,
                                          results_name = results_name)
  }

  if(!require("statgenGWAS")){
    install.packages("statgenGWAS")
  }
  library(statgenGWAS)

  #* Rename Chromosome and Position columns.
  colnames(map)[match(c("CHROM","POS"), colnames(map))] <- c("chr","pos")
  #* make the chromosome be only numbers and the position numeric
  map$chr<- as.numeric(stringr::str_remove(map$chr,"chromosome"))
  map$pos <- as.numeric(map$pos)

  #* split the phenotype dataframe into
  PhenoList <- split(x = phenotype[c("genotype", features)],
                          f = list(phenotype[,trial]),
                          drop = TRUE)


  ## Create a gData object containing map and marker information.
  gData <- createGData(geno = GWAS_mat,
                       map = map,
                       pheno = PhenoList,
                       covar = covariate)

  ##Debug:
  message("\n",Sys.time(), " - gData object created")

  gDataDedup <- tryCatch(expr = {
    message("\n\nTrying to use codeMarkers\n")
    codeMarkers(gData, impute = TRUE, verbose = TRUE)
  }, error = function(cond){
    if (grepl("cannot allocate vector", cond)){
      message(cond,"\nThe memory limit is: ~", round(memory.limit()/1000),"Gb.\nIncreasing memory limit...\n")
      s = memory.limit()
      memory.limit(s*10)
      message("\nmemory increased. Trying to use codeMarkers\n")
      codeMarkers(gData, impute = TRUE, verbose = TRUE)
    }
  }
  )

  ##Debug:
  message("\n",Sys.time(), " - impute and verbose done (cut the data)")

  if(do_save){
    message("\nSaving gData RDS file to: ", results_directory,"\n")
    Save_as_RDS(list(gData = gDataDedup,
                     mapping_info = mapping_info),
                directory = results_directory)

    results_list <- list(gData = gDataDedup,
                         directory = results_directory)
  }else{results_list <- list(gData = gDataDedup)}
  return(results_list)
}


