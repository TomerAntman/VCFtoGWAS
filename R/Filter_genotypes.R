
Filter_genotypes<- function (genotype_matrix,
                             fixed_data,
                             strains = "",
                             dir_results = getwd(),
                             results_name =  name_by_time(),
                             do_save = TRUE,
                             on_columns = TRUE,
                             filter_by = c("./.")){
  if(do_save){
    step_name = "Step1.2-Filter_genotypes"
    results_directory <- create_directory(called_from = step_name,
                                          dir_results = dir_results,
                                          results_name = results_name)
  }
  if (strains == ""){
    filtered_gt <- genotype_matrix
    fix_filt <- fixed_data
    warning(Sys.time()," - No strains were entered. No filtration required.
Output files will be identical to input files.")
  }else{
    if(!require("psych")){install.packages("psych", quiet = TRUE)}
    #* all actions assume that the strain names appear in the column names,
    #* if for some reason the entered matrix is transposed (strains in row names),
    #* It is first transposed:
    if (!on_columns){genotype_matrix <- t(genotype_matrix)}
    #* and then we continue as usual...
    #*
    #* Filter the genotypes based on the entered strains:
    filtered_gt <- genotype_matrix[, colnames(genotype_matrix) %in% strains]
    if(dim(filtered_gt)[2]==0){
      stop("no strains in the genotype matrix were found in the given strains array")}
    ###
    message(Sys.time()," - Genotype matrix filtered by strains")
    ###

    ##* Omit irrelevant SNP rows:

    if(length(filter_by)==1){
      message(Sys.time()," - Finding SNPs to keep")
      relevant_rows <- apply(filtered_gt, 1, function(x) !(all(x==filter_by |
                                                                 is.na(x))))
    }else if(length(filter_by)==2){
      message(Sys.time()," - Finding SNPs to keep")
      #* if length is 2, I expect filter_by = c("./.","0/0")
      relevant_rows <- apply(filtered_gt, 1, function(x) !(all(x==filter_by[1] |
                                                                 x==filter_by[2]|
                                                                 is.na(x))))
    }else(stop("Problem with filter_by parameter. Length isn't 1 nor 2\n
             filter_by is adivsed to be either c('./.') or c('./.','0/0')"))

    ###
    message(Sys.time()," - Update: \nOut of ",dim(filtered_gt)[1]," initial SNP rows, \n",
            length(which(relevant_rows))," were identified as relevant.\n\nFiltering data...")
    ###

    filtered_gt <- filtered_gt[relevant_rows, ]
    message(Sys.time()," - Genotype matrix filtered to have only SNP rows that contain information")
    psych::headTail(filtered_gt[,1:15])

    fixed_data <- fixed_data[relevant_rows, ]
    message(Sys.time()," - Fixed info filtered as well")
    psych::headTail(fixed_data)

    ###
    #* re-transpose the gentoype matrix if it was transposed
    if (!on_columns){filtered_gt <- t(filtered_gt)}
    ###

    if(!require("psych")){install.packages("psych", quiet = TRUE)}
    psych::headTail(fixed_data)
  }

  if(do_save){
  message("Saving files to: ", results_directory)
  Save_as_RDS(list(fix_filt = fix_filt,gt_GTonly = gt_GTonly),
              directory = results_directory)

  return_list = list(directory = results_directory,
                     fix_filt = fixed_data,
                     gt_GTonly_filt = filtered_gt)
  }else{
    return_list = list(fix_filt = fixed_data,
                       gt_GTonly_filt = filtered_gt)}

  return(return_list)

} # end function

