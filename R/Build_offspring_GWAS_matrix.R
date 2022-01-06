Build_offspring_GWAS_matrix <- function(parent_GWAS_mat,
                                        offspring_strains,
                                        mapping_info,
                                        dir_results = getwd(),
                                        do_save = TRUE,
                                        seed = 10,
                                        filter_zeros = TRUE,
                                        save_hashed_matrix = FALSE){
  # if(do_save){
  #   step_name = "Step1.5-Offspring_GWAS_matrix"
  #   #* get directory where results will be saved
  #   results_directory <- create_directory(called_from = step_name,
  #                                         dir_results = dir_results,
  #                                         results_name = results_name)
  # }
  results_directory <- dir_results

  set.seed(seed)
  names_length <- nchar(rownames(parent_GWAS_mat)[1])
  message("The length of a parent-strain name is: ", names_length)

  #* If parent1 is "BNI" and parent2 is "BKE" than the offspring name is "BNIBKE".

  message("\nWe construct the final matrix based on the name of the first parent (values in `left_vals`)\n",
          "and the name of the second parent (values in `right_vals`).\n")

  combos <- expand.grid(rownames(parent_GWAS_mat),rownames(parent_GWAS_mat))

  offspring_strains_options <- paste0(combos[,1],combos[,2])


  left_vals <- parent_GWAS_mat[substring(offspring_strains_options, 1, names_length), ]

  right_vals <- parent_GWAS_mat[substring(offspring_strains_options, names_length+1 ,names_length*2), ]

  message(Sys.time(), " - Left (first parent) & Right (second parent) matrices were created\n")

  rm(parent_GWAS_mat)

  message("Now it's time to combine the two matrices.\n",
          "Each element in left_vals (and right_vals) contains either 0,1,2 or NA.\n",
          "By combining the values, a form of 'hashed' matrix can be formed.\n")

  #* Simple mathematical operations are done element-wise in R.
  #* The 1st value in the 1st row in `left_vals` is added/multiplied by the 1st element in
  #* the 1st row of `right_vals` and so forth.

  offspring_GWAS_mat <- (left_vals + right_vals) + (left_vals * right_vals)

  rm (left_vals, right_vals)

  rownames(offspring_GWAS_mat) <- offspring_strains_options

  offspring_GWAS_mat <- offspring_GWAS_mat[rownames(offspring_GWAS_mat) %in% offspring_strains,]

  message(Sys.time()," - A hashed matrix was created.\n")

  if(save_hashed_matrix){
      message("Saving hashed matrix to: ", results_directory)
      Save_as_RDS(list(hashed_offspring_GWAS_mat = offspring_GWAS_mat),
                  directory = results_directory)

  }

  ##* Hashing:
  #* (0 + 0) + (0 * 0) = 0   => 0
  #* (0 + 1) + (0 * 1) = 1   => 0 or 1
  #* (0 + 2) + (0 * 2) = 2   => 1
  #* (1 + 1) + (1 * 1) = 3   => 0 or 1 or 2
  #* (1 + 2) + (1 * 2) = 5   => 1 or 2
  #* (2 + 2) + (2 * 2) = 8   => 2
  #* NA in any of them results in NA

  #* omit rows and cols with all NA
  message(Sys.time()," - ommiting rows/cols with all NA.\n")
  offspring_GWAS_mat <-
    offspring_GWAS_mat[rowSums(is.na(offspring_GWAS_mat)) != ncol(offspring_GWAS_mat),
                       colSums(is.na(offspring_GWAS_mat)) != nrow(offspring_GWAS_mat)]

  if(filter_zeros){
    message(Sys.time()," - removing rows with all zeros and updating the mapping info.\n")
    relevant_cols <- apply(offspring_GWAS_mat, 2, function(x) !(all(x==0 |is.na(x))))
    offspring_GWAS_mat <- offspring_GWAS_mat[, relevant_cols]
  }
  mapping_info <- mapping_info [rownames(mapping_info) %in% colnames(offspring_GWAS_mat), ]

  message(Sys.time()," - Un-hashing the matrix...\n")

  #* 1st index (the value is 0):
  offspring_GWAS_mat[which(offspring_GWAS_mat==0)] <- 0 # (0 + 0) + (0 * 0) = 0 -> *0*

  #* 2nd index (the value is 1):
  offspring_GWAS_mat[which(offspring_GWAS_mat==1)] <- sample(c(0,1),replace = T,length(which(offspring_GWAS_mat==1))) # (0 + 1) + (0 * 1) = 1 -> *0 or 1*

  #* 3rd index (the value is 2):
  offspring_GWAS_mat[which(offspring_GWAS_mat==2)] <- 1 # (0 + 2) + (0 * 2) = 2 => 1

  #* 4th index (the value is 3):
  offspring_GWAS_mat[which(offspring_GWAS_mat==3)] <- sample(c(0,1,2),replace = T, length(which(offspring_GWAS_mat==3)), prob = c(0.25,0.5,0.25)) # (1 + 1) + (1 * 1) = 3 -> *0 or 1 or 2*

  #* 6th index (the value is 5):
  offspring_GWAS_mat[which(offspring_GWAS_mat==5)] <- sample(c(1,2),replace = T,length(which(offspring_GWAS_mat==5))) # (1 + 2) + (1 * 2) = 5 -> *1 or 2*

  #* 9th index (the value is 8):
  offspring_GWAS_mat[which(offspring_GWAS_mat==8)] <- 2 # (2 + 2) + (2 * 2) = 8 -> *2*


  if(do_save){
    message("Saving offspring GWAS matrix to: ", results_directory)
    Save_as_RDS(list(offspring_GWAS_mat = offspring_GWAS_mat,
                     mapping_info_offspring = mapping_info),
                     directory = results_directory)

    results_list <- list(offspring_GWAS_mat = offspring_GWAS_mat,
                         mapping_info_offspring = mapping_info,
                         results_directory = results_directory)
  } else{ results_list <- list(offspring_GWAS_matrix = offspring_GWAS_mat) }
  return(results_list)



  # message(Sys.time()," - Un-hashing the matrix...\nThis might take a while...\n")# Progress bar added:")
  #
  #   if(!require("progress")){
  #     install.packages("progress", quiet = T)
  #   }
  #   library(progress, quietly = T)
  #
  #   n_iter <- dim(offspring_GWAS_mat)[1] * dim(offspring_GWAS_mat)[2]
  #
  #   pb <- progress_bar$new(format = "[:bar] :percent",
  #                          total = n_iter,
  #                          complete = "=",   # Completion bar character
  #                          incomplete = "-", # Incomplete bar character
  #                          current = ">",    # Current bar character
  #                          clear = FALSE,    # If TRUE, clears the bar when finish
  #                          width = 100)      # Width of the progress bar


  # for (strain in 1:nrow(offspring_GWAS_mat)) {
  #
  #   for (snp in 1:ncol(offspring_GWAS_mat)){
  #     #* if it's NA just skip it (it will stay NA)
  #     if(is.na(offspring_GWAS_mat[strain, snp])){next}
  #
  #     offspring_GWAS_mat[strain, snp] <- switch(
  #       offspring_GWAS_mat[strain, snp] + 1,
  #
  #       #* 1st index (the value is 0):
  #       0, # (0 + 0) + (0 * 0) = 0 -> *0*
  #
  #       #* 2nd index (the value is 1):
  #       sample(c(0,1),1), # (0 + 1) + (0 * 1) = 1 -> *0 or 1*
  #
  #       #* 3rd index (the value is 2):
  #       1, # (0 + 2) + (0 * 2) = 2 => 1
  #
  #       #* 4th index (the value is 3):
  #       sample(c(0,1,2), 1, prob = c(0.25,0.5,0.25)), # (1 + 1) + (1 * 1) = 3 -> *0 or 1 or 2*
  #
  #       #* 5th index (the value is 4):
  #       NULL, # there is no 4
  #
  #       #* 6th index (the value is 5):
  #       sample(c(1,2),1), # (1 + 2) + (1 * 2) = 5 -> *1 or 2*
  #
  #       #* 7th index (the value is 6):
  #       NULL, # there is no 6
  #
  #       #* 8th index (the value is 7):
  #       NULL, # there is no 7
  #
  #       #* 9th index (the value is 8):
  #       2 # (2 + 2) + (2 * 2) = 8 -> *2*
  #     )
  #     # pb$tick()
  #   }
  # }
} #end function



