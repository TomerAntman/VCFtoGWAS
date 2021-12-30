usethis::use_package("stringr")

Get_GWAS_matrix <- function(gt_GTonly_filt_expand,
                            fix_filt_expand,
                            indication,
                            only_SNPs = FALSE,
                            dir_results = getwd(),
                            results_name = name_by_time(),
                            do_save = TRUE){
  if(do_save){
    if(only_SNPs){step_name = "Step1.4-GWAS_Matrix-only_SNPs"
    }else{step_name = "Step1.4-GWAS_Matrix"}
    results_directory <- create_directory(called_from = step_name,
                                          dir_results = dir_results,
                                          results_name = results_name)
  }
  if (only_SNPs){
    message("Deleting all INDELs from the data (only SNPs will remain")
    SNP_rows <- which(nchar(fix_filt$ALT)==1 & nchar(fix_filt$REF)==1)
    fix_filt <- fix_filt[SNP_rows, ]
    gt_GTonly_filt <- gtGTonly[SNP_rows, ]
    indication <- indication[SNP_rows]
  }

  #* keep unique names
  rownames(gt_GTonly_expand) <- make.names(rownames(gt_GTonly_filt_expand), unique = TRUE)
  if(!require("stringr")){install.packages("stringr")}
  #* Create empty matrix of the same size with NAs

  indication_as_character <- as.character(round(indication))
  result_array <- stringr::str_count(gt_GTonly_filt_expand, indication_as_character)
  GWAS_mat <- matrix(result_array,
                     nrow = nrow(gt_GTonly_filt_expand),
                     ncol = ncol(gt_GTonly_filt_expand),
                     dimnames = list(rownames(gt_GTonly_filt_expand),
                                     colnames(gt_GTonly_filt_expand)))


  ### MIGHT BE UNNECESSARY: ##########
  #* creating the GWAS matrix: go cell by cell and check the list it contains.
  #* The cell in `GWAS_mat` get's a value by the number of instances of the `indication`
  #* in the matching cell of `gt_GTonly_expand`
  #*

  # if(!require("progress")){
  #   install.packages("progress")
  # }
  # library(progress)
  # #* progress creates a progress bar which allows to show a visual progress
  #
  # pb <- progress_bar$new(format = "(:spin) [:bar] :percent",
  #                        total = nrow(GWAS_mat),
  #                        complete = "=",   # Completion bar character
  #                        incomplete = "-", # Incomplete bar character
  #                        current = ">",    # Current bar character
  #                        clear = FALSE,    # If TRUE, clears the bar when finish
  #                        width = 100)      # Width of the progress bar
  #
  # for (snp in 1:nrow(GWAS_mat)){# rows
  #   for (strain in 1:ncol(GWAS_mat)){## columns
  #
  #     if(is.na(gt_GTonly_expand[snp,strain])){next} #skip NA (saves a lot of time)
  #     array(stringr::str_count(demo1[1:5,],as.character(round(indication))),dim=c(5,10))
  #     cell_value <-
  #       round(
  #         as.numeric(
  #           unlist(
  #
  #             strsplit(
  #               gt_GTonly_expand[snp,strain], "\\||/"
  #               #* \1 split each cell: "0/0" becomes a list:
  #               #*     [[1]]
  #               #*     [1] "0" "0"
  #               )#close strsplit
  #             #* \2 unlist the list of the cell to get:
  #             #* [1] "0" "0"
  #           )#close unlist
  #           #* \3 turn character list to numeric to get numeric array:
  #           #* [1] 0 0
  #         )#close as.numeric
  #         #* \4 round the array otherwise it doesn't work:
  #         #* [1] 0 0
  #       )#close round
  #
  #     #* cell_value contains the value of the cell ((3 3), (0 1), etc)
  #     #* so the info in cell_value came from the gt matrix
  #     #* we now compare the values in the cell_value to the `indication` values.
  #     #* The final value that will be introduced to the GWAS matrix is the number of appearances
  #     #* of the indication value of that row (`indication[snp]`) in the cell_value
  #
  #     GWAS_mat[snp, strain] <- sum(round(indication[snp]) == cell_value)
  #
  #     #* check the number of occurrences (instances) of the current indication in the current cell value.
  #     #* For homozygotes of the current indication we'll get 2.
  #     #* For hetero that contains it at least once, we'll get 1, and if it doesn't appear well get 0.
  #     #* If NA we'll get NA
  #
  #   } ##
  #   #* progress_bar that updates after strain:
  #   pb$tick()
  # } #

  ###############

  ##* Finish gwas matrix
  #* `rs` stands for 'reference SNP'
  rows <- sprintf("%s%0*d", "rs",1, 1:nrow(fix_filt_expand))
  rownames(fix_filt_expand) <- rows
  rownames(GWAS_mat) <- rows
  GWAS_mat <- t(GWAS_mat)

  if(!require("psych")){install.packages("psych", quiet = TRUE)}
  psych::headTail(GWAS_mat[,1:15])

  message(Sys.time(), " - Finished gwas matrix")

  if(do_save){
    message("Saving files to: ", results_directory)
    Save_as_RDS(list(GWAS_mat = GWAS_mat,
                     mapping_info = fix_filt),
                directory = results_directory)

    return_list = list(directory = results_directory,
                       mapping_info = fixed_data,
                       GWAS_mat = GWAS_mat)
  }else{return_list = list(mapping_info = fixed_data,
                           GWAS_mat = GWAS_mat)}
  return(return_list)
}











