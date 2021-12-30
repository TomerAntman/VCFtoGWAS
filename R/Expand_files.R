Expand_files <- function(gt_GTonly_filt,
                         fix_filt,
                         dir_results = getwd(),
                         results_name = name_by_time(),
                         do_save = TRUE){
  if(do_save){
    step_name = "Step1.3-Expand"
    results_directory <- create_directory(called_from = step_name,
                                          dir_results = dir_results,
                                          results_name = results_name)
  }
  rownames(fix_filt)<-NULL
  message("Expanding the data")
  #* Create the column ALT_options containing the amount of alterations per chromosomal location
  fix_filt$ALT_options <- lengths(strsplit(fix_filt$ALT, ",") )

  fix_filt_expand <-fix_filt[rep(seq.int(1,nrow(fix_filt)), fix_filt$ALT_options), ]
  fix_filt_expand$ALT <- unlist(strsplit(fix_filt$ALT,","))

  message("Number of variants: ",dim(fix_filt_expand)[1] ,"\nFixed information: \n")
  psych::headTail(fix_filt_expand)
  #* the result is that we have all the possible SNP's,
  #* even 2 in the same spot are now considered different
  indication <- as.numeric(rownames(fix_filt_expand))%%1 *10 +1

  #* Each cell might contain:
  #*  - NA
  #*  - 0 0
  #*  - 0 1 (or 1 0)
  #*  - 1 1
  #*  - 1 2
  #*  - 0 2 (or 2 0)
  #*  - 2 2
  #*  - ... even up to 6!
  #* What to do with info:
  #* I need to replicate `gt_GTonly_filt` to be the same size of `fix_filt_expand`.
  #* this can be done by doing the same operation I did on fix_filt.


  gt_GTonly_filt_expand <- tryCatch({
    #* the seq.int just creates a sequence and each number in the sequence is repeated based on the
    #* value in fix_filt$ALT_options.
    #* This causes a specific replication of the rows

    gt_GTonly_filt[rep(seq.int(1,nrow(gt_GTonly_filt)), fix_filt$ALT_options), ]

  }, error = function(err){
    if(grepl('vector memory exhausted',err, fixed = TRUE)){
    message("There was an issue with the available memmory...\n
             This was the message: ",err,"\nTrying to fix issue")
    #* I hope that this will work and solve the space issue but it is not guaranteed
      Sys.setenv('R_MAX_VSIZE'=320000000000)
      gt_GTonly_filt[rep(seq.int(1,nrow(gt_GTonly_filt)), fix_filt$ALT_options), ]

    }else if(grepl('cannot allocate vector size',err, fixed = TRUE)){
      message("There was an issue with memmory allocation...\n
             This was the message: ",err,"\nTrying to fix issue")
      memory.limit(9999999999)
      gt_GTonly_filt[rep(seq.int(1,nrow(gt_GTonly_filt)), fix_filt$ALT_options), ]
   }
 })

  if(do_save){
    message("Saving expanded information to: ", results_directory)
    Save_as_RDS(list(fix_filt_expand = fix_filt_expand,
                     indication = indication,
                     gt_GTonly_filt_expand = gt_GTonly_filt_expand),
                directory = results_directory)

    results_list <- list(fix_filt_expand = fix_filt_expand,
                         indication = indication,
                         gt_GTonly_filt_expand = gt_GTonly_filt_expand,
                         results_directory = results_directory)
   return(results_list)
  }
  results_list <- list(fix_filt_expand = fix_filt_expand,
                       indication = indication,
                       gt_GTonly_filt_expand = gt_GTonly_filt_expand)

  return(results_list)

}


