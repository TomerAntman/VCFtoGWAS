usethis::use_package("psych")
usethis::use_package("vcfR")
usethis::use_package("evaluate")
usethis::use_package("stringr")
usethis::use_package("readr")

Upload_vcf_to_R <- function(vcf_file,
                            dir_results = getwd(),
                            results_name = name_by_time(),
                            do_save = TRUE,
                            do_return = TRUE,
                            get_chr_info = TRUE,
                            fix_columns = c("CHROM","POS","REF","ALT","QUAL")){

  if (!do_save & !do_return){
    message("Nothing will be saved and nothing will be returned.\nThere is no point running the function")
    return()
  }
  if(do_save){
    step_name = "Step1.1-Upload_VCF"
    #* get directory where results will be saved
    results_directory <- create_directory(called_from = step_name,
                                          dir_results = dir_results,
                                          results_name = results_name)
  }


  if(!require("vcfR")){
    install.packages("vcfR", quiet = TRUE)
  }
  library(vcfR, quietly = TRUE)

  library(evaluate, quietly = TRUE)

  #* load the vcf file and expect trouble.
  #* that is why it is called in evaluate function
  x <- evaluate("vcf <- read.vcfR(vcf_file, verbose = FALSE)")

  #* if the vcf wasn't loaded:
  if (!exists("vcf")){
    #* It might be because of an issue with the 'Rcpp' libary (commun issue)
    #* then an error will rise and will be saved into `x[[2]]`.
    #* The try function will check if `x[[2]]` exists (it will exist if there was an error).
    #* If it doesn't exist then everything is fine meaning an error will rise for calling `x[[2]]`.
    #* If calling `x[[2]]` doesn't raise an error and it contains 'Rcpp',
    #* we try to solve the vcf loading issue and then "recall" the vcf upload
    call_vcf_again = FALSE
    try({
      if (grepl('Rcpp', unlist(x[[2]]), fixed = TRUE)){
        message(Sys.time()," - Trying to solve vcf loading issue by reinstalling 'Rcpp' package")
        ## solve issue
        install.packages('Rcpp', quiet = TRUE)
        library(Rcpp)
        ##
        call_vcf_again = TRUE
      }else(print(x))
    }, TRUE)

    #* if any action was taken, the call_vcf_again variable is changed to TRUE
    if (call_vcf_again){
      vcf <- read.vcfR(vcf_file, verbose = FALSE)
    }
  }
  if(exists("vcf")){message(Sys.time()," - The vcf file is loaded successfuly")}

  if (get_chr_info){
    message(Sys.time(), " - getting chromosome info")

    if(!require("stringr")){
      install.packages("stringr", quiet = TRUE)
    }

    if(!require("readr")){
      install.packages("readr", quiet = TRUE)
    }
    library(stringr, quietly = TRUE)

    #* find rows in meta data
    rows_with_chromosomes <- which(startsWith(vcf@meta, "##contig=<ID=chromosome"))
    #* extract rows from meta data. This is the original pattern: "##contig=<ID=chromosome1,length=230218>"
    chromosomes_and_length <- str_match(strwrap(vcf@meta[rows_with_chromosomes]), "<\\s*(.*?)\\s*>")[,2]
    #* based on the vcf pattern, split by comma and get numbers from each string
    list_of_chrom_and_len <- readr::parse_number(unlist(str_split(string = chromosomes_and_length, pattern = ",")))
    #* create data frame. the order is chromosome, length, chromosome, length, ...
    df_chrom_and_len <- data.frame(Chromosome = list_of_chrom_and_len[c(T,F)], Length = list_of_chrom_and_len[c(F,T)])
    #* add row of total length (overall)
    df_chrom_and_len[nrow(df_chrom_and_len) + 1,] = c("Total", sum(df_chrom_and_len$Length))
    print(df_chrom_and_len)
  } # end if get_chr_info

  #* Filter the fixed information (get rid of irrelevant columns):
  #* first call it from the vcf file
  fix_sub <- as.data.frame(vcf@fix)
  #* second change the columns:
  fix_sub <- fix_sub[,fix_columns]
  message(Sys.time(), " - fixed information (SNP info) was extracted and filtered")
  print("Filtered fixed information: ")

  if(!require("psych")){install.packages("psych", quiet = TRUE)}
  psych::headTail(fix_sub)


  ##* Extract genotypes from the vcf data:
  #* **GT**: genotype, encoded as allele values separated by either of / or |.
  #* The allele values are:
  #* - 0 for the reference allele (what is in the REF field)
  #* - 1 for the first allele listed in ALT
  #* - 2 for the second allele list in ALT and so on.
  #* For diploid calls examples could be 0/1, 1|0, or 1/2, etc.
  #* If a call cannot be made for a sample at a given locus, ‘.’ should be specified for each missing allele
  #* in the GT field (for example ‘./.’ for a diploid genotype and ‘.’ for haploid genotype).
  #* The meanings of the separators are as follows
  #* ◦ / : genotype unphased
  #* ◦ | : genotype phased
  #*
  #* Basically, the next line gives us a matrix with values: NA, "0/0", "./.", "0/1", "3/3", etc.
  #* The rownames are locations along the genome and the columns are genotypes (i.e strains)
  gt_GTonly <- extract.gt(vcf)

  message(Sys.time(), " - genotype info extracted successfuly\nThe output matrix is of size: ",dim(gt_GTonly),
          "\n and the first 10 columns look like this:\n")
  psych::headTail(gt_GTonly[,1:10])

  ###

  if(do_save){

    message("Saving files to: ", results_directory)
    Save_as_RDS(list(fix_sub = fix_sub,gt_GTonly = gt_GTonly),
                directory = results_directory)


    message("All files saved!\nTo read files, use the readRDS function")
    if (!do_return){return(results_directory)}
  }
  if(do_return){
    fix_and_gt <- list(fix_sub=fix_sub, gt_GTonly=gtGTonly)
    return(fix_and_gt)
  }

}#end function

