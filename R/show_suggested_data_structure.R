show_suggested_data_structure<-function(){
  if(!require("data.tree")){install.packages("data.tree",quiet = TRUE)}
  d1 <- "Results"
  d2 <- paste0(d1, "/Chapter1-VCF2GWAS")
  d3 <- paste0(d2, "/Step1.1-Upload_VCF")
  d4 <- paste0(d3, "/Step1.2-Filter_genotypes")
  d5 <- paste0(d4, "/Step1.3-Expand")
  d6 <- paste0(d5, "/Step1.4-GWAS_Matrix")
  path <- c(
    d1,
    d2,
    d3,
    paste0(d3, list("/fix_sub.RDS","/gt_GTonly.RDS")),
    d4,
    paste0(d4, list("/fix_filt.RDS","/gt_GTonly_filt.RDS")),
    d5,
    paste0(d5, list("/fix_filt_expand.RDS","/indication.RDS","/gt_GTonly_filt_expand.RDS")),
    d6,
    paste0(d6, list("/mapping_info.RDS","/GWAS_mat.RDS","/phenotypes.csv","/offspring_GWAS_mat.RDS","/offspring_phenotypes.csv")),
    "Results/Chapter2-Analysis",
    "Results/Chapter2-Analysis/Step2.1-Create_gData",
    "Results/Chapter2-Analysis/Step2.1-Create_gData/gData.RDS",
    "Results/Chapter2-Analysis/Step2.1-Create_gData/Step2.2-Single_Trait_GWAS",
    "Results/Chapter2-Analysis/Step2.1-Create_gData/Step2.2-Single_Trait_GWAS/GWAS_Result.RDS",
    "Results/Chapter2-Analysis/Step2.1-Create_gData/Step2.2-Single_Trait_GWAS/Significant_SNPs_information.csv",
    "Results/Chapter2-Analysis/Step2.1-Create_gData/Step2.2-Single_Trait_GWAS/complete_results_plots.pdf"
    # "Results/Chapter3-Conclusions",
    # "Results/Chapter3-Conclusions/Step3.1-Visualizations",
    # "Results/Chapter3-Conclusions/Step3.1-Visualizations/Results1.pdf"
    )
  data.tree::as.Node(data.frame(pathString = path))
}
