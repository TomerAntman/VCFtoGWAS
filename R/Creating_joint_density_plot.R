create_joint_denisty_plot <- function(){
# library("VCFtoGWAS")
library("statgenGWAS")
library(ggplot2)
# library(ggplotify)
# library(gridExtra, include.only = 'grid.arrange')
library(ggrepel)


# get original gData
gData <- readRDS("../Using VCFtoGWAS package/100_tests/gData.RDS")[[1]]
mapping_info <- readRDS("../Using VCFtoGWAS package/100_tests/mapping_info.RDS")

# Add another number to the discrete values based on bins (2 equal sized groups with fitness below zero)
gData$pheno$U$binfit[gData$pheno$U$binfit==-1] <- -2
Uphen<-gData$pheno$U
Uphen<-Uphen[Uphen$fit<0 & !is.na(Uphen$fit),]
Uphen<-split(Uphen, cut(Uphen$fit, 2))
lower_Uphen <- Uphen[[1]]
gData$pheno$U$binfit[gData$pheno$U$genotype %in% lower_Uphen$genotype] <- -1


gData$pheno$Y$binfit[gData$pheno$Y$binfit==-1] <- -2
Yphen<-gData$pheno$Y
Yphen<-Yphen[Yphen$fit<0 & !is.na(Yphen$fit),]
Yphen<-split(Yphen, cut(Yphen$fit, 2))
lower_Yphen <- Yphen[[1]]
gData$pheno$Y$binfit[gData$pheno$Y$genotype %in% lower_Yphen$genotype] <- -1


phen_curr <- gData$pheno$U
na_count <- length(which(is.na(phen_curr$fit)))
phen_curr$binfit[!is.na(phen_curr$fit)] <-
  as.numeric(cut_number(phen_curr$fit[!is.na(phen_curr$fit)],round((dim(phen_curr)[1]-na_count)/na_count)))
phen_curr$binfit[is.na(phen_curr$fit)] <- 0
gData$pheno$U<- phen_curr

phen_curr <- gData$pheno$Y
na_count <- length(which(is.na(phen_curr$fit)))
phen_curr$binfit[!is.na(phen_curr$fit)] <-
  as.numeric(cut_number(phen_curr$fit[!is.na(phen_curr$fit)],round((dim(phen_curr)[1]-na_count)/na_count)))
phen_curr$binfit[is.na(phen_curr$fit)] <- 0
gData$pheno$Y<- phen_curr

saving_dir <- "../Using VCFtoGWAS package/100_tests"


library(progress)

n_iter <- 900

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = n_iter,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar

results_list = list()
counter = 1
for (i in 1:10){
  pb$tick()
  #* Set the seed and the proportion of data taken for training
  set.seed(i)
  set.seed(170)
  test_relative_size = 0.7

  #* Create training and testing data
  gData_train <- gData
  gData_train$markers <- gData$markers[sample(rownames(gData$markers), round(test_relative_size*nrow(gData$markers))),]
  gData_train$pheno <- list(
    U = gData$pheno$U[gData$pheno$U$genotype %in% rownames(gData_train$markers), ],
    Y = gData$pheno$Y[gData$pheno$Y$genotype %in% rownames(gData_train$markers), ]
  )

  rownames(gData_train$pheno$U)<-NULL
  rownames(gData_train$pheno$Y) <- (1:nrow(gData_train$pheno$Y))+nrow(gData_train$pheno$U)


  gData_test <- gData
  gData_test$markers <- gData$markers[!(rownames(gData$markers) %in% rownames(gData_train$markers)),]
  gData_test$pheno <-  list(
    U = gData$pheno$U[gData$pheno$U$genotype %in% rownames(gData_test$markers), ],
    Y = gData$pheno$Y[gData$pheno$Y$genotype %in% rownames(gData_test$markers), ]
  )

  rownames(gData_test$pheno$U) <- NULL
  rownames(gData_test$pheno$Y) <-  (1:nrow(gData_test$pheno$Y))+nrow(gData_test$pheno$U)

  saveRDS(gData_test, file = paste0(saving_dir,"/bestU_gData_test.RDS"))
  saveRDS(gData_train, file = paste0(saving_dir,"/bestU_gData_train.RDS"))
  GWAS_train <- runSingleTraitGwas(gData = gData_train)
  saveRDS(GWAS_train, file = paste0(saving_dir,"/bestU_GWAS_train.RDS"))

  # message("\n",Sys.time())
  for (GWAS_trial_name in c("U","Y")){
    GWAS_trial_snps <- eval(parse(text = paste0("GWAS_train$signSnp$`",GWAS_trial_name,"`")))
    for (checked_trait in c('fit','binfit')){
      # message("\nSeed: ",i,"; counter: ",counter,"\nTrial: ",GWAS_trial_name,", Trait: ",checked_trait)

      sigsnp <- GWAS_train$signSnp[[GWAS_trial_name]]
      if (checked_trait %in% sigsnp$trait){
        sigsnp <-sigsnp[sigsnp$trait==checked_trait,"snp"]
        sigsnp<-unlist(sigsnp)
        names(sigsnp)<-NULL
      }else{ sigsnp <-character(0)}

      if (length(sigsnp)==0){next}

      pol_effect_Tr <- polygenic_effect_score(gData_train, GWAS_train,
                                              GWAS_trial_name, checked_trait)
      df_Tr <- pol_effect_Tr$effect_scores
      pol_effect_Ts <- polygenic_effect_score(gData_test, GWAS_train,
                                              GWAS_trial_name, checked_trait)
      df_Ts <- pol_effect_Ts$effect_scores

      # collect the R squared of the training on training and of training on testing (also the covariance)
      df <- data.frame(media = GWAS_trial_name,
                       trait = checked_trait,
                       rsq_TrTr = cor(df_Tr[,1],df_Tr[,2])^2,
                       rsq_TrTs = cor(df_Ts[,1],df_Ts[,2])^2,
                       cov_TrTr = cov(df_Tr[,1],df_Tr[,2]),
                       cov_TrTs = cov(df_Ts[,1],df_Ts[,2]),
                       seed = i,
                       num_SNPs = length(sigsnp),
                       snps = paste(sigsnp,collapse = ", "))


      results_list[[counter]] <- df
      counter = counter +1
    }
  }

if(i == 500){
  mid_results500 = do.call(rbind, results_list)
  saveRDS(mid_results500, file = paste0(saving_dir,"/mid_results500.RDS"))
  }
}

all_results10 = do.call(rbind, results_list)


saveRDS(all_results1000, file = paste0(saving_dir,"/all_results1000.RDS"))
# p<-ggplot(all_results, aes(x = rsq_Ts, y = rsq_Tr, size = snp_num, shape =  media,color = trait)) +
#   geom_point() +
#   theme_bw()+
#   theme(legend.position="bottom")
# p2<- ggMarginal(p, type="histogram")
# saveRDS(p2, file = paste0(saving_dir,"/plot_with_histogram1.RDS"))
}
