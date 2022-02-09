polygenic_effect_score<-function(gData,
                                 GWAS,
                                 GWAS_trial_name,
                                 checked_trait,
                                 haplo_sufficency = TRUE,
                                 do_log2 = TRUE){

  if(!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)
  try(expr = {if(!require("ggpmisc")) {install.packages("ggpmisc")}})

  bin =   grepl("bin",checked_trait)
  #* collect the effects of the significant SNPs

  # effect <- GWAS$signSnp[[GWAS_trial_name]]
  # rows_of_sig <- effect$trait==checked_trait
  # effect <- effect[rows_of_sig,"effect"][[1]]

  eval(parse(text = paste0(
    "effect <- GWAS$signSnp$",GWAS_trial_name,"$effect[GWAS$signSnp$",GWAS_trial_name,"$trait=='",checked_trait,"']"
    )))
  eval(parse(text = paste0(
    "sigsnp <- GWAS$signSnp$",GWAS_trial_name,"$snp[GWAS$signSnp$",GWAS_trial_name,"$trait=='",checked_trait,"']"
  )))

  #* Subset the markers (the original GWAS matrix)
  markers <- gData$markers[,colnames(gData$markers) %in% sigsnp]

  #* If we assume haplo sufficiency:
  if(haplo_sufficency){
    markers[which(markers==1)] <- 0
    markers[which(markers==2)] <- 1
  }


  #* calculate the polygenic score
  pol_score <- base::rowSums(as.data.frame(sweep(as.data.frame(markers), MARGIN=2, effect, `*`)))

  #* collect phenotypes
  phenos <- gData$pheno[[GWAS_trial_name]]
  strains <- phenos["genotype"][[1]]
  rownames(phenos) <- strains

  df<- data.frame(trait = phenos[strains, checked_trait], poly_score = as.vector(pol_score[strains]))
  rownames(df)<-strains
  x_label = "SNP effect score"
  if( do_log2){
    df$poly_score = log2((-1*df$poly_score)+1)
    x_label = "Log2[-(SNP effect score)+1]"
  }
  if (bin){
    pol_plot <- ggplot(df, aes(x = poly_score,
                               y = as.factor(trait),
                               color = as.factor(trait))) +
      geom_violin() +
      geom_jitter(aes(color = as.factor(trait)),width=0.15, alpha=0.5) +
      # scale_x_continuous(breaks = round(seq(min(df$poly_score), max(df$poly_score), by = 0.5),1)) +
      theme_bw() +
      labs(title = "Polygenic Effect Score",
           # subtitle = "Categorical fitness on Glycerol",
           subtitle = paste0("Trial: ", GWAS_trial_name,"; Trait: ",checked_trait),
           x = x_label,
           y = checked_trait,
           color = checked_trait)
  }else{
    df <- df[!is.na(df$trait),]
    # effect_score.lm <- lm(trait ~ poly_score, df)

    if (!requireNamespace("ggpmisc")){
      pol_plot <- ggplot(df, aes(x = poly_score, y = trait)) +
        geom_point() +
        geom_smooth(method = "lm") +
        theme_bw() +
        labs(title = "Polygenic Effect Score",
             # subtitle = "Fitness on Glycerol",
             subtitle = paste0("Trial: ", GWAS_trial_name,"; Trait: ",checked_trait),
             x = x_label,
             y = checked_trait)
    }else{
      pol_plot <- ggplot(df, aes(x = poly_score, y = trait)) +
        geom_point() +
        geom_smooth(method = "lm") +
        ggpmisc::stat_poly_eq(formula = y ~ x,
                     aes(label = ..rr.label..),
                     parse = TRUE) +
        theme_bw() +
        labs(title = "Polygenic Effect Score",
             # subtitle = "Fitness on Glycerol",
             subtitle = paste0("Trial: ", GWAS_trial_name,"; Trait: ",checked_trait),
             x = x_label,
             y = checked_trait)
    }



  }
  return(list(effect_scores = df, effect_plot = pol_plot))

}

