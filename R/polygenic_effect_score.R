polygenic_effect_score<-function(GWAS,
                                 gData,
                                 GWAS_trial_name,
                                 checked_trait,
                                 sigsnp,
                                 bin){

  if(!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)

  #* collect the effects of the significant SNPs
  effect <- GWAS$signSnp[[GWAS_trial_name]]
  rows_of_sig <- effect$trait==checked_trait
  effect <- effect[rows_of_sig,"effect"][[1]]



  #* Subset the markers (the original GWAS matrix)
  markers <- gData$markers[,sigsnp]

  #* calculate the polygenic score
  pol_score <- base::rowSums(as.data.frame(sweep(as.data.frame(markers), MARGIN=2, effect, `*`)))

  #* collect phenotypes
  phenos <- gData$pheno[[GWAS_trial_name]]
  strains <- phenos["genotype"][[1]]
  rownames(phenos) <- strains

  df<- data.frame(trait = phenos[strains, checked_trait], poly_score = as.vector(pol_score[strains]))
  rownames(df)<-strains
  if (bin){
    pol_plot <- ggplot(df, aes(x = poly_score,
                               y = as.factor(trait),
                               color = as.factor(trait))) +
      geom_violin() +
      geom_jitter(aes(color = as.factor(trait)),width=0.15, alpha=0.5) +
      theme_bw() +
      labs(color = checked_trait)+
      xlab("SNP effect score") +
      ylab(checked_trait)
  }else{
    df <- df[!is.na(df$trait),]
    pol_plot <- ggplot(df, aes(x = poly_score, y = trait)) +
      geom_point() +
      theme_bw() +
      xlab("SNP effect score") +
      ylab(checked_trait)
  }
  return(list(effect_scores = df, effect_plot = pol_plot))

}

