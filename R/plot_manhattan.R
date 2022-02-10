plot_manhattan <- function(GWAS,
                           GWAS_trial_name,
                           checked_trait){

  if(!require("ggrepel")) install.packages("ggrepel")
  library(ggrepel)
  if(!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)
  if(!require("qqman")) install.packages("qqman")
  library(qqman)
  if(!require("dplyr")) install.packages("dplyr")
  library(dplyr)

  #* data wrangle for plots #####
  gwas_data <- GWAS$GWAResult[[GWAS_trial_name]]
  gwas_data <- gwas_data[gwas_data$trait == checked_trait,c("snp","chr","pos","pValue","effect")]

  #* get y threshold of the p_value (whatever was set in the GWAS)
  if(GWAS$GWASInfo$thrType == 'fdr'){
    eval(parse(text = paste0(
      "ythr <- min(-log10(GWAS$signSnp$",GWAS_trial_name,"$pValue[GWAS$signSnp$",GWAS_trial_name,"$trait=='",checked_trait,"' & GWAS$signSnp$",GWAS_trial_name,"$snpStatus=='significant SNP']))"
    )))

  }else{ythr <- GWAS$thr[[GWAS_trial_name]][[checked_trait]]}

  #* change column names to match the plot code
  colnames(gwas_data)<- c("SNP","CHR","BP","P","EFF")

  eval(parse(text = paste0(
    "sigsnp <- GWAS$signSnp$",GWAS_trial_name,"$snp[GWAS$signSnp$",GWAS_trial_name,"$trait=='",checked_trait,"' & GWAS$signSnp$",GWAS_trial_name,"$snpStatus=='significant SNP']"
  )))
  eval(parse(text = paste0(
    "interesting_snp <- GWAS$signSnp$",GWAS_trial_name,"$snp[GWAS$signSnp$",GWAS_trial_name,"$trait=='",checked_trait,"']"
  )))


  # Prepare the dataset
  don <- gwas_data %>%
    filter(!is.na(P)) %>%

    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(gwas_data, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%

    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% interesting_snp , "yes", "no")) %>%
    mutate( is_annotate=ifelse(SNP %in% sigsnp, "yes", "no")) %>%
    mutate( is_pos_effect = ifelse(EFF>0, "yes","no"))

  # Prepare X axis
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  # Make the plot
  manhat_plot<- ggplot(don, aes(x=BPcum, y=-log10(P)))

  # make a dummy df to help with the color separation
  if (length(sigsnp!=0)){
    don_dummy <- don[c(1,2),]
    don_dummy$is_highlight <- "yes"
    don_dummy$`effect sign` <- relevel(as.factor(c("+","-")),"+")
    manhat_plot<-manhat_plot +
      geom_point(data=don_dummy, aes(fill = `effect sign`), color = c("green", "orange"),size=0)
  }


  # Show all points
  manhat_plot<-manhat_plot +
    geom_point(data=subset(don[seq(1, nrow(don),5),], is_highlight!="yes"), aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "black"), 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0.5) )      # remove space between plot area and x axis

  # Add highlighted points
  if (length(sigsnp!=0)){
    manhat_plot<-manhat_plot +
      # geom_point(data=subset(don, is_highlight=="yes" & is_pos_effect=="yes"), size = 1.3, color="green") +
      # geom_point(data=subset(don, is_highlight=="yes" & is_pos_effect=="no"), size = 1.3, color="orange") +
      # geom_point(data=subset(don, is_highlight=="yes"), aes(fill = is_pos_effect, size = 0)) +
      geom_point(data=subset(don, is_annotate=="yes" & is_pos_effect=="yes"), aes(size = abs(EFF)), color="green") +
      geom_point(data=subset(don, is_annotate=="yes" & is_pos_effect=="no"), aes(size = abs(EFF)), color="orange") +

      # scale_size_continuous(name = "|effect size|", range = (range(abs(don$EFF),na.rm = T)))+
      scale_size_continuous(name = "|effect size|") +

      scale_fill_manual(values = c("green", "orange"),drop = FALSE)

  }
  # med_value = as.vector(quantile(abs(don$EFF)[don$is_annotate=="yes"], probs = 0.7, na.rm=T))
  # don$is_annotate[abs(don$EFF) < 0.7*max(abs(don$EFF)[don$is_annotate=="yes"])] <- "no" # only mark the upper quantile of the effect sizes
  don$is_annotate[abs(don$EFF) < as.vector(quantile(abs(don$EFF)[don$is_annotate=="yes"], probs = 0.7, na.rm=T))] <- "no"
  manhat_plot<-manhat_plot +

    # Pvalue line
    geom_hline(yintercept=ythr, linetype="dashed", color = "red") +

    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=3, box.padding = 0.5, max.overlaps = Inf) +

    # Custom the theme:
    theme_bw() +
    theme(
      # legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  if (length(sigsnp!=0)){
    manhat_plot<-manhat_plot +
      guides(color = "none",
             fill = guide_legend(override.aes = list(colour = c("green","orange"),size = 3)),
             size = guide_legend(override.aes = list(colour = "black")))
  }else{
    manhat_plot<-manhat_plot +
      guides(color = "none")
  }
  manhat_plot<-manhat_plot +
    labs (title  = "Manhattan Plot",
          subtitle = paste0("Trial: ", GWAS_trial_name,"; Trait: ",checked_trait)) +
    ylab(bquote(-log[10](Pvalue))) +
    xlab("Chromosome")
  # ggtitle(paste0("Manhattan: ",plot_names))
  return(manhat_plot)

}
