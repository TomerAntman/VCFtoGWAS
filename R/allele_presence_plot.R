allele_presence_plot<-function(gData,
                               GWAS,
                               GWAS_trial_name,
                               checked_trait){


  if(!require("data.table")) install.packages("data.table")
  library(data.table)
  if(!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)

  bin =   grepl("bin",checked_trait)

  eval(parse(text = paste0(
    "sigsnp <- GWAS$signSnp$",GWAS_trial_name,"$snp[GWAS$signSnp$",GWAS_trial_name,"$trait=='",checked_trait,"']"
  )))
  markers <- gData$markers[,colnames(gData$markers) %in% sigsnp]
  phenos <- gData$pheno[[GWAS_trial_name]]

  temp_df <- as.data.frame(matrix(nrow = nrow(markers), ncol = 1))
  colnames(temp_df) <- c(checked_trait)
  rownames(temp_df) <- phenos["genotype"][[1]]
  temp_df[,1] <- phenos[,checked_trait]

  # strains <- phenos["genotype"][[1]]
  # temp_df[strains,] <- phenos[, c("genotype",checked_trait)]


  df_hist <- merge(temp_df,markers, by=0, all=TRUE)
  colnames(df_hist)[1] <- "genotype"

  molten.data <- melt(as.data.table(df_hist), id = c("genotype",checked_trait))
  colnames(molten.data)[c(2,4)]<- c("trait","allele_presence")
  if(!require("dplyr")) install.packages("dplyr")
  library(dplyr)
  molten.data <- molten.data %>%
    mutate(allele = recode(allele_presence,
                           `0` = '0/0',
                           `1` = '1/0',
                           `2` = '1/1')) %>%
    mutate(allele = as.factor(allele))
  colnames(molten.data)[2]<- "trait"

  # create colors (because there is some kind of issue)
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (bin){
    molten.data$count <- NA
    for (var in unique(molten.data$variable)){
      for (trait in unique(molten.data$trait)){
        for (allele in c(0,1,2)){
          relevant_rows <- (molten.data$variable == var & molten.data$trait == trait & molten.data$allele_presence == allele)
          molten.data$count[relevant_rows] <- length(which(relevant_rows)) # each row has information (count column about the number of times it appears)
          }
        }
      }
    color_hues <- gg_color_hue(length(unique(molten.data$trait)))

    pp <- ggplot(molten.data, aes(x = as.factor(trait),
                            y = allele,
                            colour = as.factor(allele_presence),
                            size = count), alpha = 0.8) +
      geom_point() +
      geom_text(aes(label = count),
                colour = "black",
                size = 3) +
      theme_bw() +
      scale_size_continuous(range = c(5, 20)) +
      theme(
            panel.background = element_blank()) +
      facet_wrap(~variable) +
      labs(title= "SNP Allele Presence in correlation to fitness",
           subtitle = paste0("Trial: ", GWAS_trial_name,"; Trait: ",checked_trait)) +
      xlab(checked_trait) +
      ylab("SNP allele presence") +
      guides(color = 'none', size = 'none')
               # guide_legend(override.aes =
               #                     list(colour = color_hues)))

    # pp<- ggplot(molten.data, aes(x = as.factor(trait), y = count, fill = allele)) +
    #   # geom_bar()
    #   geom_bar(position = "dodge", stat = "identity")+
    #   geom_text(aes(label=count),vjust = -0.5, color="black",
    #             position = position_dodge(0.9), size=3)+
    #   scale_y_continuous(expand = c(0,0), limits = c(0, max(molten.data$count)+100) ) +
    #   # geom_bar(position = "fill",stat = "identity")+
    #   # scale_y_continuous(labels = scales::percent_format()) +
    #   # geom_density(alpha = 0.4) +
    #   theme_bw() +
    #   facet_wrap(~variable) +
    #   xlab("SNP allele presence") +
    #   ylab("Count") +
    #   labs(fill = "Minor allele presence") +
    #   guides(fill = guide_legend(override.aes =
    #                                list(colour = color_hues)))
  }else{
    val_tab <- table(molten.data$allele_presence[molten.data$variable==sigsnp[1]])
    label <- paste0(paste0("0/0: ", val_tab["0"]),"\n",
                    paste0(" 1/0: ", val_tab["1"]),"\n",
                    paste0(" 1/1: ", val_tab["2"]))
    for (snp in sigsnp[2:length(sigsnp)]){
      val_tab <- table(molten.data$allele_presence[molten.data$variable==snp])
      label <- c(label,
                 paste0(paste0("0/0: ", val_tab["0"]),"\n",
                        paste0(" 1/0: ", val_tab["1"]),"\n",
                        paste0(" 1/1: ", val_tab["2"])))

    }
    molten.data <- molten.data[!is.na(molten.data$trait),]
    dat_text <- data.frame( label = label, variable = sigsnp)

    color_hues <- gg_color_hue(length(unique(molten.data$allele)))

    pp<- molten.data %>%
      ggplot(aes(x=trait)) +
      # geom_histogram(aes(x=trait, fill=allele), alpha = 0.5)+
      geom_density(aes(x=trait, fill=allele), alpha=0.2) +
      theme_bw() +
      facet_wrap(~variable, scales = "free") +
      xlab(checked_trait) +
      labs(title= "SNP Allele Presence in correlation to fitness",
           subtitle = paste0("Trial: ", GWAS_trial_name,"; Trait: ",checked_trait),
           fill = "SNP allele presence") +
      guides(fill = guide_legend(override.aes =
                                   list(colour = color_hues))) +
      geom_text(data = dat_text, aes(label = label),
                x = -Inf, y = Inf, hjust = -0.2, vjust = 1.3)
  }
  return(pp)

}
