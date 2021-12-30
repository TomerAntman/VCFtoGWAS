GWAS_conclusions <- function(GWAS,
                             mapping_info,
                             GWAS_file_name,
                             gData,
                             dir_to_save,
                             save_csv = TRUE,
                             drop_rand = TRUE){

  if(!require("ggrepel")) install.packages("ggrepel")
  library(ggrepel)
  if(!require("ggplot2")) install.packages("ggplot2")
  library(ggplot2)
  if(!require("ggplotify")) install.packages("ggplotify")
  library(ggplotify)
  if(!require("gridExtra")) install.packages("gridExtra")
  library(gridExtra, include.only = 'grid.arrange')
  if(!require("qqman")) install.packages("qqman")
  library(qqman)
  if(!require("dplyr")) install.packages("dplyr")
  library(dplyr)


  create_csv <- function(GWAS_trial_snps, checked_trait, name){
    #* This "temporary" function is meant to avoid writing this section twice.
    snp_full_data <- Create_complete_snp_data(GWAS_trial_snps,
                                              checked_trait,
                                              mapping_info,
                                              genes_df = VCFtoGWAS::SacCer_sgd_cds)

    # eval(parse(text = paste0(GWAS_trial_name,"_",checked_trait," = snp_full_data")))

    message("\nSaving csv for ",name,"\n")
    write.csv(snp_full_data,
              file = paste0(dir_to_save,"/",name,".csv"),
              row.names = FALSE)
  }



  pdf(file=paste0(dir_to_save,"/",GWAS_file_name,"-Plots_combined.pdf"), onefile = TRUE,
      width = 14, height = 10)

  #specify to save plots in 2x2 grid
  # par(mfrow = c(2,2))

  for (GWAS_trial_name in names(GWAS$GWAResult)){
    GWAS_trial_snps <- eval(parse(text = paste0("GWAS$signSnp$`",GWAS_trial_name,"`")))
    GWAS_trial <- eval(parse(text = paste0("GWAS$GWAResult$`",GWAS_trial_name,"`")))
    for (checked_trait in unique(GWAS_trial$trait)){
      message("\nTrial: ",GWAS_trial_name,", Trait: ",checked_trait)
      name <- paste0("SNPs-",GWAS_file_name,
                     "_Trial-",GWAS_trial_name,
                     "_Trait-",checked_trait)
      #* collect significant snps for this trait in this trial
      sigsnp <- GWAS$signSnp[[GWAS_trial_name]]
      if (checked_trait %in% sigsnp$trait){
        sigsnp <-sigsnp[sigsnp$trait==checked_trait,"snp"]
        sigsnp<-unlist(sigsnp)
        names(sigsnp)<-NULL
      }else{ sigsnp <-character(0)}

      # message("\nsigsnp:",sigsnp,"\n")
      # message("stop to debug")

      #* Create CSV of significant SNPs #####
      if (save_csv & !drop_rand & length(sigsnp)!=0){
        #* if dropping the randomized data isn't required just create the csv
        create_csv(GWAS_trial_snps, checked_trait, name)

      }else if(save_csv & drop_rand & length(sigsnp)!=0 & !grepl("Rand",checked_trait)){
        #* if dropping the randomized data *is* required then create the csv only if the trait name doesn't include "Rand"
        create_csv(GWAS_trial_snps, checked_trait,name)
      }

      ##########
      plot_names <- paste0("Trial: ",GWAS_trial_name,". Trait: ",checked_trait)
      #* data wrangle for plots #####
      gwas_data <- GWAS$GWAResult[[GWAS_trial_name]]
      gwas_data <- gwas_data[gwas_data$trait == checked_trait,c("snp","chr","pos","pValue","effect")]

      #* get y threshold of the p_value (whatever was set in the GWAS)
      ythr <- GWAS$thr[[GWAS_trial_name]][[checked_trait]]

      #* change column names to match the plot code
      colnames(gwas_data)<- c("SNP","CHR","BP","P","EFF")


      Headline_plot <-ggplot() +
        annotate("text", x = 1, y = 1, size = 10,
                 label = plot_names) +
        theme_void()

      plots <- list(Headline_plot = Headline_plot)

      #* 1) draw qq plot for all the results #####
      message("\nCreating QQ\n")
      # qqman::qq(gwas_data$P, main = paste0("QQ: ",plot_names), cex.main = 0.7)
      qq_plot <- as.ggplot(function()  qqman::qq(gwas_data$P)) +
        ggtitle("QQ plot")
        # ggtitle(paste0("QQ: ",plot_names))

      plots <- base::append(plots, list(qq_plot = qq_plot))

      #* 2) draw manhattan plot for all the results #####
      message("\nCreating manhattan\n")
      manhat_plot <-plot_manhattan(gwas_data, sigsnp, ythr)
      manhat_plot <- manhat_plot + ggtitle("Manhattan plot")

      plots <- base::append(plots, list(manhat_plot = manhat_plot))
        # ggtitle(paste0("Manhattan: ",plot_names))
      # suppressWarnings(print(manhat_plot))
      # print(manhat_plot)

      if (length(sigsnp)!=0 & !grepl("Rand",checked_trait)){

        #* 3) Polygenic score ######
        message("\nCreating Polygenic score\n")
        pol_effect <- polygenic_effect_score(GWAS, gData,
                                             GWAS_trial_name, checked_trait, sigsnp,
                                             bin =   grepl("bin",checked_trait))
        effect_details <- pol_effect$effect_scores
        effect_plot <- pol_effect$effect_plot +
          ggtitle("Polygenic Effect Score")
        # ggtitle(paste0("Polygenic Effect Score: ",plot_names))

        # print(effect_plot)

        plots <- base::append(plots, list(effect_plot = effect_plot))

        #* 4) histograms ######
        message("\nCreating allele VS fitness\n")
        alleleVSfitness <- allele_presence_plot(gData,
                                                GWAS_trial_name, checked_trait, sigsnp,
                                                bin =   grepl("bin",checked_trait))

        alleleVSfitness <- alleleVSfitness +
          ggtitle("allele VS fitness")
        # ggtitle(paste0("allele VS fitness: ",plot_names))

        # print(alleleVSfitness)

        plots <- base::append(plots, list(alleleVSfitness = alleleVSfitness))

      }

      ### PLOT EVERYTHING! ####
      # plots <- list(Headline_plot,qq_plot, manhat_plot, effect_plot, alleleVSfitness)
      # names(plots) <- c("Headline_plot","qq_plot", "manhat_plot", "effect_plot", "alleleVSfitness")
      for ( i in 1:length(plots)){
        message("\nplotting: ",names(plots)[i])
        grid.arrange(grobs = plots[i], top="", newpage = T)

      }
        ## PLOT HEAD TITLE PLOT
      # grid.arrange(grobs = list(qq_plot), top="", newpage = T)
      # grid.arrange(grobs = list(manhat_plot), top="", newpage = T)
      # grid.arrange(grobs = list(effect_plot), top="", newpage = T)
      # grid.arrange(grobs = list(alleleVSfitness), top="", newpage = T)
      #
      # grid.arrange(grobs = list(qq_plot, manhat_plot),
      #                               nrow = 2,
      #                               top=plot_names,
      #                               newpage = T)
      # grid.arrange(grobs = list(effect_plot, alleleVSfitness),
      #              nrow = 2,
      #              top=plot_names,
      #              newpage = T)
    } # for (checked_trait...
  } # for (GWAS_trial_name...



  #turn off PDF plotting
  dev.off()



  }
