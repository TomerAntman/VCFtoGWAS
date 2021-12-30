Create_complete_snp_data <- function(GWAS_trial_snps,
                                     trait,
                                     mapping_info,
                                     genes_df
                                     ){

  if (!require("Biostrings", quietly = TRUE)){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("Biostrings" , update = TRUE)
  }

  only_needed <- GWAS_trial_snps[GWAS_trial_snps$trait==trait,c("snp","chr","pos","effect","allFreq","pValue")]
  only_needed$rom_chr <- as.factor(as.character(as.roman(only_needed$chr)))
  only_needed <- cbind(only_needed, mapping_info[only_needed$snp, c("REF","ALT")])
  #* now we have the "interesting" info about the significant SNPs

  #* Now, as we look at the yeast cds data from sgd (genes_df), we first filter it to contain only
  #* chromosomes that appear in our data (shorter running time)
  genes_df_sub <- subset(genes_df, genes_df$chromosomes %in% levels(only_needed$rom_chr))
  genes_df_sub <- droplevels(genes_df_sub)
  genes_df_sub$start <- as.numeric(genes_df_sub$start)
  genes_df_sub$end <- as.numeric(genes_df_sub$end)

  #* create columns with NA
  {
    only_needed$gene <- NA
    only_needed$start <- NA
    only_needed$end <- NA
    only_needed$strand <- NA
    only_needed$description <- NA
    only_needed$ref_dna <- NA
    only_needed$alt_dna <- NA
    only_needed$ref_aa <- NA
    only_needed$alt_aa <- NA
  }

  #* Go to each snp and check in which gene it is located (if it's in a gene)
  #* and add info from `genes_df_sub` to `only_needed`

  for (snp in 1:nrow(only_needed)){
    for (gene in which(genes_df_sub$chromosomes == only_needed$rom_chr[snp])){
      if(findInterval(only_needed$pos[snp], c(genes_df_sub$start[gene],genes_df_sub$end[gene]+1))==1){
        only_needed$gene[snp] <- genes_df_sub$genes[gene]
        only_needed$start[snp] <- genes_df_sub$start[gene]
        only_needed$end[snp] <- genes_df_sub$end[gene]
        only_needed$strand[snp] <- as.character(genes_df_sub$strand[gene])
        only_needed$description[snp] <- genes_df_sub$description[gene]
        only_needed$ref_dna[snp] <- genes_df_sub$dna_seq[gene]
        only_needed$ref_aa[snp] <- genes_df_sub$aa_seq[gene]

        position_in_gene = only_needed$pos[snp] - genes_df_sub$start[gene] +1
        ref_length <- nchar(only_needed$REF[snp])
        new_dna_seq = genes_df_sub$dna_seq[gene]
        if (position_in_gene == 1){
          end_half <- substr(new_dna_seq, ref_length+1 , nchar())
          new_dna_seq <- paste0(only_needed$ALT[snp],end_half)
        }else{
          start_half <- substr(new_dna_seq, 1, position_in_gene-1)
          end_half <- substr(new_dna_seq, position_in_gene + ref_length,nchar(new_dna_seq) )
          new_dna_seq <- paste0(start_half,only_needed$ALT[snp],end_half)
        }
        only_needed$alt_dna[snp] <- new_dna_seq
        only_needed$alt_aa[snp] <- as.character(Biostrings::translate( Biostrings::DNAString(new_dna_seq)))


      } # if(findInterval ...
    } # for (gene in...
  } # for (snp in ...

return(only_needed)

} # function ...






