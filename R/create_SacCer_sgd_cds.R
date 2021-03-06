create_SacCer_sgd_cds <- function(){
  if (!require("biomartr", quietly = TRUE)){
    install.packages("biomartr")
  }


  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  if (!require("Biostrings", quietly = TRUE)){
    BiocManager::install("Biostrings" , update = TRUE)
  }
  if (!require("biomaRt", quietly = TRUE)){
    BiocManager::install("biomaRt" , update = TRUE)
  }



  cds_file_location <- "./SacCer cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa"
  cds_file_biostring<-biomartr::read_cds(cds_file_location)
  cds_file_dt<-biomartr::read_cds(cds_file_location, obj.type = "data.table")
  cds_file_df <- as.data.frame(cds_file_dt)
  rownames(cds_file_df) <- unlist(lapply(cds_file_dt$geneids, function(x) unlist(strsplit(x,"_"))[1]))

  cds_names <- names(cds_file_biostring)
  cds_names <- cds_names[!grepl(":Mito:",cds_names)]
  split_by_colon <- lapply(cds_names, function(x) strsplit(x,split = ":")[[1]])
  # split_by_colon <- do.call(rbind, split_by_colon)


  genes = unlist(lapply(split_by_colon, function(x) gsub("_.*","",x[1])))
  dna_seq = cds_file_df[genes,2]
  aa_seq = unlist(lapply(dna_seq, function(x) {
    as.character( Biostrings::translate( Biostrings::DNAString(x) ))
  }))



  SacCer_sgd_cds <- data.frame(genes = genes,
                         chromosomes = as.factor(unlist(lapply(split_by_colon, function(x) x[3]))),
                         start = unlist(lapply(split_by_colon, function(x) x[4])),
                         end = unlist(lapply(split_by_colon, function(x) x[5])),
                         strand = as.factor(unlist(lapply(split_by_colon, function(x)
                         {if (substr(x[6],1,1) == "-") "-" else "+"}))),
                         description = unlist(lapply(cds_names, function(x) gsub(".*description:","",x))),
                         dna_seq = dna_seq,
                         aa_seq = aa_seq
  )

  save(SacCer_sgd_cds, file = "./VCFtoGWAS/data/SacCer_sgd_cds.rda")
}
