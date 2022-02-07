\name{Create_complete_snp_data}
\alias{Create_complete_snp_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create Complete SNP data
%%  ~~function to do ... ~~
}
\description{
Given a set of results for a GWAS analysis (trial+trait) and the relevant mapping information of the SNPs in the data - the function creates a dataframe with information about the SNPs the were identified as significant (can be saved as csv.
}
\usage{
Create_complete_snp_data(GWAS_trial_snps, trait, mapping_info, genes_df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GWAS_trial_snps}{SNPs from a certain GWAS trial (for example`GWAS_trial_snps <- GWAS$signSnp[[GWAS_trial_name]]`)
%%     ~~Describe \code{GWAS_trial_snps} here~~
}
  \item{trait}{The desired trait (as was analysed in the GWAS)
%%     ~~Describe \code{trait} here~~
}
  \item{mapping_info}{The fixed information about the SNPs, linking between the unique names and the SNP's location and type.
%%     ~~Describe \code{mapping_info} here~~
}
  \item{genes_df}{
Dataframe that contains info about the genes (that the significant SNPs might be in). For yeast, there is an included file in the package data (SacCer_sgd_cds)
}
}
\details{
The returned information includes attributes such as location in the genome and if the SNP happens to be located in a gene coding sequence, the dataframe also includes information about the gene and the alternated amino acid sequence
}
\value{\item{SNP_information}{The information as metioned in the description and in the details}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Data source for CDS \href{http://rsat.sb-roscoff.fr/data/genomes/}{file extraction}

To create the final dataframe used as input, I used \file{create_SacCer_sgd_cds.R}


}
\author{
Tomer Antman
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
