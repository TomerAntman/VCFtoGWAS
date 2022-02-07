\name{Build_offspring_GWAS_matrix}
\alias{Build_offspring_GWAS_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build a GWAS matrix for the offspring
%%  ~~function to do ... ~~
}
\description{
This function serves for step 1 of chapter 1.1 of the "Using VCFtoGWAS package" markdown series.\cr
Based on the GWAS matrix of the parents (assuming the matrix was built based on diploid parents), a matrix is built for the offspring of the parents.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Build_offspring_GWAS_matrix(parent_GWAS_mat,
                            offspring_strains,
                            dir_results = getwd(),
                            results_name = name_by_time(),
                            do_save = TRUE,
                            seed = 10,
                            save_hashed_matrix = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parent_GWAS_mat}{A GWAS matrix (composed of 0,1,2,NA) where the rownames are the parent strains}
  \item{offspring_strains}{An array of the names of the offspring, must be composed of the parent names}

   \item{dir_results}{The directory in which a folder will be created and results will be saved. Make sure it exists!!!}

  \item{results_name}{The name of the folder in which the results will be saved within dir_results (default is a time stamp, see \code{\link{create_directory}} function)}

  \item{do_save}{Do you wish to save the results? (will be saved as RDS files) (Default is TRUE)}

  \item{seed}{The function includes random sampling. A seed can be set and changed.}

  \item{save_hashed_matrix}{Do you wish to also save the the "mid-product" of the function (the hashed matrix). Default is \verb{FALSE}}
%%     ~~Describe \code{x} here~~
}

\details{
Hashing the table:\cr
If parent1 is "BNI" and parent2 is "BKE" than the offspring name is "BNIBKE".\cr
A hashed value for "BNIBKE" (denoted BNIBKE^[i,j]) is given by:\cr
BNIBKE^[i,j] = (BNI[i,j] + BKE[i,j]) + (BNI[i,j] * BKE[i,j])\cr

Un-hashing the table:\cr
\itemize{
\item (0 + 0) + (0 * 0) = 0 ==> BNIBKE[i,j] = 0
\item (0 + 1) + (0 * 1) = 1 ==> BNIBKE[i,j] = 0 or 1
\item (0 + 2) + (0 * 2) = 2 ==> BNIBKE[i,j] = 1
\item (1 + 1) + (1 * 1) = 3 ==> BNIBKE[i,j] = 0 or 1 or 2
\item (1 + 2) + (1 * 2) = 5 ==> BNIBKE[i,j] = 1 or 2
\item (2 + 2) + (2 * 2) = 8 ==> BNIBKE[i,j] = 2
\item NA in any of them results in NA
}

}
\value{
The returned object is a list of length 3 or 1:
\item{offspring_GWAS_mat}{The created GWAS matrix of the offspring (always returned)}
\item{mapping_info_offspring}{The mapping info linking between the SNPs and their location (after the creation of the new matrix, some SNPs are removed)}
\item{results_directory}{If the results are saved, this is the directory where they were saved}
}

\author{
Tomer Antman
}
