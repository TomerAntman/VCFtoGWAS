\name{Filter_genotypes}
\alias{Filter_genotypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter the genotypes (strains)
}
\description{
This function serves for step 2 of chapter 1 of the "Using VCFtoGWAS package" markdown series.\cr
Based on an array of specific strains, the genotype matrix is filtered to only contain those strains.\cr
Afterwards, Rows that now contain irrelevant SNPs (that add no information) are omitted.
}
\usage{
Filter_genotypes(genotype_matrix,
                 fixed_data,
                 strains = "",
                 dir_results = getwd(),
                 results_name = name_by_time(),
                 do_save = TRUE,
                 on_columns = TRUE,
                 filter_zeros = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype_matrix}{
Gentoype matrix where rows are SNPs (variants) and columns are strains (genotypes) as was created in step 1 by \code{\link{Upload_vcf_to_R}}
}
  \item{fixed_data}{
%%     ~~Describe \code{fixed_data} here~~
}
  \item{strains}{
%%     ~~Describe \code{strains} here~~
}
  \item{dir_results}{
he directory in which a folder will be created and results will be saved. Make sure it exists!!!
}
  \item{results_name}{
The name of the folder in which the results will be saved within dir_results (default is a time stamp, see \code{\link{create_directory}})
}
  \item{do_save}{
Do you wish to save the results? (will be saved as RDS files) (Default is \code{TRUE})
}
  \item{on_columns}{
all actions assume that the strain names appear in the column names. If for some reason the entered matrix is transposed (strains in row names), than this parameter should be given as \code{FALSE} (default is \code{TRUE}).
}
  \item{filter_zeros}{TRUE by default. The function filters variants that are all NA. If this parameter is TRUE, it also filters variant rows for which the variant doesn't appear in any of the genotypes.
%%     ~~Describe \code{filter_by} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
