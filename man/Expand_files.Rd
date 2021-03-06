\name{Expand_files}
\alias{Expand_files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expand rows of files
}
\description{
This function serves for step 3 of chapter 1 of the "Using VCFtoGWAS package" markdown series.\cr
Expand the fixed info and the genotype (gt) matrix so that each row represents only one SNP or INDEL (up until now each row represented a position)
The output is the expanded fixed and gt and also an array (called `indication` that links between the two. Relevant for the next step)
}
\usage{
expand_files(gt_GTonly_filt,
             fix_filt,
             dir_results = getwd(),
             results_name = name_by_time(),
             do_save = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gt_GTonly_filt}{The filtered Genotype matrix (step2 output)}
  \item{fix_filt}{The filtered fixed information dataframe (step2 output)}
  \item{dir_results}{The directory in which a folder will be created and results will be saved. Make sure it exists!!!}
  \item{results_name}{The name of the folder in which the results will be saved within dir_results (default is a time stamp, see \code{\link{create_directory}})}
  \item{do_save}{Do you wish to save the results? (will be saved as RDS files) (Default is TRUE)}

}
\details{
For very large files it might take a while.\cr\cr
An explanation about the \code{indication} array that is returned:\cr
If for example, row number 3 in \code{fix_filt} represents some position that has two ALTs (like: [C, TG]) for a REF (like: [A]),  in \code{fix_filt_expand} it will be two rows (as 3.0, 3.1) with only one ALT in each.\cr
The \code{indication} array will get a "1" for row 3.0 and "2" for row 3.1. This array is the link between the data in the \code{fix_filt_expand} and the data in \code{gt_GTonly_filt_expand}.\cr
The \code{indication[i]} gives us the id number of the alternative sequence in a certain location.\cr
It is between 1 and some number (specifically here it's 6 but can be more).\cr
In simpler terms: after all the pre-processing, each row in \code{fix_filt_expand} represents only one alteration from the reference (only one ALT). Since the same position might have had several alterations (in \code{fix_filt}), the \code{indication} array indicates what alteration is represented in each row (in regards to the position in the genome).

}
\value{
  \item{fix_filt_expand}{A dataframe of the fixed information. Each row represents only one variant(one SNP or INDEL). It has an extra column: `ALT_options` which represents the number of alterations that exist per genomic position in the data}
  \item{gt_GTonly_filt_expand}{A matrix in which the number of times each row appears is affected by the `ALT_options` in \code{fix_filt_expand}}
  \item{indication}{This array is the link between the data in the \code{fix_filt_expand} and the data in \code{gt_GTonly_filt_expand}}

}

\author{
Tomer Antman
}
\note{
This should be run on the results of the \code{\link{Filter_genotypes}} function from step 2.
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
