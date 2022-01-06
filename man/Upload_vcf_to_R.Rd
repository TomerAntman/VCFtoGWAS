\name{Upload_vcf_to_R}
\alias{Upload_vcf_to_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Turn vcf files to usable R files
%%  ~~function to do ... ~~
}
\description{
This function serves for step 1 of chapter 1 of the "Using VCFtoGWAS package" markdown series.\cr
This script loads the vcf file that you want to work with  (either processed by the "usegalaxy" server or not). The file can be uploaded in its zipped version (.gz file extension). The vcf processing is based on the [vcfR package](https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html).
}
\usage{
Upload_vcf_to_R(vcf_file,
                dir_results = getwd(),
                results_name = name_by_time(),
                do_save = TRUE,
                do_return = TRUE,
                get_chr_info = TRUE,
                fix_columns = c("CHROM","POS","REF","ALT","QUAL"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf_file}{The location where the vcf file is saved}
  \item{dir_results}{The directory in which a subfolder will be created and results will be saved. Make sure it exists!!!}
  \item{results_name}{The name of the folder in which the results will be saved within dir_results (default is a time stamp, see \code{\link{create_directory}})}
  \item{do_save}{Do you wish to save the results? (will be saved as RDS files) (Default is TRUE)}
  \item{do_return}{Do you wish to return the results to your current workspace environment? (Default is TRUE)}
  \item{get_chr_info}{Information regarding the lengths of the chromosomes (Default is TRUE)}
  \item{fix_columns}{The column names that you wish to get from the fixed column of the vcf (the default is probably all you need so don't change it)}

}
\details{
These files are usually very large and it will take a while.\cr
The files exported are saved as .RDS files. They are lighter and very easy to read in R by calling readRDS(file = filepath).\cr\cr
Extract genotypes from the vcf data:\cr
GT: genotype, encoded as allele values separated by either of / or |.\cr
The allele values are:\itemize{
\item 0 for the reference allele (what is in the REF field)
\item 1 for the first allele listed in ALT
\item 2 for the second allele listed in ALT.
\item 3 for the third allele listed in ALT and so on.
}
For diploid calls examples could be 0/1, 1|0, or 1/2, etc. If a call cannot be made for a sample at a given locus, '.' is specified for each missing allele in the GT field (for example './.' for a diploid genotype and '.' for haploid genotype).\cr
The meanings of the separators are as follows:\itemize{
\item / : genotype unphased
\item | : genotype phased
}

}
\value{
%%  ~Describe the value returned
If {do_return} = TRUE:
\item{fix_and_gt}{ is a list of two matrices: filtered fixed information (without unnecessary columns) and corresponding genotype section of the VCF}
If {do_return} = FALSE:
\item{results_directory}{ a string with the directory where the results were saved.}
}
\references{
See \href{https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html}{vcfR package} for more information

And also see the "usegalaxy" \href{https://bit.ly/3oRbjx1}{VCFselectsamples} tool to pre-filter the data

%% ~put references to the literature/web site here ~
}
\author{
Tomer Antman
}
\note{
Make sure you enter proper file routes (vcf_route) such as:
%

 1) "somefolder/1011Matrix.gvcf.gz"\cr
 2) "Galaxy4_VCFselectsamples.vcf"


And also proper results route that exist (dir_results) such as:
%

 1) "somefolder"\cr
 2) "C:/Users/user/Documents"
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
files_directory <- Upload_vcf_to_R(
                                  vcf_file = "1011Matrix.gvcf.gz",
                                  dir_results = "C:/Users/user/Documents",
                                  do_return = FALSE
                                  )
print(files_directory)

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
