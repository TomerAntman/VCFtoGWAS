\name{create_directory}
\alias{create_directory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a directory
%%  ~~function to do ... ~~
}
\description{
Function ment to create subdirectories to save the results of the current step.\cr

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create_directory(called_from,
                 dir_results = \code{\link[base]{getwd}}(),
                 results_name = \code{\link{name_by_time}}())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{called_from}{
The function (and step in the process) from which the function was called. Is added to the name of the subdirectory that is created by the function.
}
  \item{dir_results}{
The directory in which a subfolder will be created and results will be saved. Make sure it exists!!!\cr
If entered directory doesn't exist (or no direcotry was entered), the subfolder will be created in the current working directory.
}
  \item{results_name}{
The results_name parameter dictates the folder (within dir_results folder) that will be created to save the results.\cr
If the results name isn't changed from the default (\code{\link{name_by_time}}), then a name based on time will be constructed. Step name is always added.
Examples: \itemize{
\item By default: \emph{"2021-12-15_10.58_Step1-Upload_VCF-results"}
\item Non default: \emph{"OnlyHomozygousStrains_Step1-Upload_VCF-results"}
} %% close itemize
} %% close results_name
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
## The function is currently defined as
function (called_from, dir_results = getwd(), results_name = name_by_time())
{
    results_name <- paste0(results_name, "_", deparse(substitute(called_from)),
        "-results")
    results_directory <- paste0(dir_results, "/", results_name)
    dir.create(results_directory)
    if (grepl("No such file or directory", names(last.warning))) {
        message(results_name, " folder couldn't be created because\n",
            dir_results, " directory doesn't exist.\nTrying to create ",
            results_name, " folder in\n", getwd())
        results_directory2 <- paste0(getwd(), "/", results_name)
        dir.create(results_directory2)
        if (grepl("No such file or directory", names(last.warning))) {
            stop("Couldn't create a directory to save the results")
        }
        else {
            results_directory <- results_directory2
            message("Results will be saved into: \n", results_directory)
        }
    }
    else if (grepl("already exists", names(last.warning))) {
        message("The results directory already exists, the results will be saved in:\n",
            results_directory)
    }
    else {
        message("Results directory created:\n", results_directory)
    }
    return(results_directory)
  }
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
