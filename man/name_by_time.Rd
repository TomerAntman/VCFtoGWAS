\name{name_by_time}
\alias{name_by_time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create name by time stamp
%%  ~~function to do ... ~~
}
\description{
The function creates a string from the current date and time (using \code{\link[base:Sys.time]{Sys.time()}}) that is used for default unique folder names
}
\usage{
name_by_time()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{results_name}{The time-stamp string that can be used as a file- or folder- name}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tomer Antman}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## The function is currently defined as
function ()
{
    results_name <- gsub(":", ".", Sys.time())
    results_name <- gsub(" ", "_", results_name)
    results_name <- substr(results_name, 1, nchar(results_name) -
        3)
    return(results_name)
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
