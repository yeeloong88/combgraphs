\encoding{UTF-8}
\name{combgraphs{OutbreakTools}}
\alias{combgraphs}
\alias{combgraphs.default}
\alias{combgraphs.list}
\title{Combining undirected graphs
}
\description{
  These functions are under development. Please email the author before using them for published work. 
  \cr\cr
  The function \code{combgraphs} generate individual graphs based on multiple distance data types such that pairs of individuals are connected if and only if they are distant by less than their respective given threshold distances. The individual graphs are combined through the union or intersection of multiple individual graphs. Graph algorithms and classes from the \code{\link[gengraph]{gengraph}} and \code{\link[igraph]{igraph}} packages are used. 
  \cr\cr
  \code{combgraphs} is a generic function with methods for \code{list} of \code{matrix} or \code{dist}.
  \cr\cr
}
\usage{
combgraphs(x, \dots)

\method{combgraphs}{list}(x, ngrp=NULL, plotAll=FALSE, outputAll=FALSE, 
            plot=TRUE, layout=layout.circle, method=c("union","intersect"), 
            col.pal=funky ,truenames=TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{list} of \code{matrix} or \code{dist} object. Each \code{matrix} or \code{dist} represents pairwise distances between considered individuals.}
  \item{ngrp}{an \code{integer} indicating the number of groups to looked for. A message is issued if the exact number could not be found.}
  \item{plotAll}{a \code{logical} indicating if the individual plots for each item in the \code{x} should be drawn; if TRUE an interactive selection interface allows the user to select the individual plots to be drawn; this operation can take time for large and highly-connected graphs; defaults to FALSE.}
  \item{outputAll}{a \code{logical} indicating if the results for each individual plot should be returned; defaults to FALSE.}
  \item{plot}{a \code{logical} indicating if the combined graph should be plotted; this operation can take time for large and highly-connected graphs; defaults to TRUE.}
  \item{layout}{a \code{function} that determines the placement of vertiecs for drawing a graph.}
  \item{method}{a \code{character} that determines if the "union" or "intersection" method should be used for combining graphs; defaults to "union".}
  \item{col.pal}{a color palette to define group colors.}
  \item{truenames}{a \code{logical} indicating if original labels should be used for plotting (TRUE), as opposed to indices of sequences (FALSE); defaults to TRUE.}
  \item{\dots}{further arguments to be used by other functions.}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The class combgraphs is a list with the following components:
\cr
\item{graph}{a graph of class \code{igraph}.}
\item{clust}{a list containing group information: \code{$membership}: an integer giving group memberships; \code{$csize}: the size of each cluster; \code{$no}: the number of clusters.}
\item{col}{the color used to plot each group. }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Original idea by Thibaut Jombart. Implementation by Low Yee Loong \email{?email?}.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The \code{\link[gengraph]{gengraph}} and \code{\link[igraph]{igraph}} package.
}
\examples{
\dontrun{
set.seed(24)
nrand <- 10

dat1 <- round(c(runif(nrand,10,15),runif(nrand,50,55),runif(nrand,90,95)))
dat2 <- round(c(runif(nrand,0,5),runif(nrand,30,35),runif(nrand,80,85)))
dat3 <- round(c(runif(nrand,0,1),runif(nrand,10,11),runif(nrand,199,200)))

dist1 <- dist(dat1)
dist2 <- dist(dat2)
dist3 <- dist(dat3)

newlist <- list(dist1,dist2,dist3)

par(mfrow=c(1,2))
combgraphs(newlist, ngrp=3, plotAll=F)  #method="union"
combgraphs(newlist, ngrp=3, plotAll=F, method="intersect")  #method="intersect"
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
