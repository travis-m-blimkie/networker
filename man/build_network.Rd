% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_network.R
\name{build_network}
\alias{build_network}
\title{Construct a PPI network from input genes and InnateDB's database}
\usage{
build_network(df, col, order, ppi_data)
}
\arguments{
\item{df}{Input data frame containing genes of interest}

\item{col}{Column of input genes as Ensembl IDs (character)}

\item{order}{Desired network order. Possible options are "zero" (default),
"first," "min_simple," or "min_steiner."}

\item{ppi_data}{Data frame of InnateDB PPI data; minimally should contain
rows of interactions as pairs of Ensembl gene IDs, e.g. "ensembl_gene_A"
and "ensembl_gene_B"}
}
\value{
\code{tidygraph} object for plotting
}
\description{
Construct a PPI network from input genes and InnateDB's database
}
\references{
None.
}
\seealso{
\url{https://www.github.com/travis-m-blimkie/networker}
}