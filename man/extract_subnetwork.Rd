% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_subnetwork.R
\name{extract_subnetwork}
\alias{extract_subnetwork}
\title{Extract a subnetwork based on pathway genes}
\usage{
extract_subnetwork(network, enrich_result, pathway_name)
}
\arguments{
\item{network}{Input network object; output from \code{build_network}}

\item{enrich_result}{Pathway enrichment result, output from \code{enrich_network}}

\item{pathway_name}{Name of the pathway determining what genes (nodes) are
pulled from the input network.}
}
\value{
\code{tidygraph} object for plotting or further analysis
}
\description{
Extract a subnetwork based on pathway genes
}
\details{
Uses functions from the igraph package to extract a minimally
connected subnetwork or module from the starting network, using genes from
a given pathway as the basis. To see what genes were pulled out for the
pathway, check \code{attr(x, "starters")}.
}
\references{
Code for network module (subnetwork) extraction was based off of
that used in jboktor/NetworkAnalystR on Github.
}
\seealso{
\url{https://www.github.com/travis-m-blimkie/networker}
}
