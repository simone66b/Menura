% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fossil_id.R
\name{fossil_id}
\alias{fossil_id}
\title{fossil_id}
\usage{
fossil_id(ftr)
}
\arguments{
\item{ftr}{Single evolutionary tree as an object of the "phylo" class in
the \code{ape} package.}
}
\description{
Create a vector of the location of the fossils in a tree
}
\details{
This function locates the tip values of fossils in the tree.  
It is primarily used in other functions of menura to handle 
fossils data but can be useful to locate the fossils in the tree.
}
