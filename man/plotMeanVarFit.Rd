% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMeanVarFit.R
\name{plotMeanVarFit}
\alias{plotMeanVarFit}
\title{Make a mean-variance plot}
\usage{
plotMeanVarFit(
  data,
  trans = "identity",
  invtrans = switch(trans, identity = "identity", log = "exp"),
  main = paste(switch(trans, identity = "No", log = "log"), "transformation"),
  log = switch(trans, identity = "", log = "y", ""),
  ...
)
}
\arguments{
\item{data}{a dataset or matrix with d1, d2 and effect column}

\item{trans, invtrans}{the transformation function for the variance and its inverse, possibly as strings}

\item{main}{the title of the plot}

\item{log}{log-transform of the axes, as in plot()}

\item{...}{passed on to plot()}
}
\value{
Plots the mean-variance trend
}
\description{
Make a mean-variance plot
}
\details{
This is a crucial graphical check for deciding on the
}
