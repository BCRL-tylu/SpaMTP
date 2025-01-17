% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConvertingBetweenObjects.R
\name{CardinalToSeurat}
\alias{CardinalToSeurat}
\title{Converts a Cardinal Object into a Seurat Object}
\usage{
CardinalToSeurat(data, run_name, seurat.coord = NULL)
}
\arguments{
\item{data}{A Cardinal Object that is being converted into a Seurat Object.}

\item{run_name}{A character string defining the run name of the Cardinal data to be converted to a Seurat Object}

\item{seurat.coord}{A Data.Frame containing two columns titled 'X_new' and 'Y_new' specifying the pixel coordinates of each data point. This is only required mapping Spatial Metabolomic data with a H&E image if the Cardinal Object coordinates are not the same as the H&E image coordinates. Else, set to NULL(default).}
}
\value{
A Seurat Object containing the mz count data of the supplied Cardinal Object
}
\description{
Converts a Cardinal Object into a Seurat Object
}
\examples{
# CardinalToSeurat(CardinalObj, run_name = "run_1", seurat.coord = NULL)
}
