% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeJuncData.R
\name{mergeJuncData}
\alias{mergeJuncData}
\title{Merge junction data of samples to a junction - sample data frame}
\usage{
mergeJuncData(sample_path_file, chrom, output_prefix = "")
}
\arguments{
\item{sample_path_file}{string Name of tab-delimited file which contain two columns: 1, unique sample ids; 2, paths to the junction files.}

\item{chrom}{string The junctions in this chromosome will be merged among samples.}

\item{output_prefix}{string Prefix of output file name. The junction-sample table will be written in prefix_chr*.txt.gz.}
}
\value{
A data frame of coordinates of unique junctions in four columns: chr, start, end, and strand.
}
\description{
Merge junction data of samples to a junction - sample data frame
}
