% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmBlackListRegion.R
\name{rmBlackListRegion}
\alias{rmBlackListRegion}
\title{Remove the regions that overlap with the provided blacklist regions}
\usage{
rmBlackListRegion(input_region, blacklist)
}
\arguments{
\item{input_region}{GRanges object The input regions, from which the ones overlaps with the blacklist will be removed.}

\item{blacklist}{string The name of BED file containing the blacklist regions}

\item{addition_annotation}{list or NULL Optional. Name and file name of additional junction annotation source. Default: list(snaptron = system.file("extdata/Reference/snaptron_gtex_anno.txt", package = "SpliPath")). Whether the input junctions are the given junction list will be shown in the output file (in column "in.annotation_name").}
}
\description{
Remove the regions that overlap with the provided blacklist regions
}
