
#' Metadata associated with one-population example
#'
#' A dataframe.
#'
#' @format ## `pedMeta`
#' A 'data.frame' of 2100 individuals (rows) with 9 variables (cols):
#' \describe{
#'   \item{id}{Integer individual ID}
#'   \item{population}{Population code. A, B or AB}
#'   \item{generation}{Generation of the individual}
#'   \item{mid}{dam ID}
#'   \item{fid}{sire ID}
#'   \item{gv1}{genetic value}
#'   \item{pv1}{phenotypic value}
#'   \item{gv2}{genetic value}
#'   \item{pv2}{phenotypic value}
#' }
#' @source Simulation
"pedMeta"




#' L inverse matrix of the one-population example pedigree
#'
#' An L inverse matrix generated from an AlphaSimR simulation of 20 generations.
#'
#'
#' @format ## `pedLInv`
#' Matrix object of class 'spam' of dimension 2100x2100,
#'     with 6100 (row-wise) nonzero elements.
#'     Density of the matrix is 0.138%.
#' Class 'spam' (32-bit)
#' @source Simulation
"pedLInv"

#' Metadata associated with the two-population example pedigree
#'
#' A dataframe.
#'
#' @format ## `pedMeta2`
#' A 'data.frame' of 2650 individuals (rows) with 12 variables (cols):
#' \describe{
#'   \item{id}{Integer individual ID}
#'   \item{population}{Population code. A, B or AB}
#'   \item{generation}{Generation of the individual}
#'   \item{mid}{dam ID}
#'   \item{fid}{sire ID}
#'   \item{gv1}{genetic value}
#'   \item{pv1}{phenotypic value}
#'   \item{gv2}{genetic value}
#'   \item{pv2}{phenotypic value}
#'   \item{gv}{genetic value}
#'   \item{pv}{phenotypic value}
#'   \item{generationPlotShift}{for plotting}
#' }
#' @source Simulation
"pedMeta2"




#' L inverse matrix of the two-population example pedigree
#'
#' An L inverse matrix generated from an AlphaSimR simulation of 20 generations.
#' Two diverged populations A and B. After a number of
#' generations, crossbreeding starts.
#'
#'
#' @format ## `pedLInv2`
#' Matrix object of class 'spam' of dimension 2650x2650,
#'     with 7750 (row-wise) nonzero elements.
#'     Density of the matrix is 0.11%.
#'     Class 'spam' (32-bit)
#' @source Simulation
"pedLInv2"



#' Metadata associated with the four-population example pedigree
#'
#' A dataframe.
#'
#' @format ## `pedMeta4`
#' A 'data.frame' of 4200 individuals (rows) with 9 variables (cols):
#' \describe{
#'   \item{id}{Integer individual ID}
#'   \item{population}{Population code. A, B, C, D, or ABCD}
#'   \item{generation}{Generation of the individual}
#'   \item{mid}{dam ID}
#'   \item{fid}{sire ID}
#'   \item{gv1}{genetic value}
#'   \item{pv1}{phenotypic value}
#'   \item{gv2}{genetic value}
#'   \item{pv2}{phenotypic value}
#' }
#' @source Simulation
"pedMeta4"




#' L inverse matrix of the four-population example pedigree
#'
#' An L inverse matrix generated from an AlphaSimR simulation of 20 generations.
#' One population, ABCD, is split into four, A, B, C, D.
#'
#'
#'
#' @format ## `pedLInv4`
#' Matrix object of class 'spam' of dimension 4200x4200,
#'     with 12200 (row-wise) nonzero elements.
#'     Density of the matrix is 0.0692%.
#'     Class 'spam' (32-bit)
#' @source Simulation
"pedLInv4"
