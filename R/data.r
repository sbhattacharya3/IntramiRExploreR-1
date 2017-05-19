#' Targets for the microRNA analyzed from Affy1 plaform using Pearson.
#'
#' A precomputed dataset containing the targets, scores and other 
#' attributes of 83 intragenic microRNAs using Pearson Correlation 
#' for plaform Affymetrix 1.
#' @format A data frame with 41845 rows and 8 variables:
#' \describe{
#'   \item{miRNA}{miRNA name, miRNA symbol}
#'   \item{GeneSymbol}{Gene name, in Gene Symbol}
#'   \item{FBGN}{Gene name, in FlybaseID}
#'   \item{CGID}{Gene name, in CGID}
#'   \item{Score}{Computed Score, in float}
#'   \item{GeneFunction}{Gene Functions, from Flybase}
#'   \item{experiments}{Experiments, from ArrayExpress}
#'   \item{TargetDatabases}{Target Database Name, from TargetDatabases}
#' }
"Affy1_Pearson_Final"
# > [1] 'Affy1_Pearson_Final'
#' Targets for the microRNA analyzed from Affy1 plaform using Distance.
#'
#' A precomputed dataset containing the targets, scores and other 
#' attributes of 83 intragenic microRNAs using Distance Correlation 
#' for plaform Affymetrix 1.
#' @format A data frame with 53399 rows and 8 variables:
#' \describe{
#'   \item{miRNA}{miRNA name, miRNA symbol}
#'   \item{GeneSymbol}{Gene name, in Gene Symbol}
#'   \item{FBGN}{Gene name, in FlybaseID}
#'   \item{CGID}{Gene name, in CGID}
#'   \item{Score}{Computed Score, in float}
#'   \item{GeneFunction}{Gene Functions, from Flybase}
#'   \item{experiments}{Experiments, from ArrayExpress}
#'   \item{TargetDatabases}{Target Database Name, from TargetDatabases}
#' }
"Affy1_Distance_Final"
# > [1] 'Affy1_Distance_Final'
#' Targets for the microRNA analyzed from Affy2 plaform using Pearson.
#'
#' A precomputed dataset containing the targets, scores and other 
#' attributes of 83 intragenic microRNAs using Pearson Correlation 
#' for plaform Affymetrix 1.
#' @format A data frame with 52913 rows and 8 variables:
#' \describe{
#'   \item{miRNA}{miRNA name, miRNA symbol}
#'   \item{GeneSymbol}{Gene name, in Gene Symbol}
#'   \item{FBGN}{Gene name, in FlybaseID}
#'   \item{CGID}{Gene name, in CGID}
#'   \item{Score}{Computed Score, in float}
#'   \item{GeneFunction}{Gene Functions, from Flybase}
#'   \item{experiments}{Experiments, from ArrayExpress}
#'   \item{TargetDatabases}{Target Database Name, from TargetDatabases}
#' }
"Affy2_Pearson_Final"
# > [1] 'Affy2_Pearson_Final'
#' Targets for the microRNA analyzed from Affy2 plaform using Distance.
#'
#' A precomputed dataset containing the targets, scores and other 
#' attributes of 83 intragenic microRNAs using Distance Correlation 
#' for plaform Affymetrix 1.
#' @format A data frame with 73374 rows and 8 variables:
#' \describe{
#'   \item{miRNA}{miRNA name, miRNA symbol}
#'   \item{GeneSymbol}{Gene name, in Gene Symbol}
#'   \item{FBGN}{Gene name, in FlybaseID}
#'   \item{CGID}{Gene name, in CGID}
#'   \item{Score}{Computed Score, in float}
#'   \item{GeneFunction}{Gene Functions, from Flybase}
#'   \item{experiments}{Experiments, from ArrayExpress}
#'   \item{TargetDatabases}{Target Database Name, from TargetDatabases}
#' }
"Affy2_Distance_Final"
# > [1] 'Affy2_Distance_Final'
#' Contains the summary for the intragenic miRNA.
#'
#' A dataset containing the summary for the intragenic miRNA.
#' @format A data frame with 257 rows and 6 variables:
#' \describe{
#'   \item{miRNA}{miRNA name, miRNA symbol}
#'   \item{Intragenic}{Responsee, in boolean}
#'   \item{Intergenic}{Responsee, in boolean}
#'   \item{Gene}{miRNA name, miRNA symbol}
#'   \item{Type.of.HostGene.mRNA.lncRNA.}{Type of Hostgene}
#'   \item{Notes}{Comments about the miRNA}
#' }
"miRNA_summary_DB"
# > [1] 'miRNA_summary_DB'
#' Contains the miRNA function information from Flybase database.
#'
#' A dataset containing the function for the intragenic miRNA.
#' @format A data frame with 66 rows and 4 variables:
#' \describe{
#'   \item{miRNA}{miRNA name, miRNA symbol}
#'   \item{FBGN}{target gene name, gene symbol}
#'   \item{miRNAFunction}{miRNA function, from Flybase}
#' }
#' @source \url{http://flybase.org/}
"miRNA_ID_to_Function"
# > [1] 'miRNA_ID_to_Function' 
