#' Gene ontology for Target Genes.
#' 
#' @param gene List  A String or vector containing the Gene names.
#' @param GO   A String depicting the chosen GO tool. Choices are
#'             "David" and "topGO"
#' @param term  A String depicting the chosen term. Choices are
#'               "GOTERM_BP_ALL","GOTERM_MF_ALL", "GOTERM_CC_ALL".
#' @param email Email Id to connect to David.
#' @param geneIdType Type of gene Id given as input. Default "ALIAS"
#' @param filename Name of the file to store Gene Ontology.
#' @param ontology Ontology selection for topGO. Choices are
#'                 "GO_BP","GO_MF","GO_CC".
#' @param path  String. The path where the data is stored if TEXT=TRUE.
#' @return Depending upon the ouput choice data is stored in the path 
#'         specified. Default option prints output to the console.
#' @examples
#' \dontrun{
#' miR="dme-miR-12"
#' a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),
#'    platform=c("Affy1"),thresh=100)
#' genes<-a$Target_GeneSymbol
#' GetGOS_ALL(genes,GO=c("topGO"),term=c("GO_BP"),path=tempdir(),
#'      filename="test")
#'  }
#' @import FGNet
#' @export
GetGOS_ALL <- function(gene, GO = c("DAVID", "topGO"), term = c("GOTERM_BP_ALL", 
    "GOTERM_MF_ALL", "GOTERM_CC_ALL"), geneIdType = "ALIAS", email, path = tempdir(), 
    ontology = c("GO_BP", "GO_MF", "GO_CC"), filename) {
    #Checking the gene length, term and ontology is entered correctly
    stopifnot(is.character(gene), length(gene) > 20)
    stopifnot(is.character(term), length(term) == 1)
    stopifnot(is.character(ontology), length(ontology) == 1)
    stopifnot(is.character(GO), length(GO) == 1, GO %in% c("DAVID", "topGO"))
    ### If GO chosen as DAVID do DAVID functional clustering Else do topGO
    ### functional Clustering
    if (identical(GO, "DAVID")) {
        ### Checking whether the email is correct
        stopifnot(is.character(email), length(email) > 0)
        ### If RDAVIDWebService is not installed, prompts you to install it.
        if (!requireNamespace("RDAVIDWebService", quietly = TRUE)) {
            stop("Please install RDAVIDWebService")
        } else {
            requireNamespace("RDAVIDWebService", quietly = TRUE)
        }
        ### Gene List input to David, preprocessing.
        geneList <- gene
        geneList1 <- gsub("FBGN", "FBgn", as.character(geneList))
        geneList1 <- gsub("FBGN", "fbgn", as.character(geneList))
        ### URL for DAVID
        url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/"
        ### Connecting to DAVID
        david <- RDAVIDWebService::DAVIDWebService(email, url)
        RDAVIDWebService::connect(david)
        ### Obtaining result from DAVID
        result <- RDAVIDWebService::addList(david, geneList1, idType = "FLYBASE_GENE_ID", 
            listName = "Targets", listType = "Gene")
        term <- as.character(term)
        RDAVIDWebService::setAnnotationCategories(david, term)
        termCluster <- RDAVIDWebService::getClusterReport(david, type = "Term")
        RDAVIDWebService::getClusterReportFile(david, type = "Term", fileName = paste(path, 
            "//", term, "_ClusterReport.tab", sep = ""))
        genetemp <- gProfileR::gconvert(unique(as.character(geneList1)), 
            organism = "dmelanogaster", target = "ENSG")
        fbgn <- as.character(genetemp$alias)
        gene <- as.character(genetemp$name)
        names(gene) <- as.character(fbgn)
        ### Result Input to FBGN
        filename <- paste(term, "_ClusterReport.tab", sep = "")
        fearesults <- FGNet::format_david(file.path(path, filename), jobName = "DavidAnalysis", 
            geneLabels = gene)
        FGNet::FGNet_report(fearesults, plotKeggPw = FALSE)
    } else {
        stopifnot(is.character(geneIdType), length(geneIdType) > 0)
        ### Preprocessing input to topgo
        filename <- paste(filename, "_topGO", sep = "")
        gene <- gene
        term <- as.character(term)
        geneIdType = geneIdType
        filename <- paste(term, "TOPGOAnalysis", sep = "")
        ### Input data to topGO via FBGN
        feaResults_topGO <- FGNet::fea_topGO(gene, geneIdType, geneLabels = NULL, 
            organisms = "Dm", annotations = ontology, genesUniverse = NULL, 
            refPackage = NULL, geneID2GO = NULL, nodeSize = 5, pValThr = 0.01, 
            testStat = NULL, jobName = file.path(path, filename))
        FGNet_report(feaResults_topGO, plotKeggPw = FALSE)
    }
    
}
