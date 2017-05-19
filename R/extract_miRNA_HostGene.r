#' Extract Intragenic miRNA for a given Host gene.
#'
#' @param gene  character. Gene Symbol.
#' @return miRf, a character string or vector containing
#'         Intragenic miRNA for the Host Gene.
#' @examples
#' gene="Gmap"
#' extract_intragenic_miR(gene)
#' @import utils 
#' @export

extract_intragenic_miR <- function(gene) {
    ###Checking the validity of the input
    stopifnot(is.character(gene), length(gene) > 0)
    ###Calling the database to get miRNAs in the 
    ###host genes
    data("miRNA_summary_DB", envir = environment())
    miRNA_summary_DB <- miRNA_summary_DB
    ###Extracting data
    gen <- as.character(miRNA_summary_DB$Gene)
    miR <- as.character(miRNA_summary_DB$miRNA)
    ###Preprocessing Search query
    g <- paste("^", gene, "$", sep = "")
    genf <- c()
    miRf <- c()
    if (length(g) > 1) {###If length of input >1
        ###Initiating the for loop
        for (ki in seq_along(g)) {
            ###Searching for the miRNAs for a 
            ##given host  gene
            d <- grep(g[ki], gen, ignore.case = TRUE)
            miRt <- c()
            if (length(d) > 1) { ###If multiple miRNAs in Host gene
                for (j in seq_along(d)) {
                    n <- d[j]
                    miRt <- c(miRt, miR[n])
                }
                ###Storing data
                miRf <- c(miRf, unique(as.character(miRt)))
                genf <- c(genf, rep(as.character(gene[ki]), 
                        length(unique(miRt))))
            } else { ###If single miRNAs in Host gene
                miRf <- c(miRf, as.character(miR[d]))
                genf <- c(genf, rep(as.character(gene[ki]), 1))
            }
        }
    } else if (length(g) == 1) {###If Length of query ==1 
        ###Searching for the miRNAs for a 
        ##given host  gene
        d <- grep(g, gen, ignore.case = TRUE)
        if (length(d) > 1) {###If multiple miRNAs in Host gene
            miRf <- c()
            ##Extracting data
            for (li in seq_along(d)) {
                n <- d[li]
                miRf <- c(miRf, unique(as.character(miR[n])))
            }
            
            miRf <- unique(miRf)
        } else {###If single miRNAs in Host gene
            miRf <- as.character(unique(miR[d]))
        }
    } else {### If Host gene does not have an intragenic miRNA
        stop("Host gene does not have a intragenic miRNA")
    }
    ###returns miRNA
    return(as.character(miRf))
}
#' Extract Host Gene for a given Intragenic miRNA.
#'
#' @param miRNA  A String containing the miRNA name.
#' @return genf, a character string or vector containing
#'         Host gene for the Intragenic miRNA.
#' @examples
#' miRNA="dme-miR-12"
#' extract_HostGene(miRNA)
#' @import utils 
#' @export
extract_HostGene <- function(miRNA) {
    ###Checking the validity of the input
    stopifnot(is.character(miRNA), length(miRNA) > 0, length(grep("dme-miR", 
        miRNA, ignore.case = TRUE)) > 0)
    ##Calls the database
    miRNA_summary_DB <- c()
    data("miRNA_summary_DB", envir = environment())
    miRNA_summary_DB <- miRNA_summary_DB
    gen <- as.character(miRNA_summary_DB$Gene)
    miR <- as.character(miRNA_summary_DB$miRNA)
    ###Preprocessing Search query
    m <- paste("^", miRNA, "$", sep = "")
    genf <- c()
    miRf <- c()
    if (length(m) > 1) {###If Length of query ==1 
        ###Searching for the host gene for a 
        ##given miRNA
        for (ki in seq_along(m)) {
            d <- grep(m[ki], miR, ignore.case = TRUE)
            gent <- c()
            if (length(d) > 1) {##if multiple genes
                for (j in 1:length(d)) {
                    n <- d[j]
                    gent <- c(gent, gen[n])
                }
                genf <- c(genf, unique(as.character(gent)))
                miRf <- c(miRf, rep(as.character(miR[ki]), 
                        length(unique(gent))))
            } else {
                genf <- c(genf, unique(as.character(gent[d])))
                miRf <- c(miRf, rep(as.character(miR[ki]), 1))
            }
        }
    } else if (length(m) == 1) {##if single genes
        d <- grep(m, miR, ignore.case = TRUE)
        if (length(d) > 1) {
            genef <- c()
            for (li in seq_along(d)) {
                n <- d[li]
                genef <- c(genef, unique(as.character(gen[n])))
            }
            genef <- unique(genef)
        } else {
            genef <- as.character(unique(gen[d]))
            genf <- genef
        }
    } else {
        stop("miRNA not intragenic")
    }
    ###returns Host Gene
    return(as.character(genf))
}
