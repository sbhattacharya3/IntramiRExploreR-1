#' Visualises the targetGene:miRNA network using Cytoscape and igraph .
#'
#' @param mRNA  character. gene Identifier.
#' @param mRNA_type  character. mRNA id type.
#'                   The choices are 'GeneSymbol','FBID' and 'CGID'.
#' @param method  character. Statistical Methods. Choices are
#'                'Pearson','Distance','Both'
#' @param platform  character. Affymetrix Platforms. Choices are
#'                'Affy1','Affy2'.
#' @param visualisation  character.Visualisation type.
#'                       Choices are 'igraph','Cytoscape','text'and "console"
#' @param path character. Path where data.frame is saved when visualisation
#'             is text. Default is tempdir().
#' @param layout  character. Network choices. Choices are
#'              'kamadakawai','reingold.tilford','fruchterman.reingold' 
#'                and 'interactive'.
#' @return Depending upon the ouput choice network image or dataframe 
#'         containg miRNAs that target the query gene are ouput.
#' @examples
#' mRNA='Syb'
#' Gene_Visualisation(mRNA,mRNA_type=c('GeneSymbol'),method=c('Pearson'),
#'                   platform=c('Affy1'), visualisation = "console")
#' @import igraph utils grDevices graphics 
#' @importFrom stats na.omit 
#' @export
Gene_Visualisation <- function(mRNA, mRNA_type = c("GeneSymbol", "FBGN", 
    "CGID"), method = c("Pearson", "Distance", "Both"), platform = c("Affy1", 
    "Affy2"), visualisation = c("igraph", "Cytoscape", "text", "console"), 
    path = tempdir(), layout = c("kamadakawai", "reingold.tilford", "fruchterman.reingold", 
        "interactive")) {
    ### Checking if input is valid
    stopifnot(is.character(mRNA), length(mRNA) > 0)
    stopifnot(is.character(mRNA_type), length(mRNA_type) > 0)
    stopifnot(is.character(method), length(method) > 0, method %in% c("Pearson", 
        "Distance", "Both"))
    stopifnot(is.character(platform), length(platform) > 0, platform %in% 
        c("Affy1", "Affy2"))
    stopifnot(is.character(visualisation), length(visualisation) > 0)
    ### If length >1
    if (length(mRNA) > 1) {
        
        dat <- c()
        ## Extracts and sorts the data in descending order
        for (ii in 1:length(mRNA)) {
            a <- tryCatch({
                genes_Stat(mRNA[ii], mRNA_type, method, Platform = platform)
                # return(expr)
            }, error = function(e) {
                # stop('mRNA record Absent')
                return(NA)
            })
            if (nrow(a) == 0) {
                stop("mRNA does not have any records")
            } else {
                commgenes <- as.character(a$Gene)
                FBGN <- as.character(a$Gene_FBGN)
                CGID <- as.character(a$Genes_CGID)
                sc <- as.numeric(a$Score)
                ord <- a[order(-sc), ]
                dat1 <- ord
                dat <- rbind(dat, dat1)
            }
            
            
            
        }
        ## If Visualisation is selected as Text If Visualisation is selected as
        ## Text, writes in a .csv file Else If Visualisation is selected as
        ## Cytoscape, writes in a .csv file in cytoscape format output Else If
        ## Visualisation is selected as igraph, outputs in the format the layout
        ## is selected Else If Visualisation is selected as console, outputs the
        ## dataframe to the console. Else prints error.
        if (identical(visualisation, "Text")) {
            ### If miRNA length >1 names concatenated with & else use the miR name
            if (length(mRNA) > 1) {
                mRNAs <- paste(mRNA, collapse = "&")
            } else {
                mRNAs = mRNA
            }
            filename = paste("Putative_miRNAs_for_", mRNAs, "_", platform, 
                "_and_", method, ".csv", sep = "")
            ### Check if the file name exists and if not print data in .csv format
            if (!(file.exists(file.path(path, filename)))) {
                write.csv(dat, file = file.path(path, filename))
            } else {
                stop("File Already Exists!!")
            }
            
        } else if (identical(visualisation, "Cytoscape")) {
            ## Extracting individual data and storing it in the database.
            miR <- as.character(dat$miRNA)
            ## Checks if the mRNA_Type selected to print in the output table is
            ## GeneSymbol, Flybase Id (FBID) or CGID
            if (identical(mRNA_type, "GeneSymbol")) {
                gene <- as.character(dat$Gene)
            } else if (identical(mRNA_type, "FBGN")) {
                gene <- as.character(dat$Gene_FBGN)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Genes_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            ## Storing data in a dataframe
            datC <- data.frame(Gene = gene, Interaction = rep("pp", length(miR)), 
                miRNA = miR, Score)
            if (length(mRNA) > 1) {
                mRNAs <- paste(mRNA, collapse = "&")
            } else {
                mRNAs = mRNA
            }
            filename = paste("Cytoscape_output_for_", mRNAs, "_", platform, 
                "_and_", method, ".csv", sep = "")
            ## Checking if the file exists and writing the data in .csv format
            if (!(file.exists(file.path(path, filename)))) {
                write.csv(datC, file = file.path(path, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else if (identical(visualisation, "igraph")) {
            ### Extracting the data and storing in the data.frame
            miR <- as.character(dat$miRNA)
            ## Checks if the mRNA_Type selected to print in the output table is
            ## GeneSymbol, Flybase Id (FBID) or CGID
            if (identical(mRNA_type, "GeneSymbol")) {
                gene <- as.character(dat$Gene)
            } else if (identical(mRNA_type, "FBGN")) {
                gene <- as.character(dat$Gene_FBGN)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Genes_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            datN <- data.frame(Gene = gene, miRNA = miR, Score)
            ### Creating the plot
            a <- graph.data.frame(datN, directed = FALSE)
            ### Creating the Edge weight
            E(a)$weight <- Score
            ### Selecting the layout to represent the network. The selected network
            ### is printed out in the format of a .tiff file.  If chosen layout is
            ### kamadakawai, prints as .tiff file If chosen layout is
            ### fruchterman.reingold, prints as .tiff file If chosen layout is
            ### reingold.tilford, prints as .tiff file If chosen layout is
            ### interactive, prints as a tkplot file Else Prints error message
            if (identical(layout, "kamadakawai")) {
                layout_graph <- layout.kamada.kawai(a)
                if (length(mRNA) > 1) {
                  mRNAs <- paste(mRNA, collapse = "&")
                } else {
                  mRNAs = mRNA
                }
                filename = paste("Graphical_output_for_", mRNAs, "_", platform, 
                  "_and_", method, ".tiff", sep = "")
                if (!(file.exists(file.path(path, filename)))) {
                  tiff(filename = file.path(path, filename), width = 1000, 
                    height = 1000, units = "px")
                  plot(a, layout = layout_graph, vertex.color = "red", 
                    edge.width = E(a)$weight)
                  dev.off()
                } else {
                  stop("File Already Exists!!")
                }
                
            } else if (identical(layout, "fruchterman.reingold")) {
                layout_graph <- layout.fruchterman.reingold(a)
                if (length(mRNA) > 1) {
                  mRNAs <- paste(mRNA, collapse = "&")
                } else {
                  mRNAs = mRNA
                }
                filename = paste("Graphical_output_for_", mRNAs, "_", platform, 
                  "_and_", method, ".tiff", sep = "")
                if (!(file.exists(file.path(path, filename)))) {
                  tiff(filename = file.path(path, filename), width = 1000, 
                    height = 1000, units = "px")
                  plot(a, layout = layout_graph, vertex.color = "red", 
                    edge.width = E(a)$weight)
                  dev.off()
                } else {
                  stop("File Already Exists!!")
                }
                
            } else if (identical(layout, "reingold.tilford")) {
                layout_graph <- layout.reingold.tilford(a)
                if (length(mRNA) > 1) {
                  mRNAs <- paste(mRNA, collapse = "&")
                } else {
                  mRNAs = mRNA
                }
                filename = paste("Graphical_output_for_", mRNAs, "_", platform, 
                  "_and_", method, ".tiff", sep = "")
                if (!(file.exists(file.path(path, filename)))) {
                  tiff(filename = file.path(path, filename), width = 1000, 
                    height = 1000, units = "px")
                  plot(a, layout = layout_graph, vertex.color = "red", 
                    edge.width = E(a)$weight)
                  dev.off()
                } else {
                  stop("File Already Exists!!")
                }
            } else if (identical(layout, "interactive")) {
                tkplot(a, vertex.color = "red", edge.width = E(a)$weight)
            } else {
                stop("Invalid Layout!!!")
            }
        } else if (identical(visualisation, "console")) {
            ## Return data to the console
            return(dat)
        } else {
            stop("Visualisation method absent")
        }
    } else {
        ## Extracts and sorts the data in descending order
        dat <- c()
        a <- tryCatch({
            genes_Stat(mRNA, mRNA_type, method, Platform = platform)
            # return(expr)
        }, error = function(e) {
            # stop('mRNA record Absent')
            return(NA)
        })
        if (nrow(a) == 0) {
            stop("mRNA does not have any records")
        } else {
            commgenes <- as.character(a$Gene)
            FBGN <- as.character(a$Gene_FBGN)
            CGID <- as.character(a$Genes_CGID)
            sc <- as.numeric(a$Score)
            ord <- a[order(-sc), ]
            dat1 <- ord
            dat <- dat1
            # print(dat[1, ])
        }
        ## If Visualisation is selected as Text If Visualisation is selected as
        ## Text, writes in a .csv file Else If Visualisation is selected as
        ## Cytoscape, writes in a .csv file in cytoscape format output Else If
        ## Visualisation is selected as igraph, outputs in the format the layout
        ## is selected Else If Visualisation is selected as console, outputs the
        ## dataframe to the console. Else prints error.
        if (identical(visualisation, "Text")) {
            
            mRNAs = mRNA
            
            filename = paste("Putative_mRNAs_for_", mRNAs, "_", platform, 
                "_and_", method, ".csv", sep = "")
            ### Check if the file name exists and if not print data in .csv format
            if (!(file.exists(file.path(path, filename)))) {
                write.csv(dat, file = file.path(path, filename))
            } else {
                stop("File Already Exists!!")
            }
            
        } else if (identical(visualisation, "Cytoscape")) {
            ## Extracting individual data and storing it in the database.
            miR <- as.character(dat$miRNA)
            ## Checks if the mRNA_Type selected to print in the output table is
            ## GeneSymbol, Flybase Id (FBID) or CGID
            if (identical(mRNA_type, "GeneSymbol")) {
                gene <- as.character(dat$Gene)
            } else if (identical(mRNA_type, "FBGN")) {
                gene <- as.character(dat$Gene_FBGN)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Genes_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            datC <- data.frame(Gene = gene, Interaction = rep("pp", length(miR)), 
                miRNA = miR, Score)
            mRNAs = mRNA
            filename = paste("Cytoscape_output_for_", mRNAs, "_", platform, 
                "_and_", method, ".csv", sep = "")
            ## Checking if the file exists and writing the data in .csv format
            if (!(file.exists(file.path(path, filename)))) {
                write.csv(datC, file = file.path(path, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else if (identical(visualisation, "igraph")) {
            ### Extracting the data and storing in the data.frame
            miR <- as.character(dat$miRNA)
            ## Checks if the mRNA_Type selected to print in the output table is
            ## GeneSymbol, Flybase Id (FBID) or CGID
            if (identical(mRNA_type, "GeneSymbol")) {
                gene <- as.character(dat$Gene)
            } else if (identical(mRNA_type, "FBGN")) {
                gene <- as.character(dat$Gene_FBGN)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Genes_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            datN <- data.frame(Gene = gene, miRNA = miR, Score)
            ### Creating the plot
            a <- graph.data.frame(datN, directed = FALSE)
            ### Creating the Edge weight
            E(a)$weight <- Score
            ### Selecting the layout to represent the network. The selected network
            ### is printed out in the format of a .tiff file.  If chosen layout is
            ### kamadakawai, prints as .tiff file If chosenfruchterman.reingold,
            ### prints as .tiff file If chosen layout isreingold.tilford, prints as
            ### .tiff file If chosen layout isinteractive, prints as a tkplot file.
            ### Else Prints error message
            if (identical(layout, "kamadakawai")) {
                layout_graph <- layout.kamada.kawai(a)
                mRNAs = mRNA
                
                filename = paste("Graphical_output_for_", mRNAs, "_", platform, 
                  "_and_", method, ".tiff", sep = "")
                if (!(file.exists(file.path(path, filename)))) {
                  tiff(filename = file.path(path, filename), width = 1000, 
                    height = 1000, units = "px")
                  plot(a, layout = layout_graph, vertex.color = "red", 
                    edge.width = E(a)$weight)
                  dev.off()
                } else {
                  stop("File Already Exists!!")
                }
            } else if (identical(layout, "fruchterman.reingold")) {
                layout_graph <- layout.fruchterman.reingold(a)
                mRNAs = mRNA
                
                filename = paste("Graphical_output_for_", mRNAs, "_", platform, 
                  "_and_", method, ".tiff", sep = "")
                if (!(file.exists(file.path(path, filename)))) {
                  tiff(filename = file.path(path, filename), width = 1000, 
                    height = 1000, units = "px")
                  plot(a, layout = layout_graph, vertex.color = "red", 
                    edge.width = E(a)$weight)
                  dev.off()
                } else {
                  stop("File Already Exists!!")
                }
            } else if (identical(layout, "reingold.tilford")) {
                layout_graph <- layout.reingold.tilford(a)
                mRNAs = mRNA
                
                filename = paste("Graphical_output_for_", mRNAs, "_", platform, 
                  "_and_", method, ".tiff", sep = "")
                if (!(file.exists(file.path(path, filename)))) {
                  tiff(filename = file.path(path, filename), width = 1000, 
                    height = 1000, units = "px")
                  plot(a, layout = layout_graph, vertex.color = "red", 
                    edge.width = E(a)$weight)
                  dev.off()
                } else {
                  stop("File Already Exists!!")
                }
            } else if (identical(layout, "interactive")) {
                tkplot(a, vertex.color = "red", edge.width = E(a)$weight)
            } else {
                stop("Invalid Layout!!!")
            }
        } else if (identical(visualisation, "console")) {
            ## Return data to the console print(dat)
            return(dat)
        } else {
            stop("Visualisation method absent")
        }
    }
    # return(dat)
}
