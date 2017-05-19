#' Visualises the targetGene:miRNA network using Cytoscape and igraph .
#'
#' @param miRNA  character. miRNA Identifier.
#' @param mRNA_type  character. mRNA id type.
#'                   The choices are 'GeneSymbol','FBID' and 'CGID'.
#' @param method  character. Statistical Methods. Choices are
#'                'Pearson','Distance','Both'
#' @param platform  character. Affymetrix Platforms. Choices are
#'                'Affy1','Affy2'.
#' @param visualisation  character.Visualisation type.
#'                       Choices are 'igraph','Cytoscape','text' and 'console'.
#' @param path character. Path where data.frame is saved when visualisation
#'             is text. Default is tempdir().
#' @param thresh integar. Threshold depicting number of rows to show.
#' @param layout  character. Network choices. Choices are
#'                'kamadakawai','reingold.tilford','fruchterman.reingold' 
#'                and 'interactive'.
#' @return Depending upon the ouput choice network image or dataframe 
#'         containg miRNAs that target the query gene are ouput.
#' @examples
#' miRNA='dme-miR-12'
#' Visualisation(miRNA,mRNA_type=c('GeneSymbol'),method=c('Pearson'),
#' platform=c('Affy1'),visualisation=c('igraph'),layout=c('kamadakawai'),
#'            path=tempdir())
#' @import igraph utils grDevices graphics  
#' @importFrom stats na.omit
#' @export

Visualisation <- function(miRNA, mRNA_type = c("GeneSymbol", "FBID", "CGID"), 
    method = c("Pearson", "Distance", "Both"), platform = c("Affy1", "Affy2"), 
    thresh = 50, visualisation = c("igraph", "Cytoscape", "Text", "console"), 
    path = tempdir(), layout = c("kamadakawai", "reingold.tilford",
    "fruchterman.reingold", "interactive")) {
    ### Checks whether the length of query is greater than 1
    stopifnot(is.character(miRNA), length(miRNA) > 0, length(grep("dme-miR", 
        miRNA, ignore.case = TRUE)) > 0)
    stopifnot(is.character(mRNA_type), length(mRNA_type) > 0)
    stopifnot(is.character(method), length(method) > 0, method %in% c("Pearson", 
        "Distance", "Both", "BothIntersect"))
    stopifnot(is.character(platform), length(platform) > 0, platform %in% 
        c("Affy1", "Affy2"))
    stopifnot(is.character(visualisation), length(visualisation) > 0)
    ## miRNA length >1
    if (length(miRNA > 1)) {
        dat <- c()
        
        ## Extracts and sorts the data in descending order
        for (ii in seq_along(miRNA)) {
            stopifnot(length(extract_HostGene(miRNA[ii])) > 0)
            # print(miRNA[ii])
            a <- tryCatch({
                miRTargets_Stat(miRNA[ii], method, Platform = platform)
                
            }, error = function(e) {
                # print('mRNA record Absent')
                return(NA)
            })
            
            if (nrow(a) == 0) {
                stop("miRNA does not have any records")
            } else {
                commgenes <- as.character(a$Target_GeneSymbol)
                FBGN <- as.character(a$Targets_FBID)
                CGID <- as.character(a$Targets_CGID)
                sc <- as.numeric(a$Score)
                ord <- a[order(-sc), ]
                dat1 <- ord[1:thresh, ]
                dat <- rbind(dat, dat1)
            }
        }
        # If Visualisation is selected as Text, writes in a .csv
        # file Else If Visualisation is selected as Cytoscape, writes in a .csv
        # file in cytoscape format output Else If Visualisation is selected as
        # igraph, outputs in the format the layout is selected Else If
        # Visualisation is selected as console returns a data frame to the
        # console Else stops with error message
        if (identical(visualisation, "Text")) {
            ### If miRNA length >1 names concatenated with & else use the miR 
            ###name
            if (length(miRNA) > 1) {
                miRNAs <- paste(miRNA, collapse = "&")
            } else {
                miRNAs = miRNA
            }
            filename = paste("Putative_target_Score_for_", miRNAs, "_", 
                platform, "_and_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(path, filename)))) {
                write.csv(dat, file = file.path(path, filename))
            } else {
                stop("File Already Exists!!")
            }
            
        } else if (identical(visualisation, "Cytoscape")) {
            miR <- as.character(dat$miRNA)
            ## Checks if the mRNA_Type selected to print in the output table is
            ## GeneSymbol, Flybase Id (FBID) or CGID
            if (identical(mRNA_type, "GeneSymbol")) {
                gene <- as.character(dat$Target_GeneSymbol)
            } else if (identical(mRNA_type, "FBID")) {
                gene <- as.character(dat$Targets_FBID)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Targets_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            ## Storing data in a dataframe
            datC <- data.frame(miRNA = miR, Interaction = rep("pp", 
                length(miR)), Gene = gene, Score)
            if (length(miRNA) > 1) {
                miRNAs <- paste(miRNA, collapse = "&")
            } else {
                miRNAs = miRNA
            }
            filename = paste("Cytoscape_output_for_", miRNAs, "_", platform, 
                "_and_", method, ".csv", sep = "")
            
            ## Checking if the file exists
            if (!(file.exists(file.path(path, filename)))) {
                write.csv(datC, file = file.path(path, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else if (identical(visualisation, "igraph")) {
            datN <- c()
            ### Extracting the data and storing in the data.frame
            miR <- as.character(dat$miRNA)
            ## Checks if the mRNA_Type selected to print in the output table is
            ## GeneSymbol, Flybase Id (FBID) or CGID
            if (identical(mRNA_type, "GeneSymbol")) {
                gene <- as.character(dat$Target_GeneSymbol)
            } else if (identical(mRNA_type, "FBID")) {
                gene <- as.character(dat$Targets_FBID)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Targets_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            datN <- data.frame(miRNA = miR, Gene = gene, Score)
            datN <- na.omit(datN)
            ### Creating the plot
            a <- graph.data.frame(datN, directed = FALSE)
            ### Creating the Edge weight
            E(a)$weight <- Score
            ### Selecting the layout to represent the network. The selected 
            ### network is printed out in the format of a .tiff file. If chosen
            ### layout iskamadakawai, prints as .tiff file If chosen layout is
            ### fruchterman.reingold, prints as .tiff file If chosen layout is
            ### reingold.tilford, prints as .tiff file If chosen layout is
            ### interactive, prints as a tkplot file Else Prints error message
            if (identical(layout, "kamadakawai")) {
                layout_graph <- layout.kamada.kawai(a)
                if (length(miRNA) > 1) {
                    miRNAs <- paste(miRNA, collapse = "&")
                } else {
                    miRNAs = miRNA
                }
                filename = paste("Graphical_output_for_", miRNAs, "_", 
                    platform, "_and_", method, ".tiff", sep = "")
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
                if (length(miRNA) > 1) {
                    miRNAs <- paste(miRNA, collapse = "&")
                } else {
                    miRNAs = miRNA
                }
                filename = paste("Graphical_output_for_", miRNAs, "_", 
                    platform, "_and_", method, ".tiff", sep = "")
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
                if (length(miRNA) > 1) {
                    miRNAs <- paste(miRNA, collapse = "&")
                } else {
                  miRNAs = miRNA
                }
                filename = paste("Graphical_output_for_", miRNAs, "_", 
                    platform, "_and_", method, ".tiff", sep = "")
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
                print("Invalid Layout!!!")
            }
        } else if (identical(visualisation, "console")) {
            
            return(dat)
        } else {
            stop("Visualisation Invalid")
        }
        
    } else {
        ## miRNA length ==1
        stopifnot(length(extract_HostGene(miRNA)) > 0)
        # dat <- c() Extracts and sorts the data in descending order
        a <- tryCatch({
            miRTargets_Stat(miRNA, method, Platform = platform)
            
        }, error = function(e) {
            # print('mRNA record Absent')
            return(NA)
        })
        if (nrow(a) == 0) {
            stop("miRNA does not have any records")
        } else {
            commgenes <- as.character
            (a$Target_GeneSymbol)
            FBGN <- as.character(a$Targets_FBID)
            CGID <- as.character(a$Targets_CGID)
            sc <- as.numeric(a$Score)
            ord <- a[order(-sc), ]
            dat1 <- ord[1:thresh, ]
            dat <- dat1
        }
        ### If Visualisation is selected as Text, writes in a .csv file Else If
        ### Visualisation is selected as Cytoscape, writes in a .csv file in
        ### cytoscape format output Else If Visualisation is selected as igraph,
        ### outputs in the format the layout is selected Else If Visualisation 
        ### is selected as console returns a data frame to the console Else stops
        ### with error message
        if (identical(visualisation, "Text")) {
            
            miRNAs = miRNA
            
            filename = paste("Putative_target_Score_for_", miRNAs, "_", 
                platform, "_and_", method, ".csv", sep = "")
            ## Checks whether the filename exists; prints if it does not previously
            ## exist
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
                gene <- as.character(dat$Target_GeneSymbol)
            } else if (identical(mRNA_type, "FBID")) {
                gene <- as.character(dat$Targets_FBID)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Targets_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            ## Storing data in a dataframe
            datC <- data.frame(miRNA = miR, Interaction = rep("pp", 
                    length(miR)), Gene = gene, Score)
            miRNAs = miRNA
            ## Checks whether the filename exists; prints if it does not
            ##  previously exist
            filename = paste("Cytoscape_output_for_", miRNAs, "_", platform, 
                "_and_", method, ".csv", sep = "")
            if (!(file.exists(file.path(path, filename)))) {
                write.csv(datC, file = file.path(path, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else if (identical(visualisation, "igraph")) {
            datN <- c()
            ## Extracting individual data and storing it in the database.
            miR <- as.character(dat$miRNA)
            ## Checks if the mRNA_Type selected to print in the output table is
            ## GeneSymbol, Flybase Id (FBID) or CGID
            if (identical(mRNA_type, "GeneSymbol")) {
                gene <- as.character(dat$Target_GeneSymbol)
            } else if (identical(mRNA_type, "FBID")) {
                gene <- as.character(dat$Targets_FBID)
            } else if (identical(mRNA_type, "CGID")) {
                gene <- as.character(dat$Targets_CGID)
            } else {
                stop("mRNA_Type Invalid!!!")
            }
            Score <- as.numeric(dat$Score)
            ## Storing data in the dataframe
            datN <- data.frame(miRNA = miR, Gene = gene, Score)
            datN <- na.omit(datN)
            ## Creating the network
            a <- graph.data.frame(datN, directed = FALSE)
            ### Score assigned as edge weight
            E(a)$weight <- Score
            ### Selecting the layout to represent the network. The selected 
            ### network is printed out in the format of a .tiff file. If chosen
            ### layout iskamadakawai, prints as .tiff file If chosen layout is
            ### fruchterman.reingold, prints as .tiff file If chosen layout is
            ### reingold.tilford, prints as .tiff file If chosen layout is
            ### interactive, prints as a tkplot file Else Prints error message
            if (identical(layout, "kamadakawai")) {
                layout_graph <- layout.kamada.kawai(a)
                miRNAs = miRNA
                
                filename = paste("Graphical_output_for_", miRNAs, "_", 
                    platform, "_and_", method, ".tiff", sep = "")
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
                miRNAs = miRNA
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
                miRNAs = miRNA
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
            stop("Visualisation Invalid")
        }
        
    }
}
