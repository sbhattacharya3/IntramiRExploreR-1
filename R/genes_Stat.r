#' Extracting miRNAs that target a query gene.
#'
#' @param gene  character. gene Identifier.
#' @param geneIDType  character. GeneIDtype choices are 'GeneSymbol', 
#'                    'FBGN', 'CGID'                     
#' @param method  character. Choices are 'Pearson','Distance','Both' and
#'                 'BothIntersected'
#' @param Platform  character. Choices are
#'                'Affy1','Affy2'.
#' @param Text  logical . To choose between storing the data as text file.
#'              Default is FALSE.
#' @param outpath  character. The path where the data is stored if TEXT=TRUE.
#'                 Default is 
#' @return Outputs the miRNA information, Target Prediction Score, miRNA 
#'         miRNA function and Target Database that predicts the interaction in
#'         a dataframe.
#'         Depending upon the ouput choice data is stored in the
#'         path specified. Default option prints output to the console. 
#' @examples
#' gene="Syb"
#' genes_Stat(gene,geneIDType="GeneSymbol",method=c("Pearson"),
#'            Platform=c("Affy1"),Text=FALSE)
#' @import utils grDevices graphics 
#' @importFrom stats na.omit
#' @export

genes_Stat <- function(gene, geneIDType = c("GeneSymbol", "FBGN", "CGID"), 
    method = c("Pearson", "Distance", "Both", "BothIntersect"), Platform = c("Affy1", 
        "Affy2"), Text = FALSE, outpath = tempdir()) {
    ## Initialising and loading the databases
    Affy2_Distance_Final <- c()
    Affy1_Distance_Final <- c()
    Affy2_Pearson_Final <- c()
    Affy1_Pearson_Final <- c()
    data("Affy1_Pearson_Final", envir = environment())
    Affy1_Pearson_Final <- Affy1_Pearson_Final
    data("Affy1_Distance_Final", envir = environment())
    Affy1_Distance_Final <- Affy1_Distance_Final
    data("Affy2_Pearson_Final", envir = environment())
    Affy2_Pearson_Final <- Affy2_Pearson_Final
    data("Affy2_Distance_Final", envir = environment())
    Affy2_Distance_Final <- Affy2_Distance_Final
    data("miRNA_ID_to_Function", envir = environment())
    miRNA_ID_to_Function <- miRNA_ID_to_Function
    ## Checking for the input variables
    stopifnot(is.character(gene), length(gene) > 0)
    stopifnot(is.character(geneIDType), length(geneIDType) > 0)
    stopifnot(is.character(method), length(method) > 0, method %in% c("Pearson", 
        "Distance", "Both", "BothIntersect"))
    stopifnot(is.character(Platform), length(Platform) > 0, Platform %in% 
        c("Affy1", "Affy2"))
    ## Checks the method, the platform and geneIDType and extracts the data
    ## for each individual query gene.  First Checks whether the method is
    ## Pearson and platform is Affy1 Else Checks If the method is is Pearson
    ## and platform is Affy2 Else Checks If the method is is Distance and
    ## platform is Affy1 Else Checks If the method is is Distance and
    ## platform is Affy2 Else Checks If the method is is Both and platform
    ## is Affy1 Else Checks If the method is is Both and platform is Affy2
    ## Else Checks If the method is is BothIntersect and platform is Affy2
    ## Else Checks If the method is is BothIntersect and platform is Affy2
    ## Else prints the error message that Method/platform is invalid
    if (identical(method, "Pearson") & identical(Platform, "Affy1")) {
        ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
        ## prints error message
        if (identical(geneIDType, "GeneSymbol")) {
            ### Extracts data when geneIDType is gene symbol and Platform is AFFY1
            dat <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$GeneSymbol) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 2], dat[, 1], dat[, 3:5], 
                  dat[, 7:8])
                miR <- finaldat1[, 2]
                mfunct <- c()
                # Extracts the miRNA function
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct <- c(mfunct,
                        as.character 
                        (miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ### Storing extracted data in a dataframe
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                  nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                  "Score", "Experiments", "TargetDatabase", "miRNAFunction", 
                  "Method")
            }
            
        } else if (identical(geneIDType, "FBGN")) {
            ### Extracts data when geneIDType is Flybase Gene ID and Platform is
            ### AFFY1
            dat <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$FBGN) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 3], dat[, 1], dat[, 2], dat[, 
                    4:5], dat[, 7:8])
                miR <- finaldat1[, 2]
                # Extracts the miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Storing the extracted data in a data frame
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("FlybaseID", "miRNA", "GeneSymbol", 
                    "Genes_CGID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
            
        } else if (identical(geneIDType, "CGID")) {
            ### Extracts data when geneIDType is Flybase Gene ID and Platform is
            ### AFFY1
            dat <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$CGID)
                %in% tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 4], dat[, 1], dat[, 2:3], 
                dat[, 5], dat[, 7:8])
                miR <- finaldat1[, 2]
                ## Extracting miRNAfunction
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Storing the extracted data in a data frame
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene_CGID", "miRNA", "GeneSymbol", 
                    "Genes_FBID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
            
        } else {
            stop("GeneIdtype Invalid")
        }
        ## If text is selected file written to outpath
        if (identical(Text, TRUE)) {
            if (length(gene) > 1) {
                Genes <- paste(gene, collapse = "&")
            } else {
                Genes = gene
            }
            ### Checks whether the file exists; write it in the path specified if
            ### file does not exist.
            filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                "_", method, ".csv", sep = "")
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
            
            
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Pearson") & identical(Platform, "Affy2")) {
        ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
        ## prints error message
        if (identical(geneIDType, "GeneSymbol")) {
            dat <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$GeneSymbol) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 2], dat[, 1], dat[, 3:5], 
                dat[, 7:8])
                miR <- finaldat1[, 2]
                ### Extracts miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Storing the data in data.frame
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                    "Score", "Experiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
            }
            
        } else if (identical(geneIDType, "FBGN")) {
            ### Extracts data when geneIDType is FBGN and Platform is AFFY2
            dat <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$FBGN) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 3], dat[, 1], dat[, 2], 
                dat[,4:5], dat[, 7:8])
                miR <- finaldat1[, 2]
                ### Extracts miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                  pa <- paste("^", miR[l], "$", sep = "")
                  dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Stores data in a data.frame
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("FlybaseID", "miRNA", "GeneSymbol", 
                    "Genes_CGID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
            
        } else if (identical(geneIDType, "CGID")) {
            
            dat <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$CGID) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 4], dat[, 1], dat[, 2:3], 
                dat[, 5], dat[, 7:8])
                miR <- finaldat1[, 2]
                ## Extracts miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                  pa <- paste("^", miR[l], "$", sep = "")
                  dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                    as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Stores data in a data.frame
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene_CGID", "miRNA", "GeneSymbol", 
                    "Genes_FBID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
        } else {
            stop("GeneIdtype Invalid")
        }
        ## If text is chosen as TRUE
        if (identical(Text, TRUE)) {
            if (length(gene) > 1) {
                Genes <- paste(gene, collapse = "&")
            } else {
                Genes = gene
            }
            filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                "_", method, ".csv", sep = "")
            ## Checks files exists or not; if it does not exist write the .csv
            ## files.
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
            
            
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Distance") & identical(Platform, "Affy1")) {
        if (identical(geneIDType, "GeneSymbol")) {
            ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
            ## prints error message
            dat <- Affy1_Distance_Final[which(tolower
                (Affy1_Distance_Final$GeneSymbol) %in% tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 2], dat[, 1], dat[, 3:5], 
                    dat[, 7:8])
                miR <- finaldat1[, 2]
                ### Extracts miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ### Stores data in a data.frame
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                     nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene", "miRNA", "Gene_FBID", 
                    "Genes_CGID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
        } else if (identical(geneIDType, "FBGN")) {
            ### Extracts data when geneIDType is Flybase ID
            dat <- Affy1_Distance_Final[which(tolower(Affy1_Distance_Final$FBGN) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 3], dat[, 1], dat[, 2], 
                dat[,4:5], dat[, 7:8])
                miR <- finaldat1[, 2]
                ## Extracting miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Storing Data in data.frame
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("FlybaseID", "miRNA", "GeneSymbol", 
                    "Genes_CGID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
        } else if (identical(geneIDType, "CGID")) {
            ### Extracts data when geneIDType is CGID
            dat <- Affy1_Distance_Final[which(tolower(Affy1_Distance_Final$CGID) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 4], dat[, 1], dat[, 2:3], 
                  dat[, 5], dat[, 7:8])
                miR <- finaldat1[, 2]
                # Extracting the miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                    as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Storing data
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene_CGID", "miRNA", "GeneSymbol", 
                    "Genes_FBID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
        } else {
            stop("GeneIdtype Invalid")
        }
        if (identical(Text, TRUE)) {
            ## If text is chosen as TRUE; if not returns value to the console Checks
            ## if gene length>1; if yes adds & after each gene name
            if (length(gene) > 1) {
                Genes <- paste(gene, collapse = "&")
            } else {
                Genes = gene
            }
            filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                "_", method, ".csv", sep = "")
            
            if (!(file.exists(file.path(outpath, filename)))) {
                ## Checks if file exists; if not write the data in .csv
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
            
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Distance") & identical(Platform, "Affy2")) {
        ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
        ## prints error message
        if (identical(geneIDType, "GeneSymbol")) {
            dat <- Affy2_Distance_Final[which
                    (tolower(Affy2_Distance_Final$GeneSymbol)
                    %in% tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 2], dat[, 1], dat[, 3:5], 
                dat[, 7:8])
                miR <- finaldat1[, 2]
                ## Extracting microRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                  pa <- paste("^", miR[l], "$", sep = "")
                  dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ### Storing Data
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                  nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene", "miRNA", "Gene_FBID", 
                    "Genes_CGID", "Score", "Experiments", "TargetDatabase",
                    "miRNAFunction", "Method")
            }
        } else if (identical(geneIDType, "FBGN")) {
            ### Extracts data when geneIDType is FLybase ID
            dat <- Affy2_Distance_Final[which
                    (tolower(Affy2_Distance_Final$FBGN) %in% tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 3], dat[, 1], dat[, 2], 
                    dat[, 4:5], dat[, 7:8])
                miR <- finaldat1[, 2]
                ## Extracting miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                  pa <- paste("^", miR[l], "$", sep = "")
                  dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                    ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                        mfunct <- c(mfunct, "NA")
                  }
                }
                ## Storing data in the data.frame.
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("FlybaseID", "miRNA", "GeneSymbol", 
                    "Genes_CGID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
        } else if (identical(geneIDType, "CGID")) {
            ### Extracts data when geneIDType is CGID
            dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$CGID) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(dat[, 4], dat[, 1], dat[, 2:3], 
                    dat[, 5], dat[, 7:8])
                miR <- finaldat1[, 2]
                ## Extracting miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                            ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                ## Storing data in the data.frame.
                finaldat1 <- cbind(finaldat1, mfunct, method = rep("Pearson", 
                    nrow(finaldat1)), row.names = NULL)
                names(finaldat1) <- c("Gene_CGID", "miRNA", "GeneSymbol", 
                    "Genes_FBID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
        } else {
            stop("GeneIdtype Invalid")
        }
        if (identical(Text, TRUE)) {
            ## If text is chosen as TRUE prints data in a .csv file else return the
            ## result to the console
            if (length(gene) > 1) {
                Genes <- paste(gene, collapse = "&")
            } else {
                Genes = gene
            }
            filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                "_", method, ".csv", sep = "")
            ## Check if file name exists. If not write it in the .csv file
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
            
            
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Both") & identical(Platform, "Affy1")) {
        ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
        ## prints error message
        if (identical(geneIDType, "GeneSymbol")) {
            # Extract data from Pearson and Distance respectively
            dat <- Affy1_Distance_Final[which(tolower(Affy1_Distance_Final$GeneSymbol) %in% 
                tolower(gene)), ]
            dat1 <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$GeneSymbol) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0 | nrow(dat1) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(Gene = dat[, 2], miRNA = dat[, 1],
                    dat[, 3:5], dat[, 7:8])
                finaldat2 <- data.frame(Gene = dat1[, 2], miRNA = dat1[, 1],
                    dat1[, 3:5], dat1[, 7:8])
                ## Extract miRNA function for Pearson and Distance respectively.
                miR <- finaldat1[, 2]
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                  if (length(dd) == 1) {
                        mfunct <- c(mfunct, 
                            as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                miR <- finaldat2[, 2]
                mfunct1 <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct1 <- c(mfunct1, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct1 <- c(mfunct1, "NA")
                  }
                }
                ## Store data in data.frame
                finaldat2 <- cbind(finaldat2, mfunct = mfunct1, 
                    method = rep("Pearson", nrow(finaldat2)), row.names = NULL)
                finaldat1 <- cbind(finaldat1, mfunct = mfunct, method = rep("Distance", 
                    nrow(finaldat1)), row.names = NULL)
                finaldat <- rbind(finaldat1, finaldat2)
                    names(finaldat) <- c("Gene", "miRNA", "Gene_FBID",
                    "Genes_CGID", "Score", "Experiments", "TargetDatabase", 
                    "miRNAFunction", "Method")
            }
        } else if (identical(geneIDType, "FBGN")) {
            ### Extracts data when geneIDType is FBGN Extract data from Pearson and
            ### Distance respectively
            dat <- Affy1_Distance_Final[which(tolower(Affy1_Distance_Final$FBGN)
                    %in% tolower(gene)), ]
            dat1 <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$FBGN)
                    %in% tolower(gene)), ]
            finaldat1 <- data.frame(Gene = dat[, 2], miRNA = dat[, 1], 
                    dat[, 3:5], dat[, 7:8])
            if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    stop("Records of the gene does not exist")
            } else {
                finaldat2 <- data.frame(Gene = dat1[, 2], miRNA = dat1[, 1],
                    dat1[, 3:5], dat1[, 7:8])
                miR <- finaldat1[, 2]
                ## miRNA function extraction
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                if (length(dd) == 1) {
                    mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                } else {
                    mfunct <- c(mfunct, "NA")
                }
                }
                miR <- finaldat2[, 2]
                mfunct1 <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                  if (length(dd) == 1) {
                        mfunct1 <- c(mfunct1, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                        mfunct1 <- c(mfunct1, "NA")
                    }
                }
                ### Storing data in form of data.frame
                finaldat2 <- cbind(finaldat2, mfunct = mfunct1,  
                        method = rep("Pearson",nrow(finaldat2)), 
                        row.names = NULL)
                finaldat1 <- cbind(finaldat1, mfunct = mfunct,  
                    method = rep("Distance", nrow(finaldat1)), 
                        row.names = NULL)
                finaldat <- rbind(finaldat1, finaldat2)
                names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                    "Score", "Experiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
            }
        } else if (identical(geneIDType, "CGID")) {
            ### Extracts data when geneIDType is CGID Extract data from Pearson and
            ### Distance respectively
            dat <- Affy1_Distance_Final[which(tolower(Affy1_Distance_Final$CGID)
                        %in% tolower(gene)), ]
            dat1 <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$CGID)
                        %in% tolower(gene)), ]
            if (nrow(dat) == 0 | nrow(dat1) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(Gene = dat[, 2], miRNA = dat[,1], 
                    dat[, 3:5], dat[, 7:8])
                finaldat2 <- data.frame(Gene = dat1[, 2], miRNA = dat1[, 1],
                    dat1[, 3:5], dat1[, 7:8])
                ## miRNA function extraction
                miR <- finaldat1[, 2]
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                        mfunct <- c(mfunct, "NA")
                  }
                }
                miR <- finaldat2[, 2]
                mfunct1 <- c()
                for (l in seq_along(miR)) {
                  pa <- paste("^", miR[l], "$", sep = "")
                  dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                            ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct1 <- c(mfunct1, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                        mfunct1 <- c(mfunct1, "NA")
                  }
                }
                ## Storing data in data.frame
                finaldat2 <- cbind(finaldat2, mfunct = mfunct1, 
                  method = rep("Pearson", nrow(finaldat2)), row.names = NULL)
                finaldat1 <- cbind(finaldat1, mfunct = mfunct,  
                  method = rep("Distance",nrow(finaldat1)), row.names = NULL)
                finaldat <- rbind(finaldat1, finaldat2)
                names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                  "Score", "Experiments", "TargetDatabase", "miRNAFunction", 
                  "Method")
            }
        } else {
            stop("GeneIdtype Invalid")
        }
        ## If text is selected as TRUE saved as .csv; else return to the console
        ## If length of gene is greater than 1
        if (identical(Text, TRUE)) {
            
            if (length(gene) > 1) {
                Genes <- paste(gene, collapse = "&")
            } else {
                Genes = gene
            }
            filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                "_", method, ".csv", sep = "")
            
            
            if (!(file.exists(file.path(outpath, filename)))) {
                ## Check if file name exists and print data.frame in .csv
                write.csv(finaldat, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
            
            
        } else {
            return(finaldat)
        }
    } else if (identical(method, "Both") & identical(Platform, "Affy2")) {
        ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
        ## prints error message
        if (identical(geneIDType, "GeneSymbol")) {
            # Extract data from Pearson and Distance respectively
            dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$GeneSymbol) %in% 
                tolower(gene)), ]
            dat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$GeneSymbol) %in% 
                tolower(gene)), ]
            if (nrow(dat) == 0 | nrow(dat1) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(Gene = dat[, 2], miRNA = dat[, 1],
                    dat[, 3:5], dat[, 7:8])
                finaldat2 <- data.frame(Gene = dat1[, 2], miRNA = dat1[, 1],
                    dat1[, 3:5], dat1[, 7:8])
                miR <- finaldat1[, 2]
                ## Extract miRNA function for mPearson and Distance respectively.
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                  if (length(dd) == 1) {
                        mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                } else {
                    mfunct <- c(mfunct, "NA")
                  }
                }
                miR <- finaldat2[, 2]
                mfunct1 <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct1 <- c(mfunct1, 
                            as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                } else {
                    mfunct1 <- c(mfunct1, "NA")
                }
                }
                ## Storing data in data.frame
                finaldat2 <- cbind(finaldat2, mfunct = mfunct1,  
                        method = rep("Pearson", nrow(finaldat2)), 
                        row.names = NULL)
                finaldat1 <- cbind(finaldat1, mfunct = mfunct, 
                        method = rep("Distance", nrow(finaldat1)), 
                        row.names = NULL)
                finaldat <- rbind(finaldat1, finaldat2)
                names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                    "Score", "Experiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
            }
        } else if (identical(geneIDType, "FBGN")) {
            ### Extracts data when geneIDType is gene symbol
            dat <- Affy2_Distance_Final[which
                    (tolower(Affy2_Distance_Final$FBGN) %in% tolower(gene)), ]
            dat1 <- Affy2_Pearson_Final[which
                    (tolower(Affy2_Pearson_Final$FBGN) %in% 
                    tolower(gene)), ]
            if (nrow(dat) == 0 | nrow(dat1) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(Gene = dat[, 2], miRNA = dat[, 1],
                    dat[, 3:5], dat[, 7:8])
                finaldat2 <- data.frame(Gene = dat1[, 2], miRNA = dat1[,1],
                    dat1[, 3:5], dat1[, 7:8])
                miR <- finaldat1[, 2]
                ## Extracts miRNA data
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                if (length(dd) == 1) {
                    mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                } else {
                    mfunct <- c(mfunct, "NA")
                }
                }
                miR <- finaldat2[, 2]
                mfunct1 <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                if (length(dd) == 1) {
                    mfunct1 <- c(mfunct1, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                } else {
                    mfunct1 <- c(mfunct1, "NA")
                }
                }
                ### Storing data in data.frame
                finaldat2 <- cbind(finaldat2, mfunct = mfunct1,  
                    method = rep("Pearson", nrow(finaldat2)), row.names = NULL)
                finaldat1 <- cbind(finaldat1, mfunct = mfunct,  
                    method = rep("Distance",nrow(finaldat1)), row.names = NULL)
                finaldat <- rbind(finaldat1, finaldat2)
                names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                    "Score", "Experiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
            }
        } else if (identical(geneIDType, "CGID")) {
            ### Extracts data when geneIDType is CGID Extracts data for Pearson and
            ### Distance respectively.
            dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$CGID)
                %in% tolower(gene)), ]
            dat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$CGID)
                %in% tolower(gene)), ]
            if (nrow(dat) == 0 | nrow(dat1) == 0) {
                stop("Records of the gene does not exist")
            } else {
                finaldat1 <- data.frame(Gene = dat[, 2], miRNA = dat[, 1], 
                    dat[, 3:5], dat[, 7:8])
                finaldat2 <- data.frame(Gene = dat1[, 2], miRNA = dat1[, 1],
                    dat1[, 3:5], dat1[, 7:8])
                miR <- finaldat1[, 2]
                ### Extract data for miRNA function
                mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                if (length(dd) == 1) {
                    mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                } else {
                    mfunct <- c(mfunct, "NA")
                }
                }
                miR <- finaldat2[, 2]
                mfunct1 <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                  if (length(dd) == 1) {
                    mfunct1 <- c(mfunct1, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                  } else {
                    mfunct1 <- c(mfunct1, "NA")
                  }
                }
                ### Storing data in data.frame
                finaldat2 <- cbind(finaldat2, mfunct = mfunct1,  
                        method = rep("Pearson",nrow(finaldat2)), 
                        row.names = NULL)
                finaldat1 <- cbind(finaldat1, mfunct = mfunct,  
                        method = rep("Distance", nrow(finaldat1)), 
                        row.names = NULL)
                finaldat <- rbind(finaldat1, finaldat2)
                names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", 
                        "Genes_CGID", "Score", "Experiments",
                        "TargetDatabase", "miRNAFunction", 
                        "Method")
            }
        } else {
            stop("GeneIdtype Invalid")
        }
        if (identical(Text, TRUE)) {
            
            if (length(gene) > 1) {
                ### If text selected as TRUE; print data in .csv file else 
                ### print it on the console Checks if gene length>1; if yes
                ### adds & after each gene name
                Genes <- paste(gene, collapse = "&")
            } else {
                Genes = gene
            }
            filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Checking if File name exists and write data into .csv files
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
            
            
        } else {
            return(finaldat)
        }
        
    } else if (identical(method, "BothIntersect") & 
                identical(Platform, "Affy1")) {
        
        if (length(gene) > 1) {
            ## Checks whether the method is BothIntersect and platform 
            ## is Affy1 If length > 1
            finaldat <- data.frame()
            ## Checks if the chosen geneIDType is GeneSymbol, 
            ## FBGN or CGID elseprints error message
            if (identical(geneIDType, "GeneSymbol")) {
                for (ii in seq_along(gene)) {
                    dat <- Affy1_Distance_Final[which
                        (tolower(Affy1_Distance_Final$GeneSymbol) %in% 
                        tolower(gene[ii])), ]
                    dat1 <- Affy1_Pearson_Final[which
                        (tolower(Affy1_Pearson_Final$GeneSymbol) %in% 
                        tolower(gene[ii])), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    stop("Records of the gene does not exist")
                } else {
                    int <- intersect(as.character(dat$miRNA), 
                        as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), ]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) %in% 
                        tolower(int)), ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 1],
                        fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 1], 
                        fdat1[, 3:5], fdat1[, 7:8])
                    ## Extracting the miRNA function
                    miR <- int
                    mfunct <- c()
                    for (l in seq_along(miR)) {
                        pa <- paste("^", miR[l], "$", sep = "")
                        dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                            ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                        mfunct <- c(mfunct, "NA")
                    }
                    }
                    finaldat1 <- cbind(Gene = dat_P[, 2], miRNA = dat_P[, 1], 
                        dat_P[, 3:5], DistanceScore = dat_D[, 5], 
                        Pexp = dat_P[, 6], Dexp = dat_D[, 6], 
                        TargetDB = dat_P[, 7], Genefunc = mfunct)
                    finaldat <- rbind(finaldat, finaldat1)
                    }
                }
                ## Storing the data in the data.frame
                finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                        nrow(finaldat)))
                names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                    "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                    "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
            } else if (identical(geneIDType, "FBGN")) {
                ### If geneIDtype selected as FBGN Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                for (ii in seq_along(gene)) {
                    dat <- Affy1_Distance_Final[which(tolower(Affy1_Distance_Final$FBGN) %in% 
                        tolower(gene[ii])), ]
                    dat1 <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$FBGN) %in% 
                        tolower(gene[ii])), ]
                    if (nrow(dat) == 0 | nrow(dat1) == 0) {
                        stop("Records of the gene does not exist")
                    } else {
                    int <- intersect(as.character(dat$miRNA), 
                        as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in%  
                        tolower(int)),]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) %in%  
                        tolower(int)),]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 1],
                        fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 1],
                        fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                    ### miRNA function extraction
                    mfunct <- c()
                    for (l in seq_along(miR)) {
                        pa <- paste("^", miR[l], "$", sep = "")
                        dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                            ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                        mfunct <- c(mfunct, "NA")
                        }
                    }
                    ### Storing data in data.frame
                    finaldat1 <- cbind(Gene = dat_P[, 3], miRNA = dat_P[,1], 
                        GS = dat_P[, 2], dat_P[, 4:5],  
                        DistanceScore = dat_D[,5], Pexp = dat_P[, 6], 
                        Dexp = dat_D[, 6], TargetDB = dat_P[, 7], 
                        Genefunc = mfunct)
                    finaldat <- rbind(finaldat, finaldat1)
                  }
                }
                finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                    nrow(finaldat)))
                names(finaldat) <- c("FBID", "miRNA", "Gene_Symbol", "Genes_CGID", 
                    "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                    "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
            } else if (identical(geneIDType, "CGID")) {
                for (ii in seq_along(gene)) {
                    dat <- Affy1_Distance_Final[ 
                    which(tolower(Affy1_Distance_Final$CGID) %in%
                    tolower(gene[ii])), ]
                    dat1 <- Affy1_Pearson_Final[which
                    (tolower(Affy1_Pearson_Final$CGID) %in% 
                    tolower(gene[ii])), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    stop("Records of the gene does not exist")
                } else {
                    int <- intersect(as.character(dat$miRNA), 
                        as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA)
                        %in% tolower(int)),]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) 
                        %in% tolower(int)), ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 1], 
                        fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 1], 
                        fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                    ### miRNA function extraction
                    mfunct <- c()
                    for (l in seq_along(miR)) {
                      pa <- paste("^", miR[l], "$", sep = "")
                      dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                        mfunct <- c(mfunct, "NA")
                        }
                    }
                    ### Storing data in data.frame
                    finaldat1 <- cbind(Gene = dat_P[, 4], miRNA = dat_P[, 1],
                        dat_P[, 2:3], PearsonScore = dat_P[, 5],  
                        DistanceScore = dat_D[, 5], Genefunc = mfunct, 
                        Pexp = dat_P[, 6], Dexp = dat_D[, 6], 
                        TargetDB = dat_P[, 7])
                    finaldat <- rbind(finaldat, finaldat1)
                    }
                }
                finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                    nrow(finaldat)))
                names(finaldat) <- c("CGID", "miRNA", "Gene_Symbol", "Genes_FBID", 
                    "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                    "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
            } else {
                stop("GeneIdtype Invalid")
                }
        } else {
            ### If length of miRNA == 1 Checks if the chosen geneIDType is
            ### GeneSymbol, FBGN or CGID else prints error message
            
            if (identical(geneIDType, "GeneSymbol")) {
                ### If geneIDtype selected as GeneSymbol Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                dat <- Affy1_Distance_Final[
                    which(tolower(Affy1_Distance_Final$GeneSymbol) %in% 
                    tolower(gene)), ]
                dat1 <- Affy1_Pearson_Final[
                    which(tolower(Affy1_Pearson_Final$GeneSymbol) %in%
                    tolower(gene)), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    stop("Records of the gene does not exist")
                } else {
                    int <- intersect(as.character(dat$miRNA), 
                        as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), ]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) %in% tolower(int)), ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 1],
                        fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 1],
                        fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                    ## Extracting miRNA function
                    mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                      ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                      mfunct <- c(mfunct, "NA")
                        }
                }
                  
                  ### Storing data in data.frame
                    finaldat <- cbind(Gene = dat_P[, 2], miRNA = dat_P[,1], 
                    dat_P[, 3:5], DistanceScore = dat_D[, 5], 
                    Genefunc = mfunct, Pexp = dat_P[, 6], Dexp = dat_D[, 6], 
                    TargetDB = dat_P[, 7])
                    finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                    nrow(finaldat)))
                    names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                    "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                    "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
                }
            } else if (identical(geneIDType, "FBGN")) {
                ### If geneIDtype selected as FBGN Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                dat <- Affy1_Distance_Final[which(tolower
                    (Affy1_Distance_Final$FBGN) %in% 
                    tolower(gene)), ]
                dat1 <- Affy1_Pearson_Final[which
                    (tolower(Affy1_Pearson_Final$FBGN) %in% 
                    tolower(gene)), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    stop("Records of the gene does not exist")
                } else {
                    int <- intersect(as.character(dat$miRNA),
                        as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), ]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) %in% 
                    tolower(int)), ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 1],
                    fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 1],
                    fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                  ## miRNA extraction for function
                  mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                      mfunct <- c(mfunct, "NA")
                    }
                }
                  ## Storing final data
                finaldat <- cbind(Gene = dat_P[, 3], miRNA = dat_P[, 1],
                GS = dat_P[, 2], dat_P[, 4:5], DistanceScore = dat_D[, 5], 
                Genefunc = mfunct, Pexp = dat_P[, 6], Dexp = dat_D[, 6],
                TargetDB = dat_P[, 7])
                finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                nrow(finaldat)))
                names(finaldat) <- c("FBID", "miRNA", "Gene_Symbol", 
                "Genes_CGID", "Pearson_Score", "Distance_Score", 
                "PearsonExperiments", "DistanceExperiments", "TargetDatabase",
                "miRNAFunction", "Method")
                }
            } else if (identical(geneIDType, "CGID")) {
                ### If geneIDtype selected as CGID Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                dat <- Affy1_Distance_Final[which
                    (tolower(Affy1_Distance_Final$CGID) %in% tolower(gene)), ]
                dat1 <- Affy1_Pearson_Final[
                    which(tolower(Affy1_Pearson_Final$CGID) 
                    %in% tolower(gene)), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    stop("Records of the gene does not exist")
                } else {
                    int <- intersect(as.character(dat$miRNA),
                        as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), ]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) %in% 
                    tolower(int)), ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 1], 
                        fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 1],
                    fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                    ## Extracting miRNA function
                    mfunct <- c()
                for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                    if (length(dd) == 1) {
                        mfunct <- c(mfunct, 
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                        mfunct <- c(mfunct, "NA")
                    }
                }
                    ### Storing data in data.frame
                    finaldat <- cbind(Gene = dat_P[, 4], miRNA = dat_P[, 1], 
                    dat_P[, 2:3], PearsonScore = dat_P[, 5],  
                    DistanceScore = dat_D[,5], Genefunc = mfunct, 
                    Pexp = dat_P[, 6], Dexp = dat_D[, 6],
                    TargetDB = dat_P[, 7], Genefunc = mfunct)
                    finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                    nrow(finaldat)))
                    names(finaldat) <- c("CGID", "miRNA", "Gene_Symbol", 
                    "Genes_FBID", "Pearson_Score", "Distance_Score", 
                    "PearsonExperiments", "DistanceExperiments", 
                    "TargetDatabase", "miRNAFunction", "Method")
                  
                }
            } else {
                stop("GeneIdtype Invalid")
            }
            ### If text selected as TRUE; print data in .csv file else print it on
            ### the console Checks if gene length>1; if yes adds & after each gene
            ### name
            if (identical(Text, TRUE)) {
                
                if (length(gene) > 1) {
                    Genes <- paste(gene, collapse = "&")
                } else {
                    Genes = gene
                }
                filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                    "_", method, ".csv", sep = "")
                ### Checking if File name exists and write data into .csv files
                if (!(file.exists(file.path(outpath, filename)))) {
                    write.csv(finaldat, file = file.path(outpath, filename))
                } else {
                    stop("File Already Exists!!")
                }
                
                
            } else {
                return(finaldat)
            }
        }
    } else if (identical(method, "BothIntersect") & identical(Platform, "Affy2")) {
        
        ## length of gene > 1
        if (length(gene) > 1) {
            finaldat <- data.frame()
            ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
            ## prints error message
            if (identical(geneIDType, "GeneSymbol")) {
                ## Extracting data from Pearson and Distance and finding the common
                ## microRNAs between them and extracting data individually for the
                ## intersected miRNAs.
                for (ii in seq_along(gene)) {
                    dat <- Affy2_Distance_Final[which
                        (tolower(Affy2_Distance_Final$GeneSymbol) %in% 
                        tolower(gene[ii])), ]
                    dat1 <- Affy2_Pearson_Final[which
                        (tolower(Affy2_Pearson_Final$GeneSymbol) %in% 
                        tolower(gene[ii])), ]
                    if (nrow(dat) == 0 | nrow(dat1) == 0) {
                        stop("Records of the gene does not exist")
                } else {
                    int <- intersect(as.character(dat$miRNA), 
                        as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), ]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) 
                        %in% tolower(int)), ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 1],
                        fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 1],
                        fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                    ### Extracting miRNA function
                    mfunct <- c()
                    for (l in seq_along(miR)) {
                        pa <- paste("^", miR[l], "$", sep = "")
                        dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                    if (length(dd) == 1) {
                         mfunct <- c(mfunct,
                        as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                        mfunct <- c(mfunct, "NA")
                      }
                    }
                    # Storing data in data.frame
                    finaldat1 <- cbind(Gene = dat_P[, 2], miRNA = dat_P[, 
                      1], dat_P[, 3:5], DistanceScore = dat_D[, 5], Genefunc = mfunct, 
                      Pexp = dat_P[, 6], Dexp = dat_D[, 6], TargetDB = dat_P[, 
                        7], Genefunc = mfunct)
                    finaldat <- rbind(finaldat, finaldat1)
                  }
                }
                
                finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                  nrow(finaldat)))
                names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                  "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                  "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                  "Method")
            } else if (identical(geneIDType, "FBGN")) {
                ### If geneIDtype selected as FBGN Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                for (ii in seq_along(gene)) {
                  dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$FBGN) %in% 
                    tolower(gene[ii])), ]
                  dat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$FBGN) %in% 
                    tolower(gene[ii])), ]
                  if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    print("Records of the gene does not exist")
                  } else {
                    int <- intersect(as.character(dat$miRNA), as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), 
                      ]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) %in% tolower(int)), 
                      ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 
                      1], fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 
                      1], fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                    ## Extracting miRNA function
                    mfunct <- c()
                    for (l in seq_along(miR)) {
                      pa <- paste("^", miR[l], "$", sep = "")
                      dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                      if (length(dd) == 1) {
                        mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                      } else {
                        mfunct <- c(mfunct, "NA")
                      }
                    }
                    ### Storing data in data.frame
                    finaldat1 <- cbind(Gene = dat_P[, 3], miRNA = dat_P[, 
                      1], GS = dat_P[, 2], dat_P[, 4:5], DistanceScore = dat_D[, 
                      5], Pexp = dat_P[, 6], Dexp = dat_D[, 6], TargetDB = dat_P[, 
                      7], Genefunc = mfunct)
                    finaldat <- rbind(finaldat, finaldat1)
                  }
                }
                finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                  nrow(finaldat)))
                names(finaldat) <- c("FBID", "miRNA", "Gene_Symbol", "Genes_CGID", 
                  "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                  "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                  "Method")
            } else if (identical(geneIDType, "CGID")) {
                ### If geneIDtype selected as CGID Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                for (ii in seq_along(gene)) {
                  dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$CGID) %in% 
                    tolower(gene[ii])), ]
                  dat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$CGID) %in% 
                    tolower(gene[ii])), ]
                  if (nrow(dat) == 0 | nrow(dat1) == 0) {
                    print("Records of the gene does not exist")
                  } else {
                    int <- intersect(as.character(dat$miRNA), as.character(dat1$miRNA))
                    fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), 
                      ]
                    fdat1 <- dat1[which(tolower(dat1$miRNA) %in% tolower(int)), 
                      ]
                    dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 
                      1], fdat[, 3:5], fdat[, 7:8])
                    dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 
                      1], fdat1[, 3:5], fdat1[, 7:8])
                    miR <- int
                    ### miRNA function extraction
                    mfunct <- c()
                    for (l in seq_along(miR)) {
                      pa <- paste("^", miR[l], "$", sep = "")
                      dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                        ignore.case = TRUE)
                      if (length(dd) == 1) {
                        mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                      } else {
                        mfunct <- c(mfunct, "NA")
                      }
                    }
                    ### Data stored in final data
                    finaldat1 <- cbind(Gene = dat_P[, 4], miRNA = dat_P[, 
                      1], dat_P[, 2:3], PearsonScore = dat_P[, 5], DistanceScore = dat_D[, 
                      5], Pexp = dat_P[, 6], Dexp = dat_D[, 6], TargetDB = dat_P[, 
                      7], Genefunc = mfunct)
                    finaldat <- rbind(finaldat, finaldat1)
                  }
                }
                finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                  nrow(finaldat)))
                names(finaldat) <- c("CGID", "miRNA", "Gene_Symbol", "Genes_FBID", 
                  "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                  "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                  "Method")
            } else {
                stop("File Already Exists!!")
            }
        } else {
            ## Checks if the chosen geneIDType is GeneSymbol, FBGN or CGID else
            ## prints error message
            if (identical(geneIDType, "GeneSymbol")) {
                ### If geneIDtype selected as GeneSymbol Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$GeneSymbol) %in% 
                  tolower(gene)), ]
                dat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$GeneSymbol) %in% 
                  tolower(gene)), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                  print("Records of the gene does not exist")
                } else {
                  int <- intersect(as.character(dat$miRNA), as.character(dat1$miRNA))
                  fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), 
                    ]
                  fdat1 <- dat1[which(tolower(dat1$miRNA) %in% tolower(int)), 
                    ]
                  dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 
                    1], fdat[, 3:5], fdat[, 7:8])
                  dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 
                    1], fdat1[, 3:5], fdat1[, 7:8])
                  miR <- int
                  ### Extracting miRNA function
                  mfunct <- c()
                  for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                      ignore.case = TRUE)
                    if (length(dd) == 1) {
                      mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                      mfunct <- c(mfunct, "NA")
                    }
                  }
                  ### Storing data in the data.frame
                  finaldat <- cbind(Gene = dat_P[, 2], miRNA = dat_P[, 
                    1], dat_P[, 3:5], DistanceScore = dat_D[, 5], Genefunc = mfunct, 
                    Pexp = dat_P[, 6], Dexp = dat_D[, 6], TargetDB = dat_P[, 
                      7], Genefunc = mfunct)
                  finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                    nrow(finaldat)))
                  names(finaldat) <- c("Gene", "miRNA", "Gene_FBID", "Genes_CGID", 
                    "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                    "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
                }
            } else if (identical(geneIDType, "FBGN")) {
                ### If geneIDtype selected as FBGN Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$FBGN) %in% 
                  tolower(gene)), ]
                dat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$FBGN) %in% 
                  tolower(gene)), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                  print("Records of the gene does not exist")
                } else {
                  int <- intersect(as.character(dat$miRNA), as.character(dat1$miRNA))
                  fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), 
                    ]
                  fdat1 <- dat1[which(tolower(dat1$miRNA) %in% tolower(int)), 
                    ]
                  dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 
                    1], fdat[, 3:5], fdat[, 7:8])
                  dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 
                    1], fdat1[, 3:5], fdat1[, 7:8])
                  miR <- int
                  ### Extracting microRNA function
                  mfunct <- c()
                  for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                      ignore.case = TRUE)
                    if (length(dd) == 1) {
                      mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                      mfunct <- c(mfunct, "NA")
                    }
                  }
                  ### Storing data in data.frame
                  finaldat <- cbind(Gene = dat_P[, 3], miRNA = dat_P[, 
                    1], GS = dat_P[, 2], dat_P[, 4:5], DistanceScore = dat_D[, 
                    5], Pexp = dat_P[, 6], Dexp = dat_D[, 6], TargetDB = dat_P[, 
                    7], Genefunc = mfunct)
                  finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                    nrow(finaldat)))
                  names(finaldat) <- c("FBID", "miRNA", "Gene_Symbol", 
                    "Genes_CGID", "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                    "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
                }
            } else if (identical(geneIDType, "CGID")) {
                ### If geneIDtype selected as CGID Extracting data from Pearson and
                ### Distance and finding the common microRNAs between them and extracting
                ### data individually for the intersected miRNAs.
                dat <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$CGID) %in% 
                  tolower(gene)), ]
                dat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$CGID) %in% 
                  tolower(gene)), ]
                if (nrow(dat) == 0 | nrow(dat1) == 0) {
                  print("Records of the gene does not exist")
                } else {
                  int <- intersect(as.character(dat$miRNA), as.character(dat1$miRNA))
                  fdat <- dat[which(tolower(dat$miRNA) %in% tolower(int)), 
                    ]
                  fdat1 <- dat1[which(tolower(dat1$miRNA) %in% tolower(int)), 
                    ]
                  dat_D <- data.frame(Gene = fdat[, 2], miRNA = fdat[, 
                    1], fdat[, 3:5], fdat[, 7:8])
                  dat_P <- data.frame(Gene = fdat1[, 2], miRNA = fdat1[, 
                    1], fdat1[, 3:5], fdat1[, 7:8])
                  miR <- int
                  ### Extracting miRNA functional data
                  mfunct <- c()
                  for (l in seq_along(miR)) {
                    pa <- paste("^", miR[l], "$", sep = "")
                    dd <- grep(pa, as.character(miRNA_ID_to_Function$miRNA), 
                      ignore.case = TRUE)
                    if (length(dd) == 1) {
                      mfunct <- c(mfunct, as.character(miRNA_ID_to_Function$miRNAFunction[dd]))
                    } else {
                      mfunct <- c(mfunct, "NA")
                    }
                  }
                  ### Storing final data
                  finaldat <- cbind(Gene = dat_P[, 4], miRNA = dat_P[, 
                    1], dat_P[, 2:3], PearsonScore = dat_P[, 5], DistanceScore = dat_D[, 
                    5], Genefunc = mfunct, Pexp = dat_P[, 6], Dexp = dat_D[, 
                    6], TargetDB = dat_P[, 7], Genefunc = mfunct)
                  finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                    nrow(finaldat)))
                  names(finaldat) <- c("CGID", "miRNA", "Gene_Symbol", 
                    "Genes_FBID", "Pearson_Score", "Distance_Score", "PearsonExperiments", 
                    "DistanceExperiments", "TargetDatabase", "miRNAFunction", 
                    "Method")
                  
                }
            } else {
                stop("File Already Exists!!")
            }
            ## If text is chosen as TRUE; if not returns value to the console Checks
            ## if gene length>1; if yes adds & after each gene name
            if (identical(Text, TRUE)) {
                
                if (length(gene) > 1) {
                  Genes <- paste(gene, collapse = "&")
                } else {
                  Genes = gene
                }
                filename = paste("miRNAs_for_Genes_", Genes, "_", Platform, 
                  "_", method, ".csv", sep = "")
                ### Checking whether file name exists amd writing data in .csv file
                if (!(file.exists(file.path(outpath, filename)))) {
                  write.csv(finaldat, file = file.path(outpath, filename))
                } else {
                  stop("File Already Exists!!")
                }
                
                
            } else {
                return(finaldat)
            }
        }
    } else {
        print("Method and/or Platform does not exist")
    }
    
}
