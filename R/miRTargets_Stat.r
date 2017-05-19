#' Extracting miRNAs that target a query gene.
#'
#' @param miR  character. miRNA symbol.
#' @param method  character. Choices are "Pearson", "Distance","Both" 
#'                and "BothIntersected"
#' @param Platform  character. Choices are "Affy1","Affy2".
#'                
#' @param Text  logical . To choose between storing the data as text file.
#'              Default is FALSE.
#' @param outpath  character. The path where the data is stored if TEXT=TRUE.
#'                 Default is tempdir().
#' @return Outputs the target information, Target Prediction Score, miRNA 
#'         target function and Target Database that predicts the interaction in
#'         a dataframe.
#'         Depending upon the ouput choice data is stored in the
#'         path specified. Default option prints output to the console. 
#' @examples
#' miRNA="dme-miR-12"
#' miRTargets_Stat (miRNA,method=c ("Pearson"),Platform=c ("Affy1"),Text=FALSE)
#' @import utils grDevices graphics 
#' @importFrom stats na.omit
#' @export 

miRTargets_Stat <- function(miR, method = c("Pearson", "Distance", "Both", 
    "BothIntersect"), Platform = c("Affy1", "Affy2"), outpath = tempdir(), 
    Text = FALSE) {
    ### Loading the Datasets & Initialising Variables
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
    ### Checking if input variables are valid
    stopifnot(is.character(miR), length(miR) > 0, 
        length(extract_HostGene(miR)) > 0)
    # stopifnot(is.character(geneIDType), length(geneIDType) > 0)
    stopifnot(is.character(method), length(method) > 0, 
        method %in% c("Pearson", "Distance", "Both", "BothIntersect"))
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
        finaldat1 <- Affy1_Pearson_Final[which(tolower
            (Affy1_Pearson_Final$miRNA) %in% tolower(miR)), ]
        ### Storing data in data.frame
        if (nrow(finaldat1) == 0) {
            print("Records of the miRNA does not exist")
        } else {
            finaldat1 <- cbind(finaldat1, method = rep("Pearson", nrow(finaldat1)))
            names(finaldat1) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "Score", "GeneFunction", "Experiments", 
                "TargetDatabases", "Method")
        }
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR 
            ###name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Pearson") & identical(Platform, "Affy2")) {
        finaldat1 <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$miRNA)
                %in% tolower(miR)), ]
        ### Storing data in data.frame
        if (nrow(finaldat1) == 0) {
            message("Records of the miRNA does not exist")
        } else {
            finaldat1 <- cbind(finaldat1, method = rep("Pearson", 
                    nrow(finaldat1)))
            names(finaldat1) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "Score", "GeneFunction", "Experiments", 
                "TargetDatabases", "Method")
        }
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR 
            ### name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Distance") & identical(Platform, "Affy1")) {
        finaldat1 <- Affy1_Distance_Final[which(tolower
                    (Affy1_Distance_Final$miRNA) %in% tolower(miR)), ]
        ### Storing data in data.frame
        if (nrow(finaldat1) == 0) {
            message("Records of the miRNA does not exist")
        } else {
            finaldat1 <- cbind(finaldat1, method = rep("Distance", 
                        nrow(finaldat1)))
            names(finaldat1) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "Score", "GeneFunction", "Experiments", 
                "TargetDatabases", "Method")
        }
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR 
            ### name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Distance") & identical(Platform, "Affy2")) {
        finaldat1 <- Affy2_Distance_Final[which(tolower
                    (Affy2_Distance_Final$miRNA) %in% tolower(miR)), ]
        ### Storing data in data.frame
        if (nrow(finaldat1) == 0) {
            message("Records of the miRNA does not exist")
        } else {
            finaldat1 <- cbind(finaldat1, method = rep("Distance", 
                    nrow(finaldat1)))
            names(finaldat1) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "Score", "GeneFunction", "Experiments", 
                "TargetDatabases", "Method")
        }
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR 
            ### name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Both") & identical(Platform, "Affy1")) {
        ### Extracting from Pearson and Distance data Individually, then
        ### concatenating the by row.
        datD <- Affy1_Distance_Final[which(tolower(Affy1_Distance_Final$miRNA)
            %in% tolower(miR)), ]
        if (nrow(datD) == 0) {
            message("Records of the miRNA does not exist")
        } else {
            datD <- cbind(datD, method = rep("Distance", nrow(datD)))
        }
        datP <- Affy1_Pearson_Final[which(tolower(Affy1_Pearson_Final$miRNA)
                %in% tolower(miR)), ]
        if (nrow(datP) == 0) {
            message("Records of the miRNA does not exist")
        } else {
            datP <- cbind(datP, method = rep("Pearson", nrow(datP)))
        }
        
        ### Storing data in data.frame
        finaldat1 <- rbind(datP, datD)
        names(finaldat1) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
            "Targets_CGID", "Score", "GeneFunction", "Experiments", 
            "TargetDatabases", "Method")
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR
            ### name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "Both") & identical(Platform, "Affy2")) {
        ### Extracting from Pearson and Distance data Individually, then
        ### concatenating the by row.
        
        datD <- Affy2_Distance_Final[which(tolower(Affy2_Distance_Final$miRNA)
        %in% tolower(miR)), ]
        if (nrow(datD) == 0) {
            message("Records of the miRNA does not exist")
        } else {
            datD <- cbind(datD, method = rep("Distance", nrow(datD)))
        }
        datP <- Affy2_Pearson_Final[which(tolower(Affy2_Pearson_Final$miRNA)
            %in% tolower(miR)), ]
        if (nrow(datD) == 0) {
            message("Records of the miRNA does not exist")
        } else {
            datP <- cbind(datP, method = rep("Pearson", nrow(datP)))
        }
        ### Storing data in data.frame
        finaldat1 <- rbind(datP, datD)
        names(finaldat1) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
            "Targets_CGID", "Score", "GeneFunction", "Experiments", 
            "TargetDatabases", "Method")
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR 
            ### name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat1, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat1)
        }
    } else if (identical(method, "BothIntersect") & identical(Platform, 
                "Affy1")) {
        ### Extracting from Pearson and Distance data Individually, and then
        ### finding intersection between the miRNAs and extracting data for the
        ### intersection miRNAs appropriately. Length miRNA>1
        finaldat <- data.frame()
        if (length(miR) > 1) {
            for (ii in seq_along(miR)) {
                
                datD <- Affy1_Distance_Final[which(tolower
                    (Affy1_Distance_Final$miRNA) %in% tolower(miR[ii])), ]
                datP <- Affy1_Pearson_Final[which(tolower
                    (Affy1_Pearson_Final$miRNA) %in% tolower(miR[ii])), ]
                if (nrow(datD) == 0 | nrow(datP) == 0) {
                  message("Records of the miRNA does not exist")
                } else {
                  ### Intersection of Distance and Pearson.
                  int <- intersect(as.character(datP$FBGN), 
                    as.character(datD$FBGN))
                  ### Extracting from individual datasets
                  
                  dat_P <- datP[which(as.character(datP$FBGN)
                        %in% as.character(int)),]
                  dat_D <- datD[which(as.character(datD$FBGN) 
                        %in% as.character(int)),]
                }
                ### Storing data for individual miRNAs
                finaldat1 <- cbind(dat_P[, 1:5], DistanceScore = dat_D[,5], 
                        Genefunc = dat_P[, 6], Pexp = dat_P[, 7], 
                        Dexp = dat_D[,7], TargetDB = dat_P[, 8])
                finaldat <- rbind(finaldat, finaldat1)
            }
            ### Storing data for teh query miRNAs in data.frame
            finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                nrow(finaldat)))
            names(finaldat) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "PearsonScore", "DistanceScore", 
                "GeneFunction", "PearsonExperiments", "DistanceExperiments",
                "TargetDatabases", "Method")
        } else {
            datD <- Affy1_Distance_Final[which(tolower
                    (Affy1_Distance_Final$miRNA) %in% tolower(miR)), ]
            datP <- Affy1_Pearson_Final[which(tolower
                (Affy1_Pearson_Final$miRNA) %in% tolower(miR)), ]
            ### Intersection of Distance and Pearson.
            if (nrow(datD) == 0 | nrow(datP) == 0) {
                message("Records of the miRNA does not exist")
            } else {
                int <- intersect(as.character(datP$FBGN), 
                    as.character(datD$FBGN))
                dat_P <- datP[which(as.character(datP$FBGN)
                %in% as.character(int)),]
                dat_D <- datD[which(as.character(datD$FBGN)
                %in% as.character(int)),]
            }
            ### Extracting from individual datasets
            finaldat1 <- cbind(dat_P[, 1:5], DistanceScore = dat_D[, 5], 
                Genefunc = dat_P[, 6], Pexp = dat_P[, 7], Dexp = dat_D[,7],
                TargetDB = dat_P[, 8])
            ### Storing data
            finaldat <- cbind(finaldat1, method = rep("PearsonDistance", 
                nrow(finaldat1)))
            names(finaldat) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "PearsonScore", "DistanceScore",  
                "GeneFunction", "PearsonExperiments", "DistanceExperiments",
                "TargetDatabases", "Method")
        }
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR 
            ### name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat)
        }
    } else if (identical(method, "BothIntersect") & 
            identical(Platform, "Affy2")) {
        ### Extracting from Pearson and Distance data Individually, and then
        ### finding intersection between the miRNAs and extracting data for the
        ### intersection miRNAs appropriately. Length miRNA>1
        finaldat <- data.frame()
        if (length(miR) > 1) {
            for (ii in seq_along(miR)) {
                datD <- Affy2_Distance_Final[which
                    (tolower(Affy2_Distance_Final$miRNA) 
                    %in% tolower(miR[ii])), ]
                datP <- Affy2_Pearson_Final[which(
                    tolower(Affy2_Pearson_Final$miRNA)
                    %in% tolower(miR[ii])), ]
                if (nrow(datD) == 0 | nrow(datP) == 0) {
                    message("Records of the miRNA does not exist")
                } else {
                    ### Intersection of Distance and Pearson.
                    int <- intersect(as.character(datP$FBGN),
                        as.character(datD$FBGN))
                    ### Extracting from individual datasets
                    dat_P <- datP[which(as.character(datP$FBGN)
                        %in% as.character(int)),]
                    dat_D <- datD[which(as.character(datD$FBGN)
                    %in% as.character(int)),]
                }
                ### Storing data for individual miRNAs
                finaldat1 <- cbind(dat_P[, 1:5], dat_D[, 5], dat_P[, 6], 
                    dat_P[, 7], dat_D[, 7], dat_P[, 8])
                finaldat <- rbind(finaldat, finaldat1)
            }
            ### Storing data for teh query miRNAs in data.frame
            finaldat <- cbind(finaldat, method = rep("PearsonDistance", 
                nrow(finaldat)))
            names(finaldat) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "PearsonScore", "DistanceScore", 
                "GeneFunction", "PearsonExperiments", "DistanceExperiments",
                "TargetDatabases", "Method")
        } else {
            datD <- Affy2_Distance_Final[which(tolower
                    (Affy2_Distance_Final$miRNA) %in% tolower(miR)), ]
            datP <- Affy2_Pearson_Final[which(tolower
                (Affy2_Pearson_Final$miRNA) %in% tolower(miR)), ]
            if (nrow(datD) == 0 | nrow(datP) == 0) {
                message("Records of the miRNA does not exist")
            } else {
                ### Intersection of Distance and Pearson.
                int <- intersect(as.character(datP$FBGN), 
                    as.character(datD$FBGN))
                ### Extracting from individual datasets
                dat_P <- datP[which(as.character(datP$FBGN)
                    %in% as.character(int)),]
                dat_D <- datD[which(as.character(datD$FBGN) 
                    %in% as.character(int)),]
            }
            finaldat1 <- cbind(dat_P[, 1:5], dat_D[, 5], dat_P[, 6], dat_P[,7],
                        dat_D[, 7], dat_P[, 8])
            ### Storing data in data.frame
            finaldat <- cbind(finaldat1, method = rep("PearsonDistance", 
                nrow(finaldat1)))
            names(finaldat) <- c("miRNA", "Target_GeneSymbol", "Targets_FBID", 
                "Targets_CGID", "PearsonScore", "DistanceScore", "GeneFunction", 
                "PearsonExperiments", "DistanceExperiments", "TargetDatabases", 
                "Method")
        }
        ### If output is selected as Text, write in .csv file Else print to the
        ### console
        if (identical(Text, TRUE)) {
            ### If miRNA length >1 names concatenated with & else use the miR
            ### name
            if (length(miR) > 1) {
                miRNAs <- paste(miR, collapse = "&")
            } else {
                miRNAs = miR
            }
            filename = paste("Targets_for_miRNAs_", miRNAs, "_", Platform, 
                "_", method, ".csv", sep = "")
            ### Check if the file name exists
            if (!(file.exists(file.path(outpath, filename)))) {
                write.csv(finaldat, file = file.path(outpath, filename))
            } else {
                stop("File Already Exists!!")
            }
        } else {
            return(finaldat)
        }
    } else {
        stop("Method and/or Platform does not exist")
    }
}
