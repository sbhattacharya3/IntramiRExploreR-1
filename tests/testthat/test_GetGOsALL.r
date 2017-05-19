###Retrieves the Targets for a given query miRNA
###using the algorithm.
#Input: miRNA, Statistical method, Platform.
#Output:Dataframe containing miRNA and Target genes (in Gene Symbol,
#       FBGN and CG IDs), score, experments where the correlation happens.

library(IntramiRExploreR)
test_GetGOsALL<-function(){
###Checking a negative case, with less number of genes than threshold
miR="dme-miR-12"
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),
                platform=c("Affy1"),thresh=20)
genes<-a$Target_GeneSymbol
expect_that(GetGOS_ALL(genes,GO=c("topGO"),ontology=c("GO_BP"),path=tempdir() ,
            filename="test"),throws_error())
###Checking a negative case; with 2 ontology as input for topGO
miR="dme-miR-12"
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),
                platform=c("Affy1"),thresh=100)
genes<-a$Target_GeneSymbol
expect_that(GetGOS_ALL(genes,GO=c("topGO"),ontology=c("GO_BP","GO_MF"),path=tempdir() ,
            filename="test"),throws_error())
###Checking a negative case; with no ontology as input for topGO
miR="dme-miR-12"
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),
                platform=c("Affy1"),thresh=100)
genes<-a$Target_GeneSymbol
expect_that(GetGOS_ALL(genes,GO=c("topGO"),path=tempdir() ,
            filename="test"),throws_error())
###Checking a negative case; with no term as input for DAVID
miR="dme-miR-12"
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),
                platform=c("Affy1"),thresh=100)
genes<-a$Target_GeneSymbol
expect_that(GetGOS_ALL(genes,GO=c("DAVID"),path=tempdir() ,
            filename="test"),throws_error())
###Checking a negative case; with 2 terms as input for DAVID
miR="dme-miR-12"
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),
                platform=c("Affy1"),thresh=100)
genes<-a$Target_GeneSymbol
expect_that(GetGOS_ALL(genes,GO=c("DAVID"), term=c("GOTERM_BP_ALL", 
            "GOTERM_MF_ALL"), path=tempdir() ,filename="test"),
  			throws_error())
###Checking a negative case; with 2 incorrect GO
miR="dme-miR-12"
a<-Visualisation(miR,mRNA_type=c("GeneSymbol"),method=c("Both"),
                platform=c("Affy1"),thresh=100)
genes<-a$Target_GeneSymbol
expect_that(GetGOS_ALL(genes,GO=c("ABC"), term=c("GOTERM_BP_ALL", 
            "GOTERM_MF_ALL"), path=tempdir() ,filename="test"),
  			throws_error())
}