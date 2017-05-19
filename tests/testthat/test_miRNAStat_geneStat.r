###Retrieves the Targets for a given query miRNA
###using the algorithm.
#Input: miRNA, Statistical method, Platform.
#Output:Dataframe containing miRNA and Target genes (in Gene Symbol,
#       FBGN and CG IDs), score, experments where the correlation happens.

library(IntramiRExploreR)
test_Visualisation<-function(){
###Checking a positive case, with single miRNA
miRNA="dme-miR-12"
result<-Visualisation(miRNA,mRNA_type=c("GeneSymbol"),method=c("Pearson"),
				platform=c("Affy1"))
expect_true(is.data.frame(result))
###Checking a positive case, with multiple miRNA
miRNA=c("dme-miR-12","dme-miR-283")
result<-Visualisation(miRNA,mRNA_type=c("GeneSymbol"),method=c("Pearson"),
				platform=c("Affy1"))
expect_true(is.data.frame(result))
###Checking a negative case with a non-intragenic miRNA
miRNA="dme-miR-120"
expect_that(Visualisation(miRNA,mRNA_type=c("GeneSymbol"),
                method=c("Pearson"),platform=c("Affy1")),throws_error())
###Checking a negative case with an invalid platform
miRNA="dme-miR-12"
expect_that(Visualisation(miRNA,mRNA_type=c("GeneSymbol"),
                method=c("Pearson"),platform=c("Affy12")),throws_error())
###Checking a negative case with an invalid method
miRNA="dme-miR-12"
expect_that(Visualisation(miRNA,mRNA_type=c("GeneSymbol"),
                method=c("Pearsons"),platform=c("Affy12")),throws_error())
###Checking a negative case with an invalid method
miRNA="dme-miR-12"
expect_that(Visualisation(miRNA,mRNA_type=c("ABC"),
                method=c("Pearsons"),platform=c("Affy12")),throws_error())

}
###Retrieves the miRNAS  for a given query Target mRNA
###using the algorithm.
#Input: mRNA, Statistical method, Platform.
#Output:Dataframe containing miRNA and Target genes (in Gene Symbol,
#       FBGN and CG IDs), score, experments where the correlation happens.

test_Genes_Visualisation<-function(){
###Checking a positive case; single mRNA
mRNA="Syb"
result<-Gene_Visualisation(mRNA,mRNA_type=c("GeneSymbol"),method=c("Pearson"),
				  platform=c("Affy1"))
expect_true(is.data.frame(result))
expect_equal(ncol(result),9)
###Checking a positive case;multiple mRNAs
mRNA=c("Syb","Ank2")
result<-Gene_Visualisation(mRNA,mRNA_type=c("GeneSymbol"),method=c("Pearson"),
				  platform=c("Affy1"))
expect_true(is.data.frame(result))
expect_equal(ncol(result),9)
###Checking a negative case with a mRNA with no record
mRNA=c("PkaC1")
expect_that(Gene_Visualisation(mRNA,mRNA_type=c("GeneSymbol"),method=c("Pearson"),
				  platform=c("Affy1")),throws_error())

###Checking a negative case; with invalid method
mRNA=c("Syb")
expect_that(Gene_Visualisation(mRNA,mRNA_type=c("GeneSymbol"),method=c("Pearsson"),
				  platform=c("Affy1")),throws_error())
###Checking a negative case; with invalid platform
mRNA=c("Syb")
expect_that(Gene_Visualisation(mRNA,mRNA_type=c("GeneSymbol"),method=c("Pearson"),
				  platform=c("Affy3")),throws_error())
###Checking a negative case; with invalid mRNA_type 
mRNA=c("Syb")
expect_that(Gene_Visualisation(mRNA,mRNA_type=c("QABCCC"),method=c("Pearson"),
				  platform=c("Affy1")),throws_error())
}
