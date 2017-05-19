###Retrieves the Host Gene for a given miRNA.
#Input: miRNA
#Output:HostGene

library(IntramiRExploreR)
test_extract_HostGene<-function(){
###Checking a positive case
result<-extract_HostGene("dme-miR-12")
expect_equal(length(result), 1)
###Checking a negative case
result1<-extract_HostGene("dme-miR-120")
expect_equal(length(result1), 0)
}
###Retrieves the Intragenic miRNA for a given host gene.
#Input: Host Gene
#Output:Intragenic miRNA
test_extract_intragenic_miR<-function(){
###Checking a positive case
result<-extract_intragenic_miR("Gmap")
expect_equal(length(result), 1)
###Checking a negative case
result1<-extract_intragenic_miR("Syb")
expect_equal(length(result1), 0)
}
