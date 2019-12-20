# dependancies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("ensembldb", quietly = TRUE)) {
    BiocManager::install("ensembldb")
}
library(ensembldb)

if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
}
library(biomaRt)
if (!requireNamespace("topGO", quietly = TRUE)) {
    BiocManager::install("topGO")
}
library(topGO)


# to read in a genelist
genelist= read.table('upinAKT.txt')

# these functions are from biomaRt
listMarts() # will have a few options
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl) # will have many options for different species
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)




idtable = getBM(filters='ensembl_transcript_id_version', values=genelist, attributes=c('external_gene_name', 'ensembl_gene_id'),mart=ensembl)

# get all human genome genes
SAPIENSGO=getBM(mart = ensembl, 
   filter="biotype", 
   value="protein_coding",
   attributes = c( "external_gene_name",
        "ensembl_gene_id", 
        "go_id", 
        "name_1006", 
        "namespace_1003"))

SAPIENSGO2= SAPIENSGO[SAPIENSGO$ensembl_gene_id != '',]
geneID2GO <- by(SAPIENSGO2$go_id, SAPIENSGO2$external_gene_name, function(x) as.character(x))
geneID2GO

all.unique.genes = names(geneID2GO)
matching = factor(as.integer(all.unique.genes%in%idtable$external_gene_name))
names(matching) = all.unique.genes
matching[c('BUB1B','AURKB','BUB1','SKA1','SKA2','SKA3','SGO1','PPP2R5C',
           'PPP2R5B','PPP2R2A','PPP2R2A','PPP2R5D','PPP2R5A','PPP2R5E',
           'PPP1CA','PPP1CC','PPP1CB','NDC80')]


# A helper function to do all the statistical testing and put them together. I have been using 'elim', but comparing the 
# result to other statistics.
GOSummary<- function(GOdata) {
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFis <- getSigGroups(GOdata, test.stat)
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata, test.stat)
  test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test",cutOff = 0.01)
  resultElim <- as.numeric( getSigGroups(GOdata, test.stat) ) # Field of interest for me. as.numeric will put in NAs for strings like '< 1e-30' though.
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  l <- list(classic = resultFis, KS = resultKS, elim = resultElim,weight = resultWeight)
  return(GenTable(object=GOdata, weight=l$weight, classic=l$classic, elim=l$elim, KS=l$KS, orderBy="weight",ranksOf = "classic", topNodes = 50))
}


BP.go = new("topGOdata", ontology='BP'
, allGenes = matching
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

MF.go = new("topGOdata", ontology='MF'
, allGenes = matching
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)

CC.go = new("topGOdata", ontology='CC'
, allGenes = matching
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)


# cuts things down to top 50 terms. See 'topNodes' argument in 'GenTable' called by 'GOSummary'.
annotated_BP = GOSummary(BP.go)
annotated_MF = GOSummary(MF.go)
annotated_CC = GOSummary(CC.go)

# now I look at the results like
View(annotated_BP)
# and maybe filter via a threshold
annotated_BP_sig = annotated_BP[annotated_BP$elim <= .1]
