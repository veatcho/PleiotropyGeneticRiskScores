# Author: Olivia Veatch

library(tidyverse)
#Load clean (QCd) WES data
WESData <- readRDS('./SSCdata/AllPDVsinSSC_QCdataset_recent.rds')
WESData <- WESData %>%
  mutate(PDV=as.integer(PDV))
WESData_PDVs <- WESData %>%
  filter(PDV>0)

#Load all protein coding genes indentified in humans from BioMaRt
ensembl <- readRDS('./PublicData/ensIDswithproteins_features.rds')

#Find PDVs in protein coding genes only
WESData_PDVs <- left_join(WESData_PDVs, ensembl[c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'uniprot_gn_id')], 
                          by = c('ENSEMBL'='ensembl_gene_id'))
WESData_PDVs <- WESData_PDVs %>%
  filter(!is.na(uniprot_gn_id))

#Find PDVs in ASD genes from SFARI Gene, https://gene.sfari.org/database/human-gene
Trait1Genes <- readRDS('./PublicData/ASDGene_SFARI2021.rds')
Trait1Genes_PDVsinDataset <- left_join(Trait1Genes[c('gene.symbol')], WESData_PDVs, by = c('gene.symbol'='external_gene_name'))
Trait1Genes_PDVsinDataset <- Trait1Genes_PDVsinDataset %>%
  mutate(PHENO='Trait1') %>%
  distinct()
nTrait1GenesVarsPDVsinData <- n_distinct(Trait1Genes_PDVsinDataset$CHR, Trait1Genes_PDVsinDataset$POS, 
                                         Trait1Genes_PDVsinDataset$END, Trait1Genes_PDVsinDataset$PRED)
#Find PDVs in Sleep Duration genes from https://sleep.hugeamp.org/phenotype.html?phenotype=ChildSleepDuration
#https://pubmed.ncbi.nlm.nih.gov/27568811/
#downloaded top common variant gene-level associations for sleep duration in children with p<=0.05
#The plot and table show genes with the most significant gene-level associations for this phenotype, calculated from bottom-line meta-analyzed genetic associations using the MAGMA (Multi-marker Analysis of GenoMic Annotation) method. Filter the plot and the table by entering a gene name and/or a p-value threshold.
Trait2Genes <- read.csv('./PublicData/gene_table.csv')
Trait2Genes_PDVsinDataset <- inner_join(Trait2Genes[c('gene')], WESData_PDVs, by = c('gene'='external_gene_name'))
Trait2Genes_PDVsinDataset <- Trait2Genes_PDVsinDataset %>%
  mutate(PHENO='Trait2') %>%
  distinct()
nSDVarsPDVsinData <- n_distinct(Trait2Genes_PDVsinDataset$CHR, Trait2Genes_PDVsinDataset$POS, 
                                Trait2Genes_PDVsinDataset$END, Trait2Genes_PDVsinDataset$PRED)

#Predict protein-protein interaction network connecting two distinct conditions with PDVs in candidate genes
library(STRINGdb)
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=0, input_directory="" )
#Select entrez ids for all genes of interest with PDVs in dataset
Trait1genesnames <- Trait1Genes_PDVsinDataset  %>%
  dplyr::select(entrezgene_id) %>%
  distinct()
nTrait1GenesPDVsinData <- n_distinct(Trait1genesnames$entrezgene_id)
Trait2genesnames <- Trait2Genes_PDVsinDataset %>%
  dplyr::select(entrezgene_id) %>%
  distinct()
nTrait2GenesPDVsinData <- n_distinct(Trait2genesnames$entrezgene_id)
Combinedgeneset_PDVs <- rbind(Trait1genesnames, Trait2genesnames)
Combinedgeneset_PDVs <- distinct(Combinedgeneset_PDVs)
Combinedgeneset_mapped <- string_db$map(Combinedgeneset_PDVs, "entrezgene_id", removeUnmappedRows = TRUE )
Combinedgeneset_PPI <- string_db$get_interactions(Combinedgeneset_mapped$STRING_id)
Combinedgeneset_PPI_mediumconfidence <- Combinedgeneset_PPI %>%
  filter(combined_score>=400)
Combinedgeneset_PPI_mediumconfidence.from <- left_join(Combinedgeneset_PPI_mediumconfidence, Combinedgeneset_mapped, 
                                              by=c('from'='STRING_id'))
Combinedgeneset_PPI_mediumconfidence <- left_join(Combinedgeneset_PPI_mediumconfidence.from, Combinedgeneset_mapped, by=c('to'='STRING_id'))
Combinedgeneset_PPI_mediumconfidence <- Combinedgeneset_PPI_mediumconfidence %>%
  mutate(from_gene=entrezgene_id.x, to_gene=entrezgene_id.y) %>%
  dplyr::select(from_gene, to_gene, combined_score) %>%
  distinct()

Combinedgeneset_PPI_mediumconfidence$from_Trait1gene <- Combinedgeneset_PPI_mediumconfidence$from_gene %in% Trait1genesnames$entrezgene_id
Combinedgeneset_PPI_mediumconfidence$to_Trait2gene <- Combinedgeneset_PPI_mediumconfidence$to_gene %in% Trait2genesnames$entrezgene_id
Combinedgeneset_PPI_mediumconfidence$from_Trait2gene <- Combinedgeneset_PPI_mediumconfidence$from_gene %in% Trait2genesnames$entrezgene_id
Combinedgeneset_PPI_mediumconfidence$to_Trait1gene <- Combinedgeneset_PPI_mediumconfidence$to_gene %in% Trait1genesnames$entrezgene_id
Combinedgeneset_PPI_mediumconfidence$Trait1_Trait2connection <- (
  Combinedgeneset_PPI_mediumconfidence$from_Trait1gene=='TRUE' & Combinedgeneset_PPI_mediumconfidence$to_Trait2gene=='TRUE') | 
  (Combinedgeneset_PPI_mediumconfidence$from_Trait2gene=='TRUE' & Combinedgeneset_PPI_mediumconfidence$to_Trait1gene=='TRUE')

Combinedgeneset_PPI_mediumconfidence <- Combinedgeneset_PPI_mediumconfidence %>%
  filter(Trait1_Trait2connection=='TRUE')
saveRDS(Combinedgeneset_PPI_mediumconfidence, './Results/SD_ASD_PPI_mediumconfidence.rds')

#What are overrepresented processes for SDASD gene network
library(topGO)
Combinedgeneset_PPI_mediumconfidence1 <- Combinedgeneset_PPI_mediumconfidence %>%
  dplyr::select(Gene=from_gene) %>%
  distinct()
Combinedgeneset_PPI_mediumconfidence2 <- Combinedgeneset_PPI_mediumconfidence %>%
  dplyr::select(Gene=to_gene) %>%
  distinct()
PPIgenes <- rbind(Combinedgeneset_PPI_mediumconfidence1, Combinedgeneset_PPI_mediumconfidence2)
PPIgenes <- distinct(PPIgenes)
PPIgenestotable <- left_join(PPIgenes, ensembl, by = c('Gene'='entrezgene_id'))
write.csv(PPIgenestotable, './Results/PPIgenes.csv')
#Run ASD gene set overrepresentation analysis
all_genes <- unique(as.character(ensembl[, 'entrezgene_id']))
geneUniverse <- rep(0, length(all_genes))
names(geneUniverse) <- all_genes
PPIgenes_topGO <- PPIgenes %>%
  filter(!is.na(Gene))
genesOfInterest <- PPIgenes$Gene
geneUniverse[names(geneUniverse) %in% genesOfInterest] <- 1
GOBPdata_PPI <- new(
  "topGOdata",
  description = "testdata",
  ontology = "BP",
  allGenes = geneUniverse,
  geneSel = function(p)
    p == 1,
  annot = annFUN.org,
  ID = "entrez",
  mapping = "org.Hs.eg.db"
)

PPIgenes_BPs <-
  runTest(
    GOBPdata_PPI, 
    algorithm="weight01", 
    statistic = "fisher"
    )

#Note p-value threshold is based on benjamini-hochberg of 0.05
testednodes <- n_distinct(GOBPdata_PPI@graph@nodes)
PPIgenes.results <-
  GenTable(GOBPdata_PPI,
           p = PPIgenes_BPs,
           topNodes = testednodes)

#Add full p-values and calculate fold enrichment
Topscores <- as.data.frame(PPIgenes_BPs@score)
write.csv(Topscores, 'Topscores.csv')
Fullpvalues <- read.csv('Topscores.csv')
PPIgenes.results_full <- left_join(PPIgenes.results, Fullpvalues, by=c('GO.ID'='X'))
PPIgenes.results <- PPIgenes.results_full %>%
  mutate(FoldEnrichment=(Significant/Expected)) %>%
  dplyr::select(GO.ID, Term, Annotated, Significant, Expected, 
                p=PPIgenes_BPs.score, FoldEnrichment)
PPIgenes.results <- PPIgenes.results %>%
  mutate(fdr=p.adjust(p, method='BH'))
  
#Pull significant results fdr correcting
PPIgenes.results_sig <- PPIgenes.results %>%
  filter(fdr < 0.05)

#Find genes with PDVs in overrepresented processes
TopBPs <- PPIgenes.results_sig
i <- seq(1, nrow(TopBPs), 1)
TopBPs <- TopBPs %>%
  mutate(GO.ID=as.character(GO.ID), GONUM=paste0("GO", i))
saveRDS(TopBPs, 'TopBPs_ASDSDPPI.rds')

#Identify genes assigned to GO terms that overlap between genes with variants
Genes.in.GOBPs <-
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = TopBPs$GO.ID,
    columns = c("ALIAS", "ENSEMBL", "ENTREZID"),
    keytype = "GOALL"
  )

Genes.in.GOBPs <- Genes.in.GOBPs %>%
  mutate(GOALL=as.factor(GOALL), ENSEMBL=as.factor(ENSEMBL), ENTREZID=as.factor(ENTREZID)) %>%
  distinct(GOALL, EVIDENCEALL, ENSEMBL, .keep_all=TRUE) %>%
  filter(!is.na(ENSEMBL))
saveRDS(Genes.in.GOBPs, 'Genes.in.ASDSDPPIGOBPs.rds')

#Build DBP scores based on ASD genes and SD genes in process
#Calculate scores to determine evidence of genetic dysfunction in biological processes (BPs) of interest to ASD
#Generate weight for level of evidence gene is assigned to process
#Assign frequency weights to evidence codes
FreqEvidenceCode.in.GOBPs <- Genes.in.GOBPs %>%
  group_by(GOALL) %>%
  mutate(EVIDENCEALL, TotalEvidences.in.GOBPs=n()) %>%
  ungroup() %>%
  group_by(GOALL, EVIDENCEALL) %>%
  mutate(ENSEMBL, EvidencebyCat.in.GOBPs=n()) %>%
  ungroup() %>%
  group_by(GOALL) %>%
  mutate(EVIDENCE.FREQ=(EvidencebyCat.in.GOBPs/TotalEvidences.in.GOBPs)) %>%
  ungroup()

#Calculate cumulative score for each gene's assignment to GOBP term
Genes.in.GOBPs_BPWeights <- FreqEvidenceCode.in.GOBPs %>%
  group_by(GOALL, ENSEMBL) %>%  
  mutate(EBP_parent=sum(EVIDENCE.FREQ)) %>%
  ungroup()

#Extract weights for parent terms by process
Genes.in.GOBPs_BPWeights_ <- list()
for (i in 1:nrow(TopBPs)) {
  
  Genes.in.GOBPs_BPWeights_[[paste0("GO", i)]] <- Genes.in.GOBPs_BPWeights %>%
    filter(GOALL==TopBPs$GO.ID[i]) %>%
    dplyr::select(GOALL, ENSEMBL, EBP_parent)
  
}
Genes.in.GOBPs_BPWeights_ALL <- bind_rows(Genes.in.GOBPs_BPWeights_)

#Pull all child terms for each top biological process
library(org.Hs.eg.db)
ChildTermsAll <- as.list(GOBPCHILDREN)

ChildTerms_TopBPs <- list()

for (i in 1:nrow(TopBPs)) {
  
  ChildTerms_TopBPs[[paste0("GO", i)]] <- ChildTermsAll[[TopBPs$GO.ID[i]]]
  
}

#Count genes assigned to child terms of each selected parent process
#Note there are no child terms for GO1 or GO8
Genes.in.ChildTerm_ <- list()

for (i in 2:nrow(TopBPs)) {
  
  ChildTermkey <- ChildTerms_TopBPs[[i]][[1]]
  Genes.in.ChildTerm_[[paste0("GO", i)]] <-
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ChildTermkey,
      columns = c("ALIAS", "ENSEMBL", "ENTREZID"),
      keytype = "GOALL"
    )
  
}

Genes.in.ChildTerms_TopBPs <- bind_rows(Genes.in.ChildTerm_)

Genes.in.ChildTerms_TopBPs <- Genes.in.ChildTerms_TopBPs %>%
  filter(!is.na(ENSEMBL)) %>%
  mutate(ENSEMBL=as.factor(ENSEMBL)) %>%
  group_by(ENSEMBL) %>%
  mutate(nCTsGO=n_distinct(GOALL)) %>%
  ungroup() %>%
  mutate(CTGOGeneScore=(nCTsGO/n_distinct(GOALL))) %>%
  distinct(ENSEMBL, CTGOGeneScore, GOCHILD=GOALL)

#Combine data for GOBP term assignments of ASD candidate genes
Genes.in.GOBPs_BPWeights_ALL <- left_join(Genes.in.GOBPs_BPWeights_ALL, TopBPs, by = c('GOALL'='GO.ID'))

GONUMParents_ChildTerms <- NULL

for (i in 1:nrow(TopBPs)) {
  
  GONUMParents_ChildTerms[[paste0("GO", i)]] <- as.data.frame(ChildTerms_TopBPs[i])
  GONUMParents_ChildTerms[[paste0("GO", i)]] <- GONUMParents_ChildTerms[[paste0("GO", i)]] %>%
    mutate(GONUM=paste0("GO", i)) %>%
    dplyr::select(GOCHILD=paste0("GO", i), GONUM) %>%
    mutate(GOCHILD=as.character(GOCHILD))
  
}

GONUMParents_ChildTerms <- bind_rows(GONUMParents_ChildTerms)

Genes.in.GOBPs_BPWeights_ALL <- left_join(Genes.in.GOBPs_BPWeights_ALL, GONUMParents_ChildTerms, by = "GONUM")
Genes.in.GOBPs_BPWeights_ALL <- left_join(Genes.in.GOBPs_BPWeights_ALL, Genes.in.ChildTerms_TopBPs, 
                                             by = c('ENSEMBL', 'GOCHILD'))
Genes.in.GOBPs_BPWeights_ALL$CTGOGeneScore[is.na(Genes.in.GOBPs_BPWeights_ALL$CTGOGeneScore)] <- 0

Genes.in.GOBPs_BPWeights_final <- Genes.in.GOBPs_BPWeights_ALL %>%
  group_by(ENSEMBL, GOALL) %>%
  mutate(GO_EBPScore=sum(
    EBP_parent#, CTGOGeneScore
    )
    ) %>%
  ungroup() %>%
  dplyr::select(c(ENSEMBL, GOALL, 
                  #GONUM,
                  GO_EBPScore)) %>%
  distinct()

saveRDS(Genes.in.GOBPs_BPWeights_final, 'Genes.in.PPIGOBPs_BPWeights_final.rds')

#Calculate final dysfunctional biological process scores for each individual in the dataset
DBPscores <- mutate(WESData_PDVs, PN=as.factor(PN), CHR=as.factor(CHR), PDV=as.numeric(PDV))
DBPscores <- inner_join(DBPscores, Genes.in.GOBPs_BPWeights_final, by='ENSEMBL')
DBPscores <- DBPscores %>%
  mutate(ENSEMBL=as.factor(ENSEMBL)) %>%
  filter(!is.na(GOALL)) %>%
  distinct()
WESDataPDVs_DBPscores <- DBPscores %>%
  group_by(GOALL) %>%
  mutate(PDVbyEBP_GO=PDV*GO_EBPScore) %>%
  mutate(ngeneBP=n_distinct(ENSEMBL)) %>%
  ungroup() %>%
  group_by(PN, GOALL) %>%
  mutate(DBPnumerator=sum(PDVbyEBP_GO)) %>%
  distinct() %>%
  mutate(DBP_GO=(DBPnumerator/ngeneBP)) %>%
  distinct(PN, DBP_GO,
           #GONUM,
           GOALL) %>%
  ungroup()

#Transpose dataframe to organize DBP scores for different BPs but same individual into one row
DBPscores.t <- pivot_wider(WESDataPDVs_DBPscores, id_cols = "PN",
                                         names_from = "GOALL", values_from = "DBP_GO")
#Replace NAs with zeros for subsequent analysis steps
DBPscores.t[is.na(DBPscores.t)] <- 0
#saveRDS(DBPscores.t, './SSCdata/DBPscores.rds')
WESDatasetIndividuals <- WESData %>%
  dplyr::select(PN) %>%
  distinct()
DBPscores <- left_join(WESDatasetIndividuals, DBPscores.t, by = 'PN')
DBPscores[is.na(DBPscores)] <- 0
saveRDS(DBPscores, './Results/DBPscores_FullDataset.rds')