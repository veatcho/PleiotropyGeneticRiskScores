#Title: Calculating Genetic Risk for Dysfunction in Pleiotropic Biological Processes Using Whole Exome Sequencing Data
#Author: Olivia J. Veatch, PhD, MS
#For details of the initial application of this pipeline to autism and sleep duration, refer to: Veatch et al., Journal of Neurodevelopmental Disorders, 2022, doi: , pmid: , pmcid:
#This code assumes variants called have been previously QC'd and damaging predictions have been made with:
#1) Variant Effect Predictor (VEP); citation("ensemblVEP")
#2) ANNOVAR; https://annovar.openbioinformatics.org/en/latest/user-guide/download/

#Load necessary packages
library(biomaRt)
library(org.Hs.eg.db)
library(STRINGdb)
library(tidyverse)
library(topGO)

#Step 1) Calculate Probably Damaging Variant (PDV) scores using results from VEP and ANNOVAR

#Load clean (QCd) WES data
#to load example data, uncomment below
#WESdata <- readRDS('./DatatoSample.rds')
WESdata_calc <- WESdata %>%
  distinct() %>%
  dplyr::select(
    SampleID, Sex, Chr, Start, End, Reference, Call, zygosity, Gene_Symbol, EnsemblID, #variant information
    VEP_Max_Impact, VEP_max_consequence, #VEP results
    SIFT_pred, Polyphen2_HVAR_pred, LRT_pred, MutationTaster_pred, MutationAssessor_pred, PROVEAN_pred,
    MetaLR_pred, M.CAP_pred, fathmm.MKL_coding_pred #ANNOVAR results
  )

#Change abbreviations for results of various in silico predictions to have common abbreviation for Frequency Damaging (FD) score calcs
WESdata_calc$SIFT_pred[WESdata_calc$SIFT_pred=='T'] <- 'B'
WESdata_calc$LRT_pred[WESdata_calc$LRT_pred=='N'] <- 'B'
WESdata_calc$LRT_pred[WESdata_calc$LRT_pred=='U'] <- '.'
WESdata_calc$Polyphen2_HVAR_pred[WESdata_calc$Polyphen2_HVAR_pred=='P'] <- 'D'
WESdata_calc$MutationTaster_pred[WESdata_calc$MutationTaster_pred=='N'] <- 'B'
WESdata_calc$MutationTaster_pred[WESdata_calc$MutationTaster_pred=='P'] <- 'B'
WESdata_calc$MutationTaster_pred[WESdata_calc$MutationTaster_pred=='A'] <- 'D'
WESdata_calc$MutationAssessor_pred[WESdata_calc$MutationAssessor_pred=='L'] <- 'B'
WESdata_calc$MutationAssessor_pred[WESdata_calc$MutationAssessor_pred=='N'] <- 'B'
WESdata_calc$MutationAssessor_pred[WESdata_calc$MutationAssessor_pred=='H'] <- 'D'
WESdata_calc$MutationAssessor_pred[WESdata_calc$MutationAssessor_pred=='M'] <- 'B'
WESdata_calc$fathmm.MKL_coding_pred[WESdata_calc$fathmm.MKL_coding_pred=='N'] <- 'B'
WESdata_calc$PROVEAN_pred[WESdata_calc$PROVEAN_pred=='N'] <- 'B'
WESdata_calc$MetaLR_pred[WESdata_calc$MetaLR_pred=='T'] <- 'B'
WESdata_calc$M.CAP_pred[WESdata_calc$M.CAP_pred=='T'] <- 'B'

#Count number of predictions that were 'damaging'
WESdata_calc <- cbind(WESdata_calc, Damaging=rowSums(WESdata_calc=='D'))
#Count number of predictions that were 'benign'
WESdata_calc <- cbind(WESdata_calc, Benign=rowSums(WESdata_calc=='B'))

#Generate FD score for each variant
WESdata_calc <- WESdata_calc %>%
  mutate(TotalPreds=(Damaging+Benign))
WESdata_calc <- WESdata_calc %>%
  mutate(FD=(((Damaging-Benign)+1)/(TotalPreds+1)))

#Make negative FD values zero
WESdata_calc$FD[WESdata_calc$FD<=0] <- 0

#Calculate PDV scores
WESdata_calc$zygosity[WESdata_calc$zygosity=='het'] <- 1
WESdata_calc$zygosity[WESdata_calc$zygosity=='hom'] <- 2
WESdata_PDVs <- WESdata_calc %>%
  mutate(Z=as.numeric(zygosity)) %>%
  mutate(PDV=FD*Z) %>%
  mutate(PDV=as.integer(PDV))%>%
  filter(PDV>0)

#Step 2) Find candidate genes with PDVs that encode for proteins in humans

#Pull names for all protein coding genes known in humans that are included in Ensembl (GRCh38.p13)
humangenes <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#listFilters(humangenes)
#listAttributes(humangenes)
hsproteins <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprot_gn_id"), 
                    mart=humangenes, uniqueRows=TRUE)
hsproteins <- hsproteins %>%
  filter(!uniprot_gn_id=="") %>%
  distinct(entrezgene_id, .keep_all=TRUE)

hsproteins <- hsproteins %>%
  mutate(Gene.Name = as.factor(external_gene_name), Human.GeneID = as.factor(entrezgene_id), EnsemblID = as.factor(ensembl_gene_id), UniProtID = as.factor(uniprot_gn_id)) %>%
  dplyr::select(Gene.Name, Human.GeneID, EnsemblID, UniProtID)

#Identify where WES data Ensembl IDs with PDVs match Ensembl IDs encoding proteins in humans
WESdata_PDVs <- inner_join(WESdata_PDVs, hsproteins[c('Gene.Name', 'Human.GeneID', 'EnsemblID', 'UniProtID')], 
                           by = 'EnsemblID')

#Find PDVs in Trait1 candidate genes
#Load dataset with candidate genes associated with first trait of interest
#Initial application of this pipeline to autism used a download of SFARI Gene, https://gene.sfari.org/database/human-gene
Trait1Genes <- readRDS('./PublicData/ASDGene_SFARI2021.rds')
Trait1Genes_PDVsinDataset <- inner_join(WESdata_PDVs, Trait1Genes[c('gene.symbol')], by = c('Gene.Name'='gene.symbol'))
Trait1Genes_PDVsinDataset <- Trait1Genes_PDVsinDataset %>%
  mutate(Phenotype=as.factor('Trait1')) %>%
  distinct()

#Find PDVs in Trait2 candidate genes 
#Load dataset with candidate genes associated with second trait of interest
#Initial application of this pipeline to sleep used a download from https://sleep.hugeamp.org/phenotype.html?phenotype=ChildSleepDuration
#https://pubmed.ncbi.nlm.nih.gov/27568811/
Trait2Genes <- readRDS('./PublicData/SDGenes_ChildSDGWAS.rds')
Trait2Genes_PDVsinDataset <- inner_join(WESdata_PDVs, Trait2Genes[c('gene')], by = c('Gene.Name'='gene'))
Trait2Genes_PDVsinDataset <- Trait2Genes_PDVsinDataset %>%
  mutate(Phenotype=as.factor('Trait2')) %>%
  distinct()

#Remove data from environment that are no longer necessary
rm(humangenes, Trait1Genes, Trait2Genes, WESdata, WESdata_calc)

#Step 3) Predict protein-protein interaction network connecting two distinct conditions with PDVs in candidate genes

#Load STRING database (https://string-db.org/)
string_db <- STRINGdb$new(version="11", species=9606, score_threshold=0, input_directory="")

#Select entrez ids for all genes of interest with PDVs in dataset distinct genes
Trait1genesnames <- Trait1Genes_PDVsinDataset  %>%
  dplyr::select(Human.GeneID) %>%
  distinct()
Trait2genesnames <- Trait2Genes_PDVsinDataset %>%
  dplyr::select(Human.GeneID) %>%
  distinct()

#Combine candidate genes with PDVS for Traits1 & 2 into one dataframe for mapping interaction predictions
Combinedgeneset_PDVs <- rbind(Trait1genesnames, Trait2genesnames)
Combinedgeneset_PDVs <- distinct(Combinedgeneset_PDVs)

#Map combined gene list to 'STRING IDs' 
Combinedgeneset_mapped <- string_db$map(Combinedgeneset_PDVs, "Human.GeneID", removeUnmappedRows = TRUE )

#Map 'STRING IDs' to STRING database predicted interactions
Combinedgeneset_PPI <- string_db$get_interactions(Combinedgeneset_mapped$STRING_id)

#Remove any interactions with confidence scores less than 400
Combinedgeneset_PPI_mediumconfidence <- Combinedgeneset_PPI %>%
  filter(combined_score>=400)

#Identify direct interactions between Trait 1 candidate proteins with PDVs and Trait 2 candidate proteins with PDVs
Combinedgeneset_PPI_mediumconfidence.from <- left_join(Combinedgeneset_PPI_mediumconfidence, Combinedgeneset_mapped, 
                                                       by=c('from'='STRING_id'))
Combinedgeneset_PPI_mediumconfidence <- left_join(Combinedgeneset_PPI_mediumconfidence.from, Combinedgeneset_mapped, 
                                                  by=c('to'='STRING_id'))
Combinedgeneset_PPI_mediumconfidence <- Combinedgeneset_PPI_mediumconfidence %>%
  mutate(from_gene=as.factor(Human.GeneID.x), to_gene=as.factor(Human.GeneID.y)) %>%
  dplyr::select(from_gene, to_gene, combined_score) %>%
  distinct()
Combinedgeneset_PPI_mediumconfidence$from_Trait1gene <- 
  Combinedgeneset_PPI_mediumconfidence$from_gene %in% Trait1genesnames$Human.GeneID
Combinedgeneset_PPI_mediumconfidence$to_Trait2gene <- 
  Combinedgeneset_PPI_mediumconfidence$to_gene %in% Trait2genesnames$Human.GeneID
Combinedgeneset_PPI_mediumconfidence$from_Trait2gene <- 
  Combinedgeneset_PPI_mediumconfidence$from_gene %in% Trait2genesnames$Human.GeneID
Combinedgeneset_PPI_mediumconfidence$to_Trait1gene <- 
  Combinedgeneset_PPI_mediumconfidence$to_gene %in% Trait1genesnames$Human.GeneID
Combinedgeneset_PPI_mediumconfidence$Trait1_Trait2connection <- (
  Combinedgeneset_PPI_mediumconfidence$from_Trait1gene=='TRUE' & Combinedgeneset_PPI_mediumconfidence$to_Trait2gene=='TRUE') | 
  (Combinedgeneset_PPI_mediumconfidence$from_Trait2gene=='TRUE' & Combinedgeneset_PPI_mediumconfidence$to_Trait1gene=='TRUE')
Combinedgeneset_PPI_mediumconfidence <- Combinedgeneset_PPI_mediumconfidence %>%
  filter(Trait1_Trait2connection=='TRUE')

#Step 4) Identify overrepresented processes for Trait1/Trait2 protein interaction network

#Generate gene list for analysis via topGO; citation("topGO")
Combinedgeneset_PPI_mediumconfidence1 <- Combinedgeneset_PPI_mediumconfidence %>%
  dplyr::select(Gene=from_gene) %>%
  distinct()
Combinedgeneset_PPI_mediumconfidence2 <- Combinedgeneset_PPI_mediumconfidence %>%
  dplyr::select(Gene=to_gene) %>%
  distinct()
PPIgenes <- rbind(Combinedgeneset_PPI_mediumconfidence1, Combinedgeneset_PPI_mediumconfidence2)
PPIgenes <- distinct(PPIgenes)

#Remove data from environment that are no longer necessary
rm(Combinedgeneset_mapped, Combinedgeneset_PDVs, Combinedgeneset_PPI,
   Combinedgeneset_PPI_mediumconfidence, Combinedgeneset_PPI_mediumconfidence.from,
   Combinedgeneset_PPI_mediumconfidence1, Combinedgeneset_PPI_mediumconfidence2,
   Trait1Genes_PDVsinDataset, Trait1genesnames, Trait2Genes_PDVsinDataset, Trait2genesnames,
   string_db)

#Prepare candidate gene and comparison gene (i.e., all human protein coding genes) lists for topGO analysis
candidategenes <- PPIgenes$Gene
all_genes <- unique(as.character(hsproteins[, 'Human.GeneID']))
geneUniverse <- rep(0, length(all_genes))
names(geneUniverse) <- all_genes
geneUniverse[names(geneUniverse) %in% candidategenes] <- 1

#Generate topGO object for downstream analysis
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

#Run gene set overrepresentation analysis to find biological processes defined in Gene Ontology (http://geneontology.org/)
PPIgenes_BPs <-
  runTest(
    GOBPdata_PPI, 
    algorithm="weight01", 
    statistic = "fisher"
  )

#Generate table with results of overrepresentation analysis for all tested biological process
testednodes <- n_distinct(GOBPdata_PPI@graph@nodes)
PPIgenes.results <-
  GenTable(GOBPdata_PPI,
           p = PPIgenes_BPs,
           topNodes = testednodes)

#Add full p-values and calculate fold enrichment
Topscores <- as.data.frame(PPIgenes_BPs@score)
Topscores <- mutate(Topscores, GO.ID=row.names(Topscores), pvalue=PPIgenes_BPs@score)
attr(Topscores$pvalue, "names") <- NULL
PPIgenes.results <- left_join(PPIgenes.results, Topscores, by='GO.ID')
PPIgenes.results <- PPIgenes.results %>%
  mutate(FoldEnrichment=(Significant/Expected)) %>%
  dplyr::select(GO.ID, Term, Annotated, Significant, Expected, 
                pvalue, FoldEnrichment)

#Adjust p-values to account for multiple testing with threshold based on benjamini-hochberg of 0.05
PPIgenes.results <- PPIgenes.results %>%
  mutate(fdr=p.adjust(pvalue, method='BH'))

#Pull significant results based on correction
PPIgenes.results <- PPIgenes.results %>%
  filter(fdr < 0.05)

#Step 5) Calculate Dysfunctional Biological Process (DBP) scores based on Trait 1 genes and Trait 2 genes with PDVs assigned to top processes

#Find all genes assigned to top overrepresented biological processes
TopBPs <- PPIgenes.results
i <- seq(1, nrow(TopBPs), 1)
TopBPs <- TopBPs %>%
  mutate(GO.ID=as.character(GO.ID), GONUM=paste0("GO", i))
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

#Calculate weight for level of evidence gene is assigned to process#
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
ChildTermsAll <- as.list(GOBPCHILDREN)
ChildTerms_TopBPs <- list()
for (i in 1:nrow(TopBPs)) {
  ChildTerms_TopBPs[[paste0("GO", i)]] <- ChildTermsAll[[TopBPs$GO.ID[i]]]
}

#Count genes assigned to child terms of each selected parent process
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

#Combine data for GOBP term assignments of candidate genes
Genes.in.GOBPs_BPWeights_ALL <- left_join(Genes.in.GOBPs_BPWeights_ALL, TopBPs, 
                                          by = c('GOALL'='GO.ID'))

#Calculate number of child terms within same branch of GO hierarchy with candidate genes assigned
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
  mutate(GO_EBPScore=sum(EBP_parent)) %>%
  ungroup() %>%
  dplyr::select(c(ENSEMBL, GOALL, GO_EBPScore)) %>%
  distinct()

#Calculate DBP scores for each individual in the dataset
DBPscores <- inner_join(WESdata_PDVs, Genes.in.GOBPs_BPWeights_final, by=c('EnsemblID'='ENSEMBL'))
DBPscores <- DBPscores %>%
  filter(!is.na(GOALL)) %>%
  distinct()
DBPscores <- DBPscores %>%
  group_by(GOALL) %>%
  mutate(PDVbyEBP_GO=PDV*GO_EBPScore) %>%
  mutate(ngeneBP=n_distinct(EnsemblID)) %>%
  ungroup() %>%
  group_by(SampleID, GOALL) %>%
  mutate(DBPnumerator=sum(PDVbyEBP_GO)) %>%
  distinct() %>%
  mutate(DBP_GO=(DBPnumerator/ngeneBP)) %>%
  distinct(SampleID, DBP_GO, GOALL) %>%
  ungroup()

#Transpose dataframe to organize DBP scores for different BPs but same individual into one row
DBPscores.t <- pivot_wider(DBPscores, id_cols = "SampleID",
                           names_from = "GOALL", values_from = "DBP_GO")

#Replace NAs with zeros for subsequent analysis steps
DBPscores.t[is.na(DBPscores.t)] <- 0

#Identify samples from QCd dataset with no evidence of dysfunction
WESDatasetIndividuals <- distinct(WESdata_PDVs, SampleID)
DBPscores <- left_join(WESDatasetIndividuals, DBPscores.t, by = 'SampleID')
DBPscores[is.na(DBPscores)] <- 0

#Save final results for DBP scores calculated for each individual
saveRDS(DBPscores, './DBPscores_FullDataset.rds')

#Clear environment
rm(list = ls())