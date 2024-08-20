library(synapser)

synapser::synLogin(email='lheath', password='Q54!A!&9iCfl')


p <- synapser::synGet('syn51757645')
biospec <- read.csv(p$path)
#biospec$WGSsource <- 'Divco'

#upload biospecimen data for wgs from amp-ad 1.0
p <- synapser::synGet('syn53352733')
biospec_amp1 <- read.csv(p$path)
biospec_amp1$WGSsource <- 'AMPAD1.0'
wgs2 <- subset(biospec_amp1, select=c(individualID, WGSsource))
#remove duplicates
wgs3 <- wgs2[!duplicated(wgs2),]

wgs4 <- merge(biospec, wgs3, all=TRUE)

tmt <- subset(wgs4, wgs4$assay=='TMT quantitation')
#how many total have tmt?
df <- unique(tmt$individualID)
length(df)

#1015 people have tmt protoemics


### for tmt manuscript, do counts ####
p2 <- synapser::synGet('syn51757646')
clindata <- read.csv(p2$path)

tmt2 <- merge(tmt, clindata)

dlpfc <- subset(tmt2, tmt2$tissue=='dorsolateral prefrontal cortex')

table(dlpfc$race, dlpfc$isHispanic)
stg <- subset(tmt2, tmt2$tissue=='superior temporal gyrus')
table(stg$race, stg$dlpfc)

#delete duplicate patients, keep one row per patient
ind_dlpfc <- dlpfc[!duplicated(dlpfc$individualID),]

ind_stg <- stg[!duplicated(stg$individualID),]

table(ind_dlpfc$race, ind_dlpfc$isHispanic)

#how many individuals have wgs from amp1.0?
tmt_one <- subset(tmt, tmt$WGSsource=='AMPAD1.0')
df <- unique(tmt_one$individualID)
length(df)

#365 individuals with TMT proteomics have WGS from amp-ad 1.0

rna <- subset(wgs4, wgs4$assay=='rnaSeq')
rna_one <- subset(rna, rna$WGSsource=='AMPAD1.0')
df <- unique(rna_one$individualID)
length(df)

#124 individuals with rnaseq have WGS from amp-ad 1.0

wgsoriginal <- subset(wgs3, wgs3$assay=='wholeGenomeSeq')

#identify individuals with rnaseq data:
rnaonly <- subset(rna, select=c(individualID, assay))
names(rnaonly)[names(rnaonly) == 'assay'] <- 'has_RNAseq'
rnaonly <- rnaonly[!duplicated(rnaonly),]
#844 people total have RNAseq

#see who ONLY has TMT protoemics and NOT rnaseq as well
tmt_rna <- merge(tmt, rnaonly, all=TRUE)
#isolate those without rnaseq 
tmt_norna <- subset(tmt_rna, is.na(tmt_rna$has_RNAseq))
df <- unique(tmt_norna$individualID)
length(df)
#362 individuals with proteomics do NOT have new rnaseq

#attach a column to indicate whether someone has NEW wgs
newWGS <- subset(biospec, biospec$assay=="wholeGenomeSeq")
newWGS <- subset(newWGS, select=c(individualID, assay))
names(newWGS)[names(newWGS) == 'assay'] <- 'has_newWGS'
df <- unique(newWGS$individualID)
length(df)


tmt_rna <- merge(tmt_rna, newWGS, all=TRUE)
df <- unique(tmt_rna$individualID)
length(df)

#who has all three NEW data types?
rna_tmt_wgs <- subset(tmt_rna, tmt_rna$assay=='TMT quantitation' & tmt_rna$has_RNAseq=="rnaSeq" & tmt_rna$has_newWGS=="wholeGenomeSeq")
df <- unique(rna_tmt_wgs$individualID)
length(df)
#508 individuals have ALL THREE NEW

#who has ANY wgs plus proteomics plus rnaseq?
rna_tmt_wgs2 <- subset(tmt_rna, tmt_rna$assay=='TMT quantitation' & tmt_rna$has_RNAseq=="rnaSeq" & (tmt_rna$has_newWGS=="wholeGenomeSeq"|tmt_rna$WGSsource=='AMPAD1.0'))
df <- unique(rna_tmt_wgs2$individualID)
length(df)

#589 have new rnaseq, protoeomics, and ANY wgs

#who has new proteomics, but NO WGS at all?
rna_tmt_NOwgs <- subset(tmt_rna, tmt_rna$assay=='TMT quantitation' & tmt_rna$has_RNAseq=="rnaSeq" & is.na(tmt_rna$has_newWGS) & is.na(tmt_rna$WGSsource))
df <- unique(rna_tmt_NOwgs$individualID)
length(df)
#64 individuals have RNAseq, proteomics, but NO WGS

#how many poeple with protoemics and NO OTHER DATA?
protsonly <- subset(tmt_rna, tmt_rna$assay=='TMT quantitation' & is.na(tmt_rna$has_RNAseq) & is.na(tmt_rna$has_newWGS) & is.na(tmt_rna$WGSsource))
df <- unique(protsonly$individualID)
length(df)
#65 individuals ONLY had proteomics

#how many have proteomics and NEW wgs but no RNAseq?
tmt_wgs_NOrna <- subset(tmt_rna, tmt_rna$assay=='TMT quantitation' & is.na(tmt_rna$has_RNAseq) & tmt_rna$has_newWGS=="wholeGenomeSeq")
df <- unique(tmt_wgs_NOrna$individualID)
length(df)

#14 have prot & NEW WGS, but no RNA

#how many have protoemics and OLD WGS but no RNAseq?
tmt_OLDwgs_NOrna <- subset(tmt_rna, tmt_rna$assay=='TMT quantitation' & is.na(tmt_rna$has_RNAseq) & tmt_rna$WGSsource=="AMPAD1.0")
df <- unique(tmt_OLDwgs_NOrna$individualID)
length(df)

#284 have prot & OLD WGS, but no RNAseq

#how many people had RNAseq and new WGS only?


table(biospec$assay)
table(biospec$specimenMetadataSource)

allbiospec <- rbind(biospec, biospec_amp1)
df <- unique(allbiospec$individualID)
length(df)


table(allbiospec$assay)
table(allbiospec$specimenMetadataSource)
table(allbiospec$WGSsource)

tmt <- subset(allbiospec, allbiospec$assay=='TMT quantitation')

df <- unique(tmt$individualID)
length(df)
table(tmt$WGSsource)


rnaseq <- subset(allbiospec, allbiospec$assay == 'rnaSeq')
p1 <- synapser::synGet('syn3191087')
rosmap1 <- read.csv(p1$path)


p2 <- synapser::synGet('syn51757646')
divco <- read.csv(p2$path)

newros <- subset(divco, divco$individualMetadataLocation=='syn3191087')
newros <- subset(newros, select=individualID)

newros2 <- dplyr::left_join(newros, rosmap1)


p3 <- synapser::synGet('syn51757643')
rnaseq <- read.csv(p3$path)





alldata <- merge(divco, biospec, all=TRUE)

p4 <- synapser::synGet('syn53352733')
wgs <- read.csv(p4$path)





p <- synapser::synTableQuery("SELECT * FROM syn51728783")

clindata <- read.table(p$filepath, sep = ',', header = TRUE)
emory1 <- subset(clindata, clindata$dataContributionGroup=="Emory")


clindata2 <- merge(clindata, biospec, all=TRUE)



p5 <- synapser::synGet('syn53405422')
emory_upload <- read.csv(p5$path)


clindata3 <- merge(emory_upload2, clindata2, all=TRUE)
clindata4 <- subset(clindata3, !is.na(clindata3$IndividualID_emory))



#upload spreadsheet with Columbia data
emoryUploadObj <- "syn53167066"
synapser::synGet(emoryUploadObj, downloadLocation = "/home/lheath/files")
sheets <- (path = "/home/lheath/files/DiverseCohorts_proteomics_biospecimen_metadata_16Dec2023.xlsx")
excel_sheets(path = sheets)
emoryUpload <- read_excel("/home/lheath/files/DiverseCohorts_proteomics_biospecimen_metadata_16Dec2023.xlsx", sheet = "template")



#merge with the biospec and clindata merged file (clindata2)
#get rid of some columns first
emoryUpload2 <- subset(emoryUpload, select=c(IndividualID, specimenID, Emory_SampleID))

emoryUploadAll <- merge(emoryUpload2, clindata2, all=TRUE)
emoryUploadAll <- subset(emoryUploadAll, !is.na(emoryUploadAll$Emory_SampleID))

#who do we have proteomics for?

df <- subset(emory, select="individualID")
df$has_biospec_data <- TRUE

clindata2 <- subset(clindata, clindata$dataContributionGroup=="Emory")
clindata3 <- merge(df, clindata2, all=TRUE)

clindata3 <- subset(clindata3, select = c(individualID, has_biospec_data, dataContributionGroup, cohort))
clindata3[is.na(clindata3)] <- FALSE
write.csv(clindata3, file="EmoryPatients_with_clinicalData.csv", row.names=FALSE)
