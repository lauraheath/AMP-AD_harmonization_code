#diverse outcomes patient data table
library(dplyr)
#install.packages('UpSetR')
library(RColorBrewer)
library(scales)

p <- synapser::synGet('syn51757646')
divco <- read.csv(p$path)
#change '90+' to 90 for median age calculation
divco$ageDeath2 <- ifelse(divco$ageDeath=='90+', 90, divco$ageDeath)
divco$ageDeath2[divco$ageDeath2=='missing or unknown'] <- NA
divco$ageDeath2 <- as.numeric(divco$ageDeath2)

table(divco$dataContributionGroup)
table(divco$cohort, divco$dataContributionGroup)
table(divco$race, divco$isHispanic)

table(divco$clinicalMetadataSource)

x <- c(306, 326, 252, 24) #numbers from paper, diverse cohorts only (not counting additional NHW for proteomics)
labels <- c("Black or African American", "Latin American", "Non-Hispanic White", "Other")
piepercent <- round(x/sum(x),2)
piepercent2 <- percent(piepercent, accuracy=1)
cols <- brewer.pal(4, "Set2")
png(file = "piechart_raceEthnicity.png")
pie(x, labels = piepercent2, main = "Diverse Cohorts Population by Race and Ethnicity", col=cols, cex=1.5)
legend("topleft", c("Black or African American", "Latin American", "Non-Hispanic White", "Other"), cex = 1,
       fill=cols)
dev.off()

divco2 <- subset(divco, divco$clinicalMetadataSource=="AMP-AD_DiverseCohorts")
AMPAD1 <- subset(divco, divco$clinicalMetadataSource=='AMP-AD 1.0 Studies')

table(divco2$dataContributionGroup, divco2$cohort)
divco2 %>%
  group_by(sex) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
summary(divco2$ageDeath2)
divco2$over65 <- ifelse(divco2$ageDeath>64.999, TRUE, FALSE)

table(divco2$race, divco2$isHispanic)

columbia <- subset(divco2, divco2$dataContributionGroup=='Columbia')
emory <- subset(divco2, divco2$dataContributionGroup=='Emory')
mayo <- subset(divco2, divco2$dataContributionGroup=='Mayo')
mssm <- subset(divco2, divco2$dataContributionGroup=='MSSM')
rush <- subset(divco2, divco2$dataContributionGroup=='Rush')

columbia %>%
  group_by(sex) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(sex) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mayo %>%
  group_by(sex) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mssm %>%
  group_by(sex) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(sex) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


summary(columbia$ageDeath2)
summary(emory$ageDeath2)
summary(mayo$ageDeath2)
summary(mssm$ageDeath2)
summary(rush$ageDeath2)

columbia %>%
  group_by(race) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
table(columbia$race, columbia$isHispanic)
emory %>%
  group_by(race) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))

table(emory$race, emory$isHispanic)
mayo %>%
  group_by(race) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
table(mayo$race, mayo$isHispanic)
mssm %>%
  group_by(race) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
table(mssm$race, mssm$isHispanic)
rush %>%
  group_by(race) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
table(rush$race, rush$isHispanic)


columbia %>%
  group_by(isHispanic) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(isHispanic) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mayo %>%
  group_by(isHispanic) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mssm %>%
  group_by(isHispanic) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(isHispanic) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


columbia %>%
  group_by(apoeGenotype) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(apoeGenotype) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mayo %>%
  group_by(apoeGenotype) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mssm %>%
  group_by(apoeGenotype) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(apoeGenotype) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))

columbia %>%
  group_by(amyThal) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(amyThal) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mayo %>%
  group_by(amyThal) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(amyThal) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))

columbia %>%
  group_by(amyCerad) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(amyCerad) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mssm %>%
  group_by(amyCerad) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(amyCerad) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


columbia %>%
  group_by(Braak) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(Braak) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mayo %>%
  group_by(Braak) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mssm %>%
  group_by(Braak) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(Braak) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))



columbia %>%
  group_by(reag) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(reag) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mayo %>%
  group_by(reag) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mssm %>%
  group_by(reag) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(reag) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


columbia %>%
  group_by(ADoutcome) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
emory %>%
  group_by(ADoutcome) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mayo %>%
  group_by(ADoutcome) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
mssm %>%
  group_by(ADoutcome) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
rush %>%
  group_by(ADoutcome) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))




#upset plot of samples
p <- synapser::synGet('syn51757645')
biospec <- read.csv(p$path)
table(biospec$dataGenerationSite, biospec$tissue)
table(biospec$dataGenerationSite, biospec$assay)



#upload clinical data and match to individualIDs
p1 <- synapser::synGet('syn51757646')
metadata <- read.csv(p1$path)

meta2 <- dplyr::left_join(biospec, metadata, by='individualID')
rnaseq <- subset(meta2, meta2$assay=='rnaSeq')
wgs <- subset(meta2, meta2$assay=='wholeGenomeSeq')

table(rnaseq$tissue, rnaseq$dataContributionGroup)

#there are duplicate individualIDs, get rid of them
rnaseq2 <- rnaseq %>% distinct(rnaseq$individualID, rnaseq$tissue, .keep_all = TRUE)
wgs2 <- wgs %>% distinct(wgs$individualID, wgs$tissue, .keep_all = TRUE)

table(rnaseq2$cohort, rnaseq2$tissue)
table(rnaseq2$dataContributionGroup, rnaseq2$tissue)

### RUSH UPSET PLOTS
rnaseqRush <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Rush")

acg <- subset(rnaseqRush, rnaseqRush$tissue=='anterior cingulate cortex')
acg <- acg$individualID
dlpfc <- subset(rnaseqRush, rnaseqRush$tissue=='dorsolateral prefrontal cortex')
dlpfc <- dlpfc$individualID
stg <- subset(rnaseqRush, rnaseqRush$tissue=='superior temporal gyrus')
stg <- stg$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Rush')
wgs3 <- wgs3$individualID

setdiff(dlpfc, acg)
listInput2 <- list(ACG_RNASEQ=acg, DLPFC_RNASEQ=dlpfc, STG_RNASEQ=stg, WGS=wgs3)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))


####COLUMBIA UPSET PLOTS
rnaseqCol <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Columbia")

#check tissues available
table(rnaseqCol$tissue)
cn <- subset(rnaseqCol, rnaseqCol$tissue=='caudate nucleus')
cn <- cn$individualID
fc <- subset(rnaseqCol, rnaseqCol$tissue=='frontal cortex')
fc <- fc$individualID
tp <- subset(rnaseqCol, rnaseqCol$tissue=='temporal pole')
tp <- tp$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Columbia')
wgs3 <- wgs3$individualID


listInput2 <- list(TP_RNASEQ=tp, CN_RNASEQ=cn, FC_RNASEQ=fc, WGS=wgs3)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))


#### MAYO UPSET PLOTS
rnaseqMayo <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Mayo")

#check tissues available
table(rnaseqMayo$tissue)
cn <- subset(rnaseqMayo, rnaseqMayo$tissue=='caudate nucleus')
cn <- cn$individualID
dlpfc <- subset(rnaseqMayo, rnaseqMayo$tissue=='dorsolateral prefrontal cortex')
dlpfc <- dlpfc$individualID
tc <- subset(rnaseqMayo, rnaseqMayo$tissue=='temporal cortex')
tc <- tc$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Mayo')
wgs3 <- wgs3$individualID


listInput2 <- list(CN_RNASEQ=cn, DLPFC_RNASEQ=dlpfc, TC_RNASEQ=tc, WGS=wgs3)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))


#### EMORY UPSET PLOTS
rnaseqEm <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Emory")

#check tissues available
table(rnaseqEm$tissue)
cn <- subset(rnaseqEm, rnaseqEm$tissue=='caudate nucleus')
cn <- cn$individualID
dlpfc <- subset(rnaseqEm, rnaseqEm$tissue=='dorsolateral prefrontal cortex')
dlpfc <- dlpfc$individualID
tc <- subset(rnaseqEm, rnaseqEm$tissue=='temporal cortex')
tc <- tc$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Emory')
wgs3 <- wgs3$individualID


listInput2 <- list(CN_RNASEQ=cn, DLPFC_RNASEQ=dlpfc, TC_RNASEQ=tc, WGS=wgs3)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))



#### MSSM UPSET PLOTS
rnaseqMSSM <- subset(rnaseq2, rnaseq2$dataContributionGroup=="MSSM")

#check tissues available
table(rnaseqMSSM$tissue)
cn <- subset(rnaseqMSSM, rnaseqMSSM$tissue=='caudate nucleus')
cn <- cn$individualID
dlpfc <- subset(rnaseqMSSM, rnaseqMSSM$tissue=='dorsolateral prefrontal cortex')
dlpfc <- dlpfc$individualID
stg <- subset(rnaseqMSSM, rnaseqMSSM$tissue=='superior temporal gyrus')
stg <- stg$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='MSSM')
wgs3 <- wgs3$individualID


listInput2 <- list(CN_RNASEQ=cn, DLPFC_RNASEQ=dlpfc, STG_RNASEQ=stg, WGS=wgs3)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))




p1 <- synapser::synGet('syn51757644')
wgsAssay <- read.csv(p1$path)

p2 <- synapser::synGet('syn51757643')
rnaAssay <- read.csv(p2$path)


df <- meta2 %>%
  filter(duplicated(.[["individualID"]]))

df <- meta2 %>%
  group_by(tissue, individualID) %>%
  count() %>%
  filter(n>1) %>%
  ungroup() %>%
  select(-n)






#### any RNAseq vs WGS
### RUSH UPSET PLOTS
rnaseqRush <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Rush")

rna <- unique(rnaseqRush$individualID)

wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Rush')
wgs3 <- wgs3$individualID

#separate by white vs non-white
rush <- subset(metadata, metadata$dataContributionGroup=='Rush')
nonwhite <- subset(rush, rush$race!='White' | rush$isHispanic==TRUE)
nonwhite <- nonwhite$individualID
white <- subset(rush, rush$race=='White' & rush$isHispanic==FALSE)
white <- white$individualID


listInput2 <- list(RNAseq=rna, WGS=wgs3, NONWHITE=nonwhite, WHITE=white)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))



####COLUMBIA UPSET PLOTS
rnaseqCol <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Columbia")
rna <- unique(rnaseqCol$individualID)

#separate by white vs non-white
rush <- subset(metadata, metadata$dataContributionGroup=='Columbia')
nonwhite <- subset(rush, rush$race!='White' | rush$isHispanic==TRUE)
nonwhite <- nonwhite$individualID
white <- subset(rush, rush$race=='White' & rush$isHispanic==FALSE)
white <- white$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Columbia')
wgs3 <- wgs3$individualID


listInput2 <- list(RNAseq=rna, WGS=wgs3, NONWHITE=nonwhite, WHITE=white)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))


#### MAYO UPSET PLOTS
rnaseqMayo <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Mayo")
rna <- unique(rnaseqMayo$individualID)

#separate by white vs non-white
rush <- subset(metadata, metadata$dataContributionGroup=='Mayo')
nonwhite <- subset(rush, rush$race!='White' | rush$isHispanic==TRUE)
nonwhite <- nonwhite$individualID
white <- subset(rush, rush$race=='White' & rush$isHispanic==FALSE)
white <- white$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Mayo')
wgs3 <- wgs3$individualID


listInput2 <- list(RNAseq=rna, WGS=wgs3, NONWHITE=nonwhite, WHITE=white)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))


#### EMORY UPSET PLOTS
rnaseqEm <- subset(rnaseq2, rnaseq2$dataContributionGroup=="Emory")
rna <- unique(rnaseqEm$individualID)

#separate by white vs non-white
rush <- subset(metadata, metadata$dataContributionGroup=='Emory')
nonwhite <- subset(rush, rush$race!='White' | rush$isHispanic==TRUE)
nonwhite <- nonwhite$individualID
white <- subset(rush, rush$race=='White' & rush$isHispanic==FALSE)
white <- white$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='Emory')
wgs3 <- wgs3$individualID


listInput2 <- list(RNAseq=rna, WGS=wgs3, NONWHITE=nonwhite, WHITE=white)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))


#### MSSM UPSET PLOTS
rnaseqMSSM <- subset(rnaseq2, rnaseq2$dataContributionGroup=="MSSM")
rna <- unique(rnaseqMSSM$individualID)

#check tissues available
rush <- subset(metadata, metadata$dataContributionGroup=='MSSM')
nonwhite <- subset(rush, rush$race!='White' | rush$isHispanic==TRUE)
nonwhite <- nonwhite$individualID
white <- subset(rush, rush$race=='White' & rush$isHispanic==FALSE)
white <- white$individualID
wgs3 <- subset(wgs2, wgs2$dataContributionGroup=='MSSM')
wgs3 <- wgs3$individualID


listInput2 <- list(RNAseq=rna, WGS=wgs3, NONWHITE=nonwhite, WHITE=white)
upset(fromList(listInput2), order.by = "freq", point.size=3, text.scale=c(2, 2, 2, 2, 2.4, 4))




#find all participants missing either wgs or rnaseq
#all participant list
individuals <- metadata$individualID
wgs_true <- wgs2$individualID
missings <- as.data.frame(setdiff(individuals, wgs_true))
names(missings)[names(missings) == 'setdiff(individuals, wgs_true)'] <- 'individualID'
missings$wgs_missing = TRUE
#missings$rna_missing = NA

rna_true <- rnaseq2$individualID
setdiff(individuals, rna_true)
missings2 <- as.data.frame(setdiff(individuals, rna_true))
names(missings2)[names(missings2) == 'setdiff(individuals, rna_true)'] <- 'individualID'
missings2$rna_missing = TRUE
#missings2$wgs_missing = NA

missings3 <- merge(missings, missings2, by='individualID', all=TRUE)
missings4 <- left_join(missings3, metadata, by='individualID')
missings4 <- subset(missings4, select=c(individualID, wgs_missing, rna_missing, dataContributionGroup, cohort, race, isHispanic))
missings5 <- missings4
missings5[is.na(missings5)] <- FALSE

#save to Abby's curation miscellaneous folder in synapse:
write.csv(missings5, file="DivCo_WGS_RNAseq_Missings.csv", row.names=FALSE)
file <- synapser::File(path='DivCo_WGS_RNAseq_Missings.csv', parentId='syn51547902')
file <- synapser::synStore(file)






biospecObj <- synapser::synGet('syn51757645')
biospec <- read.csv(biospecObj$path)
biospec <- subset(biospec, biospec$assay=='wholeGenomeSeq')

wgsObj <- synapser::synGet('syn51757644')
wgs <- read.csv(wgsObj$path)

apoe <- subset(divco, divco$apoeGenotype=='missing or unknown')

apoe2 <- left_join(apoe, biospec, by='individualID')


MissingApoes <- subset(apoe2, select=c(individualID, specimenID, dataContributionGroup, race))
write.csv(MissingApoes, file='MissingAPOE.csv', row.names = FALSE)



#tabulate AMP-AD 1.0 patients for separate tables
#need to filter out the 10x multiome patients
p <- synapser::synGet('syn51757645')
biospec <- read.csv(p$path)

AMPAD1a <- left_join(AMPAD1, biospec, by="individualID")
AMPAD1b <- subset(AMPAD1a, AMPAD1a$assay=="TMT quantitation")

AMPAD1 <- AMPAD1b
AMPAD1 %>%
  group_by(cohort) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))

AMPAD1 %>%
  group_by(sex) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))

summary(AMPAD1$ageDeath2)

AMPAD1 %>%
  group_by(race) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))
table(AMPAD1$race, AMPAD1$isHispanic)

AMPAD1 %>%
  group_by(isHispanic) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


AMPAD1 %>%
  group_by(apoeGenotype) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


AMPAD1 %>%
  group_by(amyCerad) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


AMPAD1 %>%
  group_by(Braak) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))



AMPAD1 %>%
  group_by(reag) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))



AMPAD1 %>%
  group_by(ADoutcome) %>% 
  summarise(n = n(), freq = paste0(round(n/nrow(.) * 100), "%"))


