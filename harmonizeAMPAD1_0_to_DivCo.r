#Harmonize previous AMP-AD projects to match the harmonization done for Diverse Cohorts

#upload clinical metadata from rosmap, msbb, and mayo, see what needs to be amended for full harmonization with diverse cohorts

df <- synapser::synGet('syn51757646')
divco <- read.csv(df$path)

p <- synapser::synGet('syn3191087')
rosmap <- read.csv(p$path, na.strings=c("", "NA"))

p1 <- synapser::synGet('syn6101474')
mssm <- read.csv(p1$path)

p2 <- synapser::synGet('syn23277389')
mayo <- read.csv(p2$path)

##############. ROSMAP harmonization changes
#change column headings
names(rosmap)[names(rosmap) == 'Study'] <- 'cohort'
names(rosmap)[names(rosmap) == 'msex'] <- 'sex'
names(rosmap)[names(rosmap) == 'spanish'] <- 'isHispanic'
names(rosmap)[names(rosmap) == 'age_death'] <- 'ageDeath'
names(rosmap)[names(rosmap) == 'pmi'] <- 'PMI'
names(rosmap)[names(rosmap) == 'apoe_genotype'] <- 'apoeGenotype'
names(rosmap)[names(rosmap) == 'ceradsc'] <- 'amyCerad'
names(rosmap)[names(rosmap) == 'braaksc'] <- 'Braak'

#add columns
rosmap['dataContributionGroup'] <- 'Rush'
rosmap['amyThal'] <- 'missing or unknown'
rosmap['amyA'] <- 'missing or unknown'
rosmap['mayoDx'] <- 'not applicable'

#get rid of unnecessary columns
rosmap$projid<-NULL
rosmap$educ<-NULL
rosmap$age_at_visit_max<-NULL
rosmap$age_first_ad_dx<-NULL
rosmap$cts_mmse30_first_ad_dx<-NULL
rosmap$cts_mmse30_lv<-NULL
rosmap$cogdx<-NULL
rosmap$dcfdx_lv<-NULL


#Change all the variables to match the diverse cohorts harmonization scheme
harmonized_rosmap1 <- rosmap %>% 
  mutate(sex = case_when(sex == 0 ~"female",
                         sex == 1 ~ "male"),
         race = case_when(is.na(race) ~ "missing or unknown",
                          race == 1 ~ "White",
                          race == 2 ~ "Black or African American",
                          race == 3 ~ "American Indian or Alaska Native",
                          race == 4 ~ "Native Hawaiian or Other Pacific Islander",
                          race == 5 ~ "Asian",
                          race == 6 ~ "other",
                          race == 7 ~ "missing or unknown"),
         isHispanic = case_when(is.na(isHispanic) ~ "missing or unknown",
                                isHispanic == 1 ~ "TRUE",
                                isHispanic == 2 ~ "FALSE"),
         ageDeath = case_when(is.na(ageDeath) ~ "missing or unknown",
                              TRUE ~ ageDeath),
         apoeGenotype = case_when(is.na(apoeGenotype) ~ "missing or unknown",
                                  TRUE ~ as.character(apoeGenotype)),
         PMI = case_when(is.na(PMI) ~ "missing or unknown",
                         TRUE ~ as.character(PMI)),
         amyCerad = case_when(amyCerad == 1 ~ "Frequent/Definite/C3",
                              amyCerad == 2 ~ "Moderate/Probable/C2",
                              amyCerad == 3 ~ "Sparse/Possible/C1",
                              amyCerad == 4 ~ "None/No AD/C0",
                              is.na(amyCerad) ~ "missing or unknown"),
         Braak = case_when(is.na(Braak) ~ "missing or unknown",
                           Braak == 0 ~ "None",
                           Braak == 1 ~ "Stage I",
                           Braak == 2 ~ "Stage II",
                           Braak == 3 ~ "Stage III",
                           Braak == 4 ~ "Stage IV",
                           Braak == 5 ~ "Stage V",
                           Braak == 6 ~ "Stage VI")
  )





###############. MSBB harmonization changes
#change column headings
names(mssm)[names(mssm) == 'individualIdSource'] <- 'dataContributionGroup'
names(mssm)[names(mssm) == 'ethnicity'] <- 'isHispanic'
names(mssm)[names(mssm) == 'pmi'] <- 'PMI'
names(mssm)[names(mssm) == 'CERAD'] <- 'amyCerad'

#add columns
mssm['cohort'] <- 'Mt Sinai Brain Bank'
mssm['amyThal'] <- 'missing or unknown'
mssm['amyA'] <- 'missing or unknown'
mssm['mayoDx'] <- 'not applicable'

#get rid of unnecessary columns
mssm$species<-NULL
mssm$yearsEducation<-NULL
mssm$causeDeath<-NULL
mssm$mannerDeath<-NULL
mssm$pH<-NULL
mssm$brainWeight<-NULL
mssm$diagnosis<-NULL
mssm$diagnosisCriteria<-NULL
mssm$CDR<-NULL
mssm$plaqueMean<-NULL
#change PMI into minutes
mssm$PMI <- mssm$PMI/60


#Change all the variables to match the diverse cohorts harmonization scheme
harmonized_mssm1 <- mssm %>% 
  mutate(race = case_when(is.na(race) ~ "missing or unknown",
                          race == 'W' ~ "White",
                          race == 'B' ~ "Black or African American",
                          race == 'H' ~ "other",
                          race == 'U' ~ "missing or unknown"),
         isHispanic = case_when(is.na(isHispanic) ~ "missing or unknown",
                                isHispanic == 'W' | isHispanic == 'B' | isHispanic == 'U' ~ "FALSE",
                                isHispanic == 'H' ~ "TRUE"),
         apoeGenotype = case_when(is.na(apoeGenotype) ~ "missing or unknown",
                                  TRUE ~ as.character(apoeGenotype)),
         amyCerad = case_when(amyCerad == 1 ~ "None/No AD/C0",
                              amyCerad == 2 ~ "Frequent/Definite/C3",
                              amyCerad == 3 ~ "Moderate/Probable/C2",
                              amyCerad == 4 ~ "Sparse/Possible/C1",
                              is.na(amyCerad) ~ "missing or unknown"),
         Braak = case_when(is.na(Braak) ~ "missing or unknown",
                           Braak == 0 ~ "None",
                           Braak == 1 ~ "Stage I",
                           Braak == 2 ~ "Stage II",
                           Braak == 3 ~ "Stage III",
                           Braak == 4 ~ "Stage IV",
                           Braak == 5 ~ "Stage V",
                           Braak == 6 ~ "Stage VI")
  )




##############. MAYO harmonization changes
#change column headings
names(mayo)[names(mayo) == 'pmi'] <- 'PMI'
names(mayo)[names(mayo) == 'diagnosis'] <- 'mayoDx'
names(mayo)[names(mayo) == 'Thal'] <- 'amyThal'

#add columns
mayo['dataContributionGroup'] <- 'Mayo'
mayo['isHispanic'] <- 'missing or unknown'
mayo['cohort'] <- 'Mayo Clinic'
mayo['amyA'] <- NA
mayo['amyCerad'] <- 'missing or unknown'

#get rid of unnecessary columns
mayo$individualIdSource<-NULL
mayo$species<-NULL
mayo$ethnicity<-NULL
mayo$yearsEducation<-NULL
mayo$causeDeath<-NULL
mayo$mannerDeath<-NULL
mayo$pH<-NULL
mayo$brainWeight<-NULL
mayo$diagnosisCriteria<-NULL
mayo$CERAD<-NULL


#Change all the variables to match the diverse cohorts harmonization scheme
harmonized_mayo1 <- mayo %>% 
  mutate(sex = case_when(is.na(sex) ~ "missing or unknown",
                         TRUE ~ sex),
         race = case_when(is.na(race) ~ "missing or unknown",
                          TRUE ~ race),
         ageDeath = case_when(is.na(ageDeath) ~ "missing or unknown",
                              ageDeath == "90_or_over" ~ "90+",
                              TRUE ~ ageDeath),
         apoeGenotype = case_when(is.na(apoeGenotype) ~ "missing or unknown",
                                  TRUE ~ as.character(apoeGenotype)),
         PMI = case_when(is.na(PMI) ~ "missing or unknown",
                         TRUE ~ as.character(PMI)),
         amyThal = case_when(is.na(amyThal) ~ "missing or unknown",
                            amyThal == 0 ~ "None",
                            amyThal == 1 ~ "Phase 1",
                            amyThal == 2 ~ "Phase 2",
                            amyThal == 3 ~ "Phase 3",
                            amyThal == 4 ~ "Phase 4",
                            amyThal == 5 ~ "Phase 5"),
         Braak = case_when(is.na(Braak) ~ "missing or unknown",
                           Braak == 0 ~ "None",
                           Braak == 0.5 | Braak == 1 | Braak == 1.5 ~ "Stage I",
                           Braak == 2 | Braak == 2.5 ~ "Stage II",
                           Braak == 3 | Braak == 3.5 ~ "Stage III",
                           Braak == 4 | Braak == 4.5 ~ "Stage IV",
                           Braak == 5 | Braak == 5.5 ~ "Stage V",
                           Braak == 6 ~ "Stage VI"),
         mayoDx = case_when(is.na(mayoDx) ~ "missing or unknown",
                            mayoDx == "Alzheimer Disease" ~ "AD",
                            mayoDx == "control" ~ "Control",
                            mayoDx == "pathological aging"|mayoDx == "progressive supranuclear palsy" ~ "Other",
                            TRUE ~ mayoDx)
  )


#save files:
write.csv(harmonized_rosmap1, file="ROSMAP_individual_harmonized.csv", row.names=FALSE)
file <- synapser::File(path="ROSMAP_clinical_harmonized.csv", parentId='')
file <- synapser::synStore(file)

write.csv(harmonized_mayo1, file="MAYO_individual_harmonized.csv", row.names=FALSE)
file <- synapser::File(path="MAYO_clinical_harmonized.csv", parentId='')
file <- synapser::synStore(file)

write.csv(harmonized_mssm, file="MSBB_individual_harmonized.csv", row.names=FALSE)
file <- synapser::File(path="MSBB_clinical_harmonized.csv", parentId='')
file <- synapser::synStore(file)


