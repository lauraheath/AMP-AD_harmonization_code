setwd("~/divco")

#Check the extended neuropath variables uploaded by each data contributor
#combine into one file, fix typos/harmonize, output csv
#check the attached data dictionary for variable definitions (syn61089558)

#upload spreadsheets
#Emory
emoryObj <- "syn53938671"
synapser::synGet(emoryObj, downloadLocation = "files/")
path <- "files/Emory_DiverseCohorts_Neuropath_variables_03_06_2024.xlsx"
sheetnames <- excel_sheets(path)
emory <- lapply(excel_sheets(path), read_excel, path = path)
names(emory) <- sheetnames
emory_neuro <- emory$Sheet1
emory_neuro <- as.data.frame(emory_neuro)
dim(emory_neuro)
head(emory_neuro)

#mt sinai
mssmObj <- "syn53214266"
synapser::synGet(mssmObj, downloadLocation = "files/")
path <- "files/MSSM_DiverseCohorts_Neuropath.xlsx"
sheetnames <- excel_sheets(path)
mssm <- lapply(excel_sheets(path), read_excel, path = path)
names(mssm) <- sheetnames
mssm_neuro <- mssm$Sheet1
mssm_neuro <- as.data.frame(mssm_neuro)
dim(mssm_neuro)
head(mssm_neuro)

#mayo
mayoObj <- "syn53214268"
synapser::synGet(mayoObj, downloadLocation = "files/")
path <- "files/20240504_Mayo_DiverseCohorts_Neuropath_Template.xlsx"
sheetnames <- excel_sheets(path)
mayo <- lapply(excel_sheets(path), read_excel, path = path)
names(mayo) <- sheetnames
mayo_neuro <- mayo$`Neuropath Template`
mayo_neuro <- as.data.frame(mayo_neuro)
dim(mayo_neuro)
head(mayo_neuro)

#rush
RushObj <- synapser::synGet('syn53165954')
rush_neuro <- read.csv(RushObj$path)


#change values in each dataframe to all characters for easier binding, transforming
emory_neuro <- emory_neuro %>%
  mutate(across(everything(), as.character))
mayo_neuro <- mayo_neuro %>%
  mutate(across(everything(), as.character))
exneuro <- dplyr::bind_rows(emory_neuro, mayo_neuro)


#change column headings for all rush columns & reattach data contributor/cohort designations
names(rush_neuro)[names(rush_neuro) == 'individualid'] <- 'individualID'
names(rush_neuro)[names(rush_neuro) == 'lewy_full'] <- 'LEWY_full'
names(rush_neuro)[names(rush_neuro) == 'alt_lewy_full'] <- 'alt_LEWY_full'
names(rush_neuro)[names(rush_neuro) == 'lewy_grp'] <- 'LEWY_grp'
names(rush_neuro)[names(rush_neuro) == 'lewy_any'] <- 'LEWY_ANY'
names(rush_neuro)[names(rush_neuro) == 'cdlb'] <- 'CDLB'
names(rush_neuro)[names(rush_neuro) == 'infa'] <- 'INFA'
names(rush_neuro)[names(rush_neuro) == 'hemor'] <- 'HEMOR'
names(rush_neuro)[names(rush_neuro) == 'micr'] <- 'MICR'
names(rush_neuro)[names(rush_neuro) == 'wmrcal'] <- 'WMRCAL'
names(rush_neuro)[names(rush_neuro) == 'wmr'] <- 'WMR'
names(rush_neuro)[names(rush_neuro) == 'wbvd'] <- 'WBVD'
names(rush_neuro)[names(rush_neuro) == 'hs_l'] <- 'HS_L'
names(rush_neuro)[names(rush_neuro) == 'hs_s'] <- 'HS_S'
names(rush_neuro)[names(rush_neuro) == 'alt_hs_s'] <- 'alt_HS_S'
names(rush_neuro)[names(rush_neuro) == 'tdp43'] <- 'TDP43'
names(rush_neuro)[names(rush_neuro) == 'tdp3'] <- 'TDP_3'
names(rush_neuro)[names(rush_neuro) == 'tdp_5'] <- 'TDP_5'
names(rush_neuro)[names(rush_neuro) == 'cvd_s'] <- 'CVD_S'
names(rush_neuro)[names(rush_neuro) == 'cvd_ath'] <- 'CVD_ATH'
names(rush_neuro)[names(rush_neuro) == 'cvd_ath_any'] <- 'CVD_ATH_ANY'
names(rush_neuro)[names(rush_neuro) == 'cvd_art'] <- 'CVD_ART'
names(rush_neuro)[names(rush_neuro) == 'cvd_art_any'] <- 'CVD_ART_ANY'
names(rush_neuro)[names(rush_neuro) == 'cvd_caa'] <- 'CVD_CAA'
names(rush_neuro)[names(rush_neuro) == 'cvd_caa_any'] <- 'CVD_CAA_ANY'
names(rush_neuro)[names(rush_neuro) == 'brainwt'] <- 'BRAINWT'

#upload the individual metadata to reattach data group & cohort
p1 <- synapser::synGet('syn51757646')
clinical <- read.csv(p1$path)

clinical2 <- subset(clinical, select=c(individualID, dataContributionGroup, cohort))

rush_neuro2 <- dplyr::left_join(rush_neuro, clinical2)
rush_neuro2 <- rush_neuro2 %>%
  relocate(dataContributionGroup, .after=individualID)
rush_neuro2 <- rush_neuro2 %>%
  relocate(cohort, .after=dataContributionGroup)
rush_neuro2 <- rush_neuro2 %>%
  mutate(across(everything(), as.character))

exneuro <- dplyr::bind_rows(exneuro, rush_neuro2)

#fix mssm names
names(mssm_neuro)[names(mssm_neuro) == 'Subj. ID'] <- 'individualID'
names(mssm_neuro)[names(mssm_neuro) == 'LEWT_full'] <- 'LEWY_full'
names(mssm_neuro)[names(mssm_neuro) == 'INF'] <- 'INFA'
names(mssm_neuro)[names(mssm_neuro) == 'Lacunes'] <- 'Lacunes_MSSM'
names(mssm_neuro)[names(mssm_neuro) == 'TDO43'] <- 'TDP43'
names(mssm_neuro)[names(mssm_neuro) == 'CVD_ATH'] <- 'CVD_ATH_ANY'
names(mssm_neuro)[names(mssm_neuro) == 'CVD_ART'] <- 'CVD_ART_ANY'


mssm_neuro <- mssm_neuro %>%
  mutate(across(everything(), as.character))
mssm_neuro <- dplyr::left_join(mssm_neuro, clinical2)
mssm_neuro <- mssm_neuro %>%
  relocate(dataContributionGroup, .after=individualID)
mssm_neuro <- mssm_neuro %>%
  relocate(cohort, .after=dataContributionGroup)

exneuro <- dplyr::bind_rows(exneuro, mssm_neuro)

#change all 'missing or unknown' and 'NULL' to NA
#exneuro2 <- exneuro
exneuro[exneuro == 'NULL'] <- NA
exneuro[exneuro == 'missing or unknown'] <- NA
exneuro[exneuro == 'missing or unavailable'] <- NA

#check variables
table(exneuro$dataContributionGroup)
table(exneuro$cohort, exneuro$dataContributionGroup)
table(exneuro$LEWY_full)
table(exneuro$alt_LEWY_full)
table(exneuro$LEWY_grp)
table(exneuro$LEWY_ANY)
table(exneuro$CDLB)
table(exneuro$INFA)
table(exneuro$HEMOR)
table(exneuro$MICR)
table(exneuro$WMRCAL)
table(exneuro$WMR)
table(exneuro$WBVD)
table(exneuro$HS_L)
table(exneuro$HS_S)
table(exneuro$alt_HS_S)
table(exneuro$TDP43)
table(exneuro$TDP_3)
table(exneuro$TDP_5)
table(exneuro$CVD_S)
table(exneuro$CVD_ATH)
table(exneuro$CVD_ATH_ANY)
table(exneuro$CVD_ART)
table(exneuro$CVD_CAA)
table(exneuro$CVD_CAA_ANY)
summary(exneuro$BRAINWT)


#Change all the variables to match the diverse cohorts harmonization scheme

harmonized_exneuro <- exneuro %>% 
  mutate(alt_LEWY_full = case_when(alt_LEWY_full %in% "0 = None" ~ '0',
                                   alt_LEWY_full %in% "1 = Olfactory bulb only" ~ '1',
                                   alt_LEWY_full %in% "2 = Amygdala only" ~ '2',
                                   alt_LEWY_full %in% "3 = Brainstem predominant" ~ '3',
                                   alt_LEWY_full %in% "4 = Limbic (transitional)" ~ '4',
                                   alt_LEWY_full %in% "5 = Neocortical (diffuse)" ~ '5',
                                   TRUE ~ alt_LEWY_full),
         LEWY_ANY = case_when(LEWY_ANY %in% '0 = None' ~ '0',
                              LEWY_ANY %in% '1 = Brainstem, transitional, diffuse, olfactory bulb predominant, or region unspecified' ~ '1',
                              TRUE ~ LEWY_ANY),
         CDLB = case_when(CDLB %in% '1 = Low/No Clinical Signficance' ~ '1',
                          CDLB %in% '2 = Intermediate' ~ '2',
                          CDLB %in% '3 = High' ~ '3',
                          TRUE ~ CDLB),
         INFA = case_when(INFA %in% '0 = No' ~ '0',
                          INFA %in% '0- No' ~ '0',
                          INFA %in% '1 = Yes' ~ '1',
                          INFA %in% '1- Yes' ~ '1',
                          TRUE ~ INFA),
         HEMOR = case_when(HEMOR %in% '0 = No' ~ '0',
                           HEMOR %in% '0- No' ~ '0',
                           HEMOR %in% '1 = Yes' ~'1',
                           HEMOR %in% '1- Yes' ~ '1',
                           HEMOR %in% '9- Unk' ~ NA,
                           TRUE ~ HEMOR),
         MICR = case_when(MICR %in% '0 = No' ~ '0',
                          MICR %in% '1 = Yes' ~ '1',
                          TRUE ~ MICR),
         WMRCAL = case_when(WMRCAL %in% '0 = No' ~ '0',
                            WMRCAL %in% '1 = Yes' ~ '1',
                            TRUE ~ WMRCAL),
         HS_S = case_when(HS_S %in% 'No' ~ '0',
                          HS_S %in% 'Yes' ~ '1',
                          TRUE ~ HS_S),
         alt_HS_S = case_when(alt_HS_S %in% "0 = Absence: Normal neuronal population or Focal/Patchy neuronal loss in all sectors of Ammon's horn or Sommer' s sector." ~ '0',
                              alt_HS_S %in% "1 = Presence: Marked/Severe/Extensive/Moderate neuronal loss in all sectors of Ammon's horn or Sommer' s sector." ~ '1',
                              TRUE ~ alt_HS_S),
         TDP43 = case_when(TDP43 %in% '0 = No' ~ '0',
                           TDP43 %in% '1 = Yes' ~ '1',
                           TDP43 %in% 'No' ~ '0',
                           TRUE ~ TDP43),
         CVD_S = case_when(CVD_S %in% '0 = Absence' ~ '0',
                           CVD_S %in% '1 = Presence' ~ '1',
                           TRUE ~ CVD_S),
         CVD_ATH_ANY = case_when(CVD_ATH_ANY %in% '0 = None' ~ '0',
                                 CVD_ATH_ANY %in% '0- No' ~ '0',
                                 CVD_ATH_ANY %in% '0-No' ~ '0',
                                 CVD_ATH_ANY %in% '1 = Mild, moderate, or severe' ~ '1',
                                 CVD_ATH_ANY %in% '1- Yes' ~ '1',
                                 TRUE ~ CVD_ATH_ANY),
         CVD_ART_ANY = case_when(CVD_ART_ANY %in% '0 = None' ~ '0',
                                 CVD_ART_ANY %in% '0- No' ~ '0',
                                 CVD_ART_ANY %in% '1 = Mild, moderate, or severe' ~ '1',
                                 CVD_ART_ANY %in% '1- Yes' ~ '1',
                                 TRUE ~ CVD_ART_ANY),
         CVD_CAA_ANY = case_when(CVD_CAA_ANY %in% '0 = None' ~ '0',
                                 CVD_CAA_ANY %in% '1 = Mild, Moderate, or Severe' ~ '1',
                                 TRUE ~ CVD_CAA_ANY),
         BRAINWT = case_when(BRAINWT %in% '-9' ~ NA,
                             BRAINWT %in% '9999' ~ NA,
                             TRUE ~ BRAINWT),
         Lacunes_MSSM = case_when(Lacunes_MSSM %in% '0- No' ~ '0',
                                  Lacunes_MSSM %in% '1- Yes' ~ '1',
                                  Lacunes_MSSM %in% '9- Unk' ~ NA,
                                  TRUE ~ Lacunes_MSSM)
         
         )


#some variables are derived from others, calculate these when data contributors did not

#LEWY_grp: calculate using alt_LEWY_full for Mayo, LEWY_full for others
harmonized_exneuro$LEWY_grp <- ifelse(harmonized_exneuro$alt_LEWY_full == 0 & harmonized_exneuro$dataContributionGroup=="Mayo", '0',
                                      ifelse(harmonized_exneuro$alt_LEWY_full == 3 & harmonized_exneuro$dataContributionGroup=="Mayo", '1',
                                             ifelse(harmonized_exneuro$alt_LEWY_full == 1 & harmonized_exneuro$dataContributionGroup=="Mayo", '2',
                                                    ifelse(harmonized_exneuro$alt_LEWY_full == 2 & harmonized_exneuro$dataContributionGroup=="Mayo", '2',
                                                           ifelse(harmonized_exneuro$alt_LEWY_full == 4 & harmonized_exneuro$dataContributionGroup=="Mayo", '2',
                                                                  ifelse(harmonized_exneuro$alt_LEWY_full == 4 & harmonized_exneuro$dataContributionGroup=="Mayo", '2', 
                                                                         ifelse(harmonized_exneuro$alt_LEWY_full == 5 & harmonized_exneuro$dataContributionGroup=="Mayo", '2', harmonized_exneuro$LEWY_grp)))))))

harmonized_exneuro$LEWY_grp[harmonized_exneuro$LEWY_full == '0'] <- '0'
harmonized_exneuro$LEWY_grp[harmonized_exneuro$LEWY_full == '1'] <- '1'
harmonized_exneuro$LEWY_grp[harmonized_exneuro$LEWY_full == '2'] <- '2'
harmonized_exneuro$LEWY_grp[harmonized_exneuro$LEWY_full == '3'] <- '2'
harmonized_exneuro$LEWY_grp[harmonized_exneuro$LEWY_full == '4'] <- '2'

#fill in LEWY_ANY for data groups that did not (Mayo already did,  do not need to account for them)
harmonized_exneuro$LEWY_ANY[harmonized_exneuro$LEWY_full == '0'] <- '0'
harmonized_exneuro$LEWY_ANY[harmonized_exneuro$LEWY_full == '1'] <- '1'
harmonized_exneuro$LEWY_ANY[harmonized_exneuro$LEWY_full == '2'] <- '1'
harmonized_exneuro$LEWY_ANY[harmonized_exneuro$LEWY_full == '3'] <- '1'


#CVD_S is derived, see data dictionary
harmonized_exneuro$CVD_S <- ifelse(((!is.na(harmonized_exneuro$CVD_ATH_ANY)) | (!is.na(harmonized_exneuro$CVD_ART_ANY)) | (!is.na(harmonized_exneuro$CVD_CAA_ANY))), '0', NA)
harmonized_exneuro$CVD_S[harmonized_exneuro$CVD_ATH_ANY == '1'] <- '1'
harmonized_exneuro$CVD_S[harmonized_exneuro$CVD_ART_ANY == '1'] <- '1'
harmonized_exneuro$CVD_S[harmonized_exneuro$CVD_CAA_ANY == '1'] <- '1'



#save the file in extended neuropath folder:
write.csv(harmonized_exneuro, file='Harmonized_Extended_Neuropath_Diverse_Cohorts.csv', row.names=FALSE)
file <- synapser::File(path='Harmonized_Extended_Neuropath_Diverse_Cohorts.csv', parentId='syn53060678')
file <- synapser::synStore(file)
