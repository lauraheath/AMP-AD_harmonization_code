#upset plots for amp update across 1.0, 2.0, and diverse cohorts
#rnaseq, wgs, snRNA, proteomics, metabolomics

https://github.com/avanlinden/adportal-cohorts/tree/main



library(here)
library(reticulate)
library(tidyverse)
library(fastDummies)
library(ComplexUpset)

# load reticulate and python client
synapse <- reticulate::import("synapseclient")
syn <- synapse$Synapse()

syn$login()


#look at abby's all-rosmap-specimens-datatype file
p <- synapser::synGet('syn26522644')
abby_ros <- read.csv(p$path)


### get all studies that use ROSMAP specimens (satellite studies)
rosmap_clinical <- synapser::synGet("syn3191087")
rosmap_clinical$annotations$study
# which studies have their own biospecimen files and which just use ROSMAP?
rosmap_biospecimen <- synapser::synGet("syn21323366")
#these studies use rosmap individuals but have their own biospecimen files:
rosmap_satellites <- rosmap_clinical$annotations$study[!rosmap_clinical$annotations$study %in% rosmap_biospecimen$annotations$study]
### query the fileview for those biospecimen files and download
#fv_id <- "syn11346063"
#bio_metadata_files <- synapser::synTableQuery(glue::glue("SELECT * FROM {fv_id} WHERE ((\"resourceType\" = 'metadata') AND (\"metadataType\" = 'biospecimen'))"))
#bio_metadata_files <- bio_metadata_files$asDataFrame()

bio_metadata_files <- synapser::synTableQuery(glue::glue("SELECT * FROM syn11346063 WHERE ((\"resourceType\" = 'metadata') AND (\"metadataType\" = 'biospecimen'))"))

bio_metadata_files <- read.table(bio_metadata_files$filepath, sep = ",", header = TRUE)


clean_json_string <- function(json_string, remove_spaces = FALSE) {
  if (remove_spaces) {
    return(gsub("\\[|\"|\\]| ", "", json_string))
  } else {
    return(gsub("\\[|\"|\\]|", "", json_string))
  }
}

bio_metadata_files$study <- clean_json_string(bio_metadata_files$study)
rosmap_satellite_bio_files <- bio_metadata_files %>% 
  as_tibble() %>% 
  select(id, study, metadataType) %>% 
  filter(study %in% rosmap_satellites) %>% 
  mutate(study = unlist(study))

### map through the synIDs and download into a folder
rosmap_satellite_bio_files$id %>% 
  purrr::walk(~synapser::synGet(.x, downloadLocation = here("temp/satellite-biospecimen-files/")))
### create a list of dataframes from the downloaded biospecimen files
files <- list.files(here("temp/satellite-biospecimen-files/"), full.names = TRUE) 
df_list <- files %>% 
  map(~read_csv(.x, col_types = cols(.default = col_character()), id = "file")) %>% 
  set_names(basename(files))
# add a study column
df_study <- map(df_list, ~mutate(., file = basename(file),
                                 study = str_remove(file, "_biospecimen_metadata.csv")))
# add assay columns for studies missing assays
# it's all single-nucleus data
missing_assay_col <- purrr::discard(df_study, ~any(colnames(.x) == "assay")) %>% 
  names()
df_assay <- df_study %>% 
  map_if(~all(colnames(.x) != "assay"), ~mutate(.x, assay = "snrnaSeq"))

# pull out IDs, study, organ, and assay
df_colnames <- 
  map(df_assay, ~colnames(.x))
df_sub <- map(df_assay, ~dplyr::select(., individualID, specimenID, study, organ, tissue, assay))

rosmap_satellite_specimens <- df_sub %>% 
  reduce(rbind)



#upload rosmap 1.0 biospec & clinical data:
# ROSMAP clinical file
clinical <- synapser::synGet("syn3191087")$path %>% 
  read_csv(col_select = c(individualID, projid))

# ROSMAP biospecimen file
rosmap <- synapser::synGet("syn21323366")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, cellType, assay)) %>% 
  mutate(study = "ROSMAP",
         dataStatus = "received")

# satellite specimen info

# satellites <- syn$get("syn26522612")$path %>% 
#   read_csv() %>% 
#   mutate(dataStatus = "received",
#          cellType = NA_character_)
satellites <- rosmap_satellite_specimens
satellites <- satellites %>%
  mutate(dataStatus = "received",
         cellType = NA_character_)


# join all rosmap specimen info
# filter out individuals not in the clinical file
rosmap_combined <- rosmap %>% 
  select(colnames(satellites)) %>% 
  bind_rows(satellites) %>% 
  filter(individualID %in% clinical$individualID)


#fix missing organ and tissue types fields
# chipseq - brain, dlpfc
# rnaArray - brain, dlpfc
# snpArray - just gonna stay NA
# label free proteomics - brain, dlpfc

assays <- c("ChIPSeq", "rnaArray", "label free mass spectrometry")

rosmap_combined <- rosmap
rosmap_combined <- rosmap_combined %>% 
  mutate(organ = if_else(assay %in% assays & is.na(organ) & !is.na(individualID), "brain", organ),
         tissue = if_else(assay %in% assays & is.na(tissue) & !is.na(individualID), "dorsolateral prefrontal cortex", tissue))



geneExpression <- c("mirnaArray", "rnaSeq", "rnaArray","scrnaSeq", "snrnaSeq", "snMultiome")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq", "snATACSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

rosmap_all_dtype <- rosmap_combined %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_))


#save for later:
write_csv(rosmap_all_dtype, here("temp/all-rosmap-specimens-datatype.csv"))

#upload msbb 1.0 biospec:
### MSBB specimens
#look at abby's original file
p <- synapser::synGet('syn26529170')
abby_msbb <- read.csv(p$path)

# the only MSBB satellite study with it's own biospecimens is SuperAgerEpiMap -- so far

# msbb biospecimen file
# remove missing or "unknown" individualIDs

msbb <- synapser::synGet("syn21893059")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  filter(!individualID == "Unknown") %>% 
  filter(!is.na(individualID)) %>% 
  mutate(study = "MSBB",
         dataStatus = "received")

# SuperAgerEpiMap biospecimen file - 10 individuals

superager <- synapser::synGet("syn25724246")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(study = "SuperAgerEpiMap",
         dataStatus = "received") %>% 
  select(colnames(msbb))

# combine all biospecimens 
msbb_combined <- msbb %>% 
  bind_rows(superager)

# add datatypes 
geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq", "ATACseq", "ATACSeq", "HI-C")
genomicVariants <- c("snpArray", "wholeGenomeSeq", "exomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

msbb_combined_dtypes <- msbb_combined %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_)) 
#save file for later
write_csv(msbb_combined_dtypes, here("temp/all-msbb-specimens-datatype.csv"))


#upload mayo 1.0 biospec:
### MayoRNAseq specimens for upset plot
#look at abby's mayo file:
p <- synapser::synGet('syn26529014')
abby_mayo <- read.csv(p$path)
table(abby_mayo$assay)

# get specimens from Mayo Biospecimen file
mayo <- synapser::synGet("syn20827192")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(individualID = as.character(individualID),
         study = "MayoRNAseq",
         dataStatus = "received")

mayo %>% 
  group_by(assay) %>% 
  count()

# not every row has a unique specimen ID
mayo %>% 
  mutate(uSpecimenID = paste0(specimenID, "_", assay)) %>% 
  distinct(uSpecimenID) #ok this at least does it so at least every row is unique

# add projected chipseq specimens from Mariet

mayo_expected <- synapser::synGet("syn26524978")$path %>% 
  read_csv() %>% 
  mutate(specimenID = paste0(SampleID, "_ChIPSeq")) %>% 
  separate(SampleID, into = c("individualID", "tissue"), sep = "_") %>% 
  mutate(organ = "brain",
         tissue = case_when(tissue == "TCX" ~ "temporal cortex",
                            tissue == "CER" ~ "cerebellum",
                            TRUE ~ NA_character_),
         dataStatus = "expected",
         assay = "ChIPSeq",
         study = "MayoRNAseq")



# combine current and expected specimens
# remove 30 missing individualIDs - proteomics study pools

mayo_combined <- mayo_expected %>% select(colnames(mayo)) %>% 
  bind_rows(mayo) %>% 
  filter(!is.na(individualID))

# add datatypes 
geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

mayo_combined_dtypes <- mayo_combined %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_)) 

# save for later
write_csv(mayo_combined_dtypes, here("temp/all-mayo-specimens-datatype.csv"))





### Combined Mayo, MSBB, ROSMAP upset plot

# categories: brain only, WGS, bulk brain RNAseq, all proteomics, all metabolomics

# get datasets
# download de-id data
rosmap <- read.csv(file="temp/all-rosmap-specimens-datatype.csv")
mayo <- read.csv(file="temp/all-mayo-specimens-datatype.csv")
msbb <- read.csv(file="temp/all-msbb-specimens-datatype.csv")

# combine three de-id datasets

colnames(rosmap)
colnames(msbb)
colnames(mayo)


# create new columns for removal: remove blood and sorted microglia
# then remove
rosmap_2 <- rosmap %>% 
  mutate(removeRow = case_when(cellType == "microglia" & assay == "rnaSeq" ~ TRUE,
                               organ == "blood" & assay == "rnaSeq" ~ TRUE,
                               organ == "blood" & dataType == "metabolomics" ~ TRUE,
                               organ == "blood" & dataType == "proteomics" ~ TRUE,
                               TRUE ~ FALSE)) %>% 
  filter(removeRow == FALSE) %>% 
  select(-cellType) %>% 
  mutate(cohort = "ROSMAP")


# make mayo ids character
mayo <- mayo %>% 
  mutate(individualID = as.character(individualID),
         cohort = "MAYO")

msbb <- msbb %>% 
  mutate(cohort = "MSBB")

#bind rows

comb_data <- rosmap_2 %>% 
  bind_rows(msbb) %>% 
  bind_rows(mayo)

# make categories

comb_data_categories <- comb_data %>% 
  mutate(upsetCategory = case_when(assay == "wholeGenomeSeq" ~ "WGS",
                                   assay == "rnaSeq" & organ == "brain" ~ "bulk RNAseq",
                                   assay == "TMT quantitation" | assay == "label free mass spectrometry" ~ "proteomics",
                                   dataType == "metabolomics" & organ == "brain" ~ "metabolomics",
                                   TRUE ~ NA_character_),
  ) %>% 
  filter(!is.na(upsetCategory))

comb_data_categories %>% 
  group_by(upsetCategory, cohort) %>% 
  count()

# make binary data

# use fastDummies::dummy_cols to convert to binary presence/absence matrix per individual
comb_binary_categories <- comb_data_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>%
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# convert to boolean matrix

# which cols to make boolean - ignore individualID fist column
make_boolean_cols <- colnames(comb_binary_categories)[2:length(comb_binary_categories)]
# create boolean df as copy of binary df
comb_boolean_categories <- comb_binary_categories
# convert to boolean
comb_boolean_categories[make_boolean_cols] <- comb_boolean_categories[make_boolean_cols] == 1

# join boolean categories to cohort name

comb_upset_data <- comb_boolean_categories %>% 
  left_join(distinct(select(comb_data, individualID, cohort)))

### Create upset plot ======================

# sort datatype categories for metadata stripes
combUpsetCategories <- c("WGS",
                         "bulk RNAseq",
                         "metabolomics",
                         "proteomics"
)

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# no colors upset plot:
comb_upset_data %>%
  upset(
    combUpsetCategories,
    name = "",
    min_size = 1,
    min_degree = 2,
    width_ratio = 0.2,
    height_ratio = 0.7,
    #sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    stripes = c(light_stripe_color, "white"),
    matrix = (
      intersection_matrix(
        geom = geom_point(size = 2,),
        segment = geom_segment(size = 0.3,
                               color = bar_color),
        outline_color = list(active = bar_color, inactive = inactive_dot_color)
      )
    )
    + scale_color_manual(
      values = c("TRUE" = bar_color, "FALSE" = inactive_dot_color),
      breaks = NULL
    ),
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = cohort),
        text = list(size = 3),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("ROSMAP" = "#5171C0",
                                     "MSBB" = "#5BB0B5",
                                     "MAYO" = "#E566A1"
        )
        )
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.5
    )) +
      theme(axis.ticks.x = element_line(),
      ) +
      scale_fill_manual(values = c("bars_color" = bar_color), guide = "none"),
    themes = upset_modify_themes(
      list(
        'intersections_matrix' = theme(text = element_text(size = 12),
                                       panel.grid = element_blank()),
        'overall_sizes' = theme(text = element_text(size = 10),
                                panel.grid = element_blank()),
        'Intersection size' = theme(text = element_text(size = 12),
                                    panel.grid = element_blank())
      )
    )
  ) +
  labs(caption = "Combined ROSMAP, Mayo, and MSBB brain data")

# save plot and store

ggsave("plots/final/combined_cohorts_upset.pdf",
       width = 12,
       height = 6.75,
       units = "in")








#upload diverse cohorts

# download de-id data
#divco <- read.csv(file="temp/all-divco-specimens-datatype.csv")

# get specimens from DiverseCohorts Biospecimen file
divco <- synapser::synGet("syn51757645")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(individualID = as.character(individualID),
         study = "DiverseCohorts",
         dataStatus = "received",
         removeRow = FALSE,
         cohort = "DiverseCohorts")

divco %>% 
  group_by(assay) %>% 
  count()

#fill in
# add datatypes 
geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

divco_combined_dtypes <- divco %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_)) 

# save for later
write_csv(divco_combined_dtypes, here("temp/all-divco-specimens-datatype.csv"))

#filter out individuals who are already in amp1 and amp2 data
divco_filtered <- divco_combined_dtypes[!(divco_combined_dtypes$individualID %in% comb_data$individualID),]

#join to comb_data
comb_data2 <- comb_data %>% 
  bind_rows(divco_filtered)

# make categories

comb_data_categories <- comb_data2 %>% 
  mutate(upsetCategory = case_when(assay == "wholeGenomeSeq" ~ "WGS",
                                   assay == "rnaSeq" & organ == "brain" ~ "bulk RNAseq",
                                   assay == "TMT quantitation" | assay == "label free mass spectrometry" ~ "proteomics",
                                   dataType == "metabolomics" & organ == "brain" ~ "metabolomics",
                                   TRUE ~ NA_character_),
  ) %>% 
  filter(!is.na(upsetCategory))

comb_data_categories %>% 
  group_by(upsetCategory, cohort) %>% 
  count()

# make binary data

# use fastDummies::dummy_cols to convert to binary presence/absence matrix per individual
comb_binary_categories <- comb_data_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>%
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# convert to boolean matrix

# which cols to make boolean - ignore individualID fist column
make_boolean_cols <- colnames(comb_binary_categories)[2:length(comb_binary_categories)]
# create boolean df as copy of binary df
comb_boolean_categories <- comb_binary_categories
# convert to boolean
comb_boolean_categories[make_boolean_cols] <- comb_boolean_categories[make_boolean_cols] == 1

# join boolean categories to cohort name

comb_upset_data <- comb_boolean_categories %>% 
  left_join(distinct(select(comb_data2, individualID, cohort)))

### Create upset plot ======================

# sort datatype categories for metadata stripes
combUpsetCategories <- c("WGS",
                         "bulk RNAseq",
                         "metabolomics",
                         "proteomics"
)

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# no colors upset plot:
comb_upset_data %>%
  upset(
    combUpsetCategories,
    name = "",
    min_size = 1,
    min_degree = 1,
    width_ratio = 0.2,
    height_ratio = 0.2,
    #sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    stripes = c(light_stripe_color, "white"),
    matrix = (
      intersection_matrix(
        geom = geom_point(size = 2,),
        segment = geom_segment(size = 0.1,
                               color = bar_color),
        outline_color = list(active = bar_color, inactive = inactive_dot_color)
      )
    )
    + scale_color_manual(
      values = c("TRUE" = bar_color, "FALSE" = inactive_dot_color),
      breaks = NULL
    ),
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = cohort),
        text = list(size = 5),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("ROSMAP" = "#5171C0",
                                     "MSBB" = "#5BB0B5",
                                     "MAYO" = "#E566A1",
                                     "DiverseCohorts" = "#DE9A1F"
        )
        )
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.5
    )) +
      theme(axis.ticks.x = element_line(),
      ) +
      scale_fill_manual(values = c("bars_color" = bar_color), guide = "none"),
    themes = upset_modify_themes(
      list(
        'intersections_matrix' = theme(text = element_text(size = 18),
                                       panel.grid = element_blank()),
        'overall_sizes' = theme(text = element_text(size = 12),
                                panel.grid = element_blank()),
        'Intersection size' = theme(text = element_text(size = 18),
                                    panel.grid = element_blank())
      )
    )
  ) +
  labs(caption = "Combined ROSMAP, Mayo, MSBB, and Diverse Cohorts brain data")







#diverse cohorts upset plot by itself
#upload diverse cohorts

# get specimens from DiverseCohorts Biospecimen file
divco <- synapser::synGet("syn51757645")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(individualID = as.character(individualID),
         study = "DiverseCohorts",
         dataStatus = "received",
         removeRow = FALSE)
divclinical <- synapser::synGet("syn51757646")$path %>% 
  read_csv(col_select = c(individualID, dataContributionGroup))

divcoWGS <- synapser::synGet('syn53352733')$path %>%
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>%
  mutate(individualID = as.character(individualID),
         study = "AMP1.0", 
         dataStatus = "received",
         removeRow = FALSE)

divco2 <- divco %>% 
  bind_rows(divcoWGS)


divco3 <- merge(divco2, divclinical)
names(divco3)[names(divco3) == 'dataContributionGroup'] <- 'cohort'





divco3 %>% 
  group_by(assay) %>% 
  count()

divco3 %>% 
  group_by(cohort) %>% 
  count()

divco <- divco3

#delete the NAs for now
divco <- divco[!is.na(divco$cohort),]
#change MSSM to MSBB
divco["cohort"][divco["cohort"] == "MSSM"] <- "MSBB"
divco["cohort"][divco["cohort"] == "Mayo"] <- "MAYO"
divco["cohort"][divco["cohort"] == "Emory"] <- "EMORY"
divco["cohort"][divco["cohort"] == "Rush"] <- "RUSH"
divco["cohort"][divco["cohort"] == "Columbia"] <- "COLUMBIA"

divco %>% 
  group_by(cohort) %>% 
  count()

#fill in
# add datatypes 
geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

divco_combined_dtypes <- divco %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_)) 

# save for later
write_csv(divco_combined_dtypes, here("temp/all-divco-specimens-datatype.csv"))

#filter out individuals who are already in amp1 and amp2 data
#divco_filtered <- divco_combined_dtypes[!(divco_combined_dtypes$individualID %in% comb_data$individualID),]

#join to comb_data
comb_data <- divco_combined_dtypes

# make categories

comb_data_categories <- comb_data %>% 
  mutate(upsetCategory = case_when(assay == "wholeGenomeSeq" ~ "WGS",
                                   assay == "rnaSeq" ~ "bulk RNAseq",
                                   assay == "TMT quantitation"  ~ "proteomics",
                                   TRUE ~ NA_character_),
  ) %>% 
  filter(!is.na(upsetCategory))

comb_data_categories %>% 
  group_by(upsetCategory, cohort) %>% 
  count()

# make binary data

# use fastDummies::dummy_cols to convert to binary presence/absence matrix per individual
comb_binary_categories <- comb_data_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>%
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# convert to boolean matrix

# which cols to make boolean - ignore individualID fist column
make_boolean_cols <- colnames(comb_binary_categories)[2:length(comb_binary_categories)]
# create boolean df as copy of binary df
comb_boolean_categories <- comb_binary_categories
# convert to boolean
comb_boolean_categories[make_boolean_cols] <- comb_boolean_categories[make_boolean_cols] == 1

# join boolean categories to cohort name

comb_upset_data <- comb_boolean_categories %>% 
  left_join(distinct(select(comb_data, individualID, cohort)))

### Create upset plot ======================

# sort datatype categories for metadata stripes
combUpsetCategories <- c("WGS",
                         "bulk RNAseq",
                         "proteomics"
)

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# no colors upset plot:
comb_upset_data %>%
  upset(
    combUpsetCategories,
    name = "",
    min_size = 1,
    min_degree = 1,
    width_ratio = 0.2,
    height_ratio = 0.7,
    #sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    stripes = c(light_stripe_color, "white"),
    matrix = (
      intersection_matrix(
        geom = geom_point(size = 2,),
        segment = geom_segment(size = 0.3,
                               color = bar_color),
        outline_color = list(active = bar_color, inactive = inactive_dot_color)
      )
    )
    + scale_color_manual(
      values = c("TRUE" = bar_color, "FALSE" = inactive_dot_color),
      breaks = NULL
    ),
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = cohort),
        text = list(size = 3),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("RUSH" = "#5171C0",
                                     "MSBB" = "#5BB0B5",
                                     "MAYO" = "#E566A1",
                                     "COLUMBIA" = "#DE9A1F",
                                     "EMORY" = "#2C692C"
        )
        )
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.5
    )) +
      theme(axis.ticks.x = element_line(),
      ) +
      scale_fill_manual(values = c("bars_color" = bar_color), guide = "none"),
    themes = upset_modify_themes(
      list(
        'intersections_matrix' = theme(text = element_text(size = 12),
                                       panel.grid = element_blank()),
        'overall_sizes' = theme(text = element_text(size = 10),
                                panel.grid = element_blank()),
        'Intersection size' = theme(text = element_text(size = 12),
                                    panel.grid = element_blank())
      )
    )
  ) +
  labs(caption = "Diverse Cohorts brain data")











#upset plot of all non-NHW individuals


#get specimens from DiverseCohorts Biospecimen file
divco <- synapser::synGet("syn51757645")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(individualID = as.character(individualID),
         study = "DiverseCohorts",
         dataStatus = "received",
         removeRow = FALSE)
divclinical <- synapser::synGet("syn51757646")$path %>% 
  read_csv(col_select = c(individualID, dataContributionGroup, race, isHispanic, clinicalMetadataSource))

#keep just the people who are not NHW:
divclinical2 <- subset(divclinical, divclinical$clinicalMetadataSource!="AMP-AD 1.0 Studies")
divclinical2 <- subset(divclinical2, (divclinical2$race == "White" & divclinical2$isHispanic==TRUE) | divclinical2$isHispanic==TRUE | 
                         divclinical2$race=="Black or African American" | divclinical2$race=="Asian"|divclinical2$race=="other"|divclinical2$race=="American Indian or Alaska Native")
#divclinical2 <- divclinical[ divclinical$race=="White" & ]

divclinical2$race<-NULL
divclinical2$isHispanic<-NULL
divclinical2$clinicalMetadataSource<-NULL


divco2 <- merge(divco, divclinical2)
names(divco2)[names(divco2) == 'dataContributionGroup'] <- 'cohort'





divco2 %>% 
  group_by(assay) %>% 
  count()

divco2 %>% 
  group_by(cohort) %>% 
  count()

divco <- divco2

#delete the NAs for now
divco <- divco[!is.na(divco$cohort),]
#change MSSM to MSBB
divco["cohort"][divco["cohort"] == "MSSM"] <- "MSBB"
divco["cohort"][divco["cohort"] == "Mayo"] <- "MAYO"
divco["cohort"][divco["cohort"] == "Emory"] <- "EMORY"
divco["cohort"][divco["cohort"] == "Rush"] <- "RUSH"
divco["cohort"][divco["cohort"] == "Columbia"] <- "COLUMBIA"

divco %>% 
  group_by(cohort) %>% 
  count()

#fill in
# add datatypes 
geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

divco_combined_dtypes <- divco %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_)) 



#filter out individuals who are already in amp1 and amp2 data
#divco_filtered <- divco_combined_dtypes[!(divco_combined_dtypes$individualID %in% comb_data$individualID),]

#join to comb_data
comb_data <- divco_combined_dtypes

# make categories

comb_data_categories <- comb_data %>% 
  mutate(upsetCategory = case_when(assay == "wholeGenomeSeq" ~ "WGS",
                                   assay == "rnaSeq" ~ "bulk RNAseq",
                                   assay == "TMT quantitation"  ~ "proteomics",
                                   TRUE ~ NA_character_),
  ) %>% 
  filter(!is.na(upsetCategory))

comb_data_categories %>% 
  group_by(upsetCategory, cohort) %>% 
  count()

# make binary data

# use fastDummies::dummy_cols to convert to binary presence/absence matrix per individual
comb_binary_categories <- comb_data_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>%
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# convert to boolean matrix

# which cols to make boolean - ignore individualID fist column
make_boolean_cols <- colnames(comb_binary_categories)[2:length(comb_binary_categories)]
# create boolean df as copy of binary df
comb_boolean_categories <- comb_binary_categories
# convert to boolean
comb_boolean_categories[make_boolean_cols] <- comb_boolean_categories[make_boolean_cols] == 1

# join boolean categories to cohort name

comb_upset_data <- comb_boolean_categories %>% 
  left_join(distinct(select(comb_data, individualID, cohort)))

### Create upset plot ======================

# sort datatype categories for metadata stripes
combUpsetCategories <- c("WGS",
                         "bulk RNAseq",
                         "proteomics"
)

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# no colors upset plot:
comb_upset_data %>%
  upset(
    combUpsetCategories,
    name = "",
    min_size = 1,
    min_degree = 1,
    width_ratio = 0.2,
    height_ratio = 0.7,
    #sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    stripes = c(light_stripe_color, "white"),
    matrix = (
      intersection_matrix(
        geom = geom_point(size = 2,),
        segment = geom_segment(size = 0.3,
                               color = bar_color),
        outline_color = list(active = bar_color, inactive = inactive_dot_color)
      )
    )
    + scale_color_manual(
      values = c("TRUE" = bar_color, "FALSE" = inactive_dot_color),
      breaks = NULL
    ),
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = cohort),
        text = list(size = 3),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("RUSH" = "#5171C0",
                                     "MSBB" = "#5BB0B5",
                                     "MAYO" = "#E566A1",
                                     "COLUMBIA" = "#DE9A1F",
                                     "EMORY" = "#2C692C"
        )
        )
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.5
    )) +
      theme(axis.ticks.x = element_line(),
      ) +
      scale_fill_manual(values = c("bars_color" = bar_color), guide = "none"),
    themes = upset_modify_themes(
      list(
        'intersections_matrix' = theme(text = element_text(size = 12),
                                       panel.grid = element_blank()),
        'overall_sizes' = theme(text = element_text(size = 10),
                                panel.grid = element_blank()),
        'Intersection size' = theme(text = element_text(size = 12),
                                    panel.grid = element_blank())
      )
    )
  ) +
  labs(caption = "Diverse Cohorts brain data")





######## count number of individuals with single cell data

#get datasets
# download de-id data
rosmap <- read.csv(file="temp/all-rosmap-specimens-datatype.csv")
mayo <- read.csv(file="temp/all-mayo-specimens-datatype.csv")
msbb <- read.csv(file="temp/all-msbb-specimens-datatype.csv")
divco <- read.csv(file="temp/all-divco-specimens-datatype.csv")

# combine three de-id datasets

colnames(rosmap)
colnames(msbb)
colnames(mayo)
colnames(divco)


# create new columns for removal: remove blood and sorted microglia
# then remove
rosmap_2 <- rosmap %>% 
  mutate(removeRow = case_when(cellType == "microglia" & assay == "rnaSeq" ~ TRUE,
                               organ == "blood" & assay == "rnaSeq" ~ TRUE,
                               organ == "blood" & dataType == "metabolomics" ~ TRUE,
                               organ == "blood" & dataType == "proteomics" ~ TRUE,
                               TRUE ~ FALSE)) %>% 
  filter(removeRow == FALSE) %>% 
  select(-cellType) %>% 
  mutate(cohort = "ROSMAP")


# make mayo ids character
mayo <- mayo %>% 
  mutate(individualID = as.character(individualID),
         cohort = "MAYO")

msbb <- msbb %>% 
  mutate(cohort = "MSBB")

#bind rows

comb_data <- rosmap_2 %>% 
  bind_rows(msbb) %>% 
  bind_rows(mayo) %>%
  bind_rows(divco)


table(comb_data$dataType, comb_data$assay)

geneExpression <- c("mirnaArray", "rnaSeq", "rnaArray")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")
singleCell <- c("snrnaSeq", "scrnaSeq", "10x multiome")

alldatatypes <- comb_data %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              assay %in% singleCell ~ "singleCell",
                              TRUE ~ NA_character_)) 
table(alldatatypes$dataType)

#count up individuals with single cell from diverse cohorts
#upload spreadsheet with Columbia data
# singlecell_id <- "syn52136311"
# synapser::synGet(singlecell_id, downloadLocation = "/home/lheath/files")
# sheets <- (path = "/home/lheath/files/AMPAD_DiverseCohorts_Biospecimen_Metadata_Template_Menon_Diversity_snRNAseq_20230410.xlsx")
# excel_sheets(path = sheets)
# vilas <- read_excel("/home/lheath/files/AMPAD_DiverseCohorts_Biospecimen_Metadata_Template_Menon_Diversity_snRNAseq_20230410.xlsx")
# 
# df <- unique(vilas$individualID)
# length(df)

comb_data_categories <- comb_data %>% 
  mutate(upsetCategory = case_when(assay == "wholeGenomeSeq" ~ "WGS",
                                   assay == "rnaSeq" ~ "bulk RNAseq",
                                   assay == "TMT quantitation" | assay == "label free mass spectrometry" ~ "proteomics",
                                   assay == "scrnaSeq" | assay == "10x multiome" ~ "singlecell",
                                   dataType == "metabolomics" & organ == "brain" ~ "metabolomics",
                                   TRUE ~ NA_character_),
  ) %>% 
  filter(!is.na(upsetCategory))

comb_data_categories %>% 
  group_by(upsetCategory, cohort) %>% 
  count()

# make binary data

# use fastDummies::dummy_cols to convert to binary presence/absence matrix per individual
comb_binary_categories <- comb_data_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>%
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# convert to boolean matrix

# which cols to make boolean - ignore individualID fist column
make_boolean_cols <- colnames(comb_binary_categories)[2:length(comb_binary_categories)]
# create boolean df as copy of binary df
comb_boolean_categories <- comb_binary_categories
# convert to boolean
comb_boolean_categories[make_boolean_cols] <- comb_boolean_categories[make_boolean_cols] == 1

# join boolean categories to cohort name

comb_upset_data <- comb_boolean_categories %>% 
  left_join(distinct(select(comb_data, individualID, cohort)))

#count people with rnaseq & singlecell
df <- subset(comb_upset_data, comb_upset_data$`bulk RNAseq`==TRUE & comb_upset_data$singlecell==TRUE)


### Create upset plot ======================

# sort datatype categories for metadata stripes
combUpsetCategories <- c("WGS",
                         "bulk RNAseq",
                         "metabolomics",
                         "singlecell",
                         "proteomics"
)

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# colors upset plot:
comb_upset_data %>%
  upset(
    combUpsetCategories,
    name = "",
    min_size = 1,
    min_degree = 1,
    width_ratio = 0.2,
    height_ratio = 0.2,
    #sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    stripes = c(light_stripe_color, "white"),
    matrix = (
      intersection_matrix(
        geom = geom_point(size = 2,),
        segment = geom_segment(size = 0.1,
                               color = bar_color),
        outline_color = list(active = bar_color, inactive = inactive_dot_color)
      )
    )
    + scale_color_manual(
      values = c("TRUE" = bar_color, "FALSE" = inactive_dot_color),
      breaks = NULL
    ),
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = cohort),
        text = list(size = 5),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("ROSMAP" = "#5171C0",
                                     "MSBB" = "#5BB0B5",
                                     "MAYO" = "#E566A1",
                                     "DiverseCohorts" = "#DE9A1F"
        )
        )
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.5
    )) +
      theme(axis.ticks.x = element_line(),
      ) +
      scale_fill_manual(values = c("bars_color" = bar_color), guide = "none"),
    themes = upset_modify_themes(
      list(
        'intersections_matrix' = theme(text = element_text(size = 18),
                                       panel.grid = element_blank()),
        'overall_sizes' = theme(text = element_text(size = 12),
                                panel.grid = element_blank()),
        'Intersection size' = theme(text = element_text(size = 18),
                                    panel.grid = element_blank())
      )
    )
  ) +
  labs(caption = "Combined ROSMAP, Mayo, MSBB, and Diverse Cohorts brain data")



