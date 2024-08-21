##SCRIPT: Imageomics BeetlePalooza 2024 - Group 3
##Purpose: Evaluate BioClip outputs of 2018 segmented vial images using found taxa per domain as constrained 'class'
##Primary Authors:  Laura Nagel, Sydne Record, Hilmar Lapp
##Additional Group Participants: Evan Waite, Kim Landsbergen, Isa Betancourt
##Date: August 12-15, 2024

# ---- Set file paths ----
setwd("C:/Users/lnagel/Downloads/preds/") #set to location of all files to be read in

# ---- load packages ----

#load packages and set options
library(neonOS)
library(dplyr)
library(stringr)
library(ggplot2)
options(scipen = 15)

# ---- load data ----

#read in and bind all preds output files
preds_all <- do.call(rbind,
                     lapply(list.files(path = getwd(), pattern = "preds"), read.csv))

#read in metadata from vial
#revised file corrects MOAB samples, nulls 'carabidae' from genus, removes subsp. values
metadata <- read.csv("./individual_metadata_revised.csv")

#read in domainsite Df to append domainID to metadata
domainSite <- read.csv('./domain_site_df.csv')

#read in counts of BioClip training images
train_count <- read.csv('./genus_counts_inToL.csv')

#Get NEON taxon Table using neonOS
bet_taxa <- neonOS::getTaxonList(taxonType = "BEETLE",
                                 verbose = T)


# ---- data wrangling ----

## ---- wrangling, taxon table ----

accepted <- bet_taxa %>%
  filter(acceptedTaxonID == taxonID) %>% #remove synonyms from the list
  filter(taxonRank == 'species')%>% #keep only rank of species
  filter(!is.na(tribe)) #remove rows that are missing tribes - mostly hawaii sp.

#write for bioclip testing
#write.table(accepted, file = "./bet_taxonList_palooza.txt",
#sep = "\t", fileEncoding = "UTF-8", row.names = F, na = 'NA')

#slimmed down taxon list for joining to the preds and metadata files
hierarchy <- accepted %>%
  select(taxonID, scientificName, subfamily, tribe, genus, specificEpithet)


## ---- wrangling, bioclip outputs ----

preds_all <- preds_all %>%
  mutate(parentBarcode = substr(file_name, 1, 12))%>% #extract parent barcode
  rename(scientificName = classification)%>%
  dplyr::group_by(file_name)
 
#join preds with hierarchy
preds_hier <- left_join(preds_all, hierarchy, by = 'scientificName')

#relabel metadata
metadata_edit <- metadata %>%
  mutate(parentBarcode = substr(combinedID, 1, 12),
         scientificName = paste(genus, species))%>%
  select(NEON_sampleID, parentBarcode, scientificName)%>%
  mutate(siteID = substr(NEON_sampleID, 1, 4))

metadata_edit <- left_join(metadata_edit, domainSite, by = 'siteID')
  
#join metadata with hierarchy
meta_hier <- left_join(metadata_edit, hierarchy, by = 'scientificName')

#preds and metadata joined on parent Barcode value
preds_join <- left_join(preds_hier, meta_hier, by = 'parentBarcode', multiple = 'first',
                          suffix = c('.pred', ".vial"))

#Used to re-join after summarizing in evaluations sections
testedData <- preds_join %>%
  group_by(parentBarcode) %>%
  select(parentBarcode, scientificName.vial, subfamily.vial, tribe.vial, genus.vial, specificEpithet.vial, siteID, domainID) %>%
  distinct()

#Assign T/F flags for all hierarchy levels
matchFlags <- preds_join %>%
  group_by(file_name) %>%
  mutate(subfamily.match = ifelse(subfamily.pred == subfamily.vial, T, F),
         tribe.match = ifelse(tribe.pred == tribe.vial, T, F),
         genus.match = ifelse(genus.pred == genus.vial, T, F),
         species.match = ifelse(genus.pred == genus.vial & specificEpithet.pred == specificEpithet.vial, T, F))

# ---- Evaluations ----

## ---- species score ----
summary_species <- matchFlags %>%
  group_by(file_name, scientificName.pred) %>%
  summarize(cum.score = round(sum(score), 2))

#Assign a  flag for species comparison
summary_species_flag <- summary_species %>%
  group_by(file_name) %>%
  mutate(parentBarcode = substr(file_name, 1, 12))%>%
  filter(cum.score != 0.00)

match_species <- left_join(summary_species_flag, testedData, by = "parentBarcode")

match_species <- match_species %>%
  mutate(species.match = ifelse(scientificName.pred == scientificName.vial, T, F))

match_species_top <- match_species %>%
  group_by(file_name)%>%
  slice_max(cum.score, n = 1)%>%
  mutate(testType = "species")%>%
  rename(match = species.match)

##---- genus score ----
summary_genus <- matchFlags %>%
  group_by(file_name, genus.pred) %>% #file_name = 1 image
  summarize(cum.score = sum(score))

#Assign a  flag for genus comparison
summary_genus_flag <- summary_genus %>%
  group_by(file_name) %>%
  mutate(parentBarcode = substr(file_name, 1, 12))

match_genus <- left_join(summary_genus_flag, testedData, by = "parentBarcode")
  
match_genus <- match_genus %>%
  mutate(genus.match = ifelse(genus.pred == genus.vial, T, F))

match_genus_top <- match_genus %>%
  group_by(file_name) %>%
  slice_max(cum.score, n = 1) %>%
  mutate(testType = "genus") %>%
  rename(match = genus.match)

## ---- tribe score ----
summary_tribe <- matchFlags %>%
  group_by(file_name, tribe.pred) %>%
  summarize(cum.score = round(sum(score), 2))

#Assign a  flag for tribe comparison
summary_tribe_flag <- summary_tribe %>%
  group_by(file_name) %>%
  mutate(parentBarcode = substr(file_name, 1, 12))%>%
  filter(cum.score != 0.00)

match_tribe <- left_join(summary_tribe_flag, testedData, by = "parentBarcode")

match_tribe <- match_tribe %>%
  mutate(tribe.match = ifelse(tribe.pred == tribe.vial, T, F))

match_tribe_top <- match_tribe %>%
  group_by(file_name) %>%
  slice_max(cum.score, n = 1)%>%
  mutate(testType = 'tribe')%>%
  rename(match = tribe.match)

## ---- subfamily score ----
summary_subfamily <- matchFlags %>%
  group_by(file_name, subfamily.pred) %>%
  summarize(cum.score = round(sum(score), 2))

#Assign a  flag for subfamily comparison
summary_subfamily_flag <- summary_subfamily %>%
  group_by(file_name) %>%
  mutate(parentBarcode = substr(file_name, 1, 12))%>%
  filter(cum.score != 0.00)

match_subfamily <- left_join(summary_subfamily_flag, testedData, by = "parentBarcode")

match_subfamily <- match_subfamily %>%
  mutate(subfamily.match = ifelse(subfamily.pred == subfamily.vial, T, F))

match_subfamily_top <- match_subfamily %>%
  group_by(file_name) %>%
  slice_max(cum.score, n = 1)%>%
  mutate(testType = 'subfamily')%>%
  rename(match = subfamily.match)

# ---- Convert results to long format ----

#long format of all evaluations
longmatch <- rbind(match_species_top, match_genus_top, match_tribe_top, match_subfamily_top)

# ---- visualizations ----

## ---- set factors for graphics ----
#set factors for visualization ordering
structure <- c('subfamily', 'tribe', 'genus', 'species')

## ---- success by scientificName ----
#NOTE: only completed after filtering for a subet of domain(s). 
#successCheck <- longmatch %>%
#  filter(domainID == 'D13') %>%
#  group_by(scientificName.vial, testType)%>%
#  summarize(success = sum(match == T),
#            n = n_distinct(file_name))

#param to plot horizontal line
#testcount <- successCheck %>%
#  group_by(scientificName.vial)%>%
#  summarize(n = unique(n))

#barchart
#ggplot(successCheck, aes(x = factor(testType, structure), y = success, fill = testType)) +
#  geom_bar(stat = 'identity', show.legend = F) +
#  facet_wrap(~scientificName.vial, scales = 'free')+
#  labs(x = 'test type')+
#  geom_hline(data = testcount, aes(yintercept = n))


## ---- success by test type ----

success_mean_bytype <- longmatch %>%
  group_by(testType) %>%
  dplyr::summarize(mean.match = mean(match, na.rm = TRUE))

ggplot(success_mean_bytype, aes(x = factor(testType, structure), y = mean.match, fill = testType)) +
  geom_bar(stat = 'identity', show.legend = F) +
  labs(x = 'test type',
       y = "mean",
       title = "Mean of predictive success by taxonomic level",
       subtitle = "of N = 8719 images from 2018 ethanol vials")


## ---- genus success ----

success_mean_bygenus <- longmatch %>%
  group_by(genus.vial)%>%
  dplyr::summarize(mean.match = mean(match, na.rm = TRUE))

#join bioclip image training data
success_mean_bygenus <- left_join(success_mean_bygenus, train_count, by = c("genus.vial" = "genus"))

#barchart
ggplot(success_mean_bygenus, aes(x = reorder(genus.vial, num_images), y = mean.match, fill = genus.vial)) +
  geom_bar(stat = 'identity', show.legend = F)+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  geom_text(aes(label = num_images), vjust = -0.5) +  # Position the labels slightly above the bars
  labs(x = 'genus',
       y = 'mean',
       title = 'Mean of predictive success by genus',
       subtitle = 'Arranged by increasing number of BioClip training images')
  
