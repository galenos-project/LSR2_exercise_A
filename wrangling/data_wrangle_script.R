library(dplyr)
library(stringr)
library(tibble)
library(tidyverse)
library(tidyr)
library(purrr)
library(metafor)

source("wrangling/wrangling_functions.R", local = TRUE)

LSR <- 'LSR2'

# Import SyRF outcome data
LSR2_SyRFOutcomes <- read_csv("data/Quantitative_data_-_2024_02_17_-_fd256f61-62c2-4868-9ffb-6a64d694c8ed_-_Investigators_Unblinded.csv")
###### Tidying and cleaning the data ######
#clean ; from TiAb etc
LSR2_SyRFOutcomes$Title <- gsub(";", ":", LSR2_SyRFOutcomes$Title)
LSR2_SyRFOutcomes$Abstract <- gsub(";", ":", LSR2_SyRFOutcomes$Abstract)
LSR2_SyRFOutcomes$Authors <- gsub(";", ":", LSR2_SyRFOutcomes$Authors)
LSR2_SyRFOutcomes$PublicationName <- gsub(";", ":", LSR2_SyRFOutcomes$PublicationName)
LSR2_SyRFOutcomes$AlternateName <- gsub(";", ":", LSR2_SyRFOutcomes$AlternateName)
LSR2_SyRFOutcomes$Url <- gsub(";", ":", LSR2_SyRFOutcomes$Url)
LSR2_SyRFOutcomes$AuthorAddress <- gsub(";", ":", LSR2_SyRFOutcomes$AuthorAddress)
LSR2_SyRFOutcomes$Doi <- gsub(";", ":", LSR2_SyRFOutcomes$Doi)
LSR2_SyRFOutcomes$Keywords <- gsub(";", ":", LSR2_SyRFOutcomes$Keywords)
LSR2_SyRFOutcomes$CustomId <- gsub(";", ":", LSR2_SyRFOutcomes$CustomId)
LSR2_SyRFOutcomes$Year<- gsub(";", ":", LSR2_SyRFOutcomes$Year)
LSR2_SyRFOutcomes$ReferenceType <- gsub(";", ":", LSR2_SyRFOutcomes$ReferenceType)
LSR2_SyRFOutcomes$PdfRelativePath <- gsub(";", ":", LSR2_SyRFOutcomes$PdfRelativePath)
LSR2_SyRFOutcomes$`Type of outcome:` <- gsub(";", ":", LSR2_SyRFOutcomes$`Type of outcome:`)

#change categorisation of outcome types
ot_FR <- read_csv("data/ot_FR.csv")
ot_FR <- ot_FR[,c(3,5)]
ot_FR <- subset(ot_FR, !is.na(ot_FR$Changes))

merged_df <- merge(LSR2_SyRFOutcomes, ot_FR, by.x = "OutcomeLabel", by.y = "Var2", all.x = TRUE)
merged_df <- unique(merged_df)
# Update the 'Type of outcome' column based on the condition
merged_df$`Type of outcome:` <- ifelse(!is.na(merged_df$Changes), merged_df$Changes, merged_df$`Type of outcome:`)

# Filter for reconciled studies (and rename columns for consistency with shiny outcomes/remove SYRF columnns)
LSR2_reconciled <- merged_df %>% 
  filter(`is this a reconciliation?` == TRUE) %>%  #for PTSD `is this a reconciliation?`
  rename(ModelID = `DiseaseModelId(s)`, 
         ExperimentID = ExperimentId, 
         InterventionLabel = `TreatmentLabel(s)`,
         InterventionID = `TreatmentId(s)`,
         OutcomeResult = OutcomeAverage)

## first pass for removing duplicate reconciliations
  recent_reconciledID <- LSR2_reconciled %>%
  arrange(desc(DateTimeDataEntered)) %>%
  group_by(StudyId) %>%
  slice(1) %>%
  ungroup()

# Filter all records that have a StudyIdStr and AnnotatorIdStr that is in the reconciled_study_annotator_pairs, to get all reconciled data
#reconciled_records_unique <- LSR2_reconciled %>%
#  semi_join(recent_reconciledID, by = "StudyId")

## Some studies were split due to their complexity, with RoB/ Arrive only entered once
## the ARRIVE/RoB data for the second split needs overwritten with the values from the first
#if this is the case, identify the studyIds concerned, recort in ll48 and 49, and change ll52-59 accordingly

## Match split studies [Get pairs from bibliographic download on SyRF]
## Find position of the start and end of ARRIVE ROB columns to store column names between the start and end index
col_names <- names(LSR2_reconciled)
start_index <- match("Title", col_names)
end_index <- match("Is any role of the funder in the design/analysis/reporting of study described?", col_names)
ARRIVEROB_columns <- col_names[start_index:end_index] 


## 3. Further step to remove observations from an accidental dual-reconciliation: choose the most recent reconciliation (unlikely to be necessary)
# Choose the most recent reconciliation 
recent_reconciledID <- LSR2_reconciled %>%
  arrange(desc(DateTimeDataEntered)) %>%
  group_by(StudyId) %>%
  slice(1) %>%
  ungroup()

# Filter all records that have a StudyIdStr and AnnotatorIdStr that is in the reconciled_study_annotator_pairs to get final set of reconciled data
reconciled_records_unique <- LSR2_reconciled %>%
  semi_join(recent_reconciledID, by = c("StudyId", "InvestigatorId"))

# Rearrange rows to be grouped by studies and then experiment within studies. Reorder columns for readability
reconciled_records <- reconciled_records_unique %>%
  arrange(StudyId, ExperimentID, CohortId, OutcomeId, InterventionID) %>% 
  relocate(ExperimentLabel, .after = ExperimentID) %>% 
  relocate(CohortLabel, .after = CohortId) %>% 
  relocate(OutcomeLabel, .after = OutcomeId) %>% 
  relocate(c(OutcomeLabel, OutcomeId), .after = ExperimentLabel) %>% 
  relocate(c(InterventionLabel), .after = InterventionID)

savename <- paste0('reconciled_records_',Sys.Date(),'.csv')
write.csv(reconciled_records, savename)

# Split any column which contains multiple responses (separated by ';') into seperate columns
reconciled_records <- split_columns(reconciled_records)

# Label types of cohorts
# Aim = Make each experiment one comparison (with sham where present, control and intervention)
## Differentiate between positive and negative controls and TAAR1KO's

reconciled_cohort_label <- reconciled_records %>%
  mutate(
    Treatment1Type = case_when(
      grepl("fluoxetine control", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("Vehicle", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("fluoxetine", reconciled_records$`InterventionLabel[1]`) ~ "Cointervention",
      grepl("Fluoxetine", reconciled_records$`InterventionLabel[1]`) ~ "Cointervention",
      grepl("fluoxetine control", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("Vehicle", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("sham exercise", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("sed", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("Sed", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("SED", reconciled_records$`InterventionLabel[1]`) ~ "Control",
      grepl("HIIT", reconciled_records$`InterventionLabel[1]`) ~ "Intervention",
      grepl("MIIT", reconciled_records$`InterventionLabel[1]`) ~ "Intervention",
      grepl("ex", reconciled_records$`InterventionLabel[1]`) ~ "Intervention",
      grepl("exercise", reconciled_records$`InterventionLabel[1]`) ~ "Intervention",
      grepl("Exercise", reconciled_records$`InterventionLabel[1]`) ~ "Intervention",
      grepl("sham ovariectomy", reconciled_records$`InterventionLabel[1]`) ~ "Ovariectomy control",
      grepl("ovariectomy", reconciled_records$`InterventionLabel[1]`) ~ "Ovariectomy",
      grepl("treadmill", reconciled_records$`InterventionLabel[1]`) ~ "Intervention",
      TRUE ~ NA),
    Treatment2Type = case_when(
      grepl("fluoxetine control", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("Vehicle", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("fluoxetine", reconciled_records$`InterventionLabel[2]`) ~ "Cointervention",
      grepl("Fluoxetine", reconciled_records$`InterventionLabel[2]`) ~ "Cointervention",
      grepl("fluoxetine control", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("Vehicle", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("sham exercise", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("sed", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("Sed", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("SED", reconciled_records$`InterventionLabel[2]`) ~ "Control",
      grepl("HIIT", reconciled_records$`InterventionLabel[2]`) ~ "Intervention",
      grepl("MIIT", reconciled_records$`InterventionLabel[2]`) ~ "Intervention",
      grepl("ex", reconciled_records$`InterventionLabel[2]`) ~ "Intervention",
      grepl("exercise", reconciled_records$`InterventionLabel[2]`) ~ "Intervention",
      grepl("Exercise", reconciled_records$`InterventionLabel[2]`) ~ "Intervention",
      grepl("sham ovariectomy", reconciled_records$`InterventionLabel[2]`) ~ "Ovariectomy control",
      grepl("ovariectomy", reconciled_records$`InterventionLabel[2]`) ~ "Ovariectomy",
      grepl("treadmill", reconciled_records$`InterventionLabel[2]`) ~ "Intervention",
      TRUE ~ NA),
  )

## sort the blank disease models (which are equivalent to sham)
reconciled_cohort_label <- reconciled_cohort_label %>%
  mutate(
    IsDiseaseModelControl = case_when(
      IsDiseaseModelControl == 'True' ~ TRUE,
      IsDiseaseModelControl == 'False' ~ FALSE,
      is.na(IsDiseaseModelControl) ~ TRUE,
      TRUE ~ as.logical(NA)  # default case if none of the conditions are met
    )
  )



## Make CohortType column
# Combination interventions are interventions where currently licensed treatment is an intervention

reconciled_cohort_type <- reconciled_cohort_label %>%
  mutate(CohortType = case_when(
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Control"| is.na(Treatment1Type)) & (Treatment2Type == "Control"| is.na(Treatment2Type))) ~ "Negative control",
    
    IsDiseaseModelControl == FALSE & 
      (((Treatment1Type == "Control"| is.na(Treatment1Type)) & Treatment2Type == "Cointervention") |
         ((Treatment2Type == "Control"| is.na(Treatment2Type)) & Treatment1Type == "Cointervention")) ~ "Cointervention",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Intervention" & (Treatment2Type == "Control"| is.na(Treatment2Type))) |
         (Treatment2Type == "Intervention" & (Treatment1Type == "Control"| is.na(Treatment1Type)))) ~ "Simple intervention",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Intervention" & Treatment2Type == "Cointervention") | 
         (Treatment2Type == "Intervention" & Treatment1Type == "Cointervention")) ~ "Combination intervention",

    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Intervention" & Treatment2Type == "Ovariectomy control") | 
         (Treatment2Type == "Intervention" & Treatment1Type == "Ovariectomy control")) ~ "Intervention with sham ovariectomy",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Intervention" & Treatment2Type == "Ovariectomy") | 
         (Treatment2Type == "Intervention" & Treatment1Type == "Ovariectomy")) ~ "Intervention with ovariectomy",
    
    IsDiseaseModelControl == FALSE & 
      (((Treatment1Type == "Control"| is.na(Treatment1Type)) & Treatment2Type == "Ovariectomy control") |
         ((Treatment2Type == "Control"| is.na(Treatment2Type)) & Treatment1Type == "Ovariectomy control")) ~ "Control with sham ovariectomy",
    
    IsDiseaseModelControl == FALSE & 
      (((Treatment1Type == "Control"| is.na(Treatment1Type)) & Treatment2Type == "Ovariectomy") |
         ((Treatment2Type == "Control"| is.na(Treatment2Type)) & Treatment1Type == "Ovariectomy")) ~ "Control with ovariectomy",
    
    IsDiseaseModelControl == TRUE & 
      ((Treatment1Type == "Intervention" & Treatment2Type == "Cointervention") | 
         (Treatment2Type == "Intervention" & Treatment1Type == "Cointervention")) ~ "Sham with combination",
    
   IsDiseaseModelControl == TRUE & 
      ((Treatment1Type == "Intervention" & (Treatment2Type == "Control"| is.na(Treatment2Type))) |
         (Treatment2Type == "Intervention" & (Treatment1Type == "Control"| is.na(Treatment1Type)))) ~ "Sham with intervention",
 
    IsDiseaseModelControl == TRUE & 
      (((Treatment1Type == "Control"| is.na(Treatment1Type)) & Treatment2Type == "Cointervention") |
         ((Treatment2Type == "Control"| is.na(Treatment2Type)) & Treatment1Type == "Cointervention")) ~ "Sham with cointervention",
 
   IsDiseaseModelControl == TRUE & 
     (((Treatment1Type == "Control"| is.na(Treatment1Type)) & Treatment2Type == "Ovariectomy control") |
        ((Treatment2Type == "Control"| is.na(Treatment2Type)) & Treatment1Type == "Ovariectomy control")) ~ "Sham with sham ovariectomy",
   
   IsDiseaseModelControl == TRUE & 
     (((Treatment1Type == "Control"| is.na(Treatment1Type)) & Treatment2Type == "Ovariectomy") |
        ((Treatment2Type == "Control"| is.na(Treatment2Type)) & Treatment1Type == "Ovariectomy")) ~ "Sham with ovariectomy",
   
   IsDiseaseModelControl == TRUE & 
     (((Treatment1Type == "Intervention") & Treatment2Type == "Ovariectomy control") |
        ((Treatment2Type == "Intervention") & Treatment1Type == "Ovariectomy control")) ~ "Sham with intervention and sham ovariectomy",
   
   IsDiseaseModelControl == TRUE & 
     (((Treatment1Type == "Intervention") & Treatment2Type == "Ovariectomy") |
        ((Treatment2Type == "Intervention") & Treatment1Type == "Ovariectomy")) ~ "Sham with intervention and ovariectomy",
    
    IsDiseaseModelControl == TRUE ~ "Sham"
      )) %>% 
  relocate(CohortType, .after = `InterventionID[2]`) %>% 
  relocate(c(Treatment1Type, Treatment2Type), .after = `InterventionLabel[2]`) %>% 
  relocate(IsDiseaseModelControl, .after = ModelID)

reconciled_cohort_type <- reconciled_cohort_type %>%
  mutate(GroupID = interaction(StudyId, ExperimentID, OutcomeId, TimeInMinute))


## Label how positive controls are being used in positive control v TAAR1 Ag experiments
# If for an observation with the value "Positive control" in the CohortType column, there is an observation 
# with the CohortType "Combination intervention" with the same combination of values in the StudyId, ExperimentID, 
# OutcomeId, TimeInMinute variables, then the value in the CohortType column needs changed to "Control for 
# combination intervention". Similarly, if for an observation with the value "Positive control treated sham" 
# in the CohortType column, there is an observation with the CohortType "Combination intervention" with 
# the same combination of values in the StudyId, ExperimentID, OutcomeId, TimeInMinute variables, then the 
# value in the CohortType column needs changed to "Sham for combination intervention"

## Identify experiments that are just 'Simple' comparison, those which are 'Combination' comparisons and those which have 'TAAR1KO' involved
reconciled_comparison_type <- reconciled_cohort_type %>%
  mutate(ExperimentType = case_when(
    str_detect(CohortType, "combination") | str_detect(CohortType, "Combination") ~ "Exercise:Drug combination (ED v C)",
    str_detect(CohortType, "simple") | str_detect(CohortType, "Simple") ~ "Simple intervention (E v C)",
    str_detect(CohortType, "cointervention") | str_detect(CohortType, "Cointervention") ~ "Drug (D v C)"
  ))

reconciled_cohort_role <- reconciled_comparison_type %>% 
  mutate(RoleOfCohort = case_when(
    str_detect(CohortType, "Sham") | str_detect(CohortType, "sham") ~ "S", 
    str_detect(CohortType, "Control") | str_detect(CohortType, "control") ~ "C",
    str_detect(CohortType, "Intervention") | str_detect(CohortType, "intervention") ~ "I",
    str_detect(CohortType, "Cointervention") | str_detect(CohortType, "cointervention") ~ "P"
  ))

columnnmae13 <- "TimeInMinute"
condition <- (reconciled_cohort_role$StudyId == '0384b8a0-890a-49f3-80e2-ccd95416e901' & reconciled_cohort_role$TimeInMinute == -11)
reconciled_cohort_role[condition, columnnmae13] <- 1

## Wrangle wide so each observations is a single comparison
data <- reconciled_cohort_role
data <- data %>%
  mutate(GroupID = interaction(StudyId, ExperimentID, OutcomeId, TimeInMinute))

#### SyRF correction section ends


##### Remove DA knockouts and data reported from subgroups #####
data_rem <- subset(data, data$`DiseaseModelLabel(s)` == 'DAT +/-')
data <- anti_join(data, data_rem)
#remove data reported from subgroups
data_rem <- data %>% filter(str_detect(`DiseaseModelLabel(s)`, "SUBGROUP") | str_detect(CohortLabel, 'SUBGROUP'))
data <- anti_join(data, data_rem)
data <- data %>%
  mutate(`DiseaseModelLabel(s)` = ifelse(`DiseaseModelLabel(s)` == "DAT -/-", "DAT KO", `DiseaseModelLabel(s)`))
data <- data %>% mutate_all(trimws)


##### Extract treatment names and doses ####
data <- data%>%
  mutate(drugname1 = case_when(
    grepl("fluoxetine control", data$`InterventionLabel[1]`) ~ "control",
    grepl("fluoxetine", data$`InterventionLabel[1]`) ~ "fluoxetine",
    grepl("Fluoxetine", data$`InterventionLabel[1]`) ~ "fluoxetine",
    grepl("sham exercise", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("sed", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("Sed", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("SED", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("HIIT", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("MIIT", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("ex", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("exercise", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("Exercise", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("sham ovariectomy", data$`InterventionLabel[1]`) ~ "sham ovariectomy",
    grepl("ovariectomy", data$`InterventionLabel[1]`) ~ "ovariectomy",
    grepl("treadmill", data$`InterventionLabel[1]`) ~ "exercise",
    TRUE ~ "Other"
  ))
data <- data%>%
  mutate(drugname2 = case_when(
    grepl("fluoxetine control", data$`InterventionLabel[1]`) ~ "control",
    grepl("fluoxetine", data$`InterventionLabel[1]`) ~ "fluoxetine",
    grepl("Fluoxetine", data$`InterventionLabel[1]`) ~ "fluoxetine",
    grepl("sham exercise", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("sed", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("Sed", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("SED", data$`InterventionLabel[1]`) ~ "sedentary",
    grepl("HIIT", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("MIIT", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("ex", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("exercise", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("Exercise", data$`InterventionLabel[1]`) ~ "exercise",
    grepl("sham ovariectomy", data$`InterventionLabel[1]`) ~ "sham ovariectomy",
    grepl("ovariectomy", data$`InterventionLabel[1]`) ~ "ovariectomy",
    grepl("treadmill", data$`InterventionLabel[1]`) ~ "exercise",
    TRUE ~ "Other"
  ))

data$Treatment1Label <- paste0(data$drugname1, ", ", data$`Duration of intervention (total)[1]`, " ", data$`Unit of measurement (exercise duration)[1]`)
data$Treatment2Label <- paste0(data$drugname2, ", ", data$`Duration of intervention (total)[2]`, " ", data$`Unit of measurement (exercise duration)[2]`)

#diag <- data[,c('drugname1', 'Treatment1Label','InterventionLabel[1]','Dose of treatment used:[1]','Dose of positive control treatment?[1]',
#                'Measurement unit of treatment dose:[1]','drugname2','Treatment2Label','InterventionLabel[2]',
#                'Dose of treatment used:[1]','Dose of positive control treatment?[2]','Measurement unit of treatment dose:[2]', 'GroupID','CohortId')]

#dsubset_df1 <- subset(diag, !grepl('mg', diag$Treatment1Label, ignore.case = TRUE))
#dsubset_df1a <- subset(dsubset_df1, !grepl('Other', dsubset_df1$Treatment1Label, ignore.case = TRUE))

#dsubset_df2 <- subset(diag, !grepl('mg', diag$Treatment2Label, ignore.case = TRUE))
#dsubset_df2a <- subset(dsubset_df2, !grepl('Other', dsubset_df2$Treatment2Label, ignore.case = TRUE))

#dsubset_df1b <- subset(dsubset_df1a, !grepl('Other', dsubset_df1a$Treatment2Label, ignore.case = TRUE))


#column_name <- 'Measurement unit of treatment dose:[2]'
#condition <- data$GroupID %in% dsubset_df2a$GroupID & data$CohortId %in% dsubset_df2a$CohortId
#data[condition, column_name] <- 'mg/kg'


#column_name <- 'Measurement unit of treatment dose:[1]'
#condition <- data$GroupID %in% dsubset_df1a$GroupID & data$CohortId %in% dsubset_df1a$CohortId
#data[condition, column_name] <- 'mg/kg'

#condition <- data$GroupID %in% dsubset_df1b$GroupID & data$CohortId %in% dsubset_df1b$CohortId

#column_name1 <- 'drugname1'
#column_name2 <- 'drugname2'
#column_name3 <- 'Treatment1Label'
#column_name4 <- 'Treatment2Label'
#column_name5 <- 'InterventionLabel[1]'
#column_name6 <- 'InterventionLabel[2]'
#column_name7 <- 'Dose of treatment used:[1]'
#column_name8 <- 'Dose of treatment used:[2]'
#column_name9 <- 'Dose of positive control treatment?[1]'
#column_name10 <- 'Dose of positive control treatment?[2]'
#column_name11 <- 'Measurement unit of treatment dose:[1]'
#column_name12 <- 'Measurement unit of treatment dose:[2]'

#temp <- data[condition, column_name1]
#data[condition, column_name1] <- data[condition, column_name2]
#data[condition, column_name2] <- temp

#temp <- data[condition, column_name3]
#data[condition, column_name3] <- data[condition, column_name4]
#data[condition, column_name4] <- temp

#temp <- data[condition, column_name5]
#data[condition, column_name5] <- data[condition, column_name6]
#data[condition, column_name6] <- temp

#temp <- data[condition, column_name7]
#data[condition, column_name7] <- data[condition, column_name8]
#data[condition, column_name8] <- temp

#temp <- data[condition, column_name9]
#data[condition, column_name10] <- temp
#data[condition, column_name9] <- data[condition, column_name10]

#temp <- data[condition, column_name11]
#data[condition, column_name11] <- data[condition, column_name12]
#data[condition, column_name12] <- temp

columnnmae13 <- "TimeInMinute"
condition <- (data$StudyId == '0384b8a0-890a-49f3-80e2-ccd95416e901' & data$TimeInMinute == '-11')
data[condition, columnnmae13] <- '1'

#condition <- data$ExperimentID == '40a5354c-f200-41e7-868e-2f9dbcfc2424'
#data[condition, columnnmae13] <- "C57Bl/6Jx129Sv/J (mouse)"

#condition <- data$ExperimentID == '31c009a0-cb80-4457-b529-4e0c52a97b02'
#data[condition, columnnmae13] <- "Not stated (mouse)"

#condition <- data$ExperimentID == '459641de-3ba8-4ac6-93c2-b8f7805af682'
#data[condition, columnnmae13] <- "Not stated (mouse)"

#condition <- data$StudyId == 'c064173a-747d-4877-ae38-27415dddd81e'
#data[condition, columnnmae13] <- "ICR (mouse)"

#condition <- data$ExperimentID == '1c48a8d4-0a37-4256-aa4a-34108731a48b'
#data[condition, columnnmae13] <- "NMRI (mouse)"

#condition <- data$ExperimentID == '32c6ccc3-e4e2-4a41-98c2-f36e8984775e'
#data[condition, columnnmae13] <- "NMRI (mouse)"


##### Get names of each cohort as drug, dose, unit ####
data <- data %>% 
  mutate(TreatmentLabel1 = case_when(
    grepl("Intervention", data$Treatment1Type) ~ paste0(drugname1, ", ", `Duration of intervention (total)[1]`, ' ', `Unit of measurement (exercise duration)[1]`),
  ))

data <- data %>% 
  mutate(TreatmentLabel2 = case_when(
    grepl("Intervention", data$Treatment2Type) ~ paste0(drugname2, ", ", `Duration of intervention (total)[2]`, ' ', `Unit of measurement (exercise duration)[2]`),
  ))

#### for table

data <- data %>% 
  mutate(DrugLabel1 = case_when(
    grepl("Intervention", data$Treatment1Type) ~ drugname1,
  ))

data <- data %>% 
  mutate(DrugLabel2 = case_when(
    grepl("Intervention", data$Treatment2Type) ~ drugname2, 
  ))


data <- data %>% 
  mutate(
    TreatmentLabel = case_when(
      !is.na(TreatmentLabel1) & !is.na(TreatmentLabel2) ~ paste(TreatmentLabel1, "&", TreatmentLabel2),
      !is.na(TreatmentLabel1) ~ TreatmentLabel1,
      !is.na(TreatmentLabel2) ~ TreatmentLabel2,
      TRUE ~ NA_character_
    )
  )

data <- data %>% 
  mutate(
    DrugLabel = case_when(
      !is.na(DrugLabel1) & !is.na(DrugLabel2) ~ paste(DrugLabel1, "&", DrugLabel2),
      !is.na(DrugLabel1) ~ DrugLabel1,
      !is.na(DrugLabel2) ~ DrugLabel2,
      TRUE ~ NA_character_
    )
  )


###for each cohort, label sham and control groups, along with TAAR1KO -ve control Contol for combination interventions (== positive control)
### first, we need to work out the attributes of each group - wgat types of comparisons will they allow?
# Assuming "df" is your data frame with columns "groupId" and "cohortType"

##### Count occurrences for each cohort type within each group #####
group_characteristics <- data %>%
  group_by(GroupID, CohortType) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = CohortType, values_from = count, values_fill = 0) %>%
  ungroup()
group_characteristics$total <- rowSums(select(group_characteristics, -GroupID))
group_characteristics$GroupID <- as.character(group_characteristics$GroupID)
data$GroupID <- as.character(data$GroupID)

##### Wrangling experiment type and outcome data #####
##### 3.1.1 E v Cont  with sham #####
data_TvCs <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Negative control"]),
    Sham = unique(cohset$CohortLabel [cohset$CohortType == "Sham"])
  )
  
  if (nrow(combinations) > 0) {
    cohset <- combinations %>%
      left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%
      left_join(cohset %>% filter(CohortType == "Negative control"), by = c("Control" = "CohortLabel")) %>%
      left_join(cohset %>% filter(CohortType == "Sham"), by = c("Sham" = "CohortLabel"))
    
    if (nrow(cohset) > 0) {
      data_TvCs <- bind_rows(data_TvCs, cohset)
    }}}

data_TvCs$Label <- data_TvCs$TreatmentLabel.x
data_TvCs$SortLabel <- "TvC"
data_TvCs <- subset(data_TvCs, !is.na(data_TvCs$Intervention) & !is.na(data_TvCs$Control) & !is.na(data_TvCs$Sham))

##### 3.1.2 E v Cont  without sham ####

data_TvC <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Negative control"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Negative control"), by = c("Control" = "CohortLabel"))
  data_TvC <- bind_rows(data_TvC, cohset)
}
data_TvC$Label <- data_TvC$TreatmentLabel.x
data_TvC$SortLabel <- "TvC"
data_TvC <- subset(data_TvC, !is.na(data_TvC$Intervention) & !is.na(data_TvC$Control))

data_sub_TvC <- subset(data_TvC, !data_TvC$GroupID.x %in% data_TvCs$GroupID.x)
data_TvCc <- bind_rows(data_TvCs, data_sub_TvC)

##### 3.2.1 Cont v S  : effect of model ####

data_CvS <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Sham"]),
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Negative control"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Negative control"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Sham"), by = c("Control" = "CohortLabel"))
  data_CvS <- bind_rows(data_CvS, cohset)
}
data_CvS$Label <- data_CvS$TreatmentLabel.x
data_CvS$SortLabel <- "CvS"
data_CvS <- subset(data_CvS, !is.na(data_CvS$Intervention) & !is.na(data_CvS$Control))



##### 3.3.1 E v C with sham #####
data_AvCs <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Cointervention"]),
    Sham = unique(cohset$CohortLabel [cohset$CohortType == "Sham"])
  )
  if (nrow(combinations) > 0) {
    cohset <- combinations %>%
      left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%
      left_join(cohset %>% filter(CohortType == "Cointervention"), by = c("Control" = "CohortLabel")) %>%
      left_join(cohset %>% filter(CohortType == "Sham"), by = c("Sham" = "CohortLabel"))
    
    if (nrow(cohset) > 0) {
    data_AvCs <- bind_rows(data_AvCs, cohset)
    }}}

data_AvCs$Label <- data_AvCs$TreatmentLabel.x
data_AvCs$SortLabel <- "AvC"
data_AvCs <- subset(data_AvCs, !is.na(data_AvCs$Intervention) & !is.na(data_AvCs$Control) & !is.na(data_AvCs$Sham))

##### 3.3.2 E v C without sham #####
data_AvC <- data.frame()


for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Cointervention"])

  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%

    left_join(cohset %>% filter(CohortType == "Cointervention"), by = c("Control" = "CohortLabel"))
  
  data_AvC <- bind_rows(data_AvC, cohset)
}
data_AvC$Label <- data_AvC$TreatmentLabel.x
data_AvC$SortLabel <- "AvC"
data_AvC <- subset(data_AvC, !is.na(data_AvC$Intervention) & !is.na(data_AvC$Control))

data_sub_AvC <- subset(data_AvC, !data_AvC$GroupID.x %in% data_AvCs$GroupID.x)
data_AvCc <- bind_rows(data_AvCs, data_sub_AvC)


##### combine to a single df and merge back to main #####
data2 <- bind_rows(data_TvCc, data_AvCc, data_CvS)

data2 <- data2 %>%
  rename_at(vars(ends_with(".x")), ~str_remove(., "\\.x") %>% paste0("_I")) %>%
  rename_at(vars(ends_with(".y")), ~str_remove(., "\\.y") %>% paste0("_C"))

#####colname corrections #####
data2$`Type of outcome` <- data2$`Type of outcome:_I`
data2$StudyId <- data2$StudyId_I
data2$ErrorType <- data2$ErrorType_I
data2$GreaterIsWorse <- data2$GreaterIsWorse_I
data2$`Measurement unit of treatment dose:` <- data2$`Exercise intensity units:[1]_I`
data2$`Duration of treatment` <- data2$`Duration of intervention (total)[1]_I`
data2$`Unit of measurement for treatment duration` <- data2$`Unit of measurement (exercise duration)[1]_I`
data2$`Species of animals?` <- data2$`Species of animals?_I`
data2$`Animal strain?` <- data2$`Animal strain?_I`
data2$`Sex of animals?` <- data2$`Sex of animals?_I`
data2$`Were caregivers/investigator blinded to which intervention each animal received?` <- data2$`Were caregivers/investigator blinded to which intervention each animal received?_I`
data2$`Is any role of the funder in the design/analysis/reporting of study described?` <- data2$`Is any role of the funder in the design/analysis/reporting of study described?_I`
data2$Title <- data2$Title_I
data2$Year <- data2$Year_I
data2$ModelID <- data2$ModelID_I
data2$Authors <- data2$Authors_I
data2$VoF <- data2$`Voluntary or forced exercise?[1]_I`


##### get standard placement - Single I #####
data2$F_C_L <- data2$Control
data2$F_C_n <- as.numeric(data2$NumberOfAnimals_C)
data2$F_C_m <- as.numeric(data2$OutcomeResult_C)
data2$F_C_v <- as.numeric(data2$OutcomeError_C)

data2$F_T_L <- data2$Intervention
data2$F_T_n <- as.numeric(data2$NumberOfAnimals_I)
data2$F_T_m <- as.numeric(data2$OutcomeResult_I)
data2$F_T_v <- as.numeric(data2$OutcomeError_I)

data2$F_S_L <- data2$Sham
data2$F_S_n <- as.numeric(data2$NumberOfAnimals)
data2$F_S_m <- as.numeric(data2$OutcomeResult)
data2$F_S_v <- as.numeric(data2$OutcomeError)

data_all_F <- data2


savename_all <- paste0(LSR,'data_all_',Sys.Date(),'.csv')

write_csv(data_all_F, savename_all)

###### Calculate effect size for simple interventions #####

# Data are ready to calculate effect sizes with numerical data for one comparison (test/control) per line

# Note: current code doesn't have cohort level questions split into treatment and control annotations per line
# Not an issue currently since control and treatment cohorts have had the same characteristics (sex, strain, etc) so far, but probably the subject of future edits (e.g. using pivot_wider function)



## Read in data and edit outcome label 
#pass though better
# dataall <- read_csv("dataall.csv")
dataall <- data_all_F

dataall <- dataall %>% 
  mutate(OutcomeType = case_when(str_detect(OutcomeLabel_I, regex("elevated|epm", ignore_case = TRUE)) ~ "Elevated plus maze",
                                 str_detect(OutcomeLabel_I, "BDNF") ~ "BDNF", 
                                 TRUE ~ `Type of outcome:_I`)) %>% 
  relocate(OutcomeType, .after = OutcomeLabel)

outcomeFrequencies <- dataall %>% group_by(OutcomeType, OutcomeLabel_I, `Type of outcome:_I`) %>% count()

savename_of <- paste0(LSR,'_OutcomeFrequencies',Sys.Date(),'.csv')

write.csv(outcomeFrequencies, savename_of)
## 1. Calculate SD for all comparisons 
#F = final
#C/T/S = control/treatment/sham
#L/n/m/v = cohort label/n in cohort/mean/variance

dataall <- dataall %>% 
  mutate(F_C_v.SD = case_when(ErrorType == "IQR" ~ (F_C_v/1.35), 
                              ErrorType == "SD" ~ F_C_v, 
                              ErrorType == "SEM" ~ sqrt(F_C_n)*F_C_v)) %>%  
  relocate(F_C_v.SD, .after = F_C_v) %>% 
  mutate(F_T_v.SD = case_when(ErrorType == "IQR" ~ (F_T_v/1.35), 
                              ErrorType == "SD" ~ F_T_v, 
                              ErrorType == "SEM" ~ sqrt(F_T_n)*F_T_v)) %>%  
  relocate(F_T_v.SD, .after = F_T_v) %>% 
  mutate(F_S_v.SD = case_when(ErrorType == "IQR" ~ (F_S_v/1.35), 
                              ErrorType == "SD" ~ F_S_v, 
                              ErrorType == "SEM" ~ sqrt(F_S_n)*F_S_v)) %>%  
  relocate(F_S_v.SD, .after = F_S_v)

# Check number of comparisons where NMD can be calculated

dataall <- dataall %>%
  rowwise() %>%
  mutate(`NMD_possible` = all(!is.na(F_S_L)))

## 2. Calculate effect sizes (SMD for all, NMD where possible)

### Calculate true n for control groups (n'c)

# Step 1: Number of groups served by control group
#F_C_L_frequencies <- dataall %>%
#  group_by(StudyId, OutcomeId, F_C_L) %>%
#  summarise(Frequency_FCL = n()) %>% 
#  ungroup()

# Step 2: Join the frequencies back to the original dataframe
#dataall <- dataall %>%
#  left_join(F_C_L_frequencies)

# Step 3: Divide F_C_n by the frequency count
#dataall <- dataall %>%
#  mutate(F_C_n_true = F_C_n / Frequency_FCL)

### SMD

#Hedges g to account for small sample sizes (default for SMD when using the escalc() function) - Hedge’s g (statistically corrects for variance that may be introduced when sample sizes are small (Larry V. Hedges 1981))
# m1i = control, m2i = rx 

# Hedge's g effect size
SMD_data_all.nottrue <- escalc(
  measure = "SMD", 
  m1i = dataall$F_C_m, 
  m2i = dataall$F_T_m, 
  sd1i = dataall$F_C_v.SD, 
  sd2i = dataall$F_T_v.SD, 
  n1i = dataall$F_C_n, 
  n2i = dataall$F_T_n, 
  data = dataall) %>% 
  select(yi, vi)

dataall$SMD <- SMD_data_all.nottrue$yi
dataall$SMDv <- SMD_data_all.nottrue$vi

#escalc (m1 - m2) = (control - treatment)
# so if greater is better, then *-1

#dataall$SMD_true <- SMD_data_all.true$yi #Row 112, SMD calculated, but SMD_true not calculated
#dataall$SMDv_true <- SMD_data_all.true$vi

### NMD

# Assume that treatments are closer to shams than controls are 
# So C-S > T-S 


dataall <- dataall %>% 
  mutate(`NMD` = 100*(((F_C_m - F_S_m) - (F_T_m - F_S_m))/(F_C_m - F_S_m))) %>% 
  mutate(`NMD_SDc*` = 100*((F_C_v.SD/(F_C_m - F_S_m)))) %>% 
  mutate(`NMD_SDrx*` = 100*((F_T_v.SD/(F_C_m - F_S_m)))) %>% 
  mutate(NMDv = sqrt(((`NMD_SDc*`)^2/F_C_n) + ((`NMD_SDrx*`)^2/F_T_n))) 

dataall.direction <- dataall %>% 
  mutate(SMD = if_else((GreaterIsWorse == "FALSE"), -1*SMD, SMD))

diagnostic <- dataall.direction %>% 
  select(F_C_m, F_T_m, F_S_m, SMD, NMD, GreaterIsWorse) %>% 
  mutate(CbiggerT = if_else(F_C_m > F_T_m, "Yes", "No"))

# nicely named columns for subgroup analysis

df <- dataall.direction


### Duration of treatment (categorical)

# Create variable standardised to Weeks

df$`Duration of intervention (total)[1]_I` <- as.numeric(df$`Duration of intervention (total)[1]_I`)
df$`Duration of intervention (total)[2]_I` <- as.numeric(df$`Duration of intervention (total)[2]_I`)
df$`Number of sessions:[1]_I` <- as.numeric(df$`Number of sessions:[1]_I`)
df$`Number of sessions:[2]_I` <- as.numeric(df$`Number of sessions:[2]_I`)
df$`Exercise intensity:[1]_I` <- as.numeric(df$`Exercise intensity:[1]_I`)
df$`Exercise intensity:[2]_I` <- as.numeric(df$`Exercise intensity:[2]_I`)
df$`Exercise individual session duration:[1]_I` <- as.numeric(df$`Exercise individual session duration:[1]_I`)
df$`Exercise individual session duration:[2]_I` <- as.numeric(df$`Exercise individual session duration:[2]_I`)

df <- df %>% 
  mutate(TotalExercise1 = case_when(`Exercise intensity units:[1]_I` == "(speed) centimeters per second (cm/s)" ~ `Number of sessions:[1]_I` * `Exercise intensity:[1]_I` * 0.6 * `Exercise individual session duration:[1]_I` / 1000,
                                    `Exercise intensity units:[1]_I` == "(speed) meters per minute (m/min)"  ~ `Number of sessions:[1]_I` * `Exercise intensity:[1]_I`* `Exercise individual session duration:[1]_I`/1000))

df <- df %>% 
  mutate(TotalExercise2 = case_when(`Exercise intensity units:[2]_I` == "(speed) centimeters per second (cm/s)" ~ `Number of sessions:[2]_I` * `Exercise intensity:[2]_I` * 0.6* `Exercise individual session duration:[2]_I` / 1000,
                                    `Exercise intensity units:[2]_I` == "(speed) meters per minute (m/min)"  ~ `Number of sessions:[2]_I` * `Exercise intensity:[1]_I`* `Exercise individual session duration:[2]_I` / 1000))
                                                                      
df <- df %>% 
  mutate(DurationOfTreatmentWeeks1 = case_when(`Unit of measurement (exercise duration)[1]_I` == "Days" ~ `Duration of intervention (total)[1]_I`/7, 
                                              `Unit of measurement (exercise duration)[1]_I` == "Months" ~ `Duration of intervention (total)[1]_I`*4.345, 
                                              `Unit of measurement (exercise duration)[1]_I` == "Weeks" ~ `Duration of intervention (total)[1]_I`, 
                                              `Unit of measurement (exercise duration)[1]_I` == "Single dose" ~ 1/7))
df <- df %>% 
  mutate(DurationOfTreatmentWeeks2 = case_when(`Unit of measurement (exercise duration)[2]_I` == "Days" ~ `Duration of intervention (total)[2]_I`/7, 
                                              `Unit of measurement (exercise duration)[2]_I` == "Months" ~ `Duration of intervention (total)[2]_I`*4.345, 
                                              `Unit of measurement (exercise duration)[2]_I` == "Weeks" ~ `Duration of intervention (total)[2]_I`, 
                                              `Unit of measurement (exercise duration)[2]_I` == "Single dose" ~ 1/7))

df <- df %>% 
  mutate(EI1 = case_when(`Exercise intensity units:[1]_I` == "(speed) centimeters per second (cm/s)" ~ `Exercise intensity:[1]_I` * 0.6,
                                    `Exercise intensity units:[1]_I` == "(speed) meters per minute (m/min)"  ~ `Exercise intensity:[1]_I`))

df <- df %>% 
  mutate(EI2 = case_when(`Exercise intensity units:[2]_I` == "(speed) centimeters per second (cm/s)" ~ `Exercise intensity:[2]_I` * 0.6,
                                    `Exercise intensity units:[2]_I` == "(speed) meters per minute (m/min)"  ~ `Exercise intensity:[1]_I`))


# Create categorical variable for Duration of treatment. Grouped into up to a week, between a week and 4 weeks, more than 4 weeks
df <- df %>% 
  mutate(TreatmentDurationCategory = case_when(DurationOfTreatmentWeeks1 <= 1 ~ "Less than 1 week", 
                                               DurationOfTreatmentWeeks1 > 1 & DurationOfTreatmentWeeks1 < 4 ~ "Between 1-4 weeks", 
                                               DurationOfTreatmentWeeks1 >= 4 ~ "More than 4 weeks"))




### Prophylactic or therapeutic

df$`What was the timing of the initiation of treatment[1]_I` <- as.numeric(df$`What was the timing of the initiation of treatment[1]_I`)
df <- df %>% 
  mutate(ProphylacticOrTherapeutic = case_when (`What was the timing of the initiation of treatment[1]_I` > 0  ~ "Therapeutic", 
                                                `What was the timing of the initiation of treatment[1]_I` <= 0  ~ "Prophylactic", 
                                                TRUE ~ NA))

df$`Exercise intensity:[1]_I` <- as.numeric(df$`Exercise intensity:[1]_I`)

df <- df %>% 
  mutate(IntofEx = case_when (`Exercise intensity units:[1]_I` == '(speed) centimeters per second (cm/s)' ~ (`Exercise intensity:[1]_I`*60/100), 
                                     `Exercise intensity units:[1]_I` == '(speed) meters per minute (m/min)'  ~ `Exercise intensity:[1]_I`,
                                     TRUE ~ NA))


## Give original variables better names for analysis in new columns - this is all for simple 

df <- df %>% 
  mutate(Species = `Species of animals?`, 
         Strain = `Animal strain?`,
         Sex = `Sex of animals?`, 
         DrugName = drugname1_I)



## Categorise by outcome type - requires checking with each iteration
df <-  df %>%
  mutate(outcome_type = case_when(
    `Type of outcome` == "Neurotransmitter (NT) or neurotransmitter modulator levels (e.g. dopamine, cannabinoids)" ~ "Neurotransmitter levels",
    `Type of outcome` == "Locomotor behaviour or arousal (e.g. social interaction test)"  ~ "Locomotor",
    `Type of outcome` == "Non-NT gene or protein expression/level (e.g. BDNF in plasma or western blot)" ~ "BDNF",
    `Type of outcome` == "Fear generalisation" ~ "Fear generalisation",
    `Type of outcome` == "Freezing behaviour"  ~ "Freezing",
    `Type of outcome` == "Other behavioural (enter in outcome label and in comments box)" ~ "Other behavioural",
    `Type of outcome` == "Other behavioural (enter in outcome label and in comments box)???" ~ "Freezing",
    `Type of outcome` == "Other neurobiological (enter in outcome label and in comments box)"  ~ "Other neurobiological",
    `Type of outcome` == "Fear memory (e.g. object recognition, episodic)"  ~ "Fear memory",
    `Type of outcome` == "Stress response"  ~ "Stress response",
  )) 



#RoB - first get rid of duplicates
df <-df[,-c(186:218, 346:378)]


df <- df %>%
  rename(`(RoB) Were caregivers/investigator blinded to which intervention each animal received?_I` = `Were caregivers/investigator blinded to which intervention each animal received?_I`)

# Calculate overall RoB score

df <- df %>%
  rowwise() %>%
  mutate(RoSBScoreAny = sum(c_across(contains("RoB")) == "Yes", na.rm = TRUE)) %>%
  ungroup()

df <- df %>%
  rowwise() %>%
  mutate(RoBScore = sum(c_across(contains("RoB")) == "Yes", na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(RoBScore = paste0(RoBScore, " criteria met"))

df$RoSBScoreAny <- as.numeric(df$RoSBScoreAny)

df$RoBTF <- 'No low RoB items'
for(i in 1:nrow(df)) {
  if(df[i, "RoSBScoreAny"] > 0) {
    df[i,'RoBTF'] <- 'Some low RoB items'
  }
}
##### For reporting quality subgroup analysis #####


df <- df %>% 
  rename(`(ARRIVE) Is any role of the funder in the design/analysis/reporting of the study described?_I` = `Is any role of the funder in the design/analysis/reporting of study described?_I`)

# Calculate overall ARRIVE score
df <- df %>%
  mutate(ARRIVEScore = rowSums(across(contains("ARRIVE"), ~ (.x == "Yes") | (.x == "NA (ethical approval declared)")), na.rm = TRUE)) %>% 
  mutate(ARRIVEScoreCat = case_when(ARRIVEScore <= 3 ~ "A: < 3 criteria met",
                                    ARRIVEScore > 3 & ARRIVEScore <= 7 ~ "B: 4-7 criteria met",
                                    ARRIVEScore > 7 & ARRIVEScore <= 11 ~ "C: 8-11 criteria met",
                                    ARRIVEScore > 11 & ARRIVEScore <= 15 ~ "D: 12-15 criteria met",
                                    ARRIVEScore > 15 & ARRIVEScore <= 19 ~ "E: 16-19 criteria met",
                                    ARRIVEScore > 19 ~ "F: > 20 criteria met"))



# remove studies with 999 as N or SD
df4 <- subset(df, !df$NumberOfAnimals_I == 999)
df4 <- subset(df4, !df4$NumberOfAnimals_C == 999)
df4 <- subset(df4, !df4$OutcomeError_C == 999)
df <- subset(df4, !df4$OutcomeError_I == 999)

df <- df %>% 
  mutate_all(~replace(., . == 999, NA))

df <- df %>%
  mutate(out_cat = case_when(outcome_type == 'Locomotor' ~ "Behavioural",
                             outcome_type == 'Fear memory' ~ "Behavioural",
                             outcome_type == 'Freezing' ~ "Behavioural",
                             outcome_type == 'Other behavioural' ~ "Behavioural",
                             outcome_type == 'BDNF' ~ "Neurobiological",
                             outcome_type == 'Stress response' ~ "Neurobiological",
                             outcome_type == 'Neurotransmitter levels' ~ "Neurobiological",
                             outcome_type == 'Other neurobiological' ~ "Neurobiological"))

# SAVE FILE
savefile_output <- paste0(LSR,'_','clean_data_',Sys.Date(),'.csv')
write.csv(df, savefile_output, row.names = FALSE)

