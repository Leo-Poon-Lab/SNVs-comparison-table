library(tidyverse)
library(readxl)

samples <- read_excel("../results/data_add.xlsx") %>% select(WHP2161:WHP4738_11861) %>% colnames()
samples <- sapply(samples, function(x){
	strsplit(x, "[-_]")[[1]][1]
})
data_gisaid <- read_csv("../../2021-05-05_GISAID_upload/results/gisaid_id_case_id.csv")

samples[samples %in% data_gisaid$lab_id]
data_gisaid %>% filter(lab_id %in% samples) %>% write_csv("../results/strain_submitted.csv")

# check vm no.
data_whp_vm <- read_csv("../data/WHP_VM_list_202103.csv")
data_hk_ann <- read_excel("../data/HK case annotation_20210902.xlsx")
data_metadata <- read_xlsx("../../2021-06-24_merge_metadata/results/cleaned_metadata.xlsx")

data_whp_vm <- data_whp_vm %>% select(`Lab No.`, Case, `Specimen ID`) %>% transmute(lab_id=`Lab No.`, case_id=Case, specimen_id=`Specimen ID`)
data_hk_ann <- data_hk_ann %>% select(WHP_number, `HK case`, `VM_no.`) %>% transmute(lab_id=WHP_number, case_id=`HK case`, specimen_id=`VM_no.`)

data_vm <- bind_rows(data_whp_vm, data_hk_ann) %>% unique() %>% arrange(as.numeric(case_id))
data_vm$report_year <- sapply(data_vm$case_id, function(x){
	tmp <- lubridate::year(data_metadata$`Report date`[data_metadata$case_id == x])
	tmp[1]
})

data_vm %>% filter(lab_id %in% samples[!samples %in% data_gisaid$lab_id]) %>% mutate(strain_name=paste0("hCoV-19/Hong Kong/", specimen_id, "/", report_year)) %>% write_csv("../results/strain_lack.csv")
