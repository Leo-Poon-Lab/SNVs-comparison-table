# conda activate base beforehand
library(tidyverse)
library(readxl)
library(Biostrings)

data <- read_excel("../data/genome compare table kh.xlsx")
ref_seq <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/NGS_data_input/reference.fasta")

# make sure the position schemes are consistent
check <- sapply(seq_len(nrow(data)), function(i){
	# i=1
	pos_t <- data$Position[i]
	subseq(ref_seq, pos_t, pos_t) == data$nt_ref[i]
})
stopifnot(all(check))

# WHP4529, HK_case 11774
# WHP4738, HK_case 11861
# seq_all <- readDNAStringSet("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/latest_genome_caseid.fasta")

files_tsv <- sapply(c("WHP4738", "WHP4529"), function(x){
	tmp <- list.files("../../../2020/2020-09-01_COVID_NGS_pipeline/COVID_NGS_pipeline_results_shared/variant_caller/ivar/", x, full.names = T)
	tmp <- tmp[!grepl("iseq", tmp)]
})
files_tsv <- unlist(files_tsv)

# build custom snpeff database
system("java -Xmx4g -jar ~/softwares/snpEff/snpEff.jar download NC_045512.2")
genome <- "NC_045512.2"
ref_seq_nc <- ref_seq
names(ref_seq_nc) <- genome
writeXStringSet(ref_seq_nc, "../data/refernece.fasta")

# source("./helper/prepare_vcf_mut.r")
df_ann <- lapply(seq_along(files_tsv), function(i){
	sample_t <- names(files_tsv)[i]
	file_tsv_i <- files_tsv[i]
	outfile <- paste0("../results/", sample_t, ".vcf")
	
	system(paste0("./helper/ivar_variants_to_vcf.py ", file_tsv_i, " ", outfile))
	tmp <- readLines(outfile)
	writeLines(gsub("MN908947_3", genome, tmp), outfile)

	outfile_snpeff <- paste0(outfile, ".snpeff")
	outfile_snpeff_csq <- paste0(outfile, ".snpeff.csq")
	outfile_csv <- paste0(outfile, ".csv")
	system(paste0("java -jar ~/softwares/snpEff/snpEff.jar ", genome, " ", outfile, " > ", outfile_snpeff))
	# system(paste0("bcftools csq --force --phase a -f ../data/refernece.fasta -g ../data/Sars_cov_2.ASM985889v3.101.primary_assembly.MN908947.3.gff3 ", outfile_snpeff, " > ", outfile_snpeff_csq))
	
	system(paste0("Rscript ./helper/Parse_SnpEff.r ", outfile_snpeff,  " ", outfile_csv))
	df_tmp <- read_csv(outfile_csv)
	df_tmp$sample <- sample_t
	return(df_tmp)
})

df_ann_raw <- bind_rows(df_ann)
unique(df_ann_raw$effect)
df_ann <- df_ann_raw %>% filter(!grepl("stream", effect))
df_ann <- df_ann %>% filter(effect != "frameshift_variant")
unique(df_ann$effect)
df_ann$X4

# df_ann <- df_ann %>% filter(effect %in% c("missense_variant", "synonymous_variant"))

data_nsp_region <- read_csv("../data/ORF_SCoV2.csv")
df_ann$gene_nsp <- sapply(df_ann$X2, function(x){
	data_nsp_region$sequence[data_nsp_region$start<=x & data_nsp_region$stop>=x]
})
df_ann$gene_pos <- sapply(df_ann$X2, function(x){
	check <- data_nsp_region$start<=x & data_nsp_region$stop>=x
	floor((x-data_nsp_region$start[check])/3)+1
})

pos_adjust <- data_nsp_region$stop[data_nsp_region$sequence=="nsp12_1"] - data_nsp_region$start[data_nsp_region$sequence=="nsp12_1"]+1
df_ann$gene_pos[df_ann$gene_nsp=="nsp12_2"] <- df_ann$gene_pos[df_ann$gene_nsp=="nsp12_2"]+pos_adjust/3
df_ann %>% select(X2, gene_nsp, gene_pos)

data_new <- lapply(seq_len(nrow(df_ann)), function(i){
	print(i)
	pos_i <- df_ann$X2[i]
	gene_i <- df_ann$gene_nsp[i]
	mut_aa <- gsub("p.", "", df_ann$mut_aa[i], fixed=T)
	aa_pos_i <- df_ann$gene_pos[i]
	aa_mut_i <- gsub("\\d", "", mut_aa)
	nt_ref_i <- df_ann$X4[i]
	nt_alt_i <- df_ann$X5[i]
	aa_ref_i <- strsplit(aa_mut_i, "")[[1]][1]
	aa_alt_i <- strsplit(aa_mut_i, "")[[1]][2]

	cdna_pos <- df_ann$pos_cdna[i]
	cdna_pos <- as.numeric(strsplit(cdna_pos, "/", fixed = T)[[1]][1])
	mut_pos_i <- cdna_pos%%3
	if(mut_pos_i==0){mut_pos_i <- 3}
	
	ref_codon_i <- subseq(ref_seq, pos_i-mut_pos_i+1, pos_i-mut_pos_i+3)
	codon_ref_i <- as.character(ref_codon_i)
	check <- GENETIC_CODE[[codon_ref_i]] == aa_ref_i

	if(grepl("del",mut_aa) | grepl("ins",mut_aa)){
		value_i <- "DEL"
		df_tmp <- tibble(Position=pos_i, Gene=gene_i, AmA_pos=aa_pos_i, nt_ref=nt_ref_i, nt_alt=nt_alt_i, sample=df_ann$sample[i], value=value_i, mutate_position=mut_pos_i, codon_ref=codon_ref_i, codon_alt=value_i, AmA_ref=aa_ref_i, AmA_alt=mut_aa, `Syn/Non`="non")
		return(df_tmp)
	} else {
		stopifnot(check)
		codon_alt_i <- strsplit(codon_ref_i, "")[[1]]
		codon_alt_i[mut_pos_i] <- nt_alt_i
		codon_alt_i <- paste(codon_alt_i, collapse = "")

		check_syn <- aa_ref_i==aa_alt_i
		syn_i <- ifelse(check_syn, "syn", "non")
		value_i <- ifelse(nt_ref_i==nt_alt_i, ".", NA)
		if(is.na(value_i)){value_i <- nt_alt_i}

		tibble(Position=pos_i, Gene=gene_i, AmA_pos=aa_pos_i, nt_ref=nt_ref_i, nt_alt=nt_alt_i, sample=df_ann$sample[i], value=value_i, mutate_position=mut_pos_i, codon_ref=codon_ref_i, codon_alt=codon_alt_i, AmA_ref=aa_ref_i, AmA_alt=aa_alt_i, `Syn/Non`=syn_i)
	}
})

data_new <- bind_rows(data_new)
data_new$Gene[data_new$Gene=='nsp12_2'] <- "nsp12"
data_new$`Syn/Non`
data_new$sample <- factor(data_new$sample, levels = c("WHP4529", "WHP4738"), labels = c("WHP4529_11774", "WHP4738_11861"))

df_ann_11774 <- data_new %>% filter(sample=="WHP4529_11774") %>% mutate(`WHP4529_11774`=value) %>% select(-sample, -value)
df_ann_11861 <- data_new %>% filter(sample=="WHP4738_11861") %>% mutate(`WHP4738_11861`=value) %>% select(-sample, -value)
df_ann_11861 %>% filter(Gene=="S")

data_out <- full_join(data %>% mutate_all(as.character), df_ann_11774 %>% mutate_all(as.character))
data_out <- full_join(data_out, df_ann_11861 %>% mutate_all(as.character))

func_fill <- function(x){
	x[is.na(x)] <- "."
	return(x)
}
data_out <- data_out %>% select(Position:WHP3662, WHP4529_11774, WHP4738_11861, everything())
data_out$`NC_045512.2` <- data_out$nt_ref
data_out <- data_out %>% mutate_at(vars(`61B1-P3-200015`:WHP4738_11861), func_fill)
data_out$`Syn/Non`

data_out_dup <- data_out %>% filter(Position %in% data_out$Position[duplicated(data_out$Position)])
data_out_no_dup <- data_out %>% filter(!Position %in% data_out$Position[duplicated(data_out$Position)])

func_paste <- function(x){
	if(all(x==x[1])){return(x[1])}
	x <- x[x!="."]
	out <- paste(x, collapse = "/")
}
data_out_dup_merge <- data_out_dup %>% group_by(Position) %>% summarise_all(func_paste)

data_add <- bind_rows(data_out_no_dup, data_out_dup_merge) %>% arrange(as.numeric(Position))
writexl::write_xlsx(data_add, "../results/data_add.xlsx")
