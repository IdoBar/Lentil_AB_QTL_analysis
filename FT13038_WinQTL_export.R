devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")

# Install and load needed packages
pacman::p_load(tidyverse, vcfR) # "radiator", 

options(stringsAsFactors=FALSE)


#### Save files for QTL analysis ####
# Save genetic map as csv for R/qtl2
analysis_name <- "Lentil_FT13038_Homozyg_V1_snps_Rep2_Ido"

gmap_map <- read_csv(glue::glue("./data/qtl2_files/{analysis_name}_gmap.csv"))
genos <- bind_rows(read_csv(glue::glue("./data/qtl2_files/{analysis_name}_founder_genos.csv")), 
                   read_csv(glue::glue("./data/qtl2_files/{analysis_name}_genos.csv")))
gmap_map %>% count(chr)
clean_vcf_file <- recent_file("./data/intermediate_files", 
                        glue::glue("Hari_parents_and_samples_Homozyg_V1_snps.Marker_miss15Geno_miss15.+clean.vcf"))


##### Read and prepare phenotypic data #####
# Read the data in from the excel files
# analysis_name <- "Lentil_FT13038"
pheno_data <- readxl::read_excel("./data/Phenotyping_AB_results.xlsx",na = c("", "NA", "NG"), col_types = c("skip", "text", rep("guess", 12))) 
loal_mapping <- readxl::read_excel("./data/Phenotyping_AB_results.xlsx", sheet = "LOAL_map", na = c("", "NA", "NG"), skip = 1) %>% select(RIL_no, Sample_ID)
clean_data <- pheno_data %>% 
  filter(!is.na(Leaf_score_ratio) & !is.na(Stem_score_ratio), !grepl("-mock", Variety), DPI>7, Rep==2) %>% 
  filter_at(vars(Variety:Isolate), all_vars(!is.na(.))) %>% 
  mutate_at(vars(Leaf_score_ratio:Score_sqrt), as.numeric) %>%
  mutate_at(vars(ends_with("_ratio")), list(~sqrt(.))) %>%
  # mutate_at(vars(ends_with("_ratio")), ~.*100) %>%
  mutate(Variety=sub("-.+", "", Variety)) %>% 
  mutate(Variety=loal_mapping$Sample_ID[match(Variety, loal_mapping$RIL_no)]) %>%
  mutate(Variety=sub("_", "-", Variety)) %>% 
  mutate(DPI=factor(DPI, levels=as.character(2:4*7)))

clean_data %>% filter(is.na(Variety))
pheno_data %>% filter(!is.na(Leaf_score_ratio) & !is.na(Stem_score_ratio), !grepl("-mock", Variety), DPI>7) %>% 
  filter_at(vars(Variety:Isolate), all_vars(!is.na(.))) %>% #mutate(Variety=sub("-.+", "", Variety)) %>%
  filter(! Variety %in% loal_mapping$RIL_no) %>%
  count(Variety) 
# clean_data %>% count(Variety) %>% print(n=Inf)
# Prepare the csv file for WinQTLCart
# add Chr label row under Chromosome
winqtl_map <- gmap_map %>% # rownames_to_column("marker") %>%
  dplyr::select(MkID=marker, Position=pos, Chromosome=chr) %>%
  mutate(Position=sprintf("%.3f",as.numeric(Position)), `Chr Label`=Chromosome, Chromosome=sub("LG","", Chromosome))  # , Chr_Label=paste0("LG", Chromosome)

all_pheno <- NULL
sampling_times <- unique(clean_data$DPI)
for (dpi in sampling_times){
  # dpi=sampling_times[2]
  clean_pheno <- clean_data %>% filter(DPI==dpi) %>%
    
    group_by(Variety) %>% summarise_at(vars(ends_with("_ratio")), mean)  %>%
    dplyr::rename(id=Variety) %>%
    # write_csv(., filedate(glue::glue("{analysis_name}_pheno_{dpi}dpi"), ".csv",
    #                       "data/qtl2_files",
    #                       dateformat = FALSE)) %>%
    rename_at(vars(ends_with("_ratio")), funs(paste( ., dpi,sep = "_")))
  # Remove samples that were filtered from the genotypic data
  
  # clean_pheno <- log_data  %>% dplyr::rename(id=Sample) # %>% filter(Sample %in% genos$rowname)
  
  if (is.null(all_pheno)) {
    all_pheno <- clean_pheno
  } else {
    all_pheno <- full_join(clean_pheno, all_pheno)
  }
}

# Check normality for each trait (after summary)
sw_test_results <- all_pheno %>% 
  gather(key = "variable_name", value = "value", contains("ratio")) %>% 
  group_by(variable_name)  %>% 
  do(broom::tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  select(-method)
# Visually inspect
car::qqPlot(all_pheno$Leaf_score_ratio_21)


# raw_pheno <- cbind(Name=paste0("*", rownames(trans_pheno)), trans_pheno) 
# Create a table that includes both the linkage map and the genotype data
raw_geno <-  vcfR::extract.gt(read.vcfR(clean_vcf_file))
colnames(raw_geno)
# colnames(raw_geno) <- sub("-", "", colnames(raw_geno),  fixed = TRUE)
genotypes <- colnames(raw_geno)

excl_cols <- colnames(raw_geno)[!colnames(raw_geno) %in% all_pheno$id]
# Choose whether to incorporate all the genotypes from the VCF or from the map

inc_cols <- colnames(raw_geno)[colnames(raw_geno) %in% all_pheno$id] # from vcf
inc_cols <- inc_cols[inc_cols %in% genos$rowname] # from linkage map construction
# dimnames(raw_geno)[[2]] %in% excl_cols

clean_geno <- raw_geno[,inc_cols] %>% apply(., c(1,2),  
                    function(sam) switch(sam, "0/0"=0, "1/1"=2, -1)) %>% data.frame(.) %>% 
  rownames_to_column(var = "MkID" ) %>% inner_join(winqtl_map, .) %>% 
  mutate_if(is.double, as.integer) %>% 
   data.frame(.) %>%
  mutate(Position=as.numeric(Position)) %>%
  arrange(Chromosome, Position)


  # mutate(MkID=sub("\\*", "",X1)) %>% inner_join(winqtl_map, .)  %>% data.frame(.) %>%
  # column_to_rownames(., "X1") %>% dplyr::select(-one_of(excl_cols))


# .[, colnames(raw_geno)[sample_cols %in% clean_pheno$id]]
# colnames(clean_geno)[grepl("X\\d+", colnames(clean_geno))] <- all_pheno$id

geno_csv <- as.data.frame(t(clean_geno)) %>% cbind(id=sub(".", "-", row.names(.), fixed=TRUE), .) 

# Write as csv for WinQTLCart
winqtl_file <- filedate(glue::glue("{analysis_name}_WinQTL"), ".csv", "WinQTL_files", 
                        dateformat = FALSE)
# Add prefix content to file
markers_per_lg <- gmap_map %>% count(chr)
# Create the prefix string
prefix_string <- sprintf("#START#\nPopulations,%d\nSamples,%d\nTraits,%d\nChromosomes,%d\nMarkers,%s\nMapFunction,Kosambi\nCrossType,Rf5\nMisTraitValue,NA\nTranslationTable,2,1,0,12,10,-1\nOther Traits,%s\n#END#\n", 
                1,ncol(clean_geno)-4, ncol(all_pheno)-1, nrow(markers_per_lg), 
                paste(markers_per_lg$n, collapse = ","), 
                paste(rep(0, nrow(markers_per_lg)), collapse = ","))
# write prefix string to file
cat(prefix_string, file=winqtl_file)
# Combine phenotypic data, genetic map and genotyping data
winqtl_pheno <- all_pheno %>% filter(id %in% inc_cols) %>% dplyr::select(-id) 
add_df <- do.call("rbind", replicate(3, 1:ncol(winqtl_pheno), simplify = FALSE)) %>% as.data.frame(.) %>% setNames(., colnames(winqtl_pheno))
winqtl_pheno <- winqtl_pheno %>% rbind(colnames(.), add_df, .)
# 1:ncol(.), 1:ncol(.), 1:ncol(.), .)
winqtl_csv <- cbind(geno_csv, winqtl_pheno)
# Add population number
winqtl_csv$id[1] <- "MkID (Pop = 1)"
winqtl_csv$id[grepl("Chr.Label",winqtl_csv$id)] <- "Chr Label"

winqtl_csv[1:10,1:10]
write_csv(winqtl_csv, path = winqtl_file, col_names = FALSE, append = TRUE)
# write_csv(all_pheno, filedate("Lentil_GBS_pheno_all_dpi", ".csv", "WinQTL_files", dateformat = FALSE))

# Number of chromo

##### export to MapChart ####
# Save all information to an mct file (for MapChart)
peaks_df <- readxl::read_excel(recent_file("./QTL_results", 
                                           glue::glue("LOD_peaks_{stacks_name}_mapping.+.xlsx")))
 
mc_file <- filedate(glue::glue("Lentil_ascochyta_{stacks_name}_QTL_map"), ".mct", "QTL_results")
unlink(mc_file)
map_df <- gmap_map %>% 
  rename(Group=chr, Position=pos)  %>% group_by(Group, Position) %>% slice(1) %>% ungroup() #%>% arrange(Group, Position)
# peaks_df <- unfactor(peaks_df)
LGs <- unique(map_df$Group)

traits <- unique(peaks_df$lodcolumn)
# names(time_points) <- c(0:2,7)
time_fill <- setNames(c(0:2,7), unique(peaks_df$dpi))
# for (i in seq_along(LGs)){
for (i in LGs){
  # i=1
  # ch <- LGs[i]
  # Write group name to file
  write(paste("chrom",  i), mc_file, append = file.exists(mc_file)) # paste("chrom", ch)
  write_tsv(map_df %>% filter(Group==i) %>% select(-Group), mc_file, append = TRUE, col_names = FALSE)
  if (peaks_df %>% filter(chr==i) %>% nrow(.) > 0)  write(sprintf("\nqtls"),  mc_file, append =TRUE)
  for (t in seq_along(traits)){
    for (d in names(time_fill)){
      
      # d=14
      # t=1
      trait <- traits[t]
      # save LOD scores to files
      lod_intervals <- peaks_df %>% filter(chr==i, lodcolumn==trait, dpi==d) %>% select(ci_lo, ci_hi) #%>%
      # as.character()
      j=1
      while (j<=nrow(lod_intervals)) {
        interval_string <- paste(rep(round(lod_intervals[j,], digits=3), each=2), collapse = " ")
        # Use filenames in MCT file
        write(sprintf("%s_%s %s C%d F%d", trait, d, interval_string, t, time_fill[d]),  mc_file, append =TRUE) # C4 F5
        j <- j+1
        # write(sprintf("%s_%s QTL_results/%s_%sdpi_%s.lod C1 L1 S0\n", t, dpi, t, dpi, ch),  
      }
      
      #       mc_file, append =TRUE) 
    }
    #   
  }
  write(sprintf("\n"),  mc_file, append =TRUE)
}



