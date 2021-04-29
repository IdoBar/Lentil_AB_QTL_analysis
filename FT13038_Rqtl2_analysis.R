devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# Install R/qtl2
if (!require("qtl2")) install.packages("qtl2", repos="https://rqtl.org/qtl2cran")
# if (!require("qtl2"))  devtools::install_github("rqtl/qtl2")
# install.packages("mixtools")
# devtools::install_github("KonradZych/phenotypes2genotypes")
devtools::install_github("IdoBar/LinkageMapView@1c6b560")
# Install and load needed packages
package_list <- c("qtl", "RColorBrewer","rentrez", "vcfR", "LinkageMapView", "openxlsx", "seqinr","tidyverse") # "radiator", "qtl2"
pacman::p_load(char = package_list)





options(stringsAsFactors = FALSE)
ncores=3

#### Analysis in R/qtl2  ####
analysis_name <- "Lentil_FT13038_Homozyg_V1_snps"
base_analysis <- sub("_Rep.+", "", analysis_name)
yaml_files <- dir("./data/qtl2_files", glue::glue("{analysis_name}_\\d+dpi+\\.yaml"), 
                  full.names = TRUE)

# pheno_data <- read_csv("./data/qtl2_files/Lentil_FT13038_pheno_14dpi.csv")
# yaml_files <- yaml_files[grepl(analysis_name, yaml_files)]
time_points <- as.character((2:4)*7)
thresh_df <- NULL
LOD_scores <- setNames(vector(mode = "list", length = length(time_points)), time_points)
peaks_df <- NULL
qtl_list <- setNames(vector(mode = "list", length = length(time_points)), time_points)
lod_thres <- 2.5

for (dpi in time_points){
  # dpi=time_points[1]
  yml <- yaml_files[grepl(dpi, yaml_files)]
  LogMsg(sprintf("Processing YAML file: %s", yml))
  gbs_data <- read_cross2(yml)
  # insert pseudomarkers into the genetic map
  
  # if (!exists("gen_map")) {
    # gen_map=gbs_data$gmap
    gen_map <- insert_pseudomarkers(gbs_data$gmap, step=1)
    gen_map <- reduce_map_gaps(gen_map)
    gen_map_df <- gen_map %>% imap_dfr(~data.frame(Locus=names(.x), 
                                                     Position=.x) %>%  
                                    mutate(Group=.y))
      
      
    # }
  
  # calculate the QTL genotype probabilities
  # gen_map=gbs_data$gmap
  pr <- calc_genoprob(gbs_data, gen_map, error_prob=0.002, cores=ncores)
  # Calculate error LOD 
  # calc_errorlod()
  # Calculate kinship matrix (needed for LMM)
  kinship_loco <- calc_kinship(pr, "loco", cores=ncores)
  # perform a genome scan by LMM analysis  
  qtl_list[[dpi]] <- scan1(pr, gbs_data$pheno, kinship_loco, cores=ncores)
  LOD_scores[[dpi]] <- as.data.frame(qtl_list[[dpi]]) %>% 
    rownames_to_column("Locus") %>% # [[dpi]]
    gather(key="Trait", value="LOD", ends_with("_ratio")) %>% 
    left_join(., gen_map_df) %>% dplyr::select(Position, Group, LOD, Trait) 
  # calculate LOD threshold based on permutation step
  perm_lod <- scan1perm(pr, gbs_data$pheno, kinship_loco, n_perm = 1000,
                        cores=ncores)
  thres_sum <- data.frame(summary(perm_lod)) %>% rownames_to_column("alpha") %>%
    mutate(dpi=dpi)
  thresh_df <- bind_rows(thresh_df, thres_sum)

  # Find peaks
  peaks <- find_peaks(qtl_list[[dpi]], gen_map, threshold=lod_thres, drop=1, peakdrop=2) %>%
    mutate(dpi=dpi)
  peaks_df <- bind_rows(peaks_df, peaks) # thres_sum$Leaf_lesion_percent
  # get physical position at location - need to add physical map (pmap) to the cross - check R/qtl
  # peak_Mbp <- max(qtl_list[[dpi]], gbs_data$pmap)$pos
  
}
mean(thresh_df$Stem_score_ratio)


write_xlsx(peaks_df, 
                 filedate(glue::glue("LOD_peaks_{analysis_name}_mapping"), 
                          ".xlsx", "QTL_results"), sheet ="LOD_peaks")
write_xlsx(thresh_df, 
                 filedate(glue::glue("LOD_peaks_{analysis_name}_mapping"), 
                          ".xlsx", "QTL_results"), sheet = "Sign_peaks")
# xlsx::write.xlsx(as.data.frame(peaks_df), filedate("LOD_peaks_Chr_mapping", ".xlsx", "data"), row.names = FALSE)
# thresh_df
peaks_df %>% group_by(chr) %>% summarise(max_lod=max(lod))

#### QTL SNPs ####
# Retrieve SNPs under the QTL (any snps between tag_144852_343285964 to tag_148476_364570997)
gmap <- readr::read_csv(glue::glue("./data/qtl2_files/{analysis_name}_gmap.csv"))
peaks_df <- readxl::read_excel(recent_file("./QTL_results",
                               glue::glue("LOD_peaks_{analysis_name}_mapping.+.xlsx")))
# peaks_df <- readxl::read_excel("./QTL_results/LOD_peaks_Lentil_FT13038_Homozyg_V1_snps_Rep2_mapping_29_11_2018.xlsx")
# Read in vcf, filter to include just relevant markers
shared_vcf_file <- recent_file("./data/", "Hari_parents_and_samples.+.vcf")
vcf <- readr::read_tsv(shared_vcf_file, comment = "##") %>% set_names(., sub("-1\\.sort.+", "", colnames(.)))

# marker_map <- as.data.frame(rowRanges(snp_vcf)) %>% rownames_to_column("var_id") # %>% as.tibble()
# snp_ranges <- rowRanges(snp_vcf)
# Find QTL regions for the trait on interest and manually select range
signif_lod_qtls <- peaks_df %>% filter(lod>3)
qtl_transcripts_vec <- NULL
for (chr in signif_lod_qtls$chr) {
# for (qtl_trait in unique(signif_lod_qtls$lodcolumn)){
#   for (qtl_dpi in unique(signif_lod_qtls$dpi)){
    # qtl_trait <- "Leaf_score_ratio"
    # qtl_dpi <- 14
    # qtl_region_markers <- function(gmap, region_chr, region_start, region_end){
    #   gmap %>%
    #     dplyr::filter(chr==region_chr, pos>=region_start, pos<=region_end)
    # } 
    qtl_regions <- signif_lod_qtls %>% filter(chr==chr) %>% # lodcolumn==qtl_trait, dpi==qtl_dpi
      .[,c("chr", "ci_lo", "ci_hi")] #%>% pmap_dfr(~qtl_region_markers(gmap, ..1, ..2,..3))
    ci_min <- min(qtl_regions$ci_lo)
    ci_max <- max(qtl_regions$ci_hi)
    # for (chrom in qtl_regions$)
    # Find markers under the qtl
    qtl_markers <- gmap %>%
      dplyr::filter(chr==chr, pos>=ci_min, 
                    pos<=ci_max) %>% 
      mutate(locus=sub("_\\d+$","", marker)) #%>%
      # inner_join(vcf[c("POS", "ID", "REF", "ALT", "ILL6002", "ILW-180")], by = c("marker"="ID")) %>%
      # mutate_at(vars(starts_with("IL")), ~sub(":.+", "", .)) %>%
      # dplyr::rename(snp_pos=POS)
    qtl_transcripts_vec <- c(qtl_transcripts_vec, unique(qtl_markers$locus)) %>% unique(.)
    
    
    
  }
    
# }

# Read in transcriptome
transcriptome <- seqinr::read.fasta("./data/CASSAB_K91_Reference_FINAL.fasta")
# Extract sequences
qtl_transcripts <- transcriptome[map_lgl(transcriptome, ~sub("^>", "", attr(., "Annot")) %in% qtl_transcripts_vec)]
# Save as fasta
seqinr::write.fasta(qtl_transcripts, names(qtl_transcripts), open = "w",
                    filedate(glue::glue("./QTL_results/{analysis_name}_qtl_transcripts"), ext = ".fasta"))
# Load BLAST results
# qtl_trans_annot <- readr::read_tsv(glue::glue("./QTL_results/{analysis_name}_qtl_transcripts.blastx.txt"), col_names = c('query_acc.ver','subject_acc.ver','perc_identity','alignment_length','mismatches','gap_opens','q.start','q.end','s.start','s.end','evalue','bit_score','perc_positives','query_sbjct_frames'), comment = "#")

# qtl_trans_annot <- read_csv(recent_file("./QTL_results/", glue::glue("{analysis_name}_qtl_transcripts.+.blastx.csv")), col_names = c('query_acc.ver','subject_acc.ver','perc_identity','alignment_length','mismatches','gap_opens','q.start','q.end','s.start','s.end','evalue','bit_score','perc_positives','query_sbjct_frames', "extra_field"), comment = "#")

qtl_trans_annot <- read_tsv(recent_file("./QTL_results/", glue::glue("{analysis_name}_qtl_transcripts.+.blastx.outfmt")), col_names = c('query_acc.ver','subject_acc.ver','perc_identity','alignment_length','mismatches','gap_opens','q.start','q.end','s.start','s.end','evalue','bit_score','stitle','ssciname', "scomname"), comment = "#")


top_annot <- qtl_trans_annot %>% group_by(query_acc.ver) %>% 
  arrange(dplyr::desc(bit_score)) %>% dplyr::slice(1:10)
subject_acc <-  top_annot$subject_acc.ver

# entrez_db_searchable('protein')
accession_to_title <- function(acc, db){
  ncbi_id <- rentrez::entrez_search(db, paste0(acc, "[ACCN]"))
  # LogMsg("Searching ")
  if (ncbi_id$count>0) {
    summ <- rentrez::entrez_summary(db=db, id=ncbi_id$ids) 
    # cds <- rentrez::entrez_fetch(db, id=ncbi_id$ids, rettype = "fasta_cds_aa", retmode = "text")
    return(tibble("Accession"=acc, "Title"=summ$title, "Organism"=summ$organism)) # , 
    # "Prot_seq"= cds)
  }
  # return(tibble("Accession"=acc, "Title"=NA, "Organism"=NA))
}

# ncbi_id <- rentrez::entrez_search('protein', paste0("GAU44738.1", "[ACCN]"))
iter_num <- 20
iterations <- seq(1,nrow(top_annot), length.out = iter_num)
complete_annotation <- vector(mode="list", length = iter_num-1)
for (i in 1:(iter_num-1)){
  LogMsg(glue:glue("Processing transcripts {iterations[i]}:{iterations[i+1]}"))
  complete_annotation[i] <- purrr::map_dfr(top_annot$subject_acc.ver[iterations[i]:iterations[i+1]], 
                                        ~accession_to_title(., 'protein'))
   
}
 
# annot_table <-  complete_annotation %>% right_join(top_annot, 
#                                       by = c('Accession'='subject_acc.ver')) %>%
  
annot_table <-  top_annot %>% mutate(Title=stitle) %>% 
  mutate(strand=ifelse(s.end-s.start<0, "-", "+"),
         # Prot_seq=sub("^.+\\]", "", gsub("\n", "", Prot_seq, fixed=TRUE)),
         trans_start=if_else(strand=="+", q.start, q.end), trans_end=if_else(strand=="+", q.end, q.start))
best_matches <- annot_table %>% 
  filter(!is.na(Title), !grepl("unknown|hypothetical|uncharacterized", Title, ignore.case = TRUE)) %>% 
  group_by(query_acc.ver) %>% top_n(1, dplyr::desc(bit_score)) %>% ungroup() %>% 
  write_xlsx(., excel_file = filedate(glue::glue("{analysis_name}_qtl_transcripts_annot"), ext=".xlsx", 
                                                                            outdir = here::here("QTL_results")),
                                                   sheet = glue::glue("transcript_annot"), overwritesheet=TRUE)

# Check the effect of the SNP
check_snp_effect <- function(dna_obj, snp_pos, ref_allele, alt_allele, strand, trans_start, trans_end){
  # test_mark <- marker_table[7,]
  # snp_loc <- case_when(snp_pos<trans_start ~ "5'UTR", snp_pos>trans_end ~ "3'UTR", 
  #                      TRUE ~ "cds")
  if (snp_pos<trans_start | snp_pos>trans_end ) return(NA_character_)
  # if (!grepl("cds", snp_location)) return(NA_character_)
  # qtl_transcripts[[locus]]
  mark_seq_ref <- dna_obj -> mark_seq_alt
  mark_seq_ref[snp_pos] <- ref_allele
  mark_seq_alt[snp_pos] <- alt_allele
  sense <- if_else(strand=="+", "F", "R")
  ref_orf <- getTrans(getFrag(mark_seq_ref, trans_start, trans_end), sense)
  alt_orf <- getTrans(getFrag(mark_seq_alt, trans_start, trans_end), sense)
  if (any(ref_orf != alt_orf)) return(sprintf("Non-Synonymous (%s%d%s)", ref_orf[ref_orf!=alt_orf], which(ref_orf!=alt_orf), alt_orf[alt_orf!=ref_orf]))
  case_when(all(ref_orf == alt_orf) ~ "Synonymous", 
            grepl("\\*", c2s(ref_orf)) | grepl("\\*", c2s(alt_orf)) ~ "Nonsense", 
            TRUE ~ NA_character_)
}
all_qtl_markers <- NULL
for (qtl_trait in unique(signif_lod_qtls$lodcolumn)){
  # qtl_trait="Leaf_score_ratio"
  times <- signif_lod_qtls %>% dplyr::filter(lodcolumn==qtl_trait) %>% .$dpi %>% unique()
  for (qtl_dpi in times){
    # qtl_dpi=times[1]
    chroms <- signif_lod_qtls %>% dplyr::filter(lodcolumn==qtl_trait, dpi==qtl_dpi) %>% .$chr %>% unique()
  for (chrom in chroms){
    # chrom=chroms[1]
    qtl_regions <- signif_lod_qtls %>% filter(lodcolumn==qtl_trait, dpi==qtl_dpi, chr==chrom) %>% # 
      .[,c("chr", "ci_lo", "ci_hi")] #%>% pmap_dfr(~qtl_region_markers(gmap, ..1, ..2,..3))
    ci_min <- min(qtl_regions$ci_lo)
    ci_max <- max(qtl_regions$ci_hi)
    # for (chrom in qtl_regions$)
    # Find markers under the qtl
    qtl_markers <- gmap %>%
      dplyr::filter(chr==chrom, pos>=ci_min, 
                    pos<=ci_max) %>% 
      mutate(locus=sub("_\\d+$","", marker)) %>%
      inner_join(vcf[c("POS", "ID", "REF", "ALT", "ILL6002", "ILW-180")], by = c("marker"="ID")) %>%
      mutate_at(vars(starts_with("IL")), ~sub(":.+", "", .)) %>%
      dplyr::rename(snp_pos=POS)
    all_qtl_markers <- bind_rows(qtl_markers, all_qtl_markers)
    # annot_table %>% inner_join(qtl_markers)
    # write_xlsx(qtl_markers, excel_file = filedate(glue::glue("{analysis_name}_qtl_transcripts_annot"), ext=".xlsx", 
    #                                               outdir = here::here("QTL_results")), 
    #            sheet = glue::glue("{qtl_trait}{chrom}{qtl_dpi}dpi"),overwritesheet=TRUE)
    
    # annot_table %>% filter(!is.na(Title), !grepl("unknown|hypothetical|uncharacterized", Title, ignore.case = TRUE)) %>% 
    #   write_xlsx(., excel_file = filedate(glue::glue("{analysis_name}_qtl_transcripts_annot"), ext=".xlsx", 
    #                                       outdir = here::here("QTL_results")), 
    #              sheet = glue::glue("{qtl_trait}_{chrom}_{qtl_dpi}dpi_known"), overwritesheet=TRUE)
    
    
    marker_table <- qtl_markers %>% 
      inner_join(best_matches[c("query_acc.ver","Title", "strand", "trans_start", "trans_end")], 
                 by=c("locus"="query_acc.ver")) %>%
      mutate(snp_location=case_when(snp_pos<trans_start ~ "5'UTR", snp_pos>trans_end ~ "3'UTR", 
                                    (snp_pos - trans_start) %% 3 == 0 ~ "cds_3", 
                                    TRUE ~ sprintf("cds_%d", (snp_pos - trans_start) %% 3)),
             # snp_effect=if_else(grepl("cds", snp_location),
             snp_effect=pmap_chr(., ~with(list(...), check_snp_effect(qtl_transcripts[[locus]], snp_pos, REF, ALT, strand, trans_start, trans_end)))) %>% distinct()#, NA_character_)) 
    
    write_xlsx(marker_table, excel_file = filedate(glue::glue("{analysis_name}_qtl_transcripts_annot"), ext=".xlsx", 
                                                   outdir = here::here("QTL_results")), 
               sheet = glue::glue("{qtl_trait}_{chrom}_{qtl_dpi}dpi_eff"), overwritesheet=TRUE)
  }
    
  }
}
# gsub(", ", "','", "query_acc.ver, subject_acc.ver, perc_identity, alignment_length, mismatches, gap_opens, q.start, q.end, s.start, s.end, evalue, bit_score, perc_positives, query_sbjct_frames")
# test <- accession_to_title("XP_004501916.1", 'protein')
# subject_gb <- entrez_search('protein', "XP_004501916.1[ACCN]")
# subject_seq <- entrez_fetch(db='protein', id=subject_gb$ids, rettype = "fasta_cds_aa", retmode = "text") 
# summ <- entrez_summary(db='protein', id=subject_gb$ids) 

# check_snp_effect_vect <- Vectorize(check_snp_effect)



# Compare the LOD scores for each trait with each method
# color <- c("slateblue", t_col("violetred", percent = 0))
colors <- adjustcolor(RColorBrewer::brewer.pal(4, "Set1"), alpha.f = 0.65)

#par(mar=c(4.1, 4.1, 1.6, 1.1))
par(mfrow=c(2,1))
# lod_vars <- paste0("out_pg_", time_points)

ymx <- max(peaks_df$lod)
threshold=3
pdf(filedate(glue::glue("{analysis_name}_LOD"), ".pdf", "plots"), width = 8, height = 6)
# unique(signif_lod_qtls$lodcolumn)
for(i in seq_along(colnames(gbs_data$pheno))) {
  # i=4
  plot(qtl_list[[1]], gen_map, lodcolumn=i, col=colors[1], main=colnames(gbs_data$pheno)[i], 
       xlab="Linkage Group",
       ylim=c(0, ymx*1.5))
  if (colnames(gbs_data$pheno)[i] %in% peaks_df$lodcolumn[peaks_df$lod==ymx]){
    abline(h=threshold, col="darkblue", lty="dashed", lwd=1) #ymx-1   h=thresh_df[1,i+1], col=colors[1]
    text(x=0, y=threshold, labels = "LOD threshold", adj= c(-0.05, -0.5))
  }
  for (o in 2:length(time_points)){
    # o=3
    
    plot(qtl_list[[time_points[o]]], gen_map, lodcolumn=i, col=colors[o],
         add = TRUE)
    # abline(h=2.5, col=colors[o], lty="dashed", lwd=0.8)
  }
  
  legend("topright", lwd=2, col=colors, paste(time_points, "dpi"), bg="gray90") # "H-K", , lty=c(1)
}
dev.off()


# Plot each trait and time separately
# for(i in seq_along(colnames(gbs_data$pheno))) {
#   for (o in 1:length(time_points)){
#   # i=4
#     plot(qtl_list[[time_points[o]]], gen_map, lodcolumn=i, col=colors[o], main=colnames(gbs_data$pheno)[i],
#                             xlab='Linkage Group',ylim=c(0, ymx*1.02))
#     abline(h=thresh_df[o,i+1], col=colors[o], lty="dashed", lwd=0.8)
#    legend("topright", lwd=2, col=colors, paste(time_points, "dpi"), bg="gray90")
#     # o=3
#     
#   }
#   
#    # "H-K", , lty=c(1)
# }
par(mfrow=c(1,1))
# Focus on one chromosome and only significant trait
lod_peaks <- peaks_df %>% filter(lod>threshold)
traits <- unique(lod_peaks$lodcolumn)
x_intercept_start <- min(lod_peaks$ci_lo)

x_intercept_end <- max(lod_peaks$ci_hi)
chrom <- "LG3"
x_intercept_locus <-  names(gen_map[[chrom]][gen_map[[chrom]]==intercept_start])
y_intercept <- qtl_list[[1]] %>% as.data.frame() %>% rownames_to_column(var="locus") %>% 
  filter(locus==x_intercept_locus) %>% .[,traits[1]]
# str(qtl_list)


pdf(filedate(glue::glue("{analysis_name}_LOD_{chrom}"), ".pdf", "plots"), width = 8, height = 6)
for(i in seq_along(unique(lod_peaks$lodcolumn))) {
  # i=4
  plot(qtl_list[[1]], gen_map, lodcolumn=i, col=colors[1], main=colnames(gbs_data$pheno)[i], chr=chrom,
       xlab=sub("LG", "Linkage Group ", chrom),
       ylim=c(0, ymx*1.5))
  if (colnames(gbs_data$pheno)[i] %in% peaks_df$lodcolumn[peaks_df$lod==ymx]){
    
    segments(x0 = intercept_start, y0=y_intercept, x1 = intercept_end, y1 = y_intercept, col="darkblue", lty="dashed", lwd=1) # h=thresh_df[1,i+1], col=colors[1]
    segments(x0 = intercept_start, y0=y_intercept, x1 = intercept_start, y1 = 0, col="darkblue", lty="dashed", lwd=1) # vertical line
    segments(x0 = intercept_end, y0=y_intercept, x1 = intercept_end, y1 = 0, col="darkblue", lty="dashed", lwd=1)
    text(x=intercept_start, y=0.5, labels = "QTL\nregion", adj= -0.25)
    abline(h=threshold, col="darkblue", lty="dashed", lwd=1) #ymx-1   h=thresh_df[1,i+1], col=colors[1]
    text(x=0, y=threshold, labels = "LOD threshold", adj= c(-0.05, -0.5))
  }
  
  for (o in 2:length(time_points)){
    # o=3
    
    plot(qtl_list[[time_points[o]]], gen_map, lodcolumn=i, col=colors[o],  chr=chrom,
         add = TRUE)
    # abline(h=2.5, col=colors[o], lty="dashed", lwd=0.8)
  }
  
  legend("topright", lwd=2, col=colors, paste(time_points, "dpi"), bg="gray90") # "H-K", , lty=c(1)
}
dev.off()




##### export to MapChart ####

lod_peaks <- peaks_df %>% filter(lod>3) %>% unfactor(.)
# read map from file
gmap_map <- read_csv(glue::glue("./data/qtl2_files/{base_analysis}_gmap.csv"))
gmap_df <- gmap_map %>% dplyr::select(Locus=marker, Position=pos, Group=chr) %>%
  mutate(MapChart_inc=if_else(Locus %in% all_qtl_markers$marker, "", "O"))
# Get map from R\qtl2 (including pseudomarkers)
# raw_map_df <-  names(gbs_data$gmap) %>%
#   purrr::map(function(chr) data.frame(Locus=names(gbs_data$gmap[[chr]]),
#                                Position=gbs_data$gmap[[chr]]) %>%
#         mutate(Group=chr)) %>% map_df(bind_rows) %>% 
#   mutate(MapChart_inc=if_else(Locus %in% qtl_markers$marker, "", "O"))
gmap_df %>% filter(MapChart_inc=="")
LGs <- unique(gmap_df$Group)

traits <- unique(lod_peaks$lodcolumn)
# names(time_points) <- c(0:2,7)
# MapChart colours
mapchart_cols <- set_names(1:9, 
              c("black", "red", "green", "blue", "yellow", "magenta", "lightgreen", "brown", "lightblue"))
mapchart_fills <- set_names(0:7, 
              c("blank", "full", "vert_bars", "hor_bars", "diag_bars", "diag_bars2", "grid", "diag_grid"))

fill_patterns <- c(1,0, seq(9,2,-1) )
trait_cols <-  setNames(names(mapchart_cols[seq_along(traits)+1]), traits)#mapchart_cols[seq_along(traits)+1]
time_fill <- setNames(names(mapchart_fills[seq_along(unique(lod_peaks$dpi))+1]), unique(lod_peaks$dpi))
# Save all information to an mct file (for MapChart)
mc_file <- filedate(glue::glue("{analysis_name}_QTL_map"), ".mct", "QTL_results")
unlink(mc_file)
for (i in seq_along(LGs)){
  # i=5
  ch <- LGs[i]
  # Write group name to file
  write(paste("chrom", ch), mc_file, append = file.exists(mc_file))
  write_delim(gmap_df %>% filter(Group==ch) %>% 
              dplyr::select(-Group), mc_file, append = TRUE, col_names = FALSE, delim = " ")
  if (lod_peaks %>% filter(chr==ch) %>% nrow(.) > 0)  write(sprintf("\nqtls"),  mc_file, append =TRUE)
  for (trait in traits){
    for (d in unique(lod_peaks$dpi)){
      
      # d=1
      # t=1
      # trait <- traits[t]
      # save LOD scores to files
      lod_intervals <- lod_peaks %>% filter(chr==ch, lodcolumn==trait, dpi==d) %>% 
        dplyr::select(ci_lo, ci_hi) #%>%
      # as.character()
      if (nrow(lod_intervals)>0) {
        interval_string <- paste(rep(round(lod_intervals, digits=3), each=2), collapse = " ")
        # Use filenames in MCT file
        write(sprintf("%s_%s %s C%d F%d", trait, d, interval_string, mapchart_cols[trait_cols[trait]], mapchart_fills[time_fill[d]]),  mc_file, append =TRUE) # C4 F5
        # write(sprintf("%s_%s QTL_results/%s_%sdpi_%s.lod C1 L1 S0\n", t, dpi, t, dpi, ch),  
      }
      
      #       mc_file, append =TRUE) 
    }
    #   
  }
  write(sprintf("\n"),  mc_file, append =TRUE)
}


#### Plotting map with LinkageMapView ####
## draw tickmarks at each cM from 0 to largest position of linkage groups to be drawn

# read map from file
# gmap_map <- read_csv(glue::glue("./data/qtl2_files/{base_analysis}_gmap.csv"))
# gmap_df <- gmap_map %>% dplyr::select(Locus=marker, Position=pos, Group=chr) %>%
#   mutate(MapChart_inc=if_else(Locus %in% qtl_markers$marker, "", "O"))

maxpos <- floor(max(gmap_df$Position))
at.axis <- seq(0, maxpos)

## put labels on ruler at every 10 cM
axlab <- vector()
for (lab in 0:maxpos) {
  if (!lab %% 10) {
    axlab <- c(axlab, lab)
  }
  else {
    axlab <- c(axlab, NA)
  }
}
qtl_df <- peaks_df %>% filter(lod>3) %>% 
  select(chr, qtl=lodcolumn, so=ci_lo, eo=ci_hi) %>% 
  mutate(col=trait_cols[qtl], si=so, ei=eo) %>% 
  select(chr, qtl, so, si, ei, eo, col)
gmap_plot_df <- gmap_df %>% select(group=Group, position=Position, locus=Locus) %>% as.data.frame()
# qtl_list[[2]]
unique(gmap_plot_df$group)
# Plot map with QTLs only
outfile = filedate(glue::glue("{analysis_name}_QTL_map_view"), ".pdf", "plots")
lmv.linkage.plot(gmap_plot_df,outfile,dupnbr = TRUE, ruler=TRUE, showonly=all_qtl_markers %>% filter(chr=="LG2") %>% .$marker,
                qtldf = qtl_df, lg.col = "lightblue1",denmap = FALSE, 
                posonleft=rep(FALSE,length(unique(gmap_plot_df$group))),
                pdf.width=15, pdf.height = 8) #,
                #, pdf.width=10, pdf.height = 8, cex.axis = 1, at.axis = at.axis, labels.axis = axlab)
# Plot complete map
outfile = filedate(glue::glue("{analysis_name}_map_view"), ".pdf", "plots")
lmv.linkage.plot(gmap_plot_df,outfile,dupnbr = TRUE, ruler=TRUE)
# Plot density map
outfile = filedate(glue::glue("{analysis_name}_density_map_view"), ".pdf", "plots")
lmv.linkage.plot(gmap_plot_df,outfile,dupnbr = TRUE, denmap = TRUE, cex.dens=1,
                 pdf.width=10, pdf.height = 8)
qtl_markers$marker %in% gmap_plot_df$locus # mapthese = qtl_markers$marker, 






