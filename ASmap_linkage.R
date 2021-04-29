# install.packages(c("devtools", "pacman"))
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R") # needs Rcurl installed
# Install and load needed packages
# install.packages("ASMap")
package_list <- c("tidyverse", "qtl", "RColorBrewer", "Rsamtools", "doFuture", "foreach", "ASMap", "vcfR", "pheno2geno") # "radiator", 
pacman::p_load(char=package_list)

# Onemap requires R>=3.4 and Rhtslib and zlibbioc
# git_packs <- c("augusto-garcia/onemap")
# install.deps(git_packs, repo = "git")
# devtools::install_github(git_packs)
source("./src/vcf_filtration.R")
options(stringsAsFactors=FALSE)

# Load vcf file and load associated .RData
analysis_name <- "Lentil_FT13038_Homozyg_V1_snps"
vcf_file <- "./data/hari_filtered_V1.recode.vcf"  # ref_stacks2_populations.snps.vcf
parents_vcf_file <- "./data/filtered_Hari_parents_V1.snps.vcf"
group_markers <- "LG"
vcf_basename <- tools::file_path_sans_ext(basename(vcf_file))
  

# Filter parents file
clean_parents_vcf <- vcf_filtration(parents_vcf_file, miss_rates = seq(0.2, 0.15, -0.05), 
                                    geno_miss_rate = 0, remove_hetero = TRUE, remove_multi_allele = TRUE,
                                    check_poly_parents = "IL", poly_only = TRUE)

parents_vcf <- clean_parents_vcf %>% mutate(ID=paste(`#CHROM`, POS, sep="_"))
sample_vcf <- read_tsv(vcf_file, comment = "##") %>% mutate(ID=paste(`#CHROM`, POS, sep="_"))

clean_sample_vcf <- vcf_filtration(sample_vcf, miss_rates = seq(0.5, 0.3, -0.05), geno_miss_rate = 0.5,
                                  remove_hetero = TRUE, remove_multi_allele = TRUE, poly_only = TRUE) 

is.data.frame(sample_vcf)
# clean_sample_vcf <- clean_sample_vcf 
RIL_samples <- colnames(clean_sample_vcf)[10:ncol(clean_sample_vcf)]
parents <- colnames(parents_vcf)[10:ncol(parents_vcf)]


shared_vcf <- parents_vcf %>% inner_join(clean_sample_vcf[c("ID", RIL_samples)])
# shared_vcf %>% dplyr::mutate_at(vars(one_of(RIL_samples)),
#                                           .funs = ~ifelse(grepl("^1/0", .) | grepl("^0/1", .), "./.", .))

# fix phasing
shared_vcf[1:10,10:20]

shared_vcf_file <- filedate("Hari_parents_and_samples_Homozyg_V1_snps_filtered", ext = ".vcf", 
                            outdir = "./data") 
                            
# copy shared file to folder
# Manually grep to the new file in git bash
system2("grep", args=c("'^##'", vcf_file), stdout = shared_vcf_file, stderr = shared_vcf_file)
if (file.exists(shared_vcf_file)) readr::write_tsv(shared_vcf, shared_vcf_file, append = TRUE, col_names = TRUE)
# 
# clean_vcf <- vcf_filtration(shared_vcf_file, miss_rates = seq(0.2, 0.15, -0.05), geno_miss_rate = 0.15, remove_hetero = TRUE, remove_multi_allele = TRUE, check_parents = "IL")



#### ASMap ####
# cleaned file to read
clean_vcf <- shared_vcf
# clean_vcf_file <- "./data/intermediate_files/ref_stacks2_populations.snps.miss10.clean.vcf"
# clean_vcf_basename <- sub(".clean.vcf", "", basename(clean_vcf_file))
# stacks_name <- sub("_populations.snps", "", sub("\\.miss.+", "", clean_vcf_basename))
# Manually convert to MapMaker format
# clean_vcf <- read_tsv(clean_vcf_file, comment = "##")
# Subset just the genotypes (also remove unwanted samples, RL and genotype calling meta-information)

vcf_geno <-  clean_vcf %>% #set_names(sub("-", "", names(.))) %>% 
  mutate_at(c(parents, RIL_samples), .funs = ~sub(":.+", "", .)) 
# row.names(vcf_geno) <- clean_vcf$ID
# fix phasing
vcf_geno_phased <- vcf_geno %>% filter(ILL6002=="0/0") %>% 
  bind_rows(vcf_geno %>% filter(ILL6002=="1/1") %>% 
              mutate_at(c(parents, RIL_samples), ~chartr("01","10",.x)))
  
  
 # Create the dictionary to replace the values
dict_keys <- c("0/0", "0/1", "1/0", "1/1", "./.", ".")
# For ASMap
mst_dict <- c("A", "X", "X", "B", "U", "U")
# For MapMaker
# MapMaker_dict <- c("A", "H", "H", "B", "-")

# geno_converted <- dict_values[geno]

# Apply the conversion to the entire table
# create a custom function to perform the conversion on a vector of input strings
convert_genotypes <- function(genos, diction, keys){
  names(diction) <- keys
  # Convert the 1st object in each list (the genotype) based on our dictionary
  geno_converted <- sapply(genos, function(g) {
    diction[sub(":.+", "", g)]
  })
  return(geno_converted) # return a vector of converted genotypes
}

# Import to ASMap
raw_geno <- as.data.frame(apply(vcf_geno_phased[c(parents, RIL_samples)], 2, 
                  function(col) convert_genotypes(col, mst_dict,                                                  dict_keys)), stringsAsFactors=FALSE)
raw_geno[is.na(raw_geno)] <- "U"
# row.names(raw_geno) <- paste(clean_vcf$`#CHROM`, clean_vcf$POS, sep="_") 
row.names(raw_geno) <- vcf_geno_phased$ID
raw_geno[1:10,1:10]
# illegal_characters <- apply(raw_geno, 1, function(row) sum(!grepl("[ABUX]", row)))
# illegal_characters[illegal_characters>0]
# Measure heterozygosity
# mean(apply(raw_geno, 1, function(x) sum(x=="X")/ncol(raw_geno)))
# raw_geno[row.names(raw_geno)=="Locus_617_Contig1_2965",]

range_length <- 8
# p_value_range <- 5e-10*(5^(1:range_length))
p_value_range <- 1e-15*(2^(0:range_length-1))
prettyNum(p_value_range)
# p_value_range <- 1e-8*(2^(0:(range_length-1)))
# For Hari's non-hetero set, 5.12e-11

# Create cluster with desired number of cores
registerDoFuture()
plan(multiprocess, workers=4)

# Create maps for each of the p-values
reduced_map_list <- foreach(p = p_value_range, 
          .final = function(x) setNames(x, paste0("asmap_p", seq_along(p_value_range)))) %dopar% {
  # pacman::p_load_gh(c("IdoBar/onemap@patch-3"))
  # l=3, rf=0.2
  # dummy <- length(onemap_file)
  map <- mstmap(raw_geno, pop.type = "RIL5", miss.thresh=0.25, p.value = p, bychr = FALSE)
  # Remove linkage group with small number of markers (unlinked)
  # pheno2geno::removeTooSmallChromosomes(map, minNrOfMarkers=10, 
  #                                                 verbose=FALSE)          
  # p.value = 5e-2 , objective.fun = "ML"
}


# reduced_map_list <- lapply(asmap_thres_list, 
#                     function(l) pheno2geno::removeTooSmallChromosomes(l, 
#                                                 minNrOfMarkers=10, verbose=FALSE))
# Check number of markers per groups
sapply(reduced_map_list, function(l) qtl::nmar(l))  # %>% purrr::map(~sum(.))


# process_asmap <- function(asmap){
#   # asmap <- asmap_thres_list[1]
#   tibble(map_name=names(asmap), linkage_group=names(qtl::nmar(asmap[[1]])), 
#              marker_num=qtl::nmar(asmap[[1]]))
# }

asmap_thres_df <- imap_dfr(reduced_map_list,
                    ~tibble(map_name = .y, linkage_group = names(qtl::nmar(.x)),
                                                     marker_num = qtl::nmar(.x)))

asmap_thres_df %>% filter(marker_num>150) %>% count(map_name) %>% filter(n>5)

# selected p7 (from range)
selected_p <- 9 #3
asmap <- reduced_map_list[[selected_p]]
asmap <- pheno2geno::removeTooSmallChromosomes(asmap, minNrOfMarkers=5, 
                                               verbose=TRUE)
LogMsg(glue::glue("Selected p.value is: {p_value_range[selected_p]}"))
LogMsg(glue::glue("Number of markers in each linkage group: {paste(qtl::nmar(asmap), collapse=', ')} (Total={sum(qtl::nmar(asmap))})"))
LogMsg(glue::glue("Length of each linkage group: {paste(qtl::chrlen(asmap), collapse=', ')}")) 


# use function quickEst() to re-estimate the genetic distances 
qtl::plotMissing(asmap)
sg <- statGen(asmap, bychr = FALSE, stat.type = "miss")
geno_miss_thresh <- 0.35
mapBC1 <- subset(asmap, ind = sg$miss < sum(qtl::nmar(asmap))*geno_miss_thresh)
gc <- genClones(mapBC1, tol = 0.975)
if ("cgd" %in% names(gc)) {
  mapBC2 <- fixClones(mapBC1, gc$cgd, consensus = FALSE)
} else {mapBC2 <- mapBC1}
mapBC3 <- jittermap(mapBC2)
# Check 
profileMark(mapBC3, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
# Remove highly distorted markers
dist_prop <- 0.35
mm <- statMark(mapBC3, stat.type = "marker")$marker %>% 
  rownames_to_column("marker") %>% 
  filter(AA<dist_prop  | BB<dist_prop)
LogMsg(glue::glue("Number of distorted markers:: {nrow(mm)}"))
# dm <- markernames(mapBC3)[(mm > 0.7) | (mm < 0.3)]
mapBC4 <- drop.markers(mapBC3, mm$marker)

profileMark(mapBC4, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
# Pull markers that are co-located
# mapDHs <- pullCross(mapBC4, type = "co.located")
mapDHs <- mapBC4
# Pull markers that have significant segregation distortion with a p-value less than 0.02 and have a missing
# value proportion greater than 0.03.
# # mapDHs <- pullCross(mapDHs, type = "seg.distortion", pars =
#                       list(seg.ratio = "45:20:35"))
mapDHs <- pullCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.25))
mapDHs <- pullCross(mapDHs, type = "co.located")
# names(mapDHs)
ncol(mapDHs$seg.distortion$data)
ncol(mapDHs$missing$data)
ncol(mapDHs$co.located$data)
# profileMark(mapDHs, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1, 4), type = "l")
nmar(mapDHs)
# Visualise data
heatMap(mapDHs, lmax = 8)
# Profile the statistics for the genotypes across
pg <- profileGen(mapDHs, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id =
             "Genotype", xo.lambda = 15, layout = c(1, 3), lty = 2, cex=0.7)
miss_ind <- pg$stat$miss[pg$stat$miss>sum(nmar(mapDHs))*geno_miss_thresh]
# asmap <- pheno2geno::removeTooSmallChromosomes(asmap, minNrOfMarkers=10, 
#                                                verbose=TRUE)


if (length(miss_ind>0)) mapDHs <-  subsetCross(mapDHs, ind = pg$stat$miss[pg$stat$miss<400])
  # mapRIL5_merged <- mergeCross(mapBC5, merge = list("L3" = c("L3", "L4"), 
#                                                    "L1" = c("L1", "L2"), 
#                                                   "L8" = c("L8", "L9"),
#                                                   "L5" = c("L5", "L10"),
#                                                   "L6" = c("L6", "L11"),
#                                                   "L7" = c("L7", "L12"),
#                                                   "L15" = c("L15", "L16")
                                                           # ))
mapRIL5 <- mstmap(mapDHs, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, p.value = 2) # p_value_range[selected_p]


# calculate genetic distance
# map2 <- quickEst(mapDH, map.function = "kosambi")

# Visualise data
# lod_range <- seq(10,100, by = 5)
# 
# foreach(l=lod_range) %dopar% {
#   pdf(filedate(sprintf("%s_asmap_lod%d_test", stacks_name, l), ".pdf", 
#                glue::glue("./plots/{stacks_name}/param_estimation")), width = 8, height = 8)
#   heatMap(mapRIL5, lmax = l)
#   dev.off()
#   return(NULL)
#   
# }

profileMark(mapRIL5, stat.type = c("seg.dist", "prop", "dxo", "recomb"),
            layout = c(1, 5), type = "l")
heatMap(mapRIL5, lmax = 7.5)
nmar(mapRIL5)
chrlen(mapRIL5)
# mapRIL5$geno$L.4
# Split linkage group and merge with another (based on heatmap)
# mapRIL5_split <- breakCross(mapRIL5, split = list(`L.4` = "305093:50:+"), suffix = "alpha", sep = "")
# mapRIL5_merged <- mergeCross(mapRIL5_split, merge = list(`L.4` = c("L.4", "L.13"), 
#                              `L.1` = c("L.3","L.2", "L.1")))
# mapRIL5_merged <- mstmap(mapRIL5_merged, bychr = FALSE, dist.fun = "kosambi", trace = TRUE, p.value = p_value_range[selected_p])
# Remove markers from L.10,6,7,15
# heatMap(mapRIL5_merged, lmax = 30)
# nmar(mapRIL5_merged)
# chrlen(mapRIL5_merged)
# heatMap(mapRIL5_merged, lmax = 30, chr = "L.4")
# Push back markers with higher missing rates
mapRIL5_added <- pushCross(mapRIL5, type = "co.located", 
                           pars = list(miss.thresh = 0.1, max.rf = 0.3))
# Drop small linkage groups (<20 markers)
# mapRIL5_added <- pheno2geno::removeChromosomes(mapRIL5_added, 
#                                          chromosomesToBeRmv = c("L7"), 
#                                          verbose=TRUE)
# visualise
heatMap(mapRIL5_added, lmax = 8)
nmar(mapRIL5_added) %>% sum(.)
chrlen(mapRIL5_added)
# Drop small groups
# mapRIL5_added <- pheno2geno::removeTooSmallChromosomes(mapRIL5_added, 
#                                                        minNrOfMarkers=30, 
#                                                verbose=FALSE)

# reconstruct the map (but don't reassess linkage groups - p.value=2)
mapRIL5e <- jittermap(mstmap(mapRIL5_added, bychr = TRUE, trace = TRUE, dist.fun =
                    "kosambi", p.value = 2))
# mapRIL5e_merged <- mergeCross(mapRIL5e, merge = list("L1" = c("L1", "UL3")))
# mapRIL5e <- mstmap(mapRIL5e_merged, bychr = TRUE, trace = TRUE, dist.fun =
#                      "kosambi", p.value = 2)
# # Check difference in linkage group lengths
# round(chrlen(mapRIL5e) - chrlen(mapRIL5b), 5)
# # Check difference in number of markers 
# nmar(mapRIL5e) - nmar(mapRIL5b)


# markernames(mapRIL5_chrom_final)
mapRIL5_final <- mapRIL5e
nmar(mapRIL5_final) %>% sum(.)
# Fix linkage group names
names(mapRIL5_final$geno) <- paste0("LG", 1:length(names(mapRIL5_final$geno)))

# visualise
dopdf(filedate(glue::glue("{analysis_name}_pairwise_RF_LOD_heatmap"), ext = ".pdf", outdir = "plots"), width = 9, height = 8, cmd = heatMap(mapRIL5_final, lmax = 8))

dopdf(filedate(glue::glue("{analysis_name}_linkage_map"), ext = ".pdf", outdir = "plots"), width = 8, height=6, cmd =  plot.map(mapRIL5_final, horizontal = FALSE))
 
inds <- row.names(mapRIL5_final$geno[[1]]$data)
#### Save files for QTL analysis ####
# Export map
gmap_map <- imap_dfr(mapRIL5_final$geno, 
                     ~tibble(marker = markernames(mapRIL5_final, .y),
                                     chr = .y, pos=.x$map)) %>%
  write_csv(filedate(glue::glue("{analysis_name}_gmap"), 
                     outdir = "./data/qtl2_files", 
                     ext = ".csv", dateformat = FALSE))
# Export for IRILmap
gmap_map[c("marker", "pos")] %>% write_tsv(glue::glue("./QTL_results/{analysis_name}.IRILmap"), col_names = FALSE)
##### Map statistics #####
# Average distance between markers:
marker_dist <- gmap_map %>% group_by(chr) %>% 
  mutate(marker_dist = pos - lag(pos, default = 0)) %>% 
  summarise(mean(marker_dist)) %>% 
  write_csv(filedate(glue::glue("{analysis_name}_average_marker_dist"), 
                     outdir = "./data/intermediate_files", 
                     ext = ".csv"))
# Marker density per group
marker_dens <- gmap_map %>% mutate(locus=sub("_\\d+$","", marker), 
                                   transcript_pos=sub(".+_(\\d+)$","\\1", marker)) %>%
  group_by(chr, locus) %>% summarise(marker_dens=n()) %>% summarise(mean(marker_dens)) 

map_stats_table <- inner_join(marker_dist, marker_dens) %>% mutate(LG_length=qtl::chrlen(mapRIL5_final), marker_num=nmar(mapRIL5_final)) %>% 
  write_xlsx(., excel_file = glue::glue("./QTL_results/{analysis_name}_map_stats.xlsx"), 
             sheet = "map_stats", overwritesheet = TRUE)
# xlsx::write.xlsx(as.data.frame(map_stats_table), 
#                                      glue::glue("./QTL_results/{analysis_name}_map_stats.xlsx"))
# Save physical map for R/qtl2
# Load marker map (extracted from gstacks.fa)
# stacks_dir <- "../Analysis/ref_stacks/stacks2_population_04_12_2017"
# # Load marker map (extracted from gstacks.fa)
# marker_map <- read_tsv(file.path(stacks_dir,"tag_chrom.map")) %>% 
#   dplyr::rename(CHROM_POS=POS) 
# # left_join(marker_map, c("#CHROM" = "LOC")) %>% 
# # mutate(ID=paste0("tag",`#CHROM`, POS, sep="_"), `#CHROM` = CHROM, 
# #        POS=ifelse(STRAND=="+", CHROM_POS+POS, CHROM_POS-POS)) %>% 
# # dplyr::select(-one_of(colnames(marker_map)))
# 
# marker_map <- read_tsv(file.path(stacks_dir, "tag_chrom.map")) %>%
#   dplyr::rename(CHROM_POS=POS)
# phys_map <- geno_map %>% 
#   separate(marker, into = c("tag", "LOC", "Phys_POS"), remove=FALSE, convert=TRUE) %>% 
#   inner_join(marker_map, .) %>% mutate(pos=Phys_POS) %>%
#   dplyr::select(one_of(colnames(geno_map))) %>% arrange(chr, pos)
# write_csv(phys_map, "data/qtl2_files/Lentil_GBS_pmap.csv")

# Write genotype table (decide whether to filter based on the original filtration or what was filtered during the map construction)
# genos <- as.data.frame(t(vcfR::extract.gt(vcfR::read.vcfR(shared_vcf_file)))) %>% # clean_vcf_file 
  
  
genos <- as.data.frame(t(vcf_geno_phased[c(parents, RIL_samples)])) %>% set_names(., vcf_geno_phased$ID)  %>% 
  rownames_to_column("id") %>% filter(id %in% inds) %>% .[c("id", markernames(mapRIL5_final))] 
genos[genos=="./."] <- NA# 
genos %>% filter(!grepl("^IL", id)) %>%
  write_csv(filedate(glue::glue("{analysis_name}_genos"), 
                     outdir = "./data/qtl2_files", 
                     ext = ".csv", dateformat = FALSE))
# Write founder table
genos %>% filter(grepl("^IL", id)) %>%
  write_csv(filedate(glue::glue("{analysis_name}_founder_genos"), 
                     outdir = "./data/qtl2_files", 
                     ext = ".csv", dateformat = FALSE))

# save.image(filedate(analysis_name, ext = ".RData"))
# modify YAML files (need to have the templates in place)
yamls <- list.files("./data/qtl2_files/", "Lentil_GBS_LG.+\\d+dpi.yaml", full.names = TRUE)
for (y in yamls){
  system2("sed", args=c(glue::glue("-r 's/Lentil_GBS([A-z_]+.csv)/{analysis_name}\\1/g; s/pmap/# pmap/'"),
                        y),
          stdout=sub("Lentil_GBS_LG", analysis_name, y), stderr = sub("Lentil_GBS_LG", analysis_name, y))
}

save.image(filedate(glue::glue("{analysis_name}"), ext = ".RData"))
