# Load utilities and install missing packages
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# install.packages(c("scales", "DataExplorer"))
pacman::p_load(char=c("tidyverse", "car", "e1071", "paletteer", "janitor",
                      "ggforce", "ggbeeswarm", "ggdist", "patchwork",
                      "gghalves"))#, "DataExplorer"))

# Define plotting theme
pacman::p_load_gh("Mikata-Project/ggthemr")
pale_theme <- ggthemr::ggthemr("pale", text_size = 18, set_theme = FALSE)
# plot_theme <-  theme_grey(base_size=22) +
#   theme(axis.title.y=element_text(face="bold", vjust = 1.5, size=rel(0.8)),
#         axis.title.x=element_text(face="bold", vjust = 0.1, size=rel(0.8)),
#         legend.title=element_text(size=rel(0.8), hjust=0.5),
#         legend.text=element_text(size = rel(0.7),lineheight = 1.5),
#         #panel.grid.minor=element_blank(),
#         strip.text=element_text(size=rel(0.75)),
#         strip.background=element_rect(fill = "lightskyblue"))#, colour = "black", size = 0.6))
options(stringsAsFactors = FALSE)
##### Read and prepare data #####
# Read the data in from the excel files
analysis_name <- "Lentil_FT13038_Rep2" # _Rep2
pheno_data <- readxl::read_excel("./data/Phenotyping_AB_results.xlsx", 
           col_types = c("skip", "text", rep("guess", 12)), sheet = "REP 1", 
           na = c("", "NA", "NG")) 
loal_mapping <- readxl::read_excel("./data/Phenotyping_AB_results.xlsx", 
                      sheet = "LOAL_map", na = c("", "NA", "NG"), skip = 1) %>%
  select(RIL_no, Sample_ID)

audpc_data <- readxl::read_excel("data/AUDPC_HD.xlsx")


## View bivariate continuous distribution based on `Leaf_score_ratio`
# plot_boxplot(pheno_data, by = "Leaf_score_ratio")
# Automatically detect and remove outliers (see https://datascienceplus.com/identify-describe-plot-and-removing-the-outliers-from-the-dataset/)
# source("https://goo.gl/4mthoF")
# outlierKD(pheno_data, Leaf_score_ratio)

clean_data <- pheno_data %>% 
  filter(!grepl("-mock", Variety), DPI>7,Rep==2,
         across(ends_with("ratio"), ~ !is.na(.x) ) ) %>% # DPI>7, , Rep==2
                            # !grepl("ignore", Comment, ignore.case = TRUE)) %>% 
          filter_at(vars(Variety:Isolate), all_vars(!is.na(.))) %>% 
  mutate(across(Leaf_score_ratio:Score_sqrt, as.numeric)) %>%
  # mutate_at(vars(ends_with("_ratio")), funs(sqrt(.))) %>%
  # mutate_at(vars(ends_with("_ratio")), ~.*100) %>%
  mutate(Variety=sub("-.+", "", Variety)) %>% 
  left_join(loal_mapping, by = c('Variety'='RIL_no')) %>% 
  rename(Genotype=Sample_ID) %>%
  relocate(Genotype)  %>% 
  # mutate(Variety=sub("_", "-", Variety)) %>%
  mutate(DPI=factor(DPI, levels=as.character(1:4*7)))

obs_variety <- clean_data %>% count(Variety) 
mean(obs_variety$n)
clean_data %>% count(DPI) 
clean_data %>% filter(is.na(Variety))
# pheno_data %>% filter(!is.na(Leaf_score_ratio) & !is.na(Stem_score_ratio), !grepl("-mock", Variety), DPI>7) %>% 
#   filter_at(vars(Variety:Isolate), all_vars(!is.na(.))) %>% #mutate(Variety=sub("-.+", "", Variety)) %>%
#   filter(! Variety %in% loal_mapping$RIL_no) %>% count(Variety)

##### Test data for normal distribution #####
all_pheno <- clean_data %>% #filter(DPI==dpi) %>%
  
  group_by(Variety, DPI) %>% summarise(across(ends_with("_ratio"), mean)) %>% 
  pivot_wider(names_from = DPI, values_from = ends_with("_ratio"), names_sep = "_")


# all_pheno <- NULL
# 
# 
# sampling_times <- unique(clean_data$DPI)
# for (dpi in sampling_times){
#   # dpi=sampling_times[2]
#   clean_pheno <- clean_data %>% filter(DPI==dpi) %>%
#     
#     group_by(Variety) %>% summarise_at(vars(ends_with("_ratio")), mean)  %>%
#     dplyr::rename(id=Variety) %>%
#     # write_csv(., filedate(glue::glue("{analysis_name}_pheno_{dpi}dpi"), ".csv",
#     #                       "data/qtl2_files",
#     #                       dateformat = FALSE)) %>%
#     rename(across(ends_with("_ratio"), ~paste( ., dpi,sep = "_")))
#   # Remove samples that were filtered from the genotypic data
#   
#   # clean_pheno <- log_data  %>% dplyr::rename(id=Sample) # %>% filter(Sample %in% genos$rowname)
#   
#   if (is.null(all_pheno)) {
#     all_pheno <- clean_pheno
#   } else {
#     all_pheno <- full_join(clean_pheno, all_pheno)
#   }
# }
summary(all_pheno)
# test_transform <- powerTransform(all_pheno$Leaf_score_ratio_21)
# test_transform$roundlam
sw_test_output_file <- "./data/Lentil_FT13038_pheno_SW_test.xlsx"

# Check normality for each trait (after summary)
sw_test_results <- all_pheno %>% 
  pivot_longer(-Variety, names_to =  "variable_name") %>% 
  # gather(key = "variable_name", value = "value", contains("ratio")) %>% #count(variable_name)
  group_by(variable_name)  %>% 
  do(broom::tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  select(-method) # %>% write_csv(glue::glue("./data/intermediate_files/{analysis_name}_pheno_SW_test.csv"))

write_xlsx(sw_test_results,  excel_file = sw_test_output_file, 
                    sheet = glue::glue("{analysis_name}_SW"), 
           overwritesheet = TRUE)
# Visually inspect
car::qqPlot(all_pheno$Leaf_score_ratio_21)

all_pheno_sqrt <- all_pheno %>% mutate(across(contains("_ratio"), sqrt)) 
# Check normality for each trait (after summary)
sqrt_test_results <- all_pheno_sqrt %>% 
  pivot_longer(-Variety,names_to =  "variable_name") %>% 
  group_by(variable_name)  %>% 
  do(broom::tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  select(-method) # %>% write_csv(glue::glue("./data/intermediate_files/{analysis_name}_sqrt_pheno_SW_test.csv"))
write_xlsx(sqrt_test_results, excel_file = sw_test_output_file, 
                 sheet = glue::glue("{analysis_name}_sqrt_SW"),
           overwritesheet = TRUE)
sqrt_data <- clean_data %>%
  # mutate_at(vars(ends_with("_ratio")), funs(if_else(.==0, 0.0001, .))) %>%
  mutate(across(ends_with("_ratio"), sqrt)) 
# pheno_data %>% group_by(DPI, Sample, Pot_rep) %>% summarise(n=n()) %>% filter(n!=3)%>% print(., n=Inf)
# 
# pheno_data %>% filter(is.na(Variety))
# sort(unique(pheno_data$Variety))
# clean_data <- pheno_data %>% filter(!is.na(Leaf_lesion_percent), 
#                                     !grepl("ignore", Comments, ignore.case = TRUE)) %>%
#   mutate_at(vars(Leaf_lesion_percent:Stem_lesion_percent), funs(. + 0.001)) %>%
#   mutate(DPI=factor(DPI, levels=as.character(1:4*7)))
# 
# clean_data %>% filter(Leaf_lesion_percent>100.1)

##### Summary and plotting ######

# Summarise data by parents and sampling time
parents_data <- clean_data %>% filter(grepl("^IL", Genotype)) %>% 
  group_by(Genotype, DPI) %>%
  summarise(across(ends_with("_ratio"), .fns = list(mean=mean, sd=sd), na.rm = TRUE))
# xlsx::write.xlsx(as.data.frame(parents_data), "Corrected_data/Parents_sum_data.xlsx", col.names = TRUE, row.names = FALSE)
clean_data %>% count(Genotype) %>% summary()
# Summarise data by sample and sampling time
samples_sum_data <- clean_data %>% # filter(!grepl("C\\d.+", Sample)) %>% 
  group_by(Genotype, DPI) %>%
  summarise(across(ends_with("_ratio"), .fns = list(mean=mean, sd=sd), na.rm = TRUE))
# Export to excel
write_xlsx(samples_sum_data, "./data/Phenotyping_AB_res.xlsx", "sqrt_sum_data")

summary(samples_sum_data)
samples_sum_data %>% filter(Leaf_score_ratio_sd>0.4) %>% print(n=Inf)
  # Save to excel 
# xlsx::write.xlsx(as.data.frame(samples_sum_data), "./Corrected_data/All_samples_sum_data.xlsx", 
                   # col.names = TRUE, row.names = FALSE)

# Summarise data (parents and population mean)
sum_data <- clean_data %>% filter(!grepl("^IL", Genotype)) %>% group_by(DPI)  %>%
  summarise(across(ends_with("_ratio"), 
                   .fns = list(mean=mean, sd=sd), na.rm = TRUE)) %>%
  mutate(Genotype="RIL5") %>% bind_rows(parents_data) %>% 
  mutate(Genotype=fct_relevel(factor(Genotype), "RIL5"))
# Summarise overall population mean
mean_data <- clean_data %>% group_by(DPI)  %>%
  summarise(across(ends_with("_ratio"), 
                   .fns = list(mean=mean, sd=sd), na.rm = TRUE))

plot_data <-  clean_data %>% # filter(Rep==2, DPI!="7") %>% 
  mutate(Genotype=map_chr(Genotype, ~ifelse(!grepl("^IL", .x), "RIL5", .x)),
         Genotype=fct_relevel(factor(Genotype), "RIL5"))

sampling_points <- levels(plot_data$DPI)

# define colours ####
par_pal <- paletteer_d("colorblindr::OkabeIto")[c(8, 5,6)]
# paletteer_d("ochRe::parliament", 3)[c(3,1,2)]

plot_cols <- list(
  genos=levels(plot_data$Genotype) %>% 
    setNames(as.character(par_pal), .), #%>%  paletteer_d("ggthemes::calc", length(.)))
   #  c(., c("Population"="black")), # awtools::mpalette
  dpi=sampling_points %>% 
    setNames(as.character(paletteer_d("ggthemes::calc", length(.))), .)
)
# paletteer_d("colorblindr::OkabeIto")
#### Leaf lesion distribution plots ####

  
#  c(brewer.pal(n = length(levels(sum_data$Variety)), name = "Set1"), "black") # [c(1,3)]

# Plot density plots
# ggplot(clean_data, aes(x=Leaf_score_ratio)) + # , colour=DPI
#   geom_density(size=1, colour="darkslategrey") + 
#   scale_color_manual(values=legend_cols, 
#                      guide=guide_legend(override.aes = list(linetype="solid", size=1))) +
#   labs(x="Leaf score ratio", y="Density", colour="Genotype") +
#   # scale_x_log10() + # xlim(0,1) +
#   # geom_vline(data=mean_data, aes(xintercept=Leaf_lesion_percent),
#   #            size=0.5, linetype="solid", colour=legend_cols) + #"black"
#   geom_vline(data=sum_data, aes(xintercept=Leaf_score_ratio_mean, colour=Variety),
#              linetype="dashed", size=0.75) +
#   facet_grid(DPI ~ .) +
#    plot_theme("grey", 22)   # +
#   # theme(legend.text=element_text(lineheight=.8),
#   #       legend.key.height=unit(0.5, "cm"))
# ggsave(filedate("Leaf_lesion_coverage_density_FT15124", ".pdf", "output_plots"), width = 10, height=10)

# violin plots
# quant_med <- 0.5
# quant_width <- 0.15
# quantiles <- c(quant_med-quant_width, quant_med, quant_med+quant_width)
# ggplot(samples_sum_data, aes(x=DPI, y=Leaf_score_ratio_mean)) + # , colour=DPI
  ggplot(plot_data, aes(x=DPI, y=Leaf_score_ratio)) + #
  geom_violin(fill="#C3D6DF", color = NA, scale="count") + # scale_y_log10() + #draw_quantiles = quantiles
  geom_point(data=sum_data, aes(y=Leaf_score_ratio_mean, x=DPI, colour=Genotype), 
             size=2.5) + 
  ggbeeswarm::geom_quasirandom(mapping = aes(colour=Genotype), 
                     size = 2.5, alpha = .4, width = .2, position = ) + 
  scale_fill_manual(values=plot_cols$dpi) +
  scale_color_manual(values=plot_cols$genos) + #, 
                 #    guide=guide_legend(override.aes = list(linetype="solid", size=3))) +
  # geom_errorbar(data = sum_data, mapping = aes(x=DPI, y=Leaf_lesion_percent_mean, 
  #                                    ymin=Leaf_lesion_percent_mean - Leaf_lesion_percent_sd,
  #                                    ymax=Leaf_lesion_percent_mean + Leaf_lesion_percent_sd,
  #                                    colour=Sample), width=0.1) +
  labs(y="Disease score ratio", x="Sampling time (dpi)", colour="Genotype" +
         title = "Leaf") +
  pale_theme$theme # plot_theme("grey", 22) 
ggsave(filedate("Leaf_score_mean_ratio_violin_FT13038_Rep2" , ".pdf", "plots"), width = 10, height=7)

# raincloud plots
leaf_raincloud <- ggplot(plot_data, aes(x=DPI, y=Leaf_score_ratio)) +
  ggdist::stat_halfeye(
    fill="#C3D6DF",
  adjust = .5, 
  width = .6, 
  .width = 0, 
  justification = -.2, 
  point_colour = NA) + 
  geom_boxplot( width = .125, outlier.shape = NA) +
  geom_half_point(aes(colour=Genotype, shape=Genotype), 
                  size = 3, side = "l", range_scale = .2, alpha = .45, 
                  position = position_jitter(seed = 1, width = .025)) +
  geom_point(data = sum_data,
    mapping = aes(x=DPI, y=Leaf_score_ratio_mean, 
                  colour=Genotype, shape=Genotype),
    size = 3) +
  coord_flip() + # xlim = c(1.2, NA), ylim=c(0,1), clip = "off"
  scale_y_continuous(expand = expansion(mult = c(0.01, .05))) + 
  # scale_fill_manual(values=plot_cols$dpi) +
  scale_color_manual(values=plot_cols$genos,
                     guide=guide_legend(override.aes = list(alpha=1, size=3))) +
  labs(y="Disease score ratio", x="Sampling time (dpi)", colour="Genotype",
         title = "Leaf lesions") + 
  pale_theme$theme # plot_theme("grey", 22) 
# ggsave(filedate("FT13038_Rep2_Leaf_score_ratio_raincloud" , ".pdf", "plots"), width = 10, height=7)



# ggsave(filedate("log_Leaf_lesion_coverage_violin" , ".pdf", "output_plots"), width = 10, height=7)

#### Stem lesion distribution ####

# raincloud plots
stem_raincloud <- ggplot(plot_data, aes(x=DPI, y=Stem_score_ratio)) +
  ggdist::stat_halfeye(
    fill="#C3D6DF",
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA) + 
  geom_boxplot(
    width = .125, 
    outlier.shape = NA
  ) +
  geom_half_point(aes(colour=Genotype, shape=Genotype), 
                  size = 3, side = "l", range_scale = .2, alpha = .45, 
                  position = position_jitter(seed = 1, width = .025)) +
  geom_point(data = sum_data,
             mapping = aes(x=DPI, y=Stem_score_ratio_mean, 
                           colour=Genotype, shape=Genotype),
             size = 3) +
  coord_flip() + # xlim = c(1.2, NA), ylim=c(0,1), clip = "off"
  scale_y_continuous(expand = expansion(mult = c(0.01, .05))) + 
  # scale_fill_manual(values=plot_cols$dpi) +
  scale_color_manual(values=plot_cols$genos,
                     guide=guide_legend(override.aes = list(alpha=1, size=3))) +
  labs(y="Disease score ratio", x="Sampling time (dpi)", colour="Genotype",
       title = "Stem lesions") +
  pale_theme$theme # plot_theme("grey", 22) 

# ggsave(filedate("FT13038_Rep2_Stem_score_ratio_raincloud" , ".pdf", "plots"), width = 10, height=7)
# create a composite plot
leaf_raincloud / stem_raincloud + plot_annotation(tag_levels = 'a')


ggsave(filedate("FT13038_Rep2_Leaf_Stem_score_ratio_raincloud_patchwork" , ".pdf", "plots"), width = 10, height=10)
# plot as facets
plot_data <- clean_data %>% # filter(Rep==2, DPI!="7") %>% 
  mutate(Genotype=map_chr(Genotype, ~ifelse(!grepl("^IL", .x), "RIL5", .x)),
         Genotype=factor(Genotype)) %>% 
  pivot_longer(ends_with("ratio"), names_to = "Trait", values_to = "Score_ratio") %>% 
  mutate(Trait=sub("_.+", "", Trait))

ggplot(plot_data, aes(x=DPI, y=Score_ratio)) +
  ggdist::stat_halfeye(
    fill="#C3D6DF",
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .2, 
    outlier.shape = NA
  ) +
  geom_point(
    mapping = aes(colour=Genotype), 
    size = 2.5, alpha = .45,
    position = position_jitter(
      seed = 1, width = .075
    )
  ) + 
  coord_flip(xlim = c(1.2, NA), ylim=c(0,1), clip = "off") +
  # scale_fill_manual(values=plot_cols$dpi) +
  scale_color_manual(values=plot_cols$genos,
                     guide=guide_legend(override.aes = list(alpha=1, size=3))) +
  labs(y="Disease score ratio", x="Sampling time (dpi)", colour="Genotype") +
  facet_wrap(~Trait) + 
  pale_theme$theme 
ggsave(filedate("FT13038_Rep2_Stem_Leaf_score_ratio_raincloud" , ".pdf", "plots"), width = 10, height=7)

# Plot density plots
# ggplot(clean_data, aes(x=Stem_score_ratio)) + # , colour=DPI
#   geom_density(size=1, colour="darkslategrey") + 
#   scale_color_manual(values=legend_cols, 
#                      guide=guide_legend(override.aes = list(linetype="solid", size=1))) +
#   labs(x="Stem lesion coverage (%)", y="Density", colour="Genotype") +
#   # scale_x_log10() + # xlim(0,1) +
#   # geom_vline(data=mean_data, aes(xintercept=Stem_lesion_percent),
#   #            size=0.5, linetype="solid", colour=legend_cols) + #"black"
#   geom_vline(data=sum_data, aes(xintercept=Stem_lesion_percent_mean, colour=Sample),
#              linetype="dashed", size=0.75) +
#   facet_grid(DPI ~ .) +
#   plot_theme("grey", 22)   # +
# # theme(legend.text=element_text(lineheight=.8),
# #       legend.key.height=unit(0.5, "cm"))
# ggsave(filedate("Stem_lesion_coverage_density_FT15124", ".pdf", "output_plots"), width = 10, height=10)

# violin plots

ggplot(samples_sum_data, aes(x=DPI, y=Stem_score_ratio_mean)) + # , colour=DPI
  geom_violin() + # scale_y_log10() + #draw_quantiles = quantiles
  geom_point(data=sum_data, aes(y=Stem_score_ratio_mean, x=DPI, colour=Variety), 
             size=2.5) + 
  scale_color_manual(values=legend_cols, 
          guide=guide_legend(override.aes = list(linetype="solid", size=3))) +
  # geom_errorbar(data = sum_data, mapping = aes(x=DPI, y=Stem_lesion_percent_mean, 
  #                                    ymin=Stem_lesion_percent_mean - Stem_lesion_percent_sd,
  #                                    ymax=Stem_lesion_percent_mean + Stem_lesion_percent_sd,
  #                                    colour=Sample), width=0.1) +
  labs(y="Transformed stem score ratio", x="Sampling time (dpi)", 
       colour="Genotype") +
  plot_theme("grey", 22) 
ggsave(filedate("Stem_score_ratio_mean_sqrt_violin_FT13038" , ".pdf", "output_plots"), width = 10, height=7)


##### Perform Chi-sq analysis on progeny ####

sampling_days <- levels(clean_data$DPI)
chi_table <- NULL
segregations <- data.frame(Sus=c(1,1,1, 1), Res=c(1,1.5,2,3))
quantiles_df <- data.frame(Q1=c(0.35, 0.5, 0.65), Q2=c(0.45, 0.6, 0.75))
for (d in sampling_days){
  for (j in colnames(quantiles_df)){
    for (i in 1:nrow(segregations)){
      progeny_data <- clean_data %>% filter(grepl("C\\d.+", Sample), DPI==d) 
      Leaf_lesion_quantiles <- quantile(progeny_data$Leaf_lesion_percent, probs = quantiles_df[,j])
      progeny_data$Leaf_lesion_Resistance <- NA
      progeny_data$Leaf_lesion_Resistance[progeny_data$Leaf_lesion_percent<=Leaf_lesion_quantiles[1]] <- TRUE
      progeny_data$Leaf_lesion_Resistance[progeny_data$Leaf_lesion_percent>=Leaf_lesion_quantiles[3]] <- FALSE
      # Summarise Res/Sus
      cont_table <- table(progeny_data$Leaf_lesion_Resistance)
      seg <- as.numeric(segregations[i,])
      # Chi sqaure test
      (Xsq <- chisq.test(cont_table, p=seg, rescale.p = TRUE))
      chi_table <- tibble(Dpi=d, Obs_Sus=cont_table[1], Obs_Res=cont_table[2], 
                          Exp_seg=paste(seg, collapse = ":"),
                          Exp_Res=round(Xsq$expected[1]), Exp_Sus=round(Xsq$expected[2]), 
                          Res_thresh=sprintf("<%.3f; >%.3f", Leaf_lesion_quantiles[1], Leaf_lesion_quantiles[3]), Median_at=names(Leaf_lesion_quantiles[2]),
                          Chi_sq=Xsq$statistic, p_value=Xsq$p.value) %>%
        bind_rows(chi_table, .)
    }
    
  }
}

chi_table <- chi_table %>% mutate(Signif=p_value>0.05)
xlsx::write.xlsx(as.data.frame(chi_table), filedate("Leaf_lesion_coverage_Chi_sq_test", ".xlsx"),
                 row.names = FALSE)
# Xsq


