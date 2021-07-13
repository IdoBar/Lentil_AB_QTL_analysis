# Load utilities and install missing packages
devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# install.packages(c("scales", "DataExplorer"))
pacman::p_load(char=c("tidyverse", "car", "e1071", "paletteer",
                      "janitor",
                      "ggforce", "ggbeeswarm", "ggdist",
                      "patchwork",
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
analysis_name <- "Lentil_FT13038_AUDPC" # _Rep2
pheno_file <- "./data/Phenotyping_AB_results.xlsx"
# pheno_data <- readxl::read_excel(pheno_file,
#            col_types = c("skip", "text", rep("guess", 12)), sheet = "REP 1 & 2",
#            na = c("", "NA", "NG"))
pheno_data <- readxl::read_excel(pheno_file,
                                sheet = "AUDPC",
                                 na = c("", "NA", "NG")) %>% 
  pivot_longer(contains("score"), 
               names_to = c("Trait", "DPI"), 
               names_pattern = "(.+?)_([0-9]+)",
               values_to = "Score") %>% 
  # pivot_wider(names_from = "Trait", values_from = "Score") %>% 
  rename(Sample_ID=`Sample ID`) %>% 
  mutate(id=gsub("_", "-", Sample_ID ), Genotype=sub("ILW_", "ILWL", Sample_ID),
         RIL=sub("-inoculated", "", RIL))
loal_mapping <- readxl::read_excel("./data/Phenotyping_AB_results.xlsx", 
                      sheet = "LOAL_map", na = c("", "NA", "NG"), skip = 1) %>%
  select(RIL_no, Sample_ID)

# audpc_data <- readxl::read_excel("data/AUDPC_HD.xlsx", sheet = "audpc_sqrt") %>% 
#   rename(AUDPC_sqrt=`Overall AUDPC SQRT`) %>% 
#   mutate(id=gsub("_", "-", Sample_ID ), Genotype=Sample_ID, DPI="28")

audpc_data <- pheno_data %>% 
  filter(Trait=="Overall_AUDPC_score", Rep==2) %>% 
  mutate(AUDPC_sqrt = sqrt(Score), DPI="28") #%>% 
  # mutate(id=gsub("_", "-", Sample_ID ), Genotype=sub("ILW_", "ILWL", Sample_ID))


loal_mapping %>% left_join(pheno_data) %>% filter(RIL_no!=RIL) %>% 
  select(RIL_no, Sample_ID, RIL) %>% distinct()



##### Test data for normal distribution #####

sw_test_output_file <- "./data/Lentil_FT13038_pheno_SW_test.xlsx"

# Check normality for each trait (after summary)
  sw_test_results <- pheno_data %>% filter(Rep==2) %>% 
    mutate(Score_sqrt=sqrt(Score), Score_ratio=Score/100) %>% 
    group_by(Trait, DPI) %>% 
    nest() %>% mutate(data=setNames(data, paste("Rep2",Trait,  DPI, sep = "_")),
      test = map(data, ~shapiro.test(.x$Score_sqrt)),
    tidied = map(test, broom::tidy)) %>% unnest(tidied) 
    
    
  
  # sw_test_results <- broom::tidy(shapiro.test(audpc_data$AUDPC_sqrt))

write_xlsx(sw_test_results,  excel_file = sw_test_output_file, 
                    sheet = glue::glue("{analysis_name}_SW"), 
           overwritesheet = TRUE)
# Visually inspect
# ratio transformation
iwalk(sw_test_results$data, 
      ~ car::qqPlot(.x$Score_ratio, ylab=sub("_score.*?(_[0-9]+)", "_ratio\\1", .y)))
# sqrt transformation
iwalk(sw_test_results$data, 
      ~ car::qqPlot(.x$Score_sqrt, ylab=sub("_score.*?(_[0-9]+)", "_sqrt\\1", .y)))

car::qqPlot(audpc_data$AUDPC_sqrt)


##### Summary and plotting ######
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
    mutate(Genotype=sub("ILW_", "ILWL", Sample_ID)) %>%
    relocate(Genotype)  %>% 
    # mutate(Variety=sub("_", "-", Variety)) %>%
    mutate(DPI=factor(DPI, levels=as.character(1:4*7)))



# clean_data <- audpc_data %>% mutate(Genotype=Sample_ID, DPI=28)
# Summarise data by parents and sampling time
# parents_data <- clean_data %>% filter(grepl("^IL", Genotype)) %>% 
#   group_by(Genotype, DPI) %>%
#   summarise(across(ends_with("_ratio"), .fns = list(mean=mean, sd=sd), na.rm = TRUE)) %>% 
#   left_join(audpc_data %>% select(Genotype, DPI, AUDPC_sqrt))

# xlsx::write.xlsx(as.data.frame(parents_data), "Corrected_data/Parents_sum_data.xlsx", col.names = TRUE, row.names = FALSE)
clean_data %>% count(Genotype) %>% summary()
# Summarise data by sample and sampling time
samples_sum_data <- clean_data %>% # filter(!grepl("C\\d.+", Sample)) %>% 
  group_by(Genotype, DPI) %>%
  summarise(across(ends_with("_ratio"), 
                   .fns = list(mean=mean, sd=sd), na.rm = TRUE)) %>% 
  left_join(audpc_data %>% select(Genotype, DPI, AUDPC_sqrt))
# Export to excel
write_xlsx(samples_sum_data, "./data/Phenotyping_AB_res.xlsx", "AUDPC_sqrt_sum_data")


# Summarise data (parents and population mean)
# sum_data <- clean_data %>% filter(!grepl("^IL", Genotype)) %>%  group_by(DPI)  %>%
#   summarise(across(ends_with("_sqrt"), 
#                    .fns = list(mean=mean, sd=sd), na.rm = TRUE)) %>%
#   mutate(Genotype="RIL5") %>% bind_rows(parents_data) %>% 
#   mutate(Genotype=fct_relevel(factor(Genotype), "RIL5"))
# # Summarise overall population mean
# mean_data <- clean_data %>% group_by(DPI)  %>%
#   summarise(across(ends_with("_sqrt"), 
#                    .fns = list(mean=mean, sd=sd), na.rm = TRUE))




# paletteer_d("colorblindr::OkabeIto")



#### Trait distribution ####
plot_data <- clean_data %>% # filter(Rep==2, DPI!="7") %>% 
  mutate(Genotype=map_chr(Genotype, ~ifelse(!grepl("^IL", .x), "RIL5", .x)),
         Genotype=fct_relevel(factor(Genotype), "RIL5")) 


sum_data <- plot_data %>%  group_by(Genotype, DPI)  %>%
  summarise(across(ends_with("_ratio"), 
                   .fns = list(mean=mean, sd=sd), na.rm = TRUE)) 

sampling_points <- levels(plot_data$DPI)

# define colours ####
par_pal <- paletteer_d("colorblindr::OkabeIto")[c(8, 5,6)]
# paletteer_d("ochRe::parliament", 3)[c(3,1,2)]

plot_cols <- list(
  genos=levels(plot_data$Genotype) %>% 
    setNames(as.character(par_pal), .) #%>%  paletteer_d("ggthemes::calc", length(.)))
  #  c(., c("Population"="black")), # awtools::mpalette
  # dpi=sampling_points %>% 
  #   setNames(as.character(paletteer_d("ggthemes::calc", length(.))), .)
)

# raincloud plots
leaf_raincloud <- ggplot(plot_data, aes(x=DPI, y=Leaf_score_ratio)) +
    ggdist::stat_halfeye(
        fill="#C3D6DF",
      adjust = .5, 
      width = .6, 
      .width = 0, 
      justification = -.2, 
      point_colour = NA) + 
    geom_boxplot( width = .15, outlier.shape = NA) +
    geom_half_point(aes(colour=Genotype, shape=Genotype), 
                                       size = 3, side = "l", range_scale = .2, alpha = .45, 
                                       position = position_jitter(seed = 1, width = .025)) +
    geom_point(data = sum_data,
                    mapping = aes(x=DPI, y=Leaf_score_ratio_mean, 
                                                     colour=Genotype, shape=Genotype),
                    size = 3.5) +
   coord_flip() + # xlim = c(1.2, NA), ylim=c(0,1), clip = "off"
    scale_y_continuous(expand = expansion(mult = c(0.01, .05))) + 
    # scale_fill_manual(values=plot_cols$dpi) +
    scale_color_manual(values=plot_cols$genos,
                                             guide=guide_legend(override.aes = list(alpha=1, size=3))) +
    labs(y="Disease score (ratio)", x="Sampling time (dpi)", colour="Genotype",
                   title = "Leaf lesions") + 
    pale_theme$theme # plot_theme("grey", 22)

stem_raincloud <- ggplot(plot_data, aes(x=DPI, y=Stem_score_ratio)) +
  ggdist::stat_halfeye(
    fill="#C3D6DF",
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA) + 
  geom_boxplot(width = .15, outlier.shape = NA) +
  geom_half_point(aes(colour=Genotype, shape=Genotype), 
                  size = 3, side = "l", range_scale = .2, alpha = .45, 
                  position = position_jitter(seed = 1, width = .025)) +
  geom_point(data = sum_data,
             mapping = aes(x=DPI, y=Stem_score_ratio_mean, 
                           colour=Genotype, shape=Genotype),
             size = 3.5) +
  coord_flip() + # xlim = c(1.2, NA), ylim=c(0,1), clip = "off"
  scale_y_continuous(expand = expansion(mult = c(0.0025, .05))) + 
  # scale_fill_manual(values=plot_cols$dpi) +
  scale_color_manual(values=plot_cols$genos,
                     guide=guide_legend(override.aes = list(alpha=1, size=3))) +
  labs(y="Disease score (ratio)", x="Sampling time (dpi)", colour="Genotype",
       title = "Stem lesions") +
  pale_theme$theme



# AUDPC plot #
audpc_plot_data <-  audpc_data %>% # filter(Rep==2, DPI!="7") %>% 
  mutate(Genotype=map_chr(Genotype, ~ifelse(!grepl("^IL", .x), "RIL5", .x)),
         Genotype=fct_relevel(factor(Genotype), "RIL5"))

audpc_sum_data <- audpc_plot_data %>%  group_by(Genotype, DPI)  %>%
  summarise(across(c("Score", "AUDPC_sqrt"), 
                   .fns = list(mean=mean, sd=sd), na.rm = TRUE)) 

plot_cols <- list(
  genos=levels(audpc_plot_data$Genotype) %>% 
    setNames(as.character(par_pal), .) 
)

audpc_raincloud <- ggplot(audpc_plot_data, aes(x=DPI, y=Score)) +
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
  geom_point(data = audpc_sum_data,
             mapping = aes(x=DPI, y=Score_mean, 
                           colour=Genotype, shape=Genotype),
             size = 3.5) +
  coord_flip() + # xlim = c(1.2, NA), ylim=c(0,1), clip = "off"
  scale_y_continuous(expand = expansion(add = c(35, 100)), 
                     limits = c(0, max(audpc_plot_data$Score))) + 
  # scale_fill_manual(values=plot_cols$dpi) +
  scale_color_manual(values=plot_cols$genos,
                     guide=guide_legend(override.aes = list(alpha=1, size=3))) +
  labs(y="AUDPC", x="Sampling time (dpi)", colour="Genotype",
       title = "AUDPC disease response") +
  pale_theme$theme +
  theme(panel.grid.major.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank())
  

# ggsave(filedate("FT13038_Rep2_Stem_score_ratio_raincloud" , ".pdf", "plots"), width = 10, height=7)

ggsave(filedate("FT13038_Rep2_AUDPC_raincloud" , ".pdf", "plots"), width = 10, height=3)
# create a composite plot
leaf_raincloud / stem_raincloud / audpc_raincloud + 
  plot_layout(heights = c(2,2,1), guides = 'collect') +
  plot_annotation(tag_levels = 'a')
ggsave(filedate("FT13038_Rep2_Stem_Leaf_score_ratio_AUDPC_raincloud" , ".pdf", "plots"), 
       width = 10, height=10)


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

pheno_data_28dpi <- read_csv("./data/qtl2_files/Lentil_FT13038_pheno_28dpi.csv") %>% 
  left_join(audpc_data %>% select(id, AUDPC_sqrt)) %>% 
  write_csv("./data/qtl2_files/Lentil_FT13038_pheno_28dpi.csv")
# pheno_data_28dpi %>% filter(is.na(AUDPC_sqrt))
