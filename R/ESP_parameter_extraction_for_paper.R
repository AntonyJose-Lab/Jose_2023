## This program generates many plots for analyzing systematically simulated ESP systems.
## Date: 2022_8_31
## Import libraries you might use.
library(tidyverse) # required for ggplot2, dplyr, and others.
library(grid) # required for adding annotations.

## Set working directory so that remaining paths are easy to set.
experiment <- "2022_8_31_ESP_systems"

root <- "/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/"
setwd(paste0(root, "code/R/", experiment))
dir.create(paste0(root, "analyses/R/Figures/", experiment))
dir.create(paste0(root, "analyses/R/Tables/", experiment))
path_to_figure_outputs <- paste0(root, "analyses/R/Figures/", experiment)
path_to_table_outputs <- paste0(root, "analyses/R/Tables/", experiment)

## Import the data from behaviorspace in netlogo. Notice 'skip = 6' to skip the first 6 lines, which are header information.
imported_df <- data.frame(read.csv(paste0(root, "data/2020_5_13_ESP_systems_2-16_mols.csv"), header = TRUE, skip = 6))

grouped_df <- imported_df %>% group_by(molecule.kinds, perturb.kind, perturb.phase) %>%
              summarise (n = n())
ESP_per_mlp <- grouped_df$n[1] ## Repeated runs of ESP systems per molecule kinds (i.e., 1000).

## The number of runs with perturb.kind "none" is 5x per ESP becasue they are all essentially repeats for the 5 phases because no perturbation was done.
## However, they wont give the same values because the random number seeds are different for each of the phases.
## The different robust systems recovered gives an idea of the expected variation because of differences in the random seed alone.

### Set plot themes
my_theme <- theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_line(linetype = 1, size = 0.5, lineend="round"),
                  axis.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  axis.ticks = element_line(size = 0.5), axis.ticks.length = unit(.2, "cm"),
                  axis.text = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10))

### make a plot looking at frequency of robust systems in every mode per number of molecule kinds
### process data frame.
imported_df_a <- imported_df %>% filter(perturb.kind == "none") %>% group_by(molecule.kinds) %>%
  summarise(fr_robust = sum(robust.measure)/n())

### to change labels of the facets and to reorder the factor perturb.kind by changing levels.
# Code below is not usable if only the wild type has been culled for plotting in imported_df_a above
# labels <- c(none = "wild type", gof = "gain of function", lof = "loss of function")
# imported_df_a$perturb.kind <- factor(imported_df_a$perturb.kind, levels = c("none", "lof", "gof"))

### This plot looks at molecule kinds vs fraction robust systems generaed with 50% chance of linkage and +ve regulation
zero <- ggplot(data = imported_df_a) + geom_point(mapping = aes(x = molecule.kinds, y = fr_robust)) +
  ## facet_wrap(~ perturb.kind, labeller = labeller(perturb.kind = labels), nrow = 1) + ## Uncomment for plotting lof and gof perturbations.
  my_theme +
  xlab("kinds of molecules") + ylab("fraction lasting for 250 generations") +
  scale_x_continuous(limits=c(0, 17), breaks = c(0, 4, 8, 12, 16), expand = c(0,0)) +
  scale_y_continuous(limits=c(0, 1), breaks = c(0, 0.5, 1), expand = c(0, 0))
print(zero)
ggsave(file = paste0(path_to_figure_outputs, "/zero.eps"), height = 2.7, width = 2, device = "eps", units = "in", dpi = 300)

## Get the subset that had a system surviving for 250 generations (robust.measure > 0)
## Use 'mutate' to add a new column that just registers with a 1 vs 0 if there was TEI.
robust_only_df <- imported_df %>%
                          filter(robust.measure != 0) %>%
                          mutate(binary_TEI = if_else(system.change > 0, 1, 0))

### Write a bunch of files with characteristics you care about.
dir.create("csv_files")
## names(robust_only_df) <- NULL ### eliminates column names.
write.csv(robust_only_df, file = paste0(path_to_table_outputs, "/2020_5_14_robust_run_parameters_all.csv"), row.names = FALSE)
TEI_frequency <- data.frame()
for (i in 0:4) {
tmp <- robust_only_df %>% filter(system.change == i)
TEI_frequency <- rbind(TEI_frequency, c(i, nrow(tmp)))
write.table(tmp, file = paste0(path_to_table_outputs, "/2020_5_14_robust_run_parameters_", i, "_tei.csv"), sep = ",", row.names = FALSE)
}
colnames(TEI_frequency) <- c("TEI","Frequency")

max_turtles <- max(robust_only_df$turtles.at.end)
random_20_per_turtles <- data.frame()
for (i in 2:max_turtles) {
  tmp <- robust_only_df %>% filter(turtles.at.end == i)
  write.table(tmp, file = paste0(path_to_table_outputs, "/2020_5_14_robust_run_parameters_", i, "_turtles.csv"), sep = ",", row.names = FALSE)
  set.seed(30) ## Needed here within the loop to ensure that sample_n uses the same seeds every time.
  tmp <- robust_only_df %>% filter(turtles.at.end == i) %>% sample_n(20)
  random_20_per_turtles <- rbind(random_20_per_turtles, tmp)
}
write.table(random_20_per_turtles, file = paste0(path_to_table_outputs, "/2020_5_14_robust_random_20_per_final_turtles_seed_30.csv"), sep = ",", row.names = FALSE)

max_neg_links <- max(robust_only_df$negative.links)
for (i in 2:max_turtles) {
  for (j in 1:max_neg_links) {
  tmp <- robust_only_df %>% filter(turtles.at.end == i & negative.links == j)
  write.table(tmp, file = paste0(path_to_table_outputs, "/2020_5_14_robust_run_parameters_", i, "_turtles_", j, "_neg_links.csv"), sep = ",", row.names = FALSE)
  }
}

## To find robust systems to illustrate how epistasis could be misleading.
robust_5_turtles <- robust_only_df %>% filter(turtles.at.end == 5 & binary_TEI == 0)
write.table(robust_5_turtles, file = paste0(path_to_table_outputs, "/2020_5_14_robust_5_turtles_0_tei.csv"), sep = ",", row.names = FALSE)

## generate a processed data frame in which you have performed some basic calculations you would like to perform.
processed_df_1 <- robust_only_df %>%
                group_by(molecule.kinds, perturb.kind, perturb.phase) %>%
                summarise(fraction_robust = n()/ESP_per_mlp,
                median_turtles = median(turtles.at.end), max_turtles = max(turtles.at.end),
                median_positive = median(positive.links), max_positive = max(positive.links),
                median_negative = median(negative.links), max_negative = max(negative.links),
                median_neg_reg = median (negative.links/(negative.links + positive.links)),
                median_TEI = median(system.change), max_TEI = max(system.change), fr_TEI = sum(binary_TEI)/n())


most_max_turtles <- max(processed_df_1$max_turtles) ## needed to scale the Y-axis appropriately.

## This plot looks at fraction robust vs median size of system per 1000 different ESP systems faceted by molecule kinds.
one <- ggplot(data = processed_df_1) +
  geom_point(mapping = aes(y = median_turtles, x = fraction_robust, color = median_negative, size = median_positive)) +
  facet_wrap(~ molecule.kinds, nrow = 3) + ## Uncomment for plotting lof and gof perturbations.
  my_theme +
  scale_color_gradient(high="red", low="blue") +
  xlab("fraction robust") + ylab("median number of entities") +
  guides(size = guide_legend(reverse = TRUE)) + labs(size = "pos. reg.", color = "neg. reg.") +
  scale_x_continuous(limits=c(0,1), breaks = c(0, 1), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,most_max_turtles), breaks = c(0, 5, most_max_turtles), expand=c(0,0))
print(one)
ggsave(file = paste0(path_to_figure_outputs, "/one.eps"), height = 7 , width = 7, device = "eps", units = "in", dpi = 300)

two <- ggplot(data = processed_df_1) +
   geom_point(mapping = aes(y = max_turtles, x = fraction_robust, color = max_negative, size = max_positive)) +
  ## facet_wrap(~ molecule.kinds, nrow = 3) + ## Uncomment for plotting lof and gof perturbations.
   my_theme +
   scale_color_gradient(high="red", low="blue") +
   xlab("fraction robust") + ylab("maximum number of entities") +
   scale_x_continuous(limits=c(0,1), breaks = c(0, 1), expand=c(0,0)) +
   scale_y_continuous(limits=c(0,most_max_turtles + 1), breaks = c(0, 5, most_max_turtles), expand=c(0,0))
print(two)
ggsave(file = paste0(path_to_figure_outputs, "/two.eps"), height = 7 , width = 7, device = "eps", units = "in", dpi = 300)
### Note to get this figure with transparency, you need pdf because r doesnt support semi-transparency with eps.

### Setting phase as factor so that I can use different color-blind friendly colors (cbp below) to label the different discrete phases of perturbation.
processed_df_1$perturb.phase <- as.factor(processed_df_1$perturb.phase)
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#555555", "#999999")
### Setting the order of the data.
processed_df_1$perturb.kind <- factor(processed_df_1$perturb.kind, levels = c("none", "lof", "gof"))

processed_df_1_wt <- processed_df_1 %>% filter (perturb.kind == "none")
## This plot looks at fraction of TEI vs kind of perturbation and phase.
three_a <- ggplot(data = processed_df_1_wt) +
  geom_bar(mapping = aes(y = fr_TEI, x = molecule.kinds, fill = perturb.phase), stat="identity", width = 0.8, position = position_dodge(width = 0.8)) +
  #facet_wrap(~ molecule.kinds, nrow = 1) +
  my_theme +
  #theme(axis.text.x = element_text()) +
  theme(legend.position = c(0.1, 0.7)) +
  xlab("number of entities at gen. #1") +
  ylab("fr. showing heritable epigenetic change") +
  scale_x_continuous(expand=c(0,0), breaks = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)) +
  scale_y_continuous(limits=c(0, 0.3), breaks = c(0, 0.15, 0.3), expand = c(0,0)) +
  labs(fill = "phase") +
  scale_fill_manual(values = cbp)
print(three_a)
ggsave(file = paste0(path_to_figure_outputs, "/three_a.eps"), height = 3, width = 7, device = "eps", units = "in", dpi = 300)

processed_df_1_wt <- processed_df_1 %>% filter (perturb.kind == "lof")
## This plot looks at fraction of TEI vs kind of perturbation and phase.
three_b <- ggplot(data = processed_df_1_wt) +
  geom_bar(mapping = aes(y = fr_TEI, x = molecule.kinds, fill = perturb.phase), stat="identity", width = 0.8, position = position_dodge(width = 0.8)) +
  #facet_wrap(~ molecule.kinds, nrow = 1) +
  my_theme +
  #theme(axis.text.x = element_text()) +
  theme(legend.position = c(0.1, 0.7)) +
  xlab("number of entities at gen. #1") +
  ylab("fr. showing heritable epigenetic change") +
  scale_x_continuous(expand=c(0,0), breaks = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)) +
  scale_y_continuous(limits=c(0, 0.3), breaks = c(0, 0.15, 0.3), expand = c(0,0)) +
  labs(fill = "phase") +
  scale_fill_manual(values = cbp)
print(three_b)
ggsave(file = paste0(path_to_figure_outputs, "/three_b.eps"), height = 3, width = 7, device = "eps", units = "in", dpi = 300)

processed_df_1_wt <- processed_df_1 %>% filter (perturb.kind == "gof")
## This plot looks at fraction of TEI vs kind of perturbation and phase.
three_c <- ggplot(data = processed_df_1_wt) +
  geom_bar(mapping = aes(y = fr_TEI, x = molecule.kinds, fill = perturb.phase), stat="identity", width = 0.8, position = position_dodge(width = 0.8)) +
  #facet_wrap(~ molecule.kinds, nrow = 1) +
  my_theme +
  #theme(axis.text.x = element_text()) +
  theme(legend.position = c(0.1, 0.7)) +
  xlab("number of entities at gen. #1") +
  ylab("fr. showing heritable epigenetic change") +
  scale_x_continuous(expand=c(0,0), breaks = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)) +
  scale_y_continuous(limits=c(0, 0.3), breaks = c(0, 0.15, 0.3), expand = c(0,0)) +
  labs(fill = "phase") +
  scale_fill_manual(values = cbp)
print(three_c)
ggsave(file = paste0(path_to_figure_outputs, "/three_c.eps"), height = 3, width = 7, device = "eps", units = "in", dpi = 300)

### Setting up to process even more. Ignoring the perturbation kind and perturbation phases.
processed_df_2 <- processed_df_1 %>% group_by(molecule.kinds) %>%
                   summarise (n = n())
ESP_per_m <- processed_df_2$n[1]

processed_df_3 <- processed_df_1 %>% group_by(molecule.kinds) %>%
  summarise(all_fraction_robust = n()/ESP_per_m,
            all_median_turtles = median(median_turtles), all_max_turtles = max(max_turtles),
            all_median_positive = median(median_positive), all_max_positive = max(max_positive),
            all_median_negative = median(median_negative), all_max_negative = max(max_negative),
            median_95_positive = 1.96 * sd(median_positive)/(sqrt(n())), median_95_negative = 1.96 * sd(median_negative)/(sqrt(n())),
            max_95_positive = 1.96 * sd(max_positive)/(sqrt(n())), max_95_negative = 1.96 * sd(max_negative)/(sqrt(n())))

most_max_positive <- max(processed_df_3$all_max_positive) ## needed to scale the Y-axis appropriately.

### create labels to add as annotations grobTree collects grobs (graphical objects) together.
annotation_to_add_1 <- grobTree(textGrob("negative", x=0.6,  y=0.4, hjust=0,
                                       gp=gpar(col="blue", fontsize=10, fontface="plain")),
                              textGrob("positive", x=0.6,  y=0.9, hjust=0,
                                       gp=gpar(col="red", fontsize=10, fontface="plain")))


# ## This plot looks at molecule kinds vs positive and negative regulation.
four <- ggplot(data = processed_df_3, aes(x = molecule.kinds, y = all_max_positive)) +
   my_theme +
   geom_line(aes(y = all_max_negative), color = "blue") +
   geom_line(aes(y = all_max_positive), color = "red") +
   geom_linerange(aes(ymin = all_max_positive - max_95_positive, ymax = all_max_positive + max_95_positive), position = position_dodge(.9)) +
   geom_linerange(aes(ymin = all_max_negative - max_95_negative, ymax = all_max_negative + max_95_negative), position = position_dodge(.9)) +
   xlab("number of molecules") + ylab("maximal regulation") + annotation_custom(annotation_to_add_1) +
   scale_x_continuous(limits=c(0, 17), breaks = c(0, 4, 8, 12, 16), expand = c(0,0)) +
   scale_y_continuous(limits=c(0, most_max_positive + 5), breaks = c(0, 10, 20, 30, 40, most_max_positive), expand = c(0,0))
print(four)
ggsave(file = paste0(path_to_figure_outputs, "/four.eps"), height = 2, width = 2, device = "eps", units = "in", dpi = 300)

most_median_positive <- max(processed_df_3$all_median_positive)

### create labels to add as annotations grobTree collects grobs (graphical objects) together.
annotation_to_add_2 <- grobTree(textGrob("negative", x=0.6,  y=0.1, hjust=0,
                                         gp=gpar(col="blue", fontsize=10, fontface="plain")),
                                textGrob("positive", x=0.6,  y=0.95, hjust=0,
                                         gp=gpar(col="red", fontsize=10, fontface="plain")))

## This plot looks at molecule kinds vs positive and negative regulation.
five <- ggplot(data = processed_df_3, aes(x = molecule.kinds, y = all_median_positive)) +
  my_theme +
  geom_line(aes(y = all_median_negative), color = "blue") +
  geom_line(aes(y = all_median_positive), color = "red") +
  geom_linerange(aes(ymin = all_median_positive - median_95_positive, ymax = all_median_positive + median_95_positive), position = position_dodge(.9)) +
  geom_linerange(aes(ymin = all_median_negative - median_95_negative, ymax = all_median_negative + median_95_negative), position = position_dodge(.9)) +
  xlab("number of molecules") + ylab("median regulation") + annotation_custom(annotation_to_add_2) +
  scale_x_continuous(limits=c(0, 17), breaks = c(0, 4, 8, 12, 16), expand = c(0,0)) +
  scale_y_continuous(limits=c(-0.1, most_median_positive + 1), breaks = c(0, 2, 4, most_median_positive), expand = c(0, 0))
print(five)
ggsave(file = paste0(path_to_figure_outputs, "/five.eps"), height = 2, width = 2, device = "eps", units = "in", dpi = 300)

### process data frame again for the next plots.
processed_df_4_wt <- processed_df_1 %>% filter(perturb.kind == "none") %>% group_by(molecule.kinds) %>%
                  summarise(phase_median_turtles = median(median_turtles), phase_median_turtles_95CI = 1.96 * sd(median_positive)/(sqrt(n())))
most_median_turtles <- max(processed_df_4_wt$phase_median_turtles)

processed_df_5_wt <- processed_df_1 %>% filter(perturb.kind == "none") %>% group_by(molecule.kinds) %>%
  summarise(phase_max_turtles = max(max_turtles), phase_max_turtles_95CI = 1.96 * sd(max_turtles)/(sqrt(n())))
most_max_turtles <- max(processed_df_5_wt$phase_max_turtles)

### create labels to add as annotations grobTree collects grobs (graphical objects) together.
annotation_to_add_3 <- grobTree(textGrob("median", x=0.7,  y=0.2, hjust=0,
         gp=gpar(col="blue", fontsize=10, fontface="plain")),
         textGrob("maximum", x=0.3,  y=0.9, hjust=0,
         gp=gpar(col="black", fontsize=10, fontface="plain")))

### to change labels of the facets and to reorder the factor perturb.kind by changing levels.
# Below code not needed if only wild-type is being used.
# labels <- c(none = "wild type", gof = "gain of function", lof = "loss of function")
# processed_df_4$perturb.kind <- factor(processed_df_4$perturb.kind, levels = c("none", "lof", "gof"))
# processed_df_5$perturb.kind <- factor(processed_df_5$perturb.kind, levels = c("none", "lof", "gof"))

## This plot looks at molecule kinds vs maximum suriving turtles at the end and negative regulation.
six <- ggplot() +
  geom_pointrange(data = processed_df_5_wt, mapping = aes(molecule.kinds, phase_max_turtles, ymin = phase_max_turtles - phase_max_turtles_95CI, ymax = phase_max_turtles + phase_max_turtles_95CI), position = position_dodge(.9), color = "black") +
  geom_pointrange(data = processed_df_4_wt, mapping = aes(molecule.kinds, phase_median_turtles, ymin = phase_median_turtles - phase_median_turtles_95CI, ymax = phase_median_turtles + phase_median_turtles_95CI), position = position_dodge(.9), color = "blue") +
  ## facet_wrap(~ perturb.kind, labeller = labeller(perturb.kind = labels), nrow = 1) +
  my_theme +
  xlab("number of entities at gen. 1") + ylab("entities after 250 generations") + annotation_custom(annotation_to_add_3) +
  scale_x_continuous(limits=c(0, 17), breaks = c(0, 4, 8, 12, 16), expand = c(0,0)) +
  scale_y_continuous(limits=c(0, most_max_turtles + 1), breaks = c(0, 2, 4, 6, 8, most_max_turtles), expand = c(0, 0))
print(six)
ggsave(file = paste0(path_to_figure_outputs, "/six.eps"), height = 2.7, width = 2, device = "eps", units = "in", dpi = 300)

## Setting up to process to identify proportion of systems with negative regulation for each molecule kinds.
processed_df_6 <- robust_only_df %>% mutate(binary_neg = if_else(negative.links > 0, 1, 0)) %>%
                                           group_by(molecule.kinds) %>%
                                           summarise(fr_neg = sum(binary_neg)/n()) %>%
                                           mutate(fr_neg = round(fr_neg,2))

most_neg_fr <- max(processed_df_6$fr_neg)

## This plot looks at proportion of robust systems that have negative regulation in the end.
seven <- ggplot(data = processed_df_6, aes(x = molecule.kinds, y = fr_neg)) +
  geom_point(color = "black") +
  my_theme +
  xlab("kinds of molecules") + ylab("fraction systems with negative regulation") +
  scale_x_continuous(limits=c(0, 16.5), breaks = c(0, 4, 8, 12, 16), expand = c(0,0)) +
  scale_y_continuous(limits=c(-0.01, most_neg_fr + 0.05), breaks = c(0, 0.1, 0.2, 0.3, most_neg_fr), expand = c(0, 0))
print(seven)
ggsave(file = paste0(path_to_figure_outputs, "/seven.eps"), height = 3.5, width = 2, device = "eps", units = "in", dpi = 300)

## Generate table for plot eight that would be TEI events vs. survving ESP systems.
max_Frequency <- max(TEI_frequency$Frequency)

## This plot looks at frequency of TEI events vs. surviving ESP systems.
eight <- ggplot(data = TEI_frequency, aes(x = TEI, y = Frequency)) +
  geom_point(color = "black") +
  my_theme +
  xlab("number of TEI events") + ylab("number of surviving ESP systems") +
  scale_x_continuous(limits=c(-0.2, 4.5), breaks = c(0, 2, 4), expand = c(0,0)) +
  scale_y_continuous(limits=c(1, 100000), trans = "log10", breaks = c(1, 10, 100, 1000, 10000, 100000, max_Frequency), expand = c(0, 0))
print(eight)
ggsave(file = paste0(path_to_figure_outputs, "/eight.eps"), height = 3.5, width = 3.5, device = "eps", units = "in", dpi = 300)
