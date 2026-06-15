library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)
library(stringr)
library(tidyr)
# apply cowplot theme
theme_set(theme_cowplot())

# Prepare datasets -------------------------------------------------------------

# Understanding methylation values
#meth_bins = data.frame(bin = 0:255, min_prob = (0:255)/256, max_prob = ((0:255) + 1)/256)

metadata = read.csv('/Users/avvarua/Documents/projects/diaslab/fmr1-instability/sample-metadata.tsv', sep = '\t')
all.data = read.csv('/Users/avvarua/Documents/projects/diaslab/allSamples-readsdata.txt', sep = '\t')
all.data = merge(all.data, metadata, 'bam_id')

all.data = subset(all.data, base_qual > 10)  #filter out reads with basequal <10

#converts allele_length, methylated_bases, unmethylated_bases, median_meth to numeric
all.data$allele_length = as.numeric(all.data$allele_length)
all.data$methylated_bases = as.numeric(all.data$methylated_bases)
all.data$umethylated_bases = as.numeric(all.data$umethylated_bases)
all.data$median_meth = as.numeric(all.data$median_meth)

all.data$medianmethlevel = round(all.data$median_meth, 1)
all.data$ismeth = all.data$median_meth >= 0.5
all.data$prop_meth = all.data$methylated_bases/(all.data$methylated_bases + all.data$umethylated_bases)
all.data$prop_meth[is.nan(all.data$prop_meth)] = NA
all.data$prop_meth_level = as.factor(round(all.data$prop_meth, 1))
all.data$bin <- cut(all.data$allele_length, breaks=c(0,54,200,100000), labels=c('Normal', 'Premutation', 'Pathogenic'), include.lowest=TRUE)
all.data$meth_bin <- cut(all.data$prop_meth, breaks=c(0,0.3,0.6,1), labels=c('Weak', 'Partial', 'Full'), include.lowest = TRUE, na.rm = TRUE)
all.data <- all.data %>% mutate(condition = case_when( startsWith(Individual, "FXPM") ~ "FXPM",
                                                       startsWith(Individual, "FXS")  ~ "FXS",
                                                       startsWith(Individual, "CON")  ~ "CONTROL", TRUE ~ NA_character_))

# Prepare datasets -------------------------------------------------------------

# Prepare gender datasets ------------------------------------------------------
female.data = all.data[all.data$Sex=='female',]
female.FXPM.data = female.data[str_starts(female.data$Individual, "FXPM"),]
female.FXS.data  = female.data[str_starts(female.data$Individual, "FXS"),]
female.CON.data  = female.data[str_starts(female.data$Individual, "CON"),]

male.data = all.data[all.data$Sex=='male',]
male.FXPM.data = male.data[str_starts(male.data$Individual, "FXPM"),]
male.FXS.data  = male.data[str_starts(male.data$Individual, "FXS"),]
male.CON.data  = male.data[str_starts(male.data$Individual, "CON"),]
# Prepare gender datasets ------------------------------------------------------

# Prop methylation histograms --------------------------------------------------

nfacets = distinct(male.FXS.data[, c("Individual", "Region")]) %>% nrow()
nrows = nfacets %/% 3
if (nfacets %% 3 > 0) { nrows = nrows + 1 }
ggplot(male.FXS.data, aes(x = allele_length, fill = prop_meth_level)) + 
  geom_histogram() + 
  facet_wrap(Individual~Region, ncol=3, scales = 'free') + 
  scale_fill_manual(values = c('#50b7bc', '#71c4c8', '#94d3d6', '#b4e0e2', '#daeff1', 'white', '#e6aba0', '#e6aba0', '#dd897a', '#d1644e', '#c74024'))+
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated')  
ggsave('../plots/Male_FXS_methylation_prop.pdf', width = 14, height = nrows * 3.6)

# Prop methylation histograms --------------------------------------------------

# Inidividual prop methylation -------------------------------------------------

individuals <- unique(all.data$Individual)
for (individual in individuals){
  sub.data = subset(all.data, Individual == individual)
  # Get sex from metadata
  sex <- metadata %>% filter(Individual == individual) %>% pull(Sex)
  nrows = length(unique(sub.data[, c("Region")]))
  if (nfacets %% 2 > c0) { nrows = nrows + 1 }
  ggplot(sub.data, aes(x = bin, fill = prop_meth_level)) + 
    geom_bar() + 
    facet_grid(Region~Individual, scales = 'free') + 
    scale_fill_manual(values = c('#50b7bc', '#71c4c8', '#94d3d6', '#b4e0e2', '#daeff1', 'white', '#e6aba0', '#e6aba0', '#dd897a', '#d1644e', '#c74024'),
                      breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
    labs(x = 'Allele length', y = 'Reads', fill = 'Prop. Cs in repeat methylated', title = paste0("Individual: ", individual, "    Sex: ", sex)) 
  ggsave(paste0('../plots/', individual,'_readsPropMeth.jpg'), width=7, height = nrows * 2.6)
}

# Inidividual prop methylation -------------------------------------------------

# Inidividual methylation allele bin -------------------------------------------

for (individual in individuals) {
  sub.data <- subset(all.data, Individual == individual & !is.na(prop_meth_level))
  sub.data$prop_meth_level <- as.numeric(as.character(sub.data$prop_meth_level))
  sex <- metadata %>% filter(Individual == individual) %>% pull(Sex)
  ggplot(sub.data, aes(x = prop_meth_level, fill = bin)) +
    geom_bar() + 
    scale_x_continuous(breaks=seq(0, 1, by=0.1)) +
    scale_fill_manual(values = c("Normal" = "#00c299", "Premutation" = "#faba42", "Pathogenic" = "#fc584c")) +
    labs(x = 'Proportion of Methylated Cs', y = 'Reads', fill = 'Allele length', title = paste0("Individual: ", individual, "    Sex: ", sex)) 
  ggsave(paste0('../plots/', individual,'_methvsallelebin.jpg'), width=7)
}

# Inidividual methylation allele bin -------------------------------------------

# Inidividual methylation allele bin scatter -----------------------------------

for (individual in individuals) {
  sub.data <- subset(all.data, Individual == individual & !is.na(prop_meth_level))
  sub.data$prop_meth_level <- as.numeric(as.character(sub.data$prop_meth_level))
  sex <- metadata %>% filter(Individual == individual) %>% pull(Sex)
  ggplot(sub.data, aes(x = allele_length, y = prop_meth_level, color = bin)) +
    geom_point(size = 1.5, alpha = 0.5) + 
    scale_color_manual(values = c("Normal" = "#00c299", "Premutation" = "#faba42", "Pathogenic" = "#fc584c")) +
    labs(x = 'Allele length', y = 'Proportion of Methylated Cs', fill = 'Allele length', title = paste0("Individual: ", individual, "    Sex: ", sex)) 
  ggsave(paste0('../plots/', individual,'_methvsallelescatter.jpg'), width=7)
}

# Inidividual methylation allele bin scatter -----------------------------------

for (individual in individuals) {
  
  # Subset sample data
  df <- subset(all.data, Individual == individual & !is.na(prop_meth_level))
  
  # Convert prop_meth_level to numeric (in case it's factor/character)
  df$prop_meth_level <- as.numeric(as.character(df$prop_meth_level))
  
  # Get sex from metadata
  sex <- metadata %>%
    filter(Individual == individual) %>%   # or column name that matches "Individual"
    pull(Sex)                          # or "sex" depending on your metadata
  
  # Create plot with title
  p <- ggplot(df, aes(x = prop_meth_level, fill = bin)) + 
    geom_bar() + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_fill_manual(values = c("Normal" = "#00c299", 
                                 "Premutation" = "#faba42", 
                                 "Pathogenic" = "#fc584c")) +
    labs(
      x = 'Proportion of Methylated Cs', 
      y = 'Reads', 
      fill = 'Allele length',
      title = paste0("Individual: ", individual, " Sex: ", sex)
    )
  
  # Save plot
  ggsave(
    paste0('../plots/', individual,'_methvsallelebin.jpg'), 
    plot = p,
    width = 7
  )
}


unmeth.data <- all.data[as.numeric(all.data$prop_meth)<0.5,]

for (individual in individuals){
  # if (individual == "CON5274") {
  ggplot(subset(all.data, Individual == individual), aes(x = Region, fill = bin)) + 
    # geom_histogram(breaks=c(0,55,200)) +
    geom_bar() + 
    # facet_grid(Region~Individual, scales = 'free') + 
    scale_fill_manual(values = c("0-54" = "seagreen3", "55-200" = "orangered2")) +
    # scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
    labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Allele length bin') 
  ggsave(paste0('../plots/', individual,'_test.jpg'), width=7)
  # }
}

for (individual in unique(all.data$Individual)){
  # if (individual == "FXPM1008-20-RF") {
  #   print(individual)
  ggplot(subset(all.data, Individual == individual), aes(x = allele_length, y = prop_meth_level, color = prop_meth_level)) + 
    # geom_histogram(breaks=c(0,55,200)) +
    geom_point(size=2, alpha=0.7) +
    geom_jitter(width=0.2, height = 0.2, size=2, alpha=0.7) +
    facet_grid(Region~Individual, scales = 'free') +
    scale_color_manual(values = c('#50b7bc', '#71c4c8', '#94d3d6', '#b4e0e2', '#daeff1', 'white', '#e6aba0', '#e6aba0', '#dd897a', '#d1644e', '#c74024'),
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
    # xlim(0,200) +
    # scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
    labs(x = 'Allele length (motifs)', y = 'Reads') 
  ggsave(paste0('FXS-pureTRGT-plots/', individual,'_FMR1_PureTarget_length-methylation.jpg'), width=7)
  # }
}

for (individual in unique(all.data$Individual)){
  print(individual)
  print(subset(all.data, Individual == individual))
  ggplot(subset(all.data, Individual == individual), aes(x = allele_length, fill = prop_meth_level)) + 
    geom_histogram(bins=c(0,55,200)) + 
    facet_grid(Region~Individual, scales = 'free') + 
    scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
    labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated') 
  ggsave(paste0(individual,'_FMR1_PureTarget_methylation_prop.jpg'))
}

ggplot(subset(all.data, allele_length > 400), aes(x = allele_length, fill = prop_meth_level)) +
  geom_histogram() +
  facet_wrap(~Sex) +
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated')

ggplot(all.data, aes(x = allele_length, fill = prop_meth_level)) +
  geom_histogram() +
  facet_wrap(~Sex) +
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated')

ggplot(subset(all.data, allele_length > 400), aes(x = allele_length, fill = Sex)) +
  geom_histogram() +
  labs(x = 'Allele length (motifs)', y = 'Reads')

ggplot(all.data, aes(x = median_meth, y = prop_meth)) + geom_point()

ggplot(all.data, aes(x = base_qual)) + geom_histogram()
ggplot(all.data, aes(x = base_qual, y = allele_length)) + geom_point() + ylim(0, 100)

ggplot(subset(all.data, Individual == 'FXPM4555'), aes(x = allele_length, fill = prop_meth_level)) + 
  geom_histogram() + 
  facet_grid(Region~Individual, scales = 'free') + 
  scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
  labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Prop. Cs in repeat methylated') 


ind.data <- subset(all.data, Individual=='FXPM5006')
summary.data <- ind.data %>%
                group_by(Individual, Region, bin) %>%
                summarise(
                  total_reads = n(),
                  meth_above60 = sum(prop_meth >= 0.6, na.rm = TRUE),
                  .groups = "drop"
                )
summary.data$percent_above60 <- (summary.data$meth_above60 / summary.data$total_reads) * 100
ggplot(summary.data, aes(x = Region, y = percent_above60, fill = bin)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Normal" = "#00c299",  "Premutation" = "#faba42",  "Pathogenic" = "#fc584c")) +
  labs (x = "Region", y = "Percent reads (%)", fill = "Allele", title = "Percentage of reads with >=60% of reads methylated FXPM5006") +
  theme_minimal()
ggsave('/Users/avvarua/Documents/projects/diaslab/plots/FXPM5006-alleles-meth60.jpg', width = 7, height = 4)

# Premutation reads vs Methylation bin -----------------------------------------

individual <- "FXS5319"
ind.data <- subset(all.data, Individual==individual)
summary.data <- ind.data %>%
  filter(!is.na(meth_bin)) %>%
  group_by(Individual, Region, bin, meth_bin) %>%
  summarise(
    total_reads = n(),
    .groups = "drop"
  ) %>%
  complete(Individual, Region, meth_bin, fill = list(bin = "Pathogenic", total_reads = 0))
ggplot(summary.data, aes(x = bin, y = total_reads, fill = meth_bin)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Region, nrow=2) +
  scale_fill_manual(values = c("Weak" = "blue", 
                               "Partial" = "lightblue", 
                               "Full" = "#fc584c")) +
  labs (x = "Allele bin", y = "Number of reads", fill = "Methylation", title = individual) +
  theme_minimal()
ggsave(paste0('/Users/avvarua/Documents/projects/diaslab/plots/',individual,'-reads-methbin.jpg'), width = 7, height = 4)

# Premutation reads vs Methylation bin -----------------------------------------

# Premutation reads vs Methylation bin stack bar all samples -------------------

summary.data <- all.data %>% filter(!is.na(meth_bin)) %>% filter(bin == "Premutation") %>%
  group_by(Individual, Sex, Region, bin, meth_bin) %>%
  summarise( total_reads = n(), .groups = "drop" ) %>%
  group_by(Individual, Sex, Region, bin) %>%
  mutate(percentage = total_reads / sum(total_reads) * 100) %>%
  ungroup()
  # complete(Individual, Region, meth_bin, fill = list(bin = "Pathogenic", total_reads = 0))

summary.data$label <- paste0(summary.data$Individual," (",summary.data$Region,")")
write.table(summary.data, file = "/Users/avvarua/Documents/projects/diaslab/plots/MethylationBinPercent-Premutation-Sex.tsv", sep="\t", quote=FALSE, row.names = FALSE)

summary.data <- summary.data %>% filter(Region == "BA10" & (Individual == "FXPM1008-20-RF" | Individual == "FXPM1015-09-DK"))

ggplot(summary.data, aes(x = label, y = percentage, fill = meth_bin)) +
  geom_col() +
  facet_wrap(~Sex, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("Weak" = "blue", "Partial" = "lightblue", "Full" = "#fc584c")) +
  labs (x = "Samples - Individual (Region)", y = "Percent reads (%)", fill = "Methylation", title = "Percentage of Methylated Premutation reads") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0('/Users/avvarua/Documents/projects/diaslab/plots/percent-methbin-permutation-subset.jpg'), width = 5, height = 4)


summary.data <- all.data %>% filter(!is.na(meth_bin)) %>% filter(bin == "Premutation") %>%
  group_by(Individual, Sex, Region, bin, meth_bin) %>%
  summarise( total_reads = n(), .groups = "drop" ) %>%
  group_by(Individual, Sex, Region, bin) %>%
  mutate(percentage = total_reads / sum(total_reads) * 100) %>%
  ungroup()
# complete(Individual, Region, meth_bin, fill = list(bin = "Pathogenic", total_reads = 0))

summary.data$label <- paste0(summary.data$Individual," (",summary.data$Region,")")
ggplot(summary.data, aes(x = label, y = percentage, fill = meth_bin)) +
  geom_col() +
  facet_wrap(~Sex, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("Weak" = "blue", "Partial" = "lightblue", "Full" = "#fc584c")) +
  labs (x = "Samples - Individual (Region)", y = "Percent reads (%)", fill = "Methylation", title = "Percentage of Methylated Premutation reads") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0('/Users/avvarua/Documents/projects/diaslab/plots/percent-methbin-permutation.jpg'), width = 6, height = 8)

# Premutation reads vs Methylation bin stack bar all samples -------------------

sub.data <- all.data %>% filter(all.data$Sex == "female" & all.data$condition != "FXS")
sub.data <- sub.data %>% mutate(label = paste0(Individual, "(", Region, ")"))
ggplot(sub.data, aes(x = prop_meth, y = label, color = bin, fill = bin)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2) +
  scale_color_manual(values = c("Normal" = "#00c299",  "Premutation" = "#faba42",  "Pathogenic" = "#fc584c")) +
  scale_fill_manual(values = c("Normal" = "#00c299",  "Premutation" = "#faba42",  "Pathogenic" = "#fc584c")) +
  theme_ridges() +
  theme( axis.title.y = element_blank() ) +
  labs( x = "Methylation", title = "Methylation Distribution - Males")

# for (individual in unique(all.data$Individual)){
#   ggplot(subset(all.data, Individual == individual), aes(x = allele_length, fill = as.factor(medianmethlevel))) + 
#     geom_histogram() + 
#     facet_grid(Region~Individual, scales = 'free') + 
#     scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
#     labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Median methylation') 
#   ggsave(paste0(individual,'_FMR1_PureTarget_methylation.pdf'))
# }
# ggplot(all.data, aes(x = allele_length, fill = as.factor(medianmethlevel))) + 
#   geom_histogram() + 
#   facet_wrap(Individual~Region, scales = 'free') + 
#   scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
#   labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Median methylation') 
# ggsave('FMR1_PureTarget_methylation.pdf', width = 20, height = 20)
# 
# ggplot(all.data, aes(x = allele_length, fill = as.factor(medianmethlevel))) + 
#   geom_histogram(aes(y = after_stat(density))) + 
#   facet_grid(Region~Individual, scales = 'free') + 
#   scale_fill_manual(values = c('blue4', 'blue3', 'blue2', 'blue1', 'royalblue1', 'snow2', 'tomato', 'red1', 'red2', 'red3', 'red4')) +
#   labs(x = 'Allele length (motifs)', y = 'Reads', fill = 'Median methylation')  
# ggsave('FMR1_PureTarget_methylation_byregion.pdf', width = 20, height = 20)
