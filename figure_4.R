# Figure 4 creation

load('data/rel_abund_vs_c_use_by_treatment_and_ecosystem_draft_4.RData')
load('data/lmm_with_ecosystem_draft_4.RData')
source('figure_creation_functions.R')

# color vector for axis labels
phyla <- unique(abund_c_use_full$Phylum)
phyla <- factor(c(as.character(phyla), 'Other'),
                levels=c(levels(phyla), 'Other'))
#
sankey_colors <- rev(rev(tax_color(phyla, 3920))[1:6])
sankey_colors[4] <- hsv(.1, 1, v=.6)

label_colors <- sankey_colors
label_colors <- data.frame(Phylum=levels(abund_c_use$Phylum),
                           color=label_colors)

abund_c_use <- merge(abund_c_use, label_colors, all.x=T)
abund_c_use <- abund_c_use[order(abund_c_use$Genus),]
label_colors <- abund_c_use[!duplicated(abund_c_use$Genus), 'color']
label_colors <- as.character(label_colors)
abund_c_use$color <- NULL

# If by treatment, abbreviate a single, long genus name using abbreviate() function
genus_names <- levels(abund_c_use$Genus)
genus_names[nchar(genus_names) > 20] <- abbreviate(genus_names[nchar(genus_names) > 20], minlength=20)
abund_c_use$Genus <- factor(abund_c_use$Genus, labels=genus_names)
#
abund_c_use$treatment <- factor(abund_c_use$treatment, labels=expression(''^18*'O', 'C', 'C+N'))

# reverse factor for plotting
abund_c_use$Genus <- factor(abund_c_use$Genus, levels=rev(levels(abund_c_use$Genus)))
label_colors <- rev(label_colors)

# average across ecosystems
abund_c_use <- aggregate(relative_contribution ~ ., 
                         abund_c_use[,!names(abund_c_use) %in% c('ecosystem')], 
                         mean)

# what percent of these genera, on average, use less C than their relative abundance suggests?
with(abund_c_use_full,
     sum(((c_flux_rel / rel_abund) < 1) / nrow(abund_c_use_full)))
# 75.7 %

# re-order genus levels
genus_order <- aggregate(relative_contribution ~ ., 
                         abund_c_use[,!names(abund_c_use) %in% c('treatment')], 
                         mean)
genus_order <- genus_order[genus_order$metric=='16S rRNA gene abundance',]
genus_order <- genus_order[,c('Genus', 'relative_contribution')]
genus_order <- genus_order[order(genus_order$relative_contribution, decreasing=T),]
abund_c_use$Genus <- factor(abund_c_use$Genus, levels=as.character(genus_order$Genus))

# color vector for axis labels
label_colors <- sankey_colors
label_colors <- data.frame(Phylum=levels(abund_c_use$Phylum),
                           color=label_colors)

abund_c_use <- merge(abund_c_use, label_colors, all.x=T)
abund_c_use <- abund_c_use[order(abund_c_use$Genus),]
label_colors <- abund_c_use[abund_c_use$metric=='Carbon use', 'color']
label_colors <- as.character(label_colors)
abund_c_use$color <- NULL

# If by treatment, abbreviate a single, long genus name using abbreviate() function
genus_names <- levels(abund_c_use$Genus)
genus_names[nchar(genus_names) > 20] <- abbreviate(genus_names[nchar(genus_names) > 20], minlength=20)
abund_c_use$Genus <- factor(abund_c_use$Genus, labels=genus_names)

ggplot(abund_c_use, aes(Genus, relative_contribution, color=metric)) +
  geom_line(aes(group=metric, linetype=metric)) +
  geom_point(aes(shape=metric)) +
  ylab('Relative contribution') +
  xlab('') +
  scale_color_manual(values=hsv(c(0,.55), c(.5, 1), c(1,.35))) +
  labs(tag='a') +
  facet_grid(treatment ~ .) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, face='italic', color=label_colors),
        legend.position=c(1, 1), 
        legend.justification=c(1, 1),
        legend.title=element_blank(),
        plot.margin=unit(c(5.5, 5.5, -5.5, 5.5), units='pt'))
