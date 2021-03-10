# Figure 4 creation

load('data/rel_abund_vs_c_use_by_treatment_and_ecosystem_draft_4.RData')
load('data/lmm_with_ecosystem_draft_4.RData')

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

# need to include tick marks along every 10-fold change, use the outer function
tick_marks <- outer(1:10, 1/10^(1:14), '*') # will only plot the ticks that are necessary
label_pos <- matrix('', ncol=ncol(tick_marks), nrow=nrow(tick_marks))
label_pos[1,] <- as.character(tick_marks[1,])
# only label every other for the x-axis
label_pos_x <- label_pos[1,]
blanks <- label_pos
blanks[] <- ''
blanks[1,seq(from=2, to=14, by=3)] <- label_pos_x[seq(from=2, to=14, by=3)]
label_pos_x <- blanks
label_pos[,seq(from=1, to=14, by=2)] <- ''
# NOTE, we are TRICKING the x-axis to display apparent 10-fold changes
log_tick_marks <- log10(tick_marks)

# plot log10-transformed values so you can superimpose model fit
abund_c_use_full <- within(abund_c_use_full, {
  log_c <- log10(c_flux_rel)
  log_abund <- log10(rel_abund)
})

# make data frame of linear model coefficients for 18O
gc <- function(x) summary(x)$coefficients 
line_data <- data.frame(slope=1 + gc(lmm_optim_3)[2,1],
                        intercept=gc(lmm_optim_3)[1,1])
# add C and CN
line_data <- list(line_data, line_data, line_data)
names(line_data) <- c('18O', 'C', 'CN')
line_data <- Map(function(x,y) {x$treatment <- y; x}, 
                 line_data,
                 c('18O', 'C', 'CN'))
# modify C slopes and intercepts
line_data$C$slope[1] <- 1 + gc(lmm_optim)[5,1]  # not significantly different from 1:1 line
line_data$C$intercept[1] <- line_data$`18O`$intercept + gc(lmm_optim)[3,1]  # not significantly different from 18O
# modify CN slopes and intercepts
line_data$CN$slope[1] <- 1 + gc(lmm_optim)[6,1]
line_data$CN$intercept[1] <- line_data$`18O`$intercept + gc(lmm_optim)[4,1]
# combine
line_data <- do.call(rbind, line_data)
# xmin and max
min_max <- function(x, subset) {
  data.frame(xmin=min(x[x$treatment==subset, 'log_abund'], na.rm=T), 
             xmax=max(x[x$treatment==subset, 'log_abund'], na.rm=T),
             treatment=subset)
}
#
line_map_data <- rbind(min_max(abund_c_use_full, '18O'),
                       min_max(abund_c_use_full, 'C'),
                       min_max(abund_c_use_full, 'CN'))
line_map_data <- merge(line_map_data, line_data, all.x=T)
#
line_map_data <- within(line_map_data, {
  ymin <- xmin*slope + intercept
  ymax <- xmax*slope + intercept
})

# match line map data treatment levels with plot data
line_map_data$treatment <- factor(line_map_data$treatment, labels=expression(''^18*O, 'C', 'C+N'))
abund_c_use_full$treatment <- factor(abund_c_use_full$treatment, labels=expression(''^18*O, 'C', 'C+N'))

# main scatterplot
ggplot(abund_c_use_full, aes(log_abund, log_c, color=Phylum)) +  
  geom_abline(slope=1, color=gray(.5)) +
  geom_point(size=.9) +
  geom_segment(aes(x=xmin, xend=xmax,
                   y=ymin, yend=ymax), 
               data=line_map_data, color=hsv(.6,1,0), size=.8) +
  scale_x_continuous(breaks=c(log_tick_marks, 0),
                     labels=c(label_pos_x, 1),
                     limits=c(min(abund_c_use_full$log_c), 0)) +
  scale_y_continuous(breaks=c(log_tick_marks, 0),
                     labels=c(label_pos, 1),
                     limits=c(min(abund_c_use_full$log_c), 0)) +
  xlab('16S rRNA gene relative abundance') +
  ylab('Relativized carbon use') +
  scale_color_manual(values=sankey_colors, 
                     guide=guide_legend(ncol=1,
                                        override.aes=list(size=3))) +
  facet_wrap(~ treatment, ncol=2, labeller='label_parsed') +
  geom_text(aes(x=-Inf, 
                y=Inf, 
                label=paste0('slope == ', abs(round(slope, 3))),
                hjust=-.25,
                vjust=2), 
            data=line_map_data, inherit.aes=F, parse=T, size=3) +
  geom_text(aes(x=-Inf, y=Inf, label=signif, hjust=-14, vjust=1.75), 
            data=data.frame(signif=c('', '', '*'),
                            treatment=factor(1:3, labels=expression(''^18*O, 'C', 'C+N'))), 
            inherit.aes=F, parse=F, size=4) +
  labs(tag='b') +
  theme(legend.position=c(1, .125),
        legend.justification=c(1, 0),
        legend.title=element_blank(),
        plot.margin=unit(c(-5.5, 5.5, 5.5, 5.5), units='pt'))


# average across treatments and ecosystems
abund_c_use <- aggregate(relative_contribution ~ ., 
                         abund_c_use[,!names(abund_c_use) %in% c('ecosystem', 'treatment')], 
                         mean)

# what percent of these genera, on average, use less C than their relative abundance suggests?
with(abund_c_use_full,
     sum(((c_flux_rel / rel_abund) < 1) / nrow(abund_c_use_full)))
# 75.7 %

# re-order genus levels
genus_order <- abund_c_use[abund_c_use$metric=='16S rRNA gene abundance',]
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
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, face='italic', color=label_colors),
        legend.position=c(1, 1), 
        legend.justification=c(1, 1),
        legend.title=element_blank(),
        plot.margin=unit(c(5.5, 5.5, -5.5, 5.5), units='pt'))