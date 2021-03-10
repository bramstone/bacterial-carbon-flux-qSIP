# creation of figure 5b

library(ggplot2)

load('data/rel_abund_vs_native_c_use.RData')

# reverse factor for plotting
nat_c_use$Genus <- factor(nat_c_use$Genus, levels=rev(levels(nat_c_use$Genus)))

# plot main taxa
ggplot(nat_c_use, aes(Genus, relative_contribution, color=metric)) +
  geom_line(aes(group=metric, linetype=metric)) +
  geom_point(aes(shape=metric)) +
  ylab('Relative Contribution') +
  xlab('') +
  # scale_color_manual(values=hsv(c(0,.55), c(.5, 1), c(1,.35))) +
  scale_x_discrete(label=abbreviate) +
  facet_grid(treatment ~ ecosystem) +
  coord_flip() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, face='italic'),
        legend.position='bottom',
        legend.title=element_blank(),
        legend.direction='horizontal')


# plot all taxa

# need to include tick marks along every 10-fold change, use the outer function
plot_min_exponent <- 7
tick_marks <- outer(1:10, 1/10^(1:plot_min_exponent), '*')
label_pos <- tick_marks
label_pos[2:nrow(label_pos),] <- ''
# only label every other for the x-axis
label_pos_x <- label_pos
label_pos_x[,seq(from=3, to=plot_min_exponent, by=2)] <- ''

# plot relativized native C use vs. relativized C use
ggplot(nat_c_use_full[nat_c_use_full$treatment!='18O',], 
       aes(c_flux_rel, nat_c_flux_rel, color=Phylum)) +
  geom_abline(slope=1, color=gray(.5)) +
  geom_point(size=.9) +
  coord_trans(x='log10', y='log10') +
  scale_x_continuous(breaks=c(tick_marks, 1),
                     labels=c(label_pos_x, 1),
                     limits=c(1/10^plot_min_exponent, 1)) +
  scale_y_continuous(breaks=c(tick_marks, 1),
                     labels=c(label_pos, 1),
                     limits=c(1/10^plot_min_exponent, 1)) +
  xlab('Relativized carbon use') +
  ylab('Relativized native carbon use') +
  annotate('text', x=.5, y=.5, label='1:1', color=gray(.75), vjust=1, hjust=0) +
  facet_grid(. ~ treatment) +
  theme(legend.position=c(0, 1), 
        legend.justification=c(0, 1),
        legend.title=element_blank())
