# Figure 3 creation

# Figure 2 creation

load('data/growth_ug_c_contributions_per_asv_draft_4.RData')        # biomass growth
load('data/resp_growth_ug_c_contributions_per_asv_draft_4.RData')   # growth-respiration

# keep only ones summarized at the family level
rm(list=ls(pattern='2|4'))

# change "Abundance" column names for merging
names(growth_data5)[grep('Abund', names(growth_data5))] <- 'production'
names(resp_growth_data5)[grep('Abund', names(resp_growth_data5))] <- 'resp_growth'

# merge just growth and growth-related respiration
bar_data <- Reduce(function(x,y) merge(x,y, all.x=T), list(resp_growth_data5, growth_data5))

# make total C use from growth and resp
bar_data$total_c <- apply(bar_data[,c('production', 'resp_growth')], 1, sum, na.rm=T)
bar_data <- bar_data[bar_data$total_c > 0,]

c_flux <- merge(bar_data,
                all_data[,c('taxonID', 'ecosystem', 'treatment', 'replicate', 'rel_abund_wk1')],
                all.x=T,
                by.x=c('OTU', 'ecosystem', 'treatment', 'replicate'),
                by.y=c('taxonID', 'ecosystem', 'treatment', 'replicate'))
names(c_flux)[1] <- 'taxonID'

# use for functional diversity next
funct_div <- bar_data

# plot Pielou's evenness of C use against evenness of abundance
# Pielou's evenness is Shannon's evenness / log(richness)

c_flux$production[is.na(c_flux$production)] <- 0

# relativize carbon use, re-relativize abundances
c_flux <- split(c_flux, c_flux$Sample)
c_flux <- lapply(c_flux, function(x) {
  x$rel_prod <- x$production / sum(x$production, na.rm=T)
  x$rel_resp <- x$resp_growth / sum(x$resp_growth, na.rm=T)
  x$rel_c <- x$total_c / sum(x$total_c, na.rm=T)
  x$rel_abund <- x$rel_abund_wk1 / sum(x$rel_abund_wk1, na.rm=T)
  x})

# calculate Pielou's evenness of abundance and C use for each sample
pielou <- lapply(c_flux, function(x) {
  shannon_prod <- -sum(x$rel_prod * log(x$rel_prod), na.rm=T)
  shannon_resp <- -sum(x$rel_resp * log(x$rel_resp), na.rm=T)
  shannon_c <- -sum(x$rel_c * log(x$rel_c), na.rm=T)
  shannon_abund <- -sum(x$rel_abund * log(x$rel_abund), na.rm=T)
  #
  rich <- nrow(x)
  pielou_prod <- shannon_prod / log(rich)
  pielou_resp <- shannon_resp / log(rich)
  pielou_c <- shannon_c / log(rich)
  pielou_abund <- shannon_abund / log(rich)
  #
  x <- x[1,c('ecosystem', 'treatment', 'replicate', 'Sample')]
  x <- cbind(x, pielou_prod, pielou_resp, pielou_c, pielou_abund)
})
#
pielou <- do.call(rbind, pielou)

# transform to long-form
pielou <- reshape(pielou,
                  idvar=c('ecosystem', 'treatment', 'replicate', 'Sample'),
                  varying=list(c('pielou_prod', 'pielou_resp')),
                  times=c('Production', 'Respiration'),
                  timevar='flux_type',
                  v.names='pielou_c_use',
                  direction='long')

# plot abundance evenness with C use evenness
ggplot(pielou, aes(pielou_abund, pielou_c_use, color=treatment, shape=flux_type)) +
  geom_abline(slope=1, intercept=0, color=gray(.5)) +
  geom_point(size=2) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab('Pielou\'s evenness of relative abundances') +
  ylab('Pielou\'s evenness of relativized carbon use') +
  scale_color_manual(values=trt_cols, labels=expression(''^18*O, 'C', 'C+N')) +
  scale_shape_manual(values=c(19, 1)) +
  labs(tag='a') +
  theme(legend.position=c(0, 1), 
        legend.justification=c(0,1),
        legend.text=element_text(hjust=0))

# disregard taxonID, and simply rank them by relative C flux
c_rank <- lapply(c_flux, function(x) {
  x <- x[order(x$rel_c, decreasing=T),]
  x$tax_rank <- 1:nrow(x)
  x
})

# make C flux cumulative
c_rank <- lapply(c_rank, function(x) {
  x$cumul_c <- cumsum(x$rel_c)
  x
})

# recombine, then average by treatment and taxon rank
c_rank <- do.call(rbind, c_rank)
c_rank_mean <- aggregate(cumul_c ~ tax_rank + treatment, c_rank, mean)
c_rank_ci <- aggregate(cumul_c ~ tax_rank + treatment, c_rank, function(x) se(x) * 1.96)
names(c_rank_ci)[grep('_c', names(c_rank_ci))] <- 'cumul_c_ci'

c_rank <- merge(c_rank_mean, c_rank_ci)

# calculate conf. intervals
c_rank <- within(c_rank, {
  ci_hi <- cumul_c + cumul_c_ci
  ci_lo <- cumul_c - cumul_c_ci
})

# create tick marks to coincide with tax numbers
# need to include tick marks along every 10-fold change, use the outer function
tick_marks <- outer(1:10, 10^(0:3), '*') # will only plot the ticks that are necessary
label_pos <- matrix('', ncol=ncol(tick_marks), nrow=nrow(tick_marks))
label_pos[1,] <- as.character(tick_marks[1,])
# NOTE, we are TRICKING the x-axis to display apparent 10-fold changes
log_tick_marks <- log10(tick_marks)
max_tax <- max(c_rank$tax_rank)

# plot
gg_cumul_c_use <- ggplot(c_rank, aes(log10(tax_rank), cumul_c, color=treatment)) +
  geom_hline(yintercept=.5, color=gray(.5)) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=.035) +
  geom_point() +
  ylim(0, 1) +
  xlab('Taxon rank') +
  ylab(expression('Cumulative carbon use '%+-%' 95% CI')) +
  scale_x_continuous(breaks=c(log_tick_marks, log10(max_tax)),
                     labels=c(label_pos, max_tax),
                     limits=c(0, log10(max_tax))) +
  scale_color_manual(values=trt_cols, labels=expression(''^18*O, 'C', 'C+N')) +
  annotate('text', x=log10(100), y=.5, label='50% of total\ncarbon flux', color=gray(.5), hjust=.5) +
  labs(tag='b') +
  theme(legend.position='none')