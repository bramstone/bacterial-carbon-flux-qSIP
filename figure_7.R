
# load data
load('data/absolute_growth.RData')
load('data/growth_ug_c_contributions_per_asv_draft_4.RData')        # biomass growth
load('data/resp_growth_ug_c_contributions_per_asv_draft_4.RData')   # growth-respiration
source('figure_creation_functions.R')

trt_cols <- hsv(c(0,.6,.1), c(0,1,1), c(.5,.8,.8))


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

# data for Bruce
c_flux <- merge(bar_data,
                all_data[,c('taxonID', 'ecosystem', 'treatment', 'replicate', 'rel_abund_wk1')],
                all.x=T,
                by.x=c('OTU', 'ecosystem', 'treatment', 'replicate'),
                by.y=c('taxonID', 'ecosystem', 'treatment', 'replicate'))
names(c_flux)[1] <- 'taxonID'

c_flux$production[is.na(c_flux$production)] <- 0

# relativize carbon use, re-relativize abundances
c_flux <- split(c_flux, c_flux$Sample)
c_flux <- lapply(c_flux, function(x) {
  x$rel_prod <- x$production / sum(x$production, na.rm=T)
  x$rel_resp <- x$resp_growth / sum(x$resp_growth, na.rm=T)
  x$rel_c <- x$total_c / sum(x$total_c, na.rm=T)
  x$rel_abund <- x$rel_abund_wk1 / sum(x$rel_abund_wk1, na.rm=T)
  x})


phyla_keep <- c('Actinobacteria', 'Proteobacteria', 'Acidobacteria',
                'Verrucomicrobia', 'Firmicutes', 'Chloroflexi',
                'Gemmatimonadetes', 'Bacteroidetes', 'Planctomycetes')

# split by sample AND phylum, keeping only top 9 phyla
c_rank <- do.call(rbind, c_flux)
c_rank <- c_rank[c_rank$Phylum %in% phyla_keep,]
c_rank <- droplevels(c_rank)
#
c_rank <- split(c_rank, interaction(c_rank$Sample, c_rank$Phylum))
c_rank <- c_rank[sapply(c_rank, function(x) nrow(x) > 0)]

# disregard taxonID, re-relativize C use, and simply rank them by relative C flux
c_rank <- lapply(c_rank, function(x) {
  x$rel_c <- x$rel_c / sum(x$rel_c, na.rm=T)
  x <- x[order(x$rel_c, decreasing=T),]
  x$tax_rank <- 1:nrow(x)
  x
})

# recombine, then average by treatment and taxon rank
c_rank <- do.call(rbind, c_rank)
c_rank_mean <- aggregate(rel_c ~ tax_rank + treatment + Phylum, c_rank, mean)

# make C flux cumulative
c_rank_mean <- split(c_rank_mean, interaction(c_rank_mean$treatment, c_rank_mean$Phylum))
#
c_rank_mean <- lapply(c_rank_mean, function(x) {
  x$rel_c <- x$rel_c / sum(x$rel_c, na.rm=T)
  x$cumul_c <- cumsum(x$rel_c)
  x
})
c_rank_mean <- do.call(rbind, c_rank_mean)

c_rank_ci <- aggregate(rel_c ~ tax_rank + treatment + Phylum, c_rank, function(x) se(x) * 1.96)
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
ggplot(c_rank, aes(log10(tax_rank), cumul_c, color=treatment)) +
  geom_hline(yintercept=.5, color=gray(.5)) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=.035) +
  geom_point() +
  ylim(0, 1) +
  xlab('Taxon rank') +
  ylab(expression('Cumulative carbon flux '%+-%' 95% CI')) +
  scale_x_continuous(breaks=c(log_tick_marks),
                     labels=c(label_pos),
                     limits=c(0, log10(max_tax))) +
  scale_color_manual(values=trt_cols, labels=c('Control', 'C', 'C+N')) +
  facet_wrap(~ Phylum, ncol=3) +
  theme(strip.text=element_text(size=rel(.9)),
        legend.position=c(1, 0), 
        legend.justification=c(1, 0),
        legend.text=element_text(hjust=0))
