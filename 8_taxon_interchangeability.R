# interchangeability of taxa based off their carbon use

load('data/absolute_growth.RData')
# load biomass production, and respiration data frames
# created from per_taxon_carbon_growth_respiration.R
load('data/growth_ug_c_contributions_per_asv_draft_4.RData')        # biomass growth
load('data/resp_growth_ug_c_contributions_per_asv_draft_4.RData')   # growth-respiration

# keep only ones summarized at the family level
rm(list=ls(pattern='2|4'))


# change "Abundance" column names for merging
names(growth_data5)[grep('Abund', names(growth_data5))] <- 'production'
names(resp_growth_data5)[grep('Abund', names(resp_growth_data5))] <- 'resp_growth'

# merge just growth and growth-related respiration
c_flux <- Reduce(function(x,y) merge(x,y, all.x=T), list(resp_growth_data5, growth_data5))

# make total C use from growth and resp
c_flux$total_c <- apply(c_flux[,c('production', 'resp_growth')], 1, sum, na.rm=T)
c_flux <- c_flux[c_flux$total_c > 0,]

c_flux <- merge(c_flux,
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


# 1. what is the evenness of different taxa at each taxon rank?................

# create tax ranks per sample
c_rank <- lapply(c_flux, function(x) {
  x <- x[order(x$rel_c, decreasing=T),]
  x$tax_rank <- 1:nrow(x)
  x$prop_tax_rank <- (1:nrow(x)) / nrow(x)
  x
})

# group by rank and treatment, focus only on ranks where all samples are present
rank_even <- do.call(rbind, c_rank)
rank_even <- split(rank_even, interaction(rank_even$treatment, rank_even$tax_rank))
rank_even <- rank_even[sapply(rank_even, nrow) > 1]

# sum number of times a taxon occurs in a fraction and re-relativize
rank_even <- lapply(rank_even, function(x) {
  y <- droplevels(x)
  y <- aggregate(rel_c ~ taxonID + tax_rank + treatment, y, length)
  y$rel_c <- y$rel_c / sum(y$rel_c, na.rm=T)
  y
})

# calculate: what is the extent that a certain rank is evenly made up of several different taxa?
rank_even <- lapply(rank_even, function(x) {
  shannon_c <- -sum(x$rel_c * log(x$rel_c), na.rm=T)
  rich <- nrow(x)
  pielou_c <- shannon_c / log(rich)
  simpson_c <- sum(x$rel_c^2, na.rm=T)
  #
  x <- x[1, c('tax_rank', 'treatment')]
  x <- cbind(x, pielou_c, simpson_c, rich)
})
#
rank_even <- do.call(rbind, rank_even)


# 2. what is the distribution of rank values for each taxon?...................

# let's use proportional rank

rank_vals <- do.call(rbind, c_rank)

rank_med <- aggregate(cbind(prop_tax_rank, rel_abund) ~ taxonID + treatment, rank_vals, median)
rank_med_abs <- aggregate(tax_rank ~ taxonID + treatment, rank_vals, median)
rank_sd <- aggregate(prop_tax_rank ~ taxonID + treatment, rank_vals, sd)
rank_se <- aggregate(tax_rank ~ taxonID + treatment, rank_vals, function(x) sd(x) / length(x))
rank_num <- aggregate(prop_tax_rank ~ taxonID + treatment, rank_vals, function(x) length(unique(x)))
rank_dom <- aggregate(prop_tax_rank ~ taxonID + treatment, rank_vals, function(x) {
  y <- table(x) / sum(x)
  sum(y^2, na.rm=T)
})

# rename "tax_rank" columns
rank_vals <- list(rank_med, rank_med_abs, rank_sd, rank_se, rank_num, rank_dom)
rank_vals <- Map(function(x, y) {names(x)[grep('tax_rank$', names(x))] <- paste0('tax_rank_', y); x},
                 rank_vals,
                 c('median', 'median_absolute', 'sd', 'se', 'count', 'dom'))
# merge
rank_vals <- Reduce(function(x,y) merge(x, y, by=c('taxonID', 'treatment')), rank_vals)


# save data
save(rank_even, rank_vals, file='../Data/taxon_importance_C_use_rankings.RData')