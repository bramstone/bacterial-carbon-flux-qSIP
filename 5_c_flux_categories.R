# These data were prepared for a carbon flow Sankey diagram


load('data/respiration_18O_cue_draft_4.RData')
load('data/taxonomic_data.RData')
load('data/cue_18O_method.RData')
tax <- as.matrix(tax)
tax[is.na(tax)] <- 'unknown'
tax <- as.data.frame(tax)


# initial data prep---------------------------------------------------------------------

# use resp_data but combine with taxonomic info
sankey_data_gen <- merge(resp_data[,names(resp_data) %in% c('taxonID', 'ecosystem', 'treatment', 
                                                            'replicate', 'sampleID', 'activity',
                                                            'pop_c_ug', 'qsip_maint_resp_c_ug', 'net_status',
                                                            'qsip_growth_resp_c_ug', 'qsip_resp_c_ug')], 
                         tax, 
                         all.x=T)

# since we already calculated maintenance respiration based on pop levels, turn dormant growth to 0
# sankey_data_gen$pop_c_ug[sankey_data_gen$activity=='dormant'] <- 0
names(sankey_data_gen)[names(sankey_data_gen) %in% 'pop_c_ug'] <- 'growth_c_ug'

# # this data form will work for multivariate ordinations
# ord_data <- sankey_data_gen

# use c_summary to scale respiration and growth values
c_summary <- aggregate(cbind(growth_c_ug, qsip_maint_resp_c_ug, qsip_growth_resp_c_ug, qsip_resp_c_ug) ~ ecosystem + treatment + sampleID, 
                       sankey_data_gen, 
                       sum)

c_summary <- aggregate(cbind(growth_c_ug, qsip_maint_resp_c_ug, qsip_growth_resp_c_ug, qsip_resp_c_ug) ~ ecosystem + treatment,
                       c_summary,
                       mean)

# compare against values where ecosystem:treatment averaging has taken place
sankey_data_gen <- aggregate(cbind(growth_c_ug, qsip_maint_resp_c_ug, qsip_growth_resp_c_ug, qsip_resp_c_ug) ~ ., 
                             sankey_data_gen[!names(sankey_data_gen) %in% c('replicate', 'sampleID')], 
                             mean)
sankey_data_gen <- droplevels(sankey_data_gen)

c_lowered <- aggregate(cbind(growth_c_ug, qsip_maint_resp_c_ug, qsip_growth_resp_c_ug, qsip_resp_c_ug) ~ ecosystem + treatment, 
                       sankey_data_gen, 
                       sum)

# get multiplier to adjust rates
c_adjust <- cbind(c_summary[,1:2], c_summary[,sapply(c_summary, is.numeric)] / c_lowered[,sapply(c_lowered, is.numeric)])
c_adjust <- split(c_adjust, interaction(c_adjust$ecosystem, c_adjust$treatment))
sankey_data_gen <- split(sankey_data_gen, interaction(sankey_data_gen$ecosystem, sankey_data_gen$treatment))

# use multiplier to adjust growth and respiration
sankey_data_gen <- Map(function(x, y) {
  x$growth_c_ug <- x$growth_c_ug * y$growth_c_ug
  x$qsip_maint_resp_c_ug <- x$qsip_maint_resp_c_ug * y$qsip_maint_resp_c_ug
  x$qsip_growth_resp_c_ug <- x$qsip_growth_resp_c_ug * y$qsip_growth_resp_c_ug
  x$qsip_resp_c_ug <- x$qsip_resp_c_ug * y$qsip_resp_c_ug
  x},
  sankey_data_gen,
  c_adjust)

sankey_data_gen <- do.call(rbind, sankey_data_gen)

rm(c_adjust, c_lowered, c_summary)

# identify which phyla account for 99% of total uptake across all systems
sankey_data_gen$uptake_c_ug <- with(sankey_data_gen, growth_c_ug + qsip_resp_c_ug)
uptake <- tapply(sankey_data_gen$uptake_c_ug, sankey_data_gen$Phylum, sum, na.rm=T)
uptake <- sort(uptake, decreasing=T)
uptake <- cumsum(uptake) / sum(uptake)
to_keep <- c(names(uptake[uptake < .99]),
             names(uptake[max(which(uptake < .99)) + 1])) 
# keep the top 6

# function to refactor taxonomic levels
refactor_tax <- function(x, new_levels=c(), phylum=NULL) {
  y <- as.character(x)
  y[is.na(y)] <- 'unknown'
  # if phylum is null, only use new_levels to refactor
  if(is.null(phylum)) {
    y[!y %in% new_levels] <- 'Other'
    y <- factor(y, levels=c(new_levels, 'Other'))
    return(y)
  }
  # if phylum is provided, label anything in the "Other" phylum as "Other
  if(!is.null(phylum)) {
    y[phylum=='Other'] <- 'Other'
    y[grep('metagenome|unknown', y)] <- paste0(phylum[grep('metagenome|unknown', y)], y[grep('metagenome|unknown', y)])
    y <- factor(y, levels=c(unique(y)[unique(y)!='Other'], 'Other'))
    return(y)
  }
}

# relevel phyla (and down) not in the top C users to "Other"
sankey_data_gen <- within(sankey_data_gen, {
  Phylum <- refactor_tax(Phylum, to_keep)
  Class <- refactor_tax(Class, phylum=Phylum)
  Order <- refactor_tax(Order, phylum=Phylum)
  Family <- refactor_tax(Family, phylum=Phylum)
  Genus <- refactor_tax(Genus, phylum=Phylum)
  Species <- refactor_tax(Species, phylum=Phylum)
})

# sum up carbon budget for each genus (don't need total uptake)
sankey_data_gen <- aggregate(cbind(growth_c_ug, qsip_maint_resp_c_ug, qsip_growth_resp_c_ug) ~ .,
                             sankey_data_gen[,!names(sankey_data_gen) %in% c('taxonID', 'Species', 'uptake_c_ug', 'qsip_resp_c_ug')],
                             sum, na.rm=T)

# add in a mortality/turnover/exudate budget
# We can add turnover rates calculated from our 18O-CUE calculations
# These are formulated as gross growth (or 18O-MBC) - net growth
# this leaves turnover of the initial biomass *independent* of growth
tau_avg <- aggregate(mbc_turnover ~ ecosystem + treatment, cue_calc, mean)
sankey_data_gen <- merge(sankey_data_gen, tau_avg, all.x=T)

sankey_data_gen <- within(sankey_data_gen, {
  mortality <- growth_c_ug * mbc_turnover
  mbc_turnover <- NULL
  growth_c_ug <- growth_c_ug - mortality
})


# organize into sankey form-------------------------------------------------------------

# the idea is to separate a single initial pool into separate Phyla
# C flowing into these phyla are then partitioned into maintenance or growth
# maintenance goes straight to CO2, growth is partitioned into some going to CO2, some to biomass
# some biomass is then lost to exudate production and turnover
# the final pools are CO2 and biomass

# convert numeric C fluxes to long form
sankey_data_gen <- reshape(sankey_data_gen,
                           idvar= grep('ug|mortality', names(sankey_data_gen), invert=T),
                           varying=list(c('growth_c_ug', 'qsip_maint_resp_c_ug', 'qsip_growth_resp_c_ug', 'mortality')),
                           times=c('growth', 'maint_resp', 'growth_resp', 'mortality'),
                           timevar='c_account',
                           v.names='c_flux_ug_c',
                           direction='long')
rownames(sankey_data_gen) <- NULL

# Re-order genus factor levels so that they are grouped by phyla
genus_order <- tapply(sankey_data_gen$c_flux_ug_c, sankey_data_gen[,c('Phylum', 'Genus')], sum)
genera <- colnames(genus_order)
genus_order <- split(genus_order, rownames(genus_order))
#
genus_order <- lapply(genus_order, function(x) {names(x) <- genera; x})
genus_order <- lapply(genus_order, function(x) x[!is.na(x)])
genus_order <- genus_order[c(to_keep, 'Other')]
#
genus_order <- lapply(genus_order, sort, decreasing=T)
genus_order <- lapply(genus_order, names)
genus_order <- Reduce(c, genus_order)
#
sankey_data_gen$Genus <- factor(sankey_data_gen$Genus, levels=genus_order)


# create final_pool, maint_or_growth, and initial_source columns (maybe could add in SOM)

# function which creates factors from C account and taxonomic level
# y is a list of vectors, each length 2, first item is identifiers, second is new value
# the order of the new identifiers will be the factoring order (e.g., y[[1]][2] then y[[2]][2] etc.)
sankey_grouping <- function(x, y, tax_level=NULL, reverse=F) {
  z <- x <- as.character(x)
  for(i in 1:length(y)) {
    z[grep(y[[i]][1], x)] <- y[[i]][2]
    if(length(y[[i]])==3) z[grep(y[[i]][1], x, invert=T)] <- y[[i]][3]
  }
  new_levels <- sapply(y, function(a) a[2:length(a)])
  new_levels <- Reduce(c, new_levels)
  new_levels <- new_levels[!duplicated(new_levels)]
  if(reverse) new_levels <- rev(new_levels)
  z <- factor(z, levels=new_levels)
  if(!is.null(tax_level)) z <- interaction(tax_level, z)
  return(z)
}

# could add another term parameterizing the amount of necromass kept in soil (i.e., not respired)
sankey_data_gen <- within(sankey_data_gen, {
  final_pool <- sankey_grouping(c_account,
                                list(c('resp|mortality', 'respiration', 'biomass')))
  mortality <- sankey_grouping(c_account,
                               list(c('maint_resp', 'maint_resp'),
                                    c('growth_resp', 'growth_resp'),
                                    c('mortality', 'mortality'),
                                    c('^growth$', 'biomass')),
                               Genus)
  mortality_2 <- sankey_grouping(c_account,
                                 list(c('maint_resp', 'maint_resp'),
                                      c('growth_resp', 'growth_resp'),
                                      c('mortality', 'mortality'),
                                      c('^growth$', 'biomass')))
  net_growth <- sankey_grouping(c_account,
                                list(c('maint_resp', 'maint_resp'),
                                     c('growth_resp', 'growth_resp'),
                                     c('^growth$|mortality', 'net_growth')),
                                Genus)
  maint_or_growth <- sankey_grouping(c_account,
                                     list(c('maint_resp', 'maintenance', 'growth')),
                                     Genus)
})


# save data
save(sankey_data_gen, file='data/c_flux_categories.RData')
