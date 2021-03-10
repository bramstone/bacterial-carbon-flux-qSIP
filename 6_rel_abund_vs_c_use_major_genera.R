# note: results don't appreciably change if we keep 18O soils in or not

load('data/feature_table.RData')
load('data/taxonomic_data.RData')
load('data/absolute_growth.RData')
load('data/c_flux_categories.RData')

# match sampleID with ecosystem, treatment, and replicate-----------------------

sample_data <- unique(all_data[,c('sampleID', 'ecosystem', 'treatment', 'replicate')])

# get relative abundances of genera---------------------------------------------

# convert to dgTMatrix, which stores indexing information more intuitively
asv <- as(asv, 'dgTMatrix')

# remove week 0 data and non-qSIP samples
asv_wk0 <- asv[1:16,]
asv <- asv[17:nrow(asv),]
asv <- asv[rownames(asv) %in% sample_data$sampleID,]

# construct data frame from abundances
seq_data <- data.frame(taxonID=asv@Dimnames[[2]][asv@j + 1], 
                        sampleID=asv@Dimnames[[1]][asv@i + 1],
                        abund=asv@x,
                        stringsAsFactors=T)

# merge with taxonomic info
seq_data <- merge(seq_data, tax, all.x=T)

# merge by genus
seq_data <- aggregate(abund ~ ., seq_data[,!names(seq_data) %in% c('taxonID', 'Species')], sum)

# construct relative abundances
seq_data <- split(seq_data, seq_data$sampleID)
seq_data <- lapply(seq_data, function(x) {x$rel_abund <- x$abund / sum(x$abund); x})
seq_data <- do.call(rbind, seq_data)
row.names(seq_data) <- NULL

# match with ecosystem and treatment info
seq_data <- merge(seq_data, sample_data, all.x=T)

# average relative abundance across replicates
seq_data <- aggregate(rel_abund ~ ., 
                      seq_data[,!names(seq_data) %in% c('sampleID', 'replicate', 'abund')], 
                      mean)

# combine relative abundances with C use data-----------------------------------

c_use  <- aggregate(c_flux_ug_c ~., 
                    sankey_data_gen[!names(sankey_data_gen) %in% c('c_account', 'maint_or_growth', 'net_growth', 
                                                                   'mortality', 'mortality_2', 'final_pool')],
                    sum)

# relativize C use per treatment and ecosystem
c_use <- split(c_use, interaction(c_use$treatment, c_use$ecosystem))
c_use <- lapply(c_use, function(x) {
  x$c_flux_rel <- x$c_flux_ug_c / sum(x$c_flux_ug_c)
  x$c_flux_ug_c <- NULL
  x})
c_use <- do.call(rbind, c_use)

# combine 
abund_c_use <- merge(c_use, seq_data, all=F)
abund_c_use <- droplevels(abund_c_use)
nlevels(abund_c_use$Genus)
# there are 103 genera with a C contribution

# determine how many genera to look at.............
# look at cumulative abundance 
abund_c_use <- abund_c_use[order(abund_c_use$rel_abund, decreasing=T),]
abund_c_use <- split(abund_c_use, interaction(abund_c_use$treatment, abund_c_use$ecosystem))
abund_c_use <- lapply(abund_c_use, function(x) {x$cumul <- cumsum(x$rel_abund); x})
sapply(abund_c_use, function(x) min(which(x$cumul > .5))) # top 8-27 to get to 50% rel abundance


# look at cumulative C use
abund_c_use <- lapply(abund_c_use, function(x) x[order(x$c_flux_rel, decreasing=T),])
abund_c_use <- lapply(abund_c_use, function(x) {x$cumul <- cumsum(x$c_flux_rel); x})
sapply(abund_c_use, function(x) min(which(x$cumul > .5))) # only need 1-8 to get 50% C use
abund_c_use <- do.call(rbind, abund_c_use)

# full (103 genera)
abund_c_use_full <- abund_c_use
abund_c_use_full$cumul <- NULL

# let's just keep with top 22 genera in each treatment:ecosystem combo
abund_c_use$cumul <- NULL
abund_c_use <- abund_c_use[order(abund_c_use$rel_abund, decreasing=T),]
abund_c_use <- split(abund_c_use, interaction(abund_c_use$treatment, abund_c_use$ecosystem))
abund_c_use <- lapply(abund_c_use, '[', 1:22,)
abund_c_use <- do.call(rbind, abund_c_use)

# remove some of the lowest contributors to C flux
genus_contrib <- aggregate(c_flux_rel ~ Genus, abund_c_use, mean)
genus_contrib <- genus_contrib[order(genus_contrib$c_flux_rel, decreasing=T),]
genus_contrib$cumul <- cumsum(genus_contrib$c_flux_rel)
# keep only genera that, on average, represent 95% of C cycling
# this is the FINAL selection criteria that should be reported
abund_c_use <- abund_c_use[abund_c_use$Genus %in% genus_contrib$Genus[1:min(which(genus_contrib$cumul > .95))],]
abund_c_use <- droplevels(abund_c_use)

# re-level factors based off relative abundance order
abund_c_use$Genus <- factor(abund_c_use$Genus, levels=unique(abund_c_use$Genus))

# convert relative contributions to single column
abund_c_use <- reshape(abund_c_use,
                       idvar=names(abund_c_use)[1:8],
                       varying=list(c('c_flux_rel', 'rel_abund')),
                       times=c('Carbon use', '16S rRNA gene abundance'),
                       timevar='metric',
                       v.names='relative_contribution',
                       direction='long')
row.names(abund_c_use) <- NULL

# save data
save(abund_c_use, abund_c_use_full, file='data/rel_abund_vs_c_use_by_treatment_and_ecosystem_draft_4.RData')