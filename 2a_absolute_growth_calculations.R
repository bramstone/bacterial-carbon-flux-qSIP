# Another calculation of absolute birth and death both exponential and linear growth rates
# and using both week 0 abundances and week 1
# Using PER-SAMPLE growth rates

# load in data
load('data/phyloseq_data_silva.RData')
# ASV table from phyloseq data, abundances summed at the replicate level, in S4 Matrix format
load('data/feature_table.RData') 
load('data/birth_rates_asv_silva_per_sample.RData')
load('data/death_rates_asv_silva_per_sample.RData')
#
treatments <- read.csv('data/density_and_16S_and_soil.csv')
copy_data <- read.csv('data/asv_copy_no_cell_mass.csv')
field_replicates <- read.table('data/old_map_merged.txt', header=T, check.names=F, comment.char='')

# load useful functions
source('figure_creation_functions.R')


# add in replicate numbers for exact replicate matching for growth rates
field_replicates <- field_replicates[,c('Description', 'Replicate')]
names(field_replicates) <- c('sampleID', 'replicate')
field_replicates <- field_replicates[grep('w0|W1_', field_replicates$sampleID),]
field_replicates$sampleID <- gsub('-', '_', toupper(field_replicates$sampleID))

# convert ASVs to relative abundance while keeping in sparse Matrix
# note, abundances have already been standardized to grams of dry-soil
# if you get an error about S4 method calling, restart R and try again
abund <- Matrix::rowSums(asv)
asv <- Matrix::Diagonal(x=1 / abund) %*% asv

# select taxa desired
tax_names <- list(birth$taxonID, death$taxonID)
tax_names <- Reduce(union, tax_names)
asv <- asv[,colnames(asv) %in% tax_names]
asv <- t(as(asv, 'matrix'))

# put abundances into data.frame
abund <- data.frame(sampleID=names(abund), total_abund=abund)


# Convert to long form--------------------------------------------------------------------

# abundances and relative abundances...........

# week 0 abundances and relative abundances
week0 <- asv[,grep('W0', colnames(asv))]
week0 <- qsip_long(week0, group_name='sampleID', value='rel_abund_wk0')
week0 <- merge(week0, abund, all.x=T)
week0 <- within(week0, {
  ecosystem <- sub('W0_(\\w{2})_(\\d?)', '\\1', sampleID)
  replicate <- sub('W0_(\\w{2})_(\\d?)', '\\2', sampleID)
  sampleID <- NULL
  avg_16S_copies <- rel_abund_wk0 * total_abund
  total_abund <- NULL
})

# week 1 abundances and relative abundances
week1 <- asv[,grep('W0', colnames(asv), invert=T)]
week1 <- qsip_long(week1, group_name='sampleID', value='rel_abund_wk1')
week1 <- merge(week1, abund, all.x=T)
week1 <- within(week1, {
  avg_16S_copies_wk1 <- rel_abund_wk1 * total_abund
  total_abund <- NULL
})

# combine
# abund <- merge(week1, week0, all.x=T)
# abund <- abund[abund$avg_16S_copies > 0 | abund$avg_16S_copies_wk1 > 0,]


# birth rates..................................

# make negative birth rate values 0
birth_vals <- as.matrix(birth[,-1])
birth_vals[birth_vals < 0] <- 0
birth[,-1] <- birth_vals

birth <- qsip_long(birth, 'sampleID', 'birth_rate')
birth <- within(birth, {
  ecosystem <- sub('W1_(\\w{2})_(\\w*)', '\\1', sampleID)
})

birth <- merge(birth, field_replicates, all.x=T)


# convert death rates.........................

# regular death rate adjustments
rownames(death) <- death$taxonID
death <- death_per_site <- as.matrix(death[,-1])
adjust <- sort(death, decreasing=T)
adjust <- quantile(adjust, .95)
# death <- death - adjust
death[death > 0] <- 0
death <- data.frame(taxonID=rownames(death), death)

# PER-SITE death rate adjustments
adjusts <- list('MC', 'PP', 'PJ', 'GL')
adjusts <- lapply(adjusts, function(x) sort(death_per_site[,grep(x, colnames(death_per_site))], decreasing=T))
adjusts <- lapply(adjusts, quantile, .95)
# adjust death values
death_per_site <- lapply(list('MC', 'PP', 'PJ', 'GL'), function(x) death_per_site[,grep(x, colnames(death_per_site))])
death_per_site <- Map('-', death_per_site, adjusts)
death_per_site <- do.call(cbind, death_per_site)
#
death_per_site[death_per_site > 0] <- 0
death_per_site <- data.frame(taxonID=rownames(death_per_site), death_per_site)


death <- qsip_long(death, 'sampleID', 'death_rate')
death <- within(death, {
  ecosystem <- sub('W1_(\\w{2})_(\\w*)', '\\1', sampleID)
})

death_per_site <- qsip_long(death_per_site, 'sampleID', 'death_rate')
death_per_site <- within(death_per_site, {
  ecosystem <- sub('W1_(\\w{2})_(\\w*)', '\\1', sampleID)
})
names(death_per_site)[3] <- 'death_rate_per_site'

death <- merge(death, field_replicates, all.x=T)
death_per_site <- merge(death_per_site, field_replicates, all.x=T)


# Adjust and combine data-----------------------------------------------------------------

# map samples to treatment
sam_trt <- unique(treatments[, c('sampleID', 'isotope.treatment')])
names(sam_trt)[2] <- 'treatment'
sam_trt <- sam_trt[grep('18O', sam_trt$treatment),]
sam_trt$treatment <- factor(factor(sam_trt$treatment),
                            levels=c('18O', '12C_18O', '12C_18O_N'),
                            labels=c('18O', 'C', 'CN'))

# add sample number and treatment data to qSIP rates
qsip_rates <- list(birth, death, death_per_site)
qsip_rates <- lapply(qsip_rates, merge, sam_trt, all.x=T)
qsip_rates <- Reduce(function(x, y) merge(x, y, all=T), qsip_rates)

# combine qSIP and abundances, and remove 0 abundances in BOTH week0 and week1 cases
# when you remove missing cases from BOTH, then using either formulation of net growth will produce the same result
all_data <- Reduce(function(x, y) merge(x, y, all.x=T), list(qsip_rates, week0, week1))
all_data <- all_data[all_data$avg_16S_copies > 0 & all_data$avg_16S_copies_wk1 > 0,]


# get copy number data estimates.............................
# note: (only 14K taxon have good estimates, others should be imputed based on taxonomy)

# one taxon, e4d4fc890c5ad90695f3601f36a5c650, has multiple entries, average those
arcobacter <- copy_data[copy_data$taxonID=='e4d4fc890c5ad90695f3601f36a5c650',]
accession_no <- arcobacter$Assembly_Accession[1]
#
arcobacter <- arcobacter[,sapply(arcobacter, is.numeric)]
arcobacter <- lapply(arcobacter, mean)
arcobacter <- data.frame(taxonID='e4d4fc890c5ad90695f3601f36a5c650', Assembly_Accession=accession_no, arcobacter)
arcobacter$copy_number <- floor(arcobacter$copy_number)
arcobacter <- arcobacter[,names(copy_data)]
#
copy_data <- copy_data[copy_data$taxonID!='e4d4fc890c5ad90695f3601f36a5c650',]
copy_data <- rbind(copy_data, arcobacter)
rm(arcobacter)

# get general taxonomic info, make "uncultured" or "unknown" type identifications be NA
tax <- as(dat_silva@tax_table, 'matrix')
tax[grep('uncultured|uncultivated|unknown|metagenome', tax, ignore.case=T)] <- NA
tax <- data.frame(taxonID=rownames(tax), tax)
rownames(tax) <- NULL
for(i in 1:ncol(tax)) attr(tax[,i], 'names') <- NULL

# get taxonomic data for these 14K taxa ONLY
copy_data <- merge(copy_data, tax, all.x=T)

# get taxa represented in qSIP experiment
qsip_tax <- all_data[!duplicated(all_data$taxonID), 'taxonID', drop=F]
qsip_tax <- merge(qsip_tax, tax, all.x=T)

# identify where taxonomic info drops off, up to family
na_position <- apply(qsip_tax[,-1], 1, function(x) {y <- which(is.na(x)); ifelse(length(y)==0, 9, y)})
na_position[na_position > 5] <- 5

# must determine if best ID matches with any of Junhui's IDs, if not, must raise that ID to the next level up
na_position <- match_repeat(na_position, qsip_tax, copy_data)
qsip_tax$best_id <- factor(na_position - 1, levels=1:5, labels=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family'))

# identify those with matching information from Junhui's work, and those that must be estimated
qsip_tax_true <- merge(qsip_tax, copy_data[,c('taxonID', 'copy_number', 'genome_size_bp', 'cell_mass_g')], all=F)
qsip_tax_est <- merge(qsip_tax, copy_data[,c('taxonID', 'copy_number', 'genome_size_bp', 'cell_mass_g')], all.x=T)
qsip_tax_est <- qsip_tax_est[!qsip_tax_est$taxonID %in% unique(qsip_tax_true$taxonID),]
qsip_tax_est <- qsip_tax_est[,!names(qsip_tax_est) %in% c('copy_number', 'genome_size_bp', 'cell_mass_g')]

# split those taxa that must be estimated by lowest taxa possible
qsip_tax_est <- split(qsip_tax_est, qsip_tax_est$best_id)
qsip_tax_est <- qsip_tax_est[!names(qsip_tax_est) %in% 'Family']

# estimate copy number for different taxonomic levels (up to Order since Junhui calculated values up to Family)
kingdom_copy <- aggregate(cbind(copy_number, genome_size_bp, cell_mass_g) ~ Kingdom, copy_data, median)
phyla_copy <- aggregate(cbind(copy_number, genome_size_bp, cell_mass_g) ~ Phylum, copy_data, median)
class_copy <- aggregate(cbind(copy_number, genome_size_bp, cell_mass_g) ~ Class, copy_data, median)
order_copy <- aggregate(cbind(copy_number, genome_size_bp, cell_mass_g) ~ Order, copy_data, median)
copies <- list(kingdom_copy, phyla_copy, class_copy, order_copy)
names(copies) <- c('Kingdom', 'Phylum', 'Class', 'Order')
rm(list=ls(pattern='_copy'))

# Match estimated copy numbers, genome sizes, and cell masses with average values
qsip_tax_est <- Map(function(x, y, z) merge(x, y[,c(z, 'copy_number', 'genome_size_bp', 'cell_mass_g')], all.x=T), 
                    qsip_tax_est,
                    copies,
                    c('Kingdom', 'Phylum', 'Class', 'Order'))
qsip_tax_est <- do.call(rbind, qsip_tax_est)

# combine known taxa with estimated taxa
qsip_tax_est <- qsip_tax_est[,names(qsip_tax_true)]
qsip_copy_data <- rbind(qsip_tax_true, qsip_tax_est)
rm(list=ls(pattern='qsip_tax|copies|na'))

# combine with pop flux data
all_data <- merge(all_data, qsip_copy_data[,c('taxonID', 'copy_number', 'cell_mass_g')], all.x=T)


# Calculate population flux---------------------------------------------------------------

# If using Nt - Nt*e^-gt, set reverse=T
reverse <- T

all_data <- within(all_data, {
  # copy number and relative abundance to use
  abundance_16S <- avg_16S_copies_wk1
  rel_abund <- rel_abund_wk1
  # pop_cell <- abundance_16S / copy_number
  pop_cell <- abundance_16S
  #
  # exponential
  birth <- birth_rate
  death <- death_rate
  net <- birth %+% death
  #
  growth <- pop_flux(birth, pop_cell, reverse=reverse)
  loss <- pop_flux(death, pop_cell, reverse=reverse)
  net_calc <- pop_flux(net, pop_cell, reverse=reverse)
  net_diff <- growth %+% loss
  #
  # growth_ug <- pop_flux(birth, pop_cell_ug, reverse=reverse)
  # loss_ug <- pop_flux(death, pop_cell_ug, reverse=reverse)
  # net_calc_ug <- pop_flux(net, pop_cell_ug, reverse=reverse)
  # net_diff_ug <- growth_ug %+% loss_ug
  #
  abundance_16S <- NULL
  rel_abund <- NULL
  birth <- NULL
  death <- NULL
  net <- NULL
})

# remove NAs
# None of these are the Thaumarchaeota
all_data <- all_data[!is.na(all_data$copy_number),]


# save estimates of absolute growth
save(all_data, file='data/absolute_growth.RData')
