# modeling respiration values taking into account substrate availability and activity/dormancy

# Respir = Rg + Rm = [(r / CUE) - r] + [(m / CUE) - m]
#
# where:
#   r = observed growth rate
#   m = maintenance rate = r / 1000
#   CUE = carbon use efficiency

library(ggplot2)
library(readxl)


# load in data
load('data/absolute_growth.RData')
load('data/feature_table.RData')
load('data/taxonomic_data.RData')
#
load('data/co2_data_mbc.RData')
load('data/cue_18O_method.RData') 
copy_data <- read.csv('data/asv_copy_no_cell_mass.csv')

# load useful functions for 16S copy number matching
source('figure_creation_functions.R')

# # get general taxonomic info, rename "uncultured" or "unknown" IDs to NA
# load('../Data/phyloseq_data_silva.RData')
# tax <- as(dat_silva@tax_table, 'matrix')
# tax[grep('uncultured|uncultivated|unknown', tax, ignore.case=T)] <- NA
# tax <- data.frame(taxonID=rownames(tax), tax)
# rownames(tax) <- NULL
# for(i in 1:ncol(tax)) attr(tax[,i], 'names') <- NULL
# save(tax, file='../Data/taxonomic_data.RData')


# prepare data--------------------------------------------------------------------------------------

# convert to dgTMatrix, which stores indexing information more intuitively
asv <- as(asv, 'dgTMatrix')

# remove week 0 data
asv_wk0 <- asv[1:16,]
asv <- asv[17:nrow(asv),]

# construct data frame from relative abundances
resp_data <- data.frame(taxonID=asv@Dimnames[[2]][asv@j + 1], 
                        sampleID=asv@Dimnames[[1]][asv@i + 1],
                        abund=asv@x,
                        stringsAsFactors=T)
#
total_abund <- resp_data

# convert ASVs to relative abundance while keeping in sparse Matrix
# note, abundances have already been standardized to grams of dry-soil
abund <- Matrix::rowSums(asv)
asv <- Matrix::Diagonal(x=1 / abund) %*% asv

# construct data frame from relative abundances
rel_data <- data.frame(taxonID=asv@Dimnames[[2]][asv@j + 1], 
                       sampleID=asv@Dimnames[[1]][asv@i + 1],
                       rel_abund=asv@x,
                       stringsAsFactors=T)

# combine abundances and relative abundances (these are WEEK 1 abundances)
resp_data <- merge(resp_data, rel_data)

# remove light (unlabeled samples)
# if you calculate per group parameters (including birth rate, you could use all)
resp_data <- merge(resp_data, unique(all_data[,c('sampleID', 'replicate', 'ecosystem', 'treatment')]), all.x=T)
resp_data <- resp_data[!is.na(resp_data$treatment),]

# add in absolute growth and birth rates
resp_data <- merge(resp_data, 
                   all_data[,c('taxonID', 'sampleID', 'growth', 'birth_rate')],
                   all.x=T)

# remove rows with NA for birth rate (will have NA for other rates too)
resp_data <- resp_data[!is.na(resp_data$birth_rate),]


# imputation of 16S copy number using taxonomic information-----------------------------------------

# one taxon, e4d4fc890c5ad90695f3601f36a5c650, has multiple entries in copy_data, average those
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


# get taxonomic data for these 14K taxa ONLY
copy_data <- merge(copy_data, tax, all.x=T)

# get taxa represented in labeled samples in qSIP experiment
qsip_tax <- resp_data[!duplicated(resp_data$taxonID), 'taxonID', drop=F]
qsip_tax <- merge(qsip_tax, tax, all.x=T)

# identify where taxonomic info drops off, up to family
na_position <- apply(qsip_tax[,-1], 1, function(x) {y <- which(is.na(x)); ifelse(length(y)==0, 9, min(y))})
na_position[na_position > 5] <- 5

# must determine if best ID matches with any of Junhui's IDs, if not, must raise that ID to the next level up
best_id <- match_repeat(na_position, qsip_tax[,-1], copy_data)
qsip_tax$best_id <- factor(best_id, levels=1:5, labels=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family'))

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

# combine with respiration data
resp_data <- merge(resp_data, qsip_copy_data[,c('taxonID', 'copy_number', 'cell_mass_g', 'Kingdom', 'Family')], all.x=T)


# Calculate change in population biomass (productivity)---------------------------------------------

# activity is based on the total budget of active biomass we expect

# calculate population biomass-C in ug
resp_data <- within(resp_data, {
  cell_c_ug <- cell_mass_g * 1E6 * .2
  cell_c_ug <- cell_c_ug / 10
  pop_c_ug <- (abund / copy_number) * cell_c_ug
  growth[is.na(growth)] <- 0
})

# re-order data
resp_data <- resp_data[order(resp_data$birth_rate, resp_data$taxonID, decreasing=T),]


# model respiration rates---------------------------------------------------------------------------

# use 18O-MBC CUE as "ground-truth" of CUE for each treatment
cue_min <- aggregate(cue ~ treatment + ecosystem, cue_calc, min)
names(cue_min)[3] <- 'cue_min'
cue_range <- aggregate(cue ~ treatment + ecosystem, cue_calc, range)
cue_range$cue_range <- cue_range$cue[,2] - cue_range$cue[,1]
#
resp_data <- merge(resp_data, cue_min, all.x=T)
resp_data <- merge(resp_data, cue_range[,c('treatment', 'ecosystem', 'cue_range')], all.x=T)

# make CUE calculations
resp_data <- within(resp_data, {
  cue <- -(cue_range/.25) * (birth_rate - .5)^2 + cue_range + cue_min 
  cue[birth_rate==0 & !is.na(birth_rate)] <- 0
})

# Respir = Rg + Rm = (g/CUE  -  g)  +  (g/CUE  -  g) * 0.01

# calculate respiration rates
resp_data <- within(resp_data, {
  Rg_rate <- (birth_rate / cue) - birth_rate
  Rm_rate <- Rg_rate / 1E2
  resp_rate <- Rg_rate + Rm_rate
})


# create taxonomically insensitive model------------------------------------------------------------

# Here, organism's contribution to the C cycle should apply ONLY in proportion to it's abundance 
# (i.e., don't take into account growth)

resp_data_b <- total_abund[total_abund$sampleID %in% cue_calc$sampleID,]
#
conversions <- aggregate(abund ~ sampleID, resp_data_b, sum)
conversions <- merge(conversions, cue_calc[,c('sampleID', 'g_co2_c_g_soil_week', 'mbc_18O',
                                              'ecosystem', 'treatment')],
                     all.x=T)

# conversion between abundances against MBC and CO2 production
co2_per_16s <- coefficients(lm(g_co2_c_g_soil_week ~ abund + 0, conversions))
mbc_per_16s <- coefficients(lm(mbc_18O ~ abund + 0, conversions))

# apply conversions to per-taxon abundances
resp_data_b <- within(resp_data_b, {
  co2_per_tax <- abund * co2_per_16s
  mbc_per_tax <- abund * mbc_per_16s
  c_use_per_tax <- co2_per_tax + mbc_per_tax
})

# sum up C fluxes and abundances per sample
resp_data_b <- aggregate(cbind(abund, c_use_per_tax, mbc_per_tax, co2_per_tax) ~ sampleID, resp_data_b, sum)
resp_data_b <- merge(resp_data_b, cue_calc[,c('sampleID', 'ecosystem', 'treatment', 
                                              'replicate', 'g_co2_c_g_soil_week', 'mbc_18O')], 
                     all.x=T)


# calculate respiration flux------------------------------------------------------------------------

# respiration rate is built on growth rate which is per day, change to per week
resp_data <- within(resp_data, {
  qsip_resp_c_ug <- resp_rate * 7 * pop_c_ug
  qsip_growth_resp_c_ug <- Rg_rate * 7 * pop_c_ug
  qsip_maint_resp_c_ug <- Rm_rate * 7 * pop_c_ug
})

# change factor levels
resp_data$ecosystem <- factor(resp_data$ecosystem, levels=c('MC', 'PP', 'PJ', 'GL'))

# sum up respiration per sample
resp_summary <- aggregate(qsip_resp_c_ug ~ ecosystem + treatment + sampleID + replicate, resp_data, sum)

# compare to CO2 data
load('data/co2_data_mbc.RData')
co2 <- co2[co2$experiment=='high substrate',]
co2 <- co2[,!names(co2) %in% c('experiment', 'week')]
names(co2)[grep('replicate', names(co2))] <- 'replicate'
co2 <- within(co2, {
  treatment <- factor(treatment, levels=levels(treatment), labels=c('18O', 'C', 'CN', 'C', 'CN'))
  replicate <- factor(replicate)
})

resp_summary <- merge(resp_summary, co2, all.x=T)
#
resp_summary_per_sample <- resp_summary

resp_summary_se <- aggregate(cbind(qsip_resp_c_ug, g_co2_c_g_soil_week) ~ ecosystem + treatment, resp_summary, se)
resp_summary <- aggregate(cbind(qsip_resp_c_ug, g_co2_c_g_soil_week) ~ ecosystem + treatment, resp_summary, mean)

resp_summary <- merge(resp_summary, resp_summary_se, 
                      by=c('ecosystem', 'treatment'),
                      suffixes=c('', '_se'))


# save data
save(resp_summary_per_sample, resp_summary, resp_data, resp_data_b, file='data/respiration_18O_cue_draft_4.RData')
