# Calculation of per-capita birth/death for Dimensions data set, converting abundances to 16S copies per tube per g soil
# Using PER-SAMPLE calculations

# install devtools and BiocManager
install.packages('devtools')
install.packages('BiocManager')

# install phyloseq and qsip using the utilities on devtools and BiocManagerK
BiocManager::install('phyloseq')
devtools::install_github('bramstone/qsip', ref='master')

# load packages
library(qsip)
library(phyloseq)

# sequencing data have been organized into phyloseq format
load('data/phyloseq_data_silva.RData')
treatments <- read.csv('data/density_and_16S_and_soil.csv')
field_replicates <- read.table('data/old_map_merged.txt', header=T, check.names=F, comment.char='')

# remove samples with low sequence read counts
dat_silva <- prune_samples(colSums(otu_table(dat_silva)) > 4000, dat_silva) # 32 samples, ~2.1% of total

# fix a problematic sample where treatment was coded incorrectly
dat_silva@sam_data$treatment[dat_silva@sam_data$sampleID=='W1_PJ_24' & dat_silva@sam_data$treatment=='C_N'] <- 'C'


# convert abundances from 16S copies per uL to 16S copies per g soil----------------------

sam <- as(sample_data(dat_silva), 'data.frame')

abund <- split(sam[,c('sampleID', 'isotope', 'treatment', 
                      'ecosystem', 'avg_16S_copies', 'timepoint',
                      'dna_ng_ul', 'uL_added', 'dry_soil_wt_g')], 
               sam$sampleID)

# get abundances per replicate
abund <- lapply(abund, function(x) {y <- x[1,]; y$avg_16S_copies <- sum(x$avg_16S_copies, na.rm=T); y})
abund <- do.call(rbind, abund)

# ~1ug DNA was added to fractioned tubes, so standardize using DNA content in each tube
abund <- within(abund, {
  dna_ug_tube <- (dna_ng_ul / 1000) * 100
  avg_16S_copies <- avg_16S_copies * dna_ug_tube
  avg_16S_g_soil <- avg_16S_copies / dry_soil_wt_g
})

# adjust 16S abundances in each fraction to the total amount of DNA in the tube
sam <- split(sam, sam$sampleID)
abund <- split(abund$avg_16S_g_soil, abund$sampleID)
sam <- sam[match(names(abund), names(sam))]
sam <- Map(function(x,y) {x$avg_16S_g_soil <- y * (x$avg_16S_copies/sum(x$avg_16S_copies)); x}, sam, abund)
sam <- do.call(rbind, sam)

# add back to phyloseq object
dat_silva@sam_data$avg_16S_g_soil <- sam$avg_16S_g_soil[match(dat_silva@sam_data$SampleID, sam$SampleID)]
rm(sam, abund)

# we will convert 16S copies to per-cell when we calculate absolute growth

# add in replicate numbers for exact replicate matching for growth rates
field_replicates <- field_replicates[,c('Description', 'Replicate')]
names(field_replicates) <- c('sampleID', 'replicate')
field_replicates <- field_replicates[grep('w0|W1_', field_replicates$sampleID),]
field_replicates$sampleID <- gsub('-', '_', toupper(field_replicates$sampleID))

# add back to phyloseq object
dat_silva@sam_data$replicate <- field_replicates$replicate[match(dat_silva@sam_data$sampleID, field_replicates$sampleID)]
rm(field_replicates)

# 16O and 18O-------------------------------------------------------------------------------------------------

# select data from only 18O amendment (ignore C and C+N treatments) and keep week 0 data
dat_con <- prune_samples(dat_silva@sam_data$treatment=='control' | is.na(dat_silva@sam_data$treatment), dat_silva)

# specify qsip-related dat_cona from phyloseq object, create phylosip object
dat_con <- specify_qsip(dat_con,
                        density='Density.g.ml',
                        abund='avg_16S_g_soil',
                        rep_id='sampleID',
                        rep_num='replicate',
                        rep_group='ecosystem',
                        iso='18O',
                        iso_trt='isotope',
                        timepoint='timepoint')

# specify filtering level:
#   for each ecosystem (the grouping factor), taxa must occur in more than 2 of the 3 replicates
#   and in more than 5 of the fractions that make up a replicate
dat_con@qsip@filter_levels <- create_filters(2, 5, soft=1)

# calculate population growth and death
dat_con <- calc_pop(dat_con, separate_label=T, filter=T, correction=T, match_replicate=T)

birth <- dat_con@qsip[['birth_rate']]
death <- dat_con@qsip[['death_rate']]

par(mfrow=c(2,1), mar=c(2,2,1,1)+.1)
hist(birth); hist(death)
# Focus just on birth to relate to CO2 flux
birth_con <- birth


# Look at C addition treatment------------------------------------------------------------------------------------------

# isolate just the C treatment, 12C and 12C + 18O
dat_c <- prune_samples(dat_silva@sam_data$treatment=='C' | is.na(dat_silva@sam_data$treatment), dat_silva)
treatments_c <- treatments[match(dat_c@sam_data$sampleID, treatments$sampleID),]
dat_c@sam_data$isotope_2 <- treatments_c$isotope.treatment
dat_c <- prune_samples(dat_c@sam_data$isotope_2!='13C' & dat_c@sam_data$isotope_2!='18O'| dat_c@sam_data$timepoint==0, dat_c)

# specify qsip-related data from phyloseq object, create phylosip object
dat_c <- specify_qsip(dat_c,
                      density='Density.g.ml',
                      abund='avg_16S_g_soil',
                      rep_id='sampleID',
                      rep_num='replicate',
                      rep_group='ecosystem',
                      iso='18O',
                      iso_trt='isotope_2',
                      timepoint='timepoint')

# specify filtering level
dat_c@qsip@filter_levels <- create_filters(2, 5, soft=1)

# calculate population growth and death
dat_c <- calc_pop(dat_c, separate_label=T, filter=T, correction=T, match_replicate=T)

birth <- dat_c@qsip[['birth_rate']]
death <- dat_c@qsip[['death_rate']]

par(mfrow=c(2,1), mar=c(2,2,1,1)+.1)
hist(birth); hist(death)
birth_c <- birth


# Look at C+N addition treatment------------------------------------------------------------------------------------------

# isolate just the C+N treatment, 12C and 12C + 18O
dat_cn <- prune_samples(dat_silva@sam_data$treatment=='C_N' | is.na(dat_silva@sam_data$treatment), dat_silva)
treatments_cn <- treatments[match(dat_cn@sam_data$sampleID, treatments$sampleID),]
dat_cn@sam_data$isotope_2 <- sub('_N', '', treatments_cn$isotope.treatment)
dat_cn <- prune_samples(dat_cn@sam_data$isotope_2!='13C' & dat_cn@sam_data$isotope_2!='18O' | dat_cn@sam_data$timepoint==0, dat_cn)

# specify qsip-related data from phyloseq object, create phylosip object
dat_cn <- specify_qsip(dat_cn,
                       density='Density.g.ml',
                       abund='avg_16S_g_soil',
                       rep_id='sampleID',
                       rep_num='replicate',
                       rep_group='ecosystem',
                       iso='18O',
                       iso_trt='isotope_2',
                       timepoint='timepoint')

# specify filtering levels
dat_cn@qsip@filter_levels <- create_filters(2, 5, soft=1)

# calculate population growth and death
dat_cn <- calc_pop(dat_cn, separate_label=T, filter=T, correction=T, match_replicate=T)

birth <- dat_cn@qsip[['birth_rate']]
death <- dat_cn@qsip[['death_rate']]

par(mfrow=c(2,1), mar=c(2,2,1,1)+.1)
hist(birth); hist(death)
birth_cn <- birth


# 13C enrichment under high C treatment----------------------------------------------------------

dat_13c <- prune_samples(dat_silva@sam_data$treatment=='C' | is.na(dat_silva@sam_data$treatment), dat_silva)
dat_13c <- prune_samples(dat_13c@sam_data$isotope!='12C_18O' & dat_13c@sam_data$isotope!='18O' | dat_13c@sam_data$timepoint==0, dat_13c)

# specify qsip-related data from phyloseq object, create phylosip object
dat_13c <- specify_qsip(dat_13c,
                        density='Density.g.ml',
                        abund='avg_16S_copies',
                        rep_id='sampleID',
                        rep_group='ecosystem',
                        iso='13C',
                        iso_trt='isotope',
                        timepoint='timepoint')

# calculate atom excess fraction / atom percent excess
dat_13c@qsip@filter_levels <- create_filters(2, 5, soft=1)
dat_13c <- calc_excess(dat_13c, separate_label=T, filter=T, correction=T)

# extract excess atom fraction (EAF)
eaf_c <- dat_13c@qsip[['atom_excess']]


# 13C enrichment under high CN treatment---------------------------------------------------------

dat_13cn <- prune_samples(dat_silva@sam_data$treatment=='C_N' | is.na(dat_silva@sam_data$treatment), dat_silva)
dat_13cn <- prune_samples(dat_13cn@sam_data$isotope!='12C_18O' & dat_13cn@sam_data$isotope!='18O' | dat_13cn@sam_data$timepoint==0, dat_13cn)

# specify qsip-related data from phyloseq object, create phylosip object
dat_13cn <- specify_qsip(dat_13cn,
                         density='Density.g.ml',
                         abund='avg_16S_copies',
                         rep_id='sampleID',
                         rep_group='ecosystem',
                         iso='13C',
                         iso_trt='isotope',
                         timepoint='timepoint')

# calculate atom excess fraction / atom percent excess
dat_13cn@qsip@filter_levels <- create_filters(2, 5, soft=1)
dat_13cn <- calc_excess(dat_13cn, separate_label=T, filter=T, correction=T)

# extract EAF
eaf_cn <- dat_13cn@qsip[['atom_excess']]


# Combine growth values from all treatments--------------------------------------------------

make.data.frame <- function(x) {
  taxonID <- rownames(x)
  y <- as.data.frame(x)
  y <- data.frame(taxonID, y, stringsAsFactors=F)
  rownames(y) <- NULL
  return(y)
}

# merge into one data.frame
birth <- Reduce(function(x,y) merge(x,y,by='taxonID', all=T), 
                list(make.data.frame(birth_con),
                     make.data.frame(birth_c),
                     make.data.frame(birth_cn)))

# remove rows with all NA
all_nas <- apply(sapply(birth[,2:ncol(birth)], is.na), 1, all)
birth <- birth[!all_nas,]
rownames(birth) <- NULL


# combine death values from all treatments-------------------------------------------------------

data_names <- paste0('dat_', c('con', 'c', 'cn'))
output_names <- paste0('death_', c('con', 'c', 'cn'))

# for-loop to extract all aef values and clean rownames
for(i in 1:length(data_names)) {
  y <- get(data_names[i])@qsip[['death_rate']]
  y <- as(y, 'matrix')
  assign(output_names[i], y)
}

# merge into one data.frame
death <- Reduce(function(x,y) merge(x,y,by='taxonID', all=T), 
                list(make.data.frame(death_con),
                     make.data.frame(death_c),
                     make.data.frame(death_cn)))

# remove rows with all NA
all_nas <- apply(sapply(death[,2:ncol(death)], is.na), 1, all)
death <- death[!all_nas,]
rownames(death) <- NULL


# combine enrichment from 13C treatments---------------------------------------------------------

eaf_13c <- merge(eaf_c, eaf_cn, by='taxonID', all=T)

# remove rows with all NA
all_nas <- apply(sapply(eaf_13c[,2:ncol(eaf_13c)], is.na), 1, all)
eaf_13c <- eaf_13c[!all_nas,]
rownames(eaf_13c) <- NULL


# export birth, death, and enrichment------------------------------------------------------------


save(birth, file='data/birth_rates_asv_silva_per_sample.RData')
save(death, file='data/death_rates_asv_silva_per_sample.RData')
save(eaf_13c, file='data/enrichment_13c_silva_per_sample.RData')