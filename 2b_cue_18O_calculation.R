# calculation of 18O-CUE using qSIP enrichment

install.packages('readxl')

library(ggplot2)

load('data/absolute_growth.RData') 
load('data/original-data.RData') 
load('data/co2_data.RData')
mbc <- readxl::read_excel('data/Priming_Dimensions_Database_Archive.xlsx',
                          sheet='mbc',
                          range='B4:K132',
                          na='N/A')
mbc <- data.frame(mbc[,-c(2,6)])


# how much of the total relative abundance pool do these frequent taxa make up?
# important since we want to extrapolate enrichment in these taxa to enrichment to the ENTIRE community
tax_depth <- split(all_data[,c('taxonID', 'sampleID', 'rel_abund_wk1')], all_data$sampleID)
tax_depth <- lapply(tax_depth, function(x) {y <- sum(x$rel_abund_wk1, na.rm=T); data.frame(x$sampleID[1], y)})
tax_depth <- do.call(rbind, tax_depth)
names(tax_depth) <- c('sampleID', 'qsip_coverage')


tax_depth <- merge(tax_depth, 
                   unique(all_data[,c('sampleID', 'ecosystem', 'treatment', 'replicate')]), 
                   all.x=T)
tax_depth$ecosystem <- factor(tax_depth$ecosystem, levels=c('MC', 'PP', 'PJ', 'GL'))

ggplot(tax_depth, aes(ecosystem, qsip_coverage)) +
  geom_boxplot() +
  geom_hline(data=aggregate(qsip_coverage ~ treatment, tax_depth, mean), 
             aes(yintercept=qsip_coverage), color='red') +
  facet_grid(treatment ~ .)
# looks like 20 - 50% coverage from the qSIP-relevant taxa
# little bias by treatment

# renormalize relative abundances to 1
all_data <- merge(all_data, tax_depth, all.x=T)
all_data$rel_abund_wk1 <- all_data$rel_abund_wk1 / all_data$qsip_coverage
all_data$qsip_coverage <- NULL


# Calculate enrichment------------------------------------------------------------------

enrich <- split(all_data[,c('taxonID', 'sampleID', 'rel_abund_wk1', 'aef')], all_data$sampleID)
enrich_cumul <- enrich
enrich <- lapply(enrich, function(x) {y <- sum(x$rel_abund_wk1 * x$aef, na.rm=T); data.frame(x$sampleID[1], y)})
enrich <- do.call(rbind, enrich)
rownames(enrich) <- NULL
names(enrich) <- c('sampleID', 'enrichment')

# look at cumulative enrichment
enrich_cumul <- lapply(enrich_cumul, function(x) {
  x$enrich_contrib <- x$rel_abund_wk1 * x$aef
  x$enrich_contrib[is.na(x$enrich_contrib)] <- 0
  x})
#
enrich_cumul <- lapply(enrich_cumul, function(x) x[order(x$enrich_contrib, decreasing=T),])
enrich_cumul <- lapply(enrich_cumul, function(x) {
  x$enrich_contrib <- cumsum(x$enrich_contrib)
  x$tax_rank <- 1:nrow(x)
  x})
#
enrich_cumul <- do.call(rbind, enrich_cumul)
#
enrich_cumul <- merge(enrich_cumul, unique(all_data[,c('sampleID', 'ecosystem', 'treatment')]), all.x=T)

# plot cumulative enrichment
ggplot(enrich_cumul, aes(tax_rank, enrich_contrib, color=treatment)) +
  geom_line(aes(group=sampleID), size=1) +
  facet_wrap(~ ecosystem, ncol=2) +
  ylab(expression('Cumulative contribution to '^18*O~'enrichment')) +
  xlab('Taxon rank')
# looks like using relative abundance-weighted enrichemnt may approximate
# the main patterns of assemblage-level enrichment


# Calculate DNA per g dry soil----------------------------------------------------------

abund <- within(abund, {
  dna_ug_g_ds <- dna_ug_tube / dry_soil_wt_g
  replicate <- as.numeric(sub('(.*)_(\\d+)$', '\\2', sampleID))
})

# focus on pre-incubation data
abund_sub <- abund[abund$timepoint==0,]


# Calculate MBC ~ DNA relationship-------------------------------------------------------

# fix mbc names
names(mbc) <- gsub('[^[:alnum:]]+', '_', names(mbc)) 
names(mbc) <- gsub('^_?|_?$', '', names(mbc)) 
names(mbc) <- tolower(names(mbc))
names(mbc)[c(1,3:5)] <- c('experiment', 'ecosystem', 'replicate', 'treatment')

# fix MBC columns
mbc$treatment <- factor(mbc$treatment, levels=unique(mbc$treatment), labels=c('18O', 'C', 'CN'))
mbc$experiment <- factor(tolower(mbc$experiment))

# limit to high substrate experiment
mbc <- mbc[mbc$experiment=='high substrate',]
mbc$experiment <- NULL

# limit to pre-incubation data
mbc_sub <- mbc[mbc$week==0,]

mbc_dna <- merge(mbc_sub, abund_sub, by=c('ecosystem', 'replicate'))

mbc_dna_lm <- lm(mbc_g_g_soil ~ dna_ug_g_ds + 0, mbc_dna)
mbc_dna_lm$coefficients # linear slope is 14.95, slope in log-space is 90.61

# plot
ggplot(mbc_dna, aes(dna_ug_g_ds, mbc_g_g_soil)) +
  geom_point() +
  geom_smooth(method='lm', formula=y ~ x + 0)


# Calculate CUE for 18O soils------------------------------------------------------------

# fix CO2 data
co2$week <- NULL
co2$treatment <- factor(co2$treatment,
                        levels=levels(co2$treatment), 
                        labels=c('18O' ,'C', 'CN', 'C', 'CN'))
co2$ecosystem <- factor(co2$ecosystem, levels=c('MC', 'PP', 'PJ', 'GL'))
# limit to high substrate experiment
co2 <- co2[co2$experiment=='high substrate',]
co2$experiment <- NULL
names(co2)[2] <- 'replicate'
# get sampleID from all_data
co2 <- merge(co2, 
             unique(all_data[,c('sampleID', 'ecosystem', 'treatment', 'replicate')]), 
             by=c('ecosystem', 'replicate', 'treatment'),
             all.x=T)

# combine data
cue_calc <- Reduce(function(x, y) merge(x, y, by='sampleID', all.x=T), 
                   list(enrich, 
                        abund[,c('sampleID', 'dna_ug_g_ds')], 
                        co2[,c('sampleID',  'ecosystem', 'treatment', 'replicate', 'g_co2_c_g_soil_week')]))

# calculate new MBC produced from 18O enrichment, then CUE
cue_calc <- within(cue_calc, {
  mbc_18O <- enrichment * dna_ug_g_ds * mbc_dna_lm$coefficients  # linear transformation
  cue <- mbc_18O / (mbc_18O + g_co2_c_g_soil_week)
})

# plot CUE
ggplot(cue_calc, aes(treatment, cue, fill=ecosystem)) +
  geom_boxplot() +
  ylim(0, 1) +
  ylab('CUE calculated from 18O enrichment')

# much smaller CUE in C, C+N plots
range(cue_calc[cue_calc$treatment=='18O', 'cue'])
range(cue_calc[cue_calc$treatment!='18O', 'cue'])
# range in C, C+N treatments is: 0.03 - 0.13


# Compare CUE here with CUE from original estimate---------------------------------------

# split MBC by replicate, calculate 5 week differences in MBC, divide by 5 to get 1 week estimate of MBC FLUX
mbc <- split(mbc, mbc$replicate)
mbc <- lapply(mbc, function(x) {
  y <- x[x$week==5,]
  y$mbc_g_g_soil_week <- (x[x$week==5,'mbc_g_g_soil'] - rep(x[x$week==0,'mbc_g_g_soil'], each=3)) / 5
  y} )
mbc <- do.call(rbind, mbc)

# merge net MBC change with CUE data
cue_calc <- merge(cue_calc, 
                  mbc[,c('ecosystem', 'treatment', 'replicate', 'mbc_g_g_soil_week')],
                  by=c('ecosystem', 'treatment', 'replicate'),
                  all.x=T)
# merge initial MBC pool with CUE data
cue_calc <- merge(cue_calc, 
                  mbc_sub[,c('ecosystem', 'replicate', 'mbc_g_g_soil')],
                  by=c('ecosystem', 'replicate'),
                  all.x=T)
names(cue_calc)[grep('^mbc_g_g_soil$', names(cue_calc))] <- 'initial_mbc'

# calculate CUE from net MBC change
# use the difference between total growth (18O-MBC) and net change as a turnover rate
# the formulation depends on the defintion of turnover; we want turnover *independent* of labeling and growth
# other formulations of just 18O_MBC / init_MBC miss the effect of the loss of unlabeled biomass which allows greater EAF of MBC
cue_calc <- within(cue_calc, {
  cue_from_net_mbc <- mbc_g_g_soil_week / (mbc_g_g_soil_week + g_co2_c_g_soil_week)
  mbc_turnover <- (mbc_18O - mbc_g_g_soil_week) / initial_mbc
  # mbc_turnover <- mbc_18O / initial_mbc  # the traditional definition of turnover (= biomass-specific growth)
  mbc_turnover[mbc_turnover < 0] <- 0
})


# Compare turnover here with turnover calculated by Xiao-Jun-----------------------------

turn <- readxl::read_excel('../Data/Priming_Dimensions_Database_Archive.xlsx',
                           sheet='turnover',
                           range='B4:G68',
                           na='N/A')
turn <- data.frame(turn)
names(turn) <- c('ecosystem', 'treatment', 'experiment', 'CN_treatment', 'replicate', 'bulk_turnover')

turn <- turn[turn$experiment=='high',]
turn <- within(turn, {
  treatment <- sub('L|H', '', treatment)
  CN_treatment <- NULL
  experiment <- NULL
})

cue_calc <- merge(cue_calc, turn, all.x=T)

cue_calc$replicate <- factor(cue_calc$replicate)


# save CUE data calculated from 18O-enrichment
save(cue_calc, file='../Data/cue_18O_method.RData')