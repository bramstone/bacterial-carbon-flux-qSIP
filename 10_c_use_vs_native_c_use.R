# note: results don't appreciably change if we keep 18O soils in or not

load('data/feature_table.RData')
load('data/taxonomic_data.RData')
load('data/absolute_growth.RData')
# load('data/c_flux_categories.RData')
load('data/respiration_18O_cue_draft_4.RData')
load('data/enrichment_13c_silva_per_sample.RData')

# for 13C data, add in replicate numbers for exact replicate matching
field_replicates <- read.table('data/old_map_merged.txt', header=T, check.names=F, comment.char='')
#
field_replicates <- field_replicates[,c('Description', 'Replicate', 'Soil', 'Treatment')]
names(field_replicates) <- c('sampleID', 'replicate', 'ecosystem', 'treatment')
field_replicates <- field_replicates[grep('w0|W1_', field_replicates$sampleID),]
field_replicates$sampleID <- gsub('-', '_', toupper(field_replicates$sampleID))


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
seq_data$sampleID <- NULL

# # average relative abundance across replicates
# seq_data <- aggregate(rel_abund ~ ., 
#                       seq_data[,!names(seq_data) %in% c('sampleID', 'replicate', 'abund')], 
#                       mean)

# average 13C-EAF per genus................................

eaf_13c <- reshape(eaf_13c,
                   idvar='taxonID',
                   varying=list(names(eaf_13c)[-1]),
                   times=names(eaf_13c)[-1],
                   timevar='sampleID',
                   v.names='eaf_13c',
                   direction='long')
rownames(eaf_13c) <- NULL

eaf_13c <- eaf_13c[!is.na(eaf_13c$eaf_13c),]

# OK, need ecosystem, treatment, and replicate info
eaf_13c <- merge(eaf_13c, field_replicates, all.x=T)

# refactor treatment to match 18O data
eaf_13c <- within(eaf_13c, {
  treatment <- factor(treatment, labels=c('C', 'CN'))
})

# combine relative abundances with C use data-----------------------------------

# combine productivity with respiration to get total C use
c_use <- within(resp_data, {
  c_flux_ug_c <- pop_c_ug + qsip_resp_c_ug
  pop_c_ug <- NULL 
  qsip_resp_c_ug <- NULL
  sampleID <- NULL
})

c_use <- c_use[c_use$treatment!='18O',]
c_use$treatment <- droplevels(c_use$treatment)

# combine with tax data
c_use <- merge(c_use[,!names(c_use) %in% c('Kingdom', 'Family')], tax, all.x=T)

# combine with 13 enrichment (ignore sample-ID names since they won't match)
c_use <- merge(c_use, eaf_13c[,!names(eaf_13c) %in% 'sampleID'], all.x=T)

# average C use and 13C enrichment across genera
c_use <- aggregate(cbind(c_flux_ug_c, eaf_13c) ~ ecosystem + 
                     treatment + replicate + Kingdom + Phylum + Class + Order +
                     Family + Genus, c_use, mean, na.rm=T)

# taxa without detectable 13C enrichment are assumed to have 0
c_use$eaf_13c[is.na(c_use$eaf_13c)] <- 0

# relativice total and NATIVE C use per treatment and ecosystem
c_use <- split(c_use, interaction(c_use$treatment, c_use$ecosystem))
c_use <- lapply(c_use, function(x) {
  native_c_use <- x$c_flux_ug_c * (1 - x$eaf_13c)
  glucose_c_use <- x$c_flux_ug_c * x$eaf_13c
  #
  x$c_flux_rel <- x$c_flux_ug_c / sum(x$c_flux_ug_c)
  x$nat_c_flux_rel <- native_c_use / sum(native_c_use, na.rm=T)
  x$glucose_c_flux_rel <- glucose_c_use / sum(glucose_c_use)
  x$glucose_c_flux_rel[is.nan(x$glucose_c_flux_rel) | x$glucose_c_flux_rel <= 0] <- NA
  x})

c_use <- do.call(rbind, c_use)

# combine 
abund_c_use <- merge(c_use, seq_data, all=F)
abund_c_use <- droplevels(abund_c_use)
nlevels(abund_c_use$Genus)
# there are 95 genera with a C contribution

# # function to refactor taxonomic levels
# refactor_tax <- function(x, new_levels=c(), phylum=NULL) {
#   y <- as.character(x)
#   y[is.na(y)] <- 'unknown'
#   # if phylum is null, only use new_levels to refactor
#   if(is.null(phylum)) {
#     y[!y %in% new_levels] <- 'Other'
#     y <- factor(y, levels=c(new_levels, 'Other'))
#     return(y)
#   }
#   # if phylum is provided, label anything in the "Other" phylum as "Other
#   if(!is.null(phylum)) {
#     y[phylum=='Other'] <- 'Other'
#     y[grep('metagenome|unknown', y)] <- paste0(phylum[grep('metagenome|unknown', y)], y[grep('metagenome|unknown', y)])
#     y <- factor(y, levels=c(unique(y)[unique(y)!='Other'], 'Other'))
#     return(y)
#   }
# }
# 
# top_phyla <- c('Actinobacteria', 'Proteobacteria', 'Acidobacteria',
#                'Gemmatimonadetes', 'Planctomycetes', 'Bacteroidetes')
# 
# # relevel phyla (and down) not in the top C users to "Other"
# c_use <- within(c_use, {
#   Phylum <- refactor_tax(Phylum, levels(top_phyla))
#   Class <- refactor_tax(Class, phylum=Phylum)
#   Order <- refactor_tax(Order, phylum=Phylum)
#   Family <- refactor_tax(Family, phylum=Phylum)
#   Genus <- refactor_tax(Genus, phylum=Phylum)
# })


# determine how many genera to look at.............
# look at cumulative abundance 
abund_c_use <- abund_c_use[order(abund_c_use$rel_abund, decreasing=T),]
abund_c_use <- split(abund_c_use, interaction(abund_c_use$treatment, abund_c_use$ecosystem))
abund_c_use <- lapply(abund_c_use, function(x) {x$cumul <- cumsum(x$rel_abund); x})
sapply(abund_c_use, function(x) min(which(x$cumul > .5))) # top 1-7 to get to 50% rel abundance


# look at cumulative C use
abund_c_use <- lapply(abund_c_use, function(x) x[order(x$c_flux_rel, decreasing=T),])
abund_c_use <- lapply(abund_c_use, function(x) {x$cumul <- cumsum(x$c_flux_rel); x})
sapply(abund_c_use, function(x) min(which(x$cumul > .5))) # only need 1-3 to get 50% C use
abund_c_use <- do.call(rbind, abund_c_use)

# full (95 genera)
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
                       varying=list(c('c_flux_rel', 'rel_abund', 
                                      'nat_c_flux_rel', 'glucose_c_flux_rel')),
                       times=c('Carbon use', '16S rRNA gene abundance',
                               'Native carbon use', 'Glucose carbon use'),
                       drop='eaf_13c',
                       timevar='metric',
                       v.names='relative_contribution',
                       direction='long')
row.names(abund_c_use) <- NULL

nat_c_use <- abund_c_use
nat_c_use_full <- abund_c_use_full

# save(nat_c_use, nat_c_use_full, file='data/rel_abund_vs_native_c_use.RData')


# statistical tests-------------------------------------------------------------

# Question: Is there more variation between total C use and native C use in the CN treatment?
# Said another way, how does nutrient addition affect variation about the trend line?

# create log10-transformed values for linear model fitting
test_data <- within(abund_c_use_full, {
  log_c <- log10(c_flux_rel)
  log_nat_c <- log10(nat_c_flux_rel)
  log_gluc <- log10(glucose_c_flux_rel)
  log_abund <- log10(rel_abund)
})

# test_data <- test_data[test_data$treatment!='18O',]
# test_data$treatment <- droplevels(test_data$treatment)

# calculate residuals from trend-line (which is essentially the 1:1 line)
trend_lm <- lm(log_nat_c ~ log_gluc, test_data)
data_resid <- rep(NA, nrow(test_data))
data_resid[!is.na(test_data$glucose_c_flux_rel)] <- residuals(trend_lm)
test_data$residuals <- data_resid

# difference in residuals across treatments?
test_data$eco_rep <- paste(test_data$ecosystem, 
                           test_data$replicate, 
                           sep='_')
var_lm <- lm(residuals ~ eco_rep:treatment, test_data)
car::leveneTest(var_lm)
# F[22, 770] = 3.53  p < 0.001

# correlation between natural C and glucose C?
cor.test(test_data$log_nat_c, test_data$log_gluc, use='pairwise.complete.obs')

# get variance in each treatment by summing squared residuals
aggregate(residuals ~ treatment, test_data, function(x) sum(x^2))

# SO SIGNIFICANTLY MORE VARIATION AROUND THE 1:1 LINE IN C+N SOILS
# I.E., IN C SOILS, TOTAL C UTILIZATION AND NATIVE C UTILIZATION WERE GENERALLY IN PROPORTION
# BUT IN C+N SOILS, NATIVE C UTILIZATION DIFFERS MORE FROM TOTAL C UTILIZATION, SO AN ORGANISM'S
# USE OF NATIVE SOIL C DOES NOT AS EASILY FOLLOW ITS TOTAL USE