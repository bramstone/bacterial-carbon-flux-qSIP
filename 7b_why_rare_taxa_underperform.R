# Why do rare taxa underperform in terms of C use in C+N treatment?

library(ggplot2)
library(lme4)
library(lmerTest)

load('data/rel_abund_vs_c_use_by_treatment_and_ecosystem_draft_4.RData')
load('data/lmm_with_ecosystem_draft_4.RData')
load('data/taxonomic_data.RData')
load('data/absolute_growth.RData')


# get residuals and data from linear mixed model with 1:1 offest
data <- lmm_optim_3@frame
data$log_c_resid <- residuals(lmm_optim_3)

# now combine residual difference from 1:1 lines with growth rate, CUE, cell mass, 16S count
# first aggregate those data to the genus level
all_data <- merge(all_data, tax, all.x=T)
keep <- c('ecosystem', 'treatment', 'birth_rate', 
          'death_rate', 'copy_number', 'cell_mass_g',
          'Phylum', 'Class', 'Order', 'Family', 'Genus')
all_data <- aggregate(cbind(cell_mass_g, copy_number, birth_rate, death_rate) ~ ., data=all_data[,keep], mean)

data <- merge(data, all_data, all.x=T)

# check for multicollinearity
cor(data[,c('log_c_resid', 'cell_mass_g', 'copy_number', 'birth_rate')])
# probably the highest may be between copy number and cell mass
# which to remove? copy number correlates less strongly to residual C use, so remove that

lm_data <- data[,c('log_c_resid', 'cell_mass_g', 'birth_rate')]

# compare models with and without primary interactions
AIC(lm(log_c_resid ~ .^2, lm_data),
    lm(log_c_resid ~ ., lm_data))
# marginally better without the interactions
# if using alternative data, then it is better with the interactions

# for alternative data: any interactions to be deleted?
drop1(lm(log_c_resid ~ .^2, lm_data))
# better to keep interactions

# any single terms that can be deleted?
drop1(lm(log_c_resid ~ ., lm_data))
# cell mass could be removed
# none for alternative data

# so what is a significant factor in lower C use values by rare taxa?
summary(lm(log_c_resid ~ ., lm_data))

car::Anova(lm(log_c_resid ~ ., lm_data))


# plot residuals by birth rate for each treatment
ggplot(data, aes(birth_rate, log_c_resid)) +
  geom_hline(yintercept=0, color=gray(.5)) +
  geom_point() +
  xlab('per-capita birth rate') +
  ylab('Residual variation in relativized C use ~ relative abundance') +
  facet_grid(treatment ~ .)