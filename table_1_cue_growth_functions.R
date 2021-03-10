# compare across diversity of CUE ~ growth functions proposed

library(ggplot2)

load('data/respiration_18O_cue_draft_4.RData')  # ignore respiration values already calculated
load('data/co2_data_mbc.RData')
load('data/cue_18O_method.RData')
rm(resp_data_b, resp_summary, resp_summary_per_sample)

# standard error function
se <- function(x) sd(x, na.rm=T) / sqrt(length(x[!is.na(x)]))

# treatment colors for plotting
trt_cols <- hsv(c(0,.6,.1), c(0,1,1), c(.5,.8,.8))

resp_data$ecosystem <- factor(resp_data$ecosystem, levels=c('MC', 'PP', 'PJ', 'GL'))

# prepare CO2 data for comparison against qSIP-based respiration
co2 <- co2[co2$experiment=='high substrate',]
co2 <- co2[,!names(co2) %in% c('experiment', 'week')]
names(co2)[grep('replicate', names(co2))] <- 'replicate'
co2 <- within(co2, {
  treatment <- factor(treatment, levels=levels(treatment), labels=c('18O', 'C', 'CN', 'C', 'CN'))
  replicate <- factor(replicate)
})

# use 18O-MBC CUE as "ground-truth" of CUE for each treatment and ecosystem
cue_min <- aggregate(cue ~ treatment + ecosystem, cue_calc, min)
names(cue_min)[3] <- 'cue_min'
cue_range <- aggregate(cue ~ treatment + ecosystem, cue_calc, range)
cue_range$cue_range <- cue_range$cue[,2] - cue_range$cue[,1]
#
resp_data <- merge(resp_data, cue_min, all.x=T)
resp_data <- merge(resp_data, cue_range[,c('treatment', 'ecosystem', 'cue_range')], all.x=T)


# function which applies per-taxon CUE and to scale against some value----------

#   CO2 related values are: qsip_resp_c_ug, g_co2_c_g_soil_week
#   CUE related values are: cue_contrib, total_cue
#   scale=T/F indicates whether or not to z-transform the variables
ecosystem_scale <- function(x, y, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week', scale=T) {
  x <- within(x, {
    cue[birth_rate==0 & !is.na(birth_rate)] <- 0
    #
    Rg_rate <- (birth_rate / cue) - birth_rate
    Rm_rate <- Rg_rate / 1E2
    resp_rate <- Rg_rate + Rm_rate
    #
    qsip_resp_c_ug <- resp_rate * 7 * pop_c_ug
    qsip_growth_resp_c_ug <- Rg_rate * 7 * pop_c_ug
    qsip_maint_resp_c_ug <- Rm_rate * 7 * pop_c_ug
  })
  #
  x <- aggregate(get(x_val) ~ ecosystem + treatment + sampleID + replicate, x, sum)
  names(x)[grep('get\\(x_val\\)|V1', names(x))] <- x_val
  z <- merge(x, y, all.x=T)
  # z-transform the variables
  if(scale) {
    z[,x_val] <- scale(z[,x_val])
    z[,y_val] <- scale(z[,y_val])
  }
  #
  z_se <- aggregate(cbind(get(x_val), get(y_val)) ~ ecosystem + treatment, z, se)
  z <- aggregate(cbind(get(x_val), get(y_val)) ~ ecosystem + treatment, z, mean)
  z <- merge(z, z_se, 
             by=c('ecosystem', 'treatment'),
             suffixes=c('', '_se'))
  names(z) <- gsub('V1', x_val, names(z))
  names(z) <- gsub('V2', y_val, names(z))
  return(z)
}


#-------------------------------------------------------------------------------
# Unconstrained methods---------------------------------------------------------
#-------------------------------------------------------------------------------

# linear negative..........................

# # identify the shape of the function
# load('unconstrained_linear_neg_cue_growth.RData')
# best_cue$best_A <- -1 * best_cue$best_A
# best_cue[best_cue$slope_close_to_1=='close',]
# # use y-intercept of 0.14 and slope of -1.6

# apply shape to calculate CUE
resp_data <- within(resp_data, {
  cue <- -(1.6 * birth_rate) + 0.14
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
lin_neg <- merge(scale_co2, scale_cue)
#
u_lin_neg_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, lin_neg)
u_lin_neg_lm_cue <- lm(cue ~ cue_contrib, lin_neg)

ggplot(resp_data, aes(birth_rate, cue)) +
  geom_line() +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Unconstrained linear negative')

# linear positive..........................

# # identify the shape of the function
# load('unconstrained_linear_pos_cue_growth.RData')
# best_cue[best_cue$slope_close_to_1=='close',]
# # use y-intercept of 0.13 and slope of 0.1

# apply shape to calculate CUE
resp_data <- within(resp_data, {
  cue <- (.1 * birth_rate) + .13
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
lin_pos <- merge(scale_co2, scale_cue)
#
u_lin_pos_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, lin_pos)
u_lin_pos_lm_cue <- lm(cue ~ cue_contrib, lin_pos)

ggplot(resp_data, aes(birth_rate, cue)) +
  geom_line() +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Unconstrained linear positive')

# unimodal 0.5-centered....................

# # identify the shape of the function
# load('unconstrained_unimodal.RData')
# best_cue[best_cue$slope_close_to_1=='close',]
# # use lower-bound CUE of 0.13 and range of 0.01

# apply shape to calculate CUE
resp_data <- within(resp_data, {
  cue <- -(.01/.25) * (birth_rate - .5)^2 + (.01 + .13)
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
unimod <- merge(scale_co2, scale_cue)
#
u_unimod_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, unimod)
u_unimod_lm_cue <- lm(cue ~ cue_contrib, unimod)

ggplot(resp_data, aes(birth_rate, cue)) +
  geom_line() +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Unconstrained unimodal 0.5-centered')

# unimodal med. growth rate-centered.......

# # identify the shape of the function
# load('unconstrained_unimodal_med_growth.RData')
# best_cue[best_cue$slope_close_to_1=='close',]
# # no good choice: use lower bound CUE of 0.05, and range of 0.65

# apply shape to calculate CUE
resp_data <- within(resp_data, {
  cue <- dgamma(birth_rate*10, 1.052632, 1.052632) * .65 + .05
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
unimod_med <- merge(scale_co2, scale_cue)
#
u_unimod_med_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, unimod_med)
u_unimod_med_lm_cue <- lm(cue ~ cue_contrib, unimod_med)

ggplot(resp_data, aes(birth_rate, cue)) +
  geom_line() +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Unconstrained unimodal\ncentered on median growth rate')

# exponential decline......................

# # identify the shape of the function
# load('unconstrained_exponential_decline.RData')
# best_cue[best_cue$slope_close_to_1=='close',]
# # use max CUE of 0.21, exponential decrease rate of 14.7

# apply shape to calculate CUE
resp_data <- within(resp_data, {
  cue <- exp(-birth_rate * 14.7) * .21
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
expntl <- merge(scale_co2, scale_cue)
#
u_expntl_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, expntl)
u_expntl_lm_cue <- lm(cue ~ cue_contrib, expntl)

ggplot(resp_data, aes(birth_rate, cue)) +
  geom_line() +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Unconstrained exponential decline')


#-------------------------------------------------------------------------------
# Constrained functions---------------------------------------------------------
#-------------------------------------------------------------------------------

# linear negative relationship......................

# calculate CUE
resp_data <- within(resp_data, {
  cue <- -(cue_range * birth_rate) + (cue_min + cue_range)
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
lin_neg <- merge(scale_co2, scale_cue)
#
lin_neg_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, lin_neg)
lin_neg_lm_cue <- lm(cue ~ cue_contrib, lin_neg)

ggplot(resp_data, aes(birth_rate, cue, color=treatment)) +
  geom_line(aes(group=sampleID)) +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Constrained linear negative') +
  scale_color_manual(values=trt_cols) +
  theme(legend.position=c(1, 1), 
        legend.justification=c(1,1))


# linear positive relationship......................

# calculate CUE
resp_data <- within(resp_data, {
  cue <- (cue_range * birth_rate) + cue_min
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
lin_pos <- merge(scale_co2, scale_cue)
#
lin_pos_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, lin_pos)
lin_pos_lm_cue <- lm(cue ~ cue_contrib, lin_pos)

ggplot(resp_data, aes(birth_rate, cue, color=treatment)) +
  geom_line(aes(group=sampleID)) +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Constrained linear positive') +
  scale_color_manual(values=trt_cols) +
  theme(legend.position=c(1, 1), 
        legend.justification=c(1,1))


# unimodal centered on global median relationship...

# calculate CUE
resp_data <- within(resp_data, {
  cue <- dgamma(birth_rate*10, 1.052632, 1.052632) * cue_range + cue_min
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
unimod_med <- merge(scale_co2, scale_cue)
#
unimod_med_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, unimod_med)
unimod_med_lm_cue <- lm(cue ~ cue_contrib, unimod_med)

ggplot(resp_data, aes(birth_rate, cue, color=treatment)) +
  geom_line(aes(group=sampleID)) +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Constrained unimodal\ncentered on median growth rate') +
  scale_color_manual(values=trt_cols) +
  theme(legend.position=c(1, 1), 
        legend.justification=c(1,1))


# unimodal centered on 0.5..........................

# calculate CUE
resp_data <- within(resp_data, {
  cue <- -(cue_range/.25) * (birth_rate - .5)^2 + (cue_range + cue_min)
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
unimod <- merge(scale_co2, scale_cue)
#
unimod_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, unimod)
unimod_lm_cue <- lm(cue ~ cue_contrib, unimod)

ggplot(resp_data, aes(birth_rate, cue, color=treatment)) +
  geom_line(aes(group=sampleID)) +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Constrained unimodal 0.5-centered') +
  scale_color_manual(values=trt_cols) +
  theme(legend.position=c(1, 1), 
        legend.justification=c(1,1))


# exponential decline...............................

# calculate CUE
resp_data <- within(resp_data, {
  cue <- exp(-birth_rate * 1) * (cue_range + cue_min)
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
expntl <- merge(scale_co2, scale_cue)
#
expntl_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, expntl)
expntl_lm_cue <- lm(cue ~ cue_contrib, expntl)

ggplot(resp_data, aes(birth_rate, cue, color=treatment)) +
  geom_line(aes(group=sampleID)) +
  labs(x='qSIP-measured growth rates', y='estimated CUE',
       title='Constrained exponential decline') +
  scale_color_manual(values=trt_cols) +
  theme(legend.position=c(1, 1), 
        legend.justification=c(1,1))


# difference in growth from 18O.....................

# Change in growth rates across treatments
growth_change <- reshape(resp_data[c('taxonID', 'ecosystem', 'treatment', 
                                     'replicate', 'birth_rate')],
                         idvar=c('taxonID', 'ecosystem', 'replicate'),
                         timevar='treatment',
                         direction='wide')
#
growth_change <- reshape(growth_change,
                         idvar=c('taxonID', 'ecosystem', 'replicate'),
                         varying=list(c('birth_rate.C', 'birth_rate.CN')),
                         timevar='treatment',
                         v.names='birth_rate_trt',
                         times=c('C', 'CN'),
                         direction='long')
growth_change <- growth_change[!is.na(growth_change$birth_rate_trt),]
rownames(growth_change) <- NULL

# - Proportional change from this rate equates to change in CUE
growth_change <- within(growth_change, {
  prop_growth_change <- birth_rate_trt / birth_rate.18O
  prop_growth_change[is.infinite(prop_growth_change) | is.nan(prop_growth_change)] <- NA
})
#
resp_data$prop_growth_change <- NULL
resp_data <- merge(resp_data,
                   growth_change[,!names(growth_change) %in% c('birth_rate.18O', 'birth_rate_trt')],
                   all.x=T)

# calculate CUE
resp_data <- within(resp_data, {
  prop_growth_change[treatment=='18O'] <- 1
  cue <- .19 * .19^(abs(log(prop_growth_change))*1)
  cue[is.na(prop_growth_change)] <- .19
  cue_contrib <- rel_abund * cue
})

# scale against CO2 and CUE
scale_co2 <- ecosystem_scale(resp_data, co2, x_val='qsip_resp_c_ug', y_val='g_co2_c_g_soil_week')
scale_cue <- ecosystem_scale(resp_data, cue_calc, x_val='cue_contrib', y_val='cue')
diff_18O <- merge(scale_co2, scale_cue)
#
diff_18O_lm <- lm(g_co2_c_g_soil_week ~ qsip_resp_c_ug, diff_18O)
diff_18O_lm_cue <- lm(cue ~ cue_contrib, diff_18O)

ggplot(resp_data[resp_data$treatment!='18O',], 
       aes(prop_growth_change, cue, color=treatment)) +
  geom_line(aes(group=sampleID)) +
  labs(x='Proportional change in growth rate\nfrom18O to treatment', y='estimated CUE',
       title='Difference in growth rate\nfrom 18O to treatment') +
  scale_color_manual(values=trt_cols[2:3]) +
  theme(legend.position=c(1, 1), 
        legend.justification=c(1,1))


#-------------------------------------------------------------------------------
# Compare models----------------------------------------------------------------
#-------------------------------------------------------------------------------

# CO2 scaling..............................
aic_compare <- AIC(lin_neg_lm, lin_pos_lm,
                   unimod_med_lm, unimod_lm,
                   expntl_lm, diff_18O_lm,
                   u_lin_neg_lm, u_lin_pos_lm,
                   u_unimod_med_lm, u_unimod_lm,
                   u_expntl_lm)

# CUE scaling..................
aic_compare_cue <- AIC(lin_neg_lm_cue, lin_pos_lm_cue,
                       unimod_med_lm_cue, unimod_lm_cue,
                       expntl_lm_cue, diff_18O_lm_cue,
                       u_lin_neg_lm_cue, u_lin_pos_lm_cue,
                       u_unimod_med_lm_cue, u_unimod_lm_cue,
                       u_expntl_lm_cue)

#
aic_compare$AIC_cue <- aic_compare_cue$AIC
aic_compare$AIC_combn <- 2*(aic_compare$AIC) + aic_compare$AIC_cue

# calculate delta_AIC values and Akaike weights
delt_aic <- function(x) {
  best_aic <- min(x)
  x - best_aic
}

aic_weight <- function(x) {
  d <- delt_aic(x)
  exp(-.5 * d) / sum(exp(-.5 * d))
}

delt_aic_compare <- within(aic_compare, {
  AIC <- delt_aic(AIC)
  AIC_cue <- delt_aic(AIC_cue)
  AIC_combn <- 2*(AIC) + AIC_cue
})

aic_weight_compare <- within(aic_compare, {
  w <- aic_weight(AIC)
  w_cue <- aic_weight(AIC_cue)
  w_combn <- 2*(w) + w_cue
  AIC <- NULL
  AIC_cue <- NULL
  AIC_combn <- NULL
})

# report Delta-AIC values
delt_aic_compare$shape <- rownames(delt_aic_compare)
delt_aic_compare$df <- NULL

delt_aic_compare <- within(delt_aic_compare, {
  AIC <- round(AIC, 2)
  AIC_cue <- round(AIC_cue, 2)
  AIC_combn <- round(AIC_combn, 2)
  df <- NULL
})

delt_aic_compare <- delt_aic_compare[,c('shape', 'AIC', 'AIC_cue', 'AIC_combn')]
delt_aic_compare <- delt_aic_compare[order(delt_aic_compare$AIC_combn, decreasing=F),]

names(delt_aic_compare) <- sub('AIC', 'delt_AIC', names(delt_aic_compare))

# print out data - Table 1 in manuscript
delt_aic_compare