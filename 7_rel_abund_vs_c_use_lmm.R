# linear mixed model of relative abundance vs. relative C use
# Question: What is the relationship between C use and abundance, and how does treatment affect this?
#   taking into account variance due to taxonomy and ecosystem, if need be

library(lme4)
library(lmerTest)
library(ggplot2)

load('data/rel_abund_vs_c_use_by_treatment_and_ecosystem_draft_4.RData')

# get data for plotting
data <- abund_c_use_full
data <- within(data, {
  log_c <- log10(c_flux_rel)
  log_abund <- log10(rel_abund)
})

# first test of linear model to see if it fits asssumptions
lm_1 <- lm(log_c ~ log_abund*treatment, data)
plot(lm_1, which=1)
# Looks like homogeneity isn't an issue! Will want to remove pseudo-replication


# 1. Find optimal random structure....................................

lm_r1 <- lmer(log_c ~ log_abund*treatment + (1|Genus),
              data=data, REML=F)

lm_r2 <- lmer(log_c ~ log_abund*treatment + (log_abund|Genus),
              data=data, REML=F) # failed to converge

lm_r3 <- lmer(log_c ~ log_abund*treatment + (1|ecosystem),
              data=data, REML=F)

lm_r4 <- lmer(log_c ~ log_abund*treatment + (log_abund|ecosystem),
              data=data, REML=F)

lm_r5 <- lmer(log_c ~ log_abund*treatment + (1|ecosystem) + (1|Genus),
              data=data, REML=F)

lm_r6 <- lmer(log_c ~ log_abund*treatment + (1|ecosystem) + (1|Phylum/Genus),
              data=data, REML=F) # failed to converge

lm_r7 <- lmer(log_c ~ log_abund*treatment + (log_abund|ecosystem) + (log_abund|Genus),
              data=data, REML=F) # failed to converge

lm_r8 <- lmer(log_c ~ log_abund*treatment + (log_abund|ecosystem) + (log_abund|Phylum/Genus),
              data=data, REML=F) # failed to converge

# compare/test, keeping lm objet at the end
anova(lm_r1, lm_r3, lm_r4, lm_r5, lm_1)
rand_aic <- AIC(lm_r1, lm_r3, lm_r4, lm_r5, lm_1)
rand_aic[order(rand_aic$AIC),]
# lowest AIC score is lm_r2 and lm_r8 (and r5)
# these models were significantly better than initial lm
# r2 includes random slopes/intercepts for Genus,
# r8 includes random slopes & intercepts for ecosystem and genus,
# r5 includes random intercepts for ecosystem and genus

anova(lm_r2, lm_r5) # r2 better
anova(lm_r2, lm_r8) # no diff.
anova(lm_r5, lm_r8) # r8 better
# there's no sig. difference between r2 and r8. These are best
# HOWEVER, model fails to converge with REML=T unless we use r5

# 2. Validate optimal variance structure..............................

lmm_full <- lmer(log_c ~ log_abund*treatment + 
                   (1|ecosystem) + 
                   (1|Genus),
                 data=data, REML=T)

# plot residuals
resid_data <- cbind(model.frame(lmm_full), 
                    resid=resid(lmm_full, scaled=T),
                    fitted=fitted(lmm_full))

# standardized (scaled, or normalized) residuals vs fitted
ggplot(resid_data, aes(fitted, resid)) + 
  geom_point(pch=1) +
  geom_hline(yintercept=0, color='red') +
  xlab('Fitted values') +
  ylab('Normalized residuals')
# generally pretty good, having treatment as fixed effect removes some patterns in residuals

# residuals against relative abundance
ggplot(resid_data, aes(log_abund, resid)) + 
  geom_point(pch=1) +
  geom_hline(yintercept=0, color='red') +
  xlab(expression('log'[10]*'(Relative abundance)')) +
  ylab('Normalized residuals')
# looks good


# 1. Find optimal fixed structure.....................................

# specify full formula
form_full <- formula(log_c ~ log_abund + treatment + log_abund:treatment +                    
                       (1|ecosystem) + 
                       (1|Genus))

# remove interaction
lm_f1 <- lmer(update(form_full, . ~ . - log_abund:treatment), data=data, REML=F)
anova(lmm_full, lm_f1)
# keep pop:treatment interaction

# can't remove any other fixed variables

set.seed(11020)
lmm_optim <- lmer(log_c ~ log_abund*treatment + 
                    (1|ecosystem) + 
                    (1|Genus),
                  data=data, REML=T)

# IF WE WANT TO TEST IF SLOPES ARE SIGNIFICANTLY DIFFERENT THAN 1, USE THE OFFSET FUNCTION
lmm_optim_3 <- lmer(log_c ~ log_abund*treatment + offset(log_abund) +
                      (1|ecosystem) + 
                      (1|Genus),
                    data=data, REML=T)

# get final model slopes and intercepts
summary(lmm_optim)
anova(lmm_optim)

# slopes significantly different than 1?
summary(lmm_optim_3)
# 18O and C are not, CN is significantly higher

save(lmm_optim, lmm_optim_3, file='../Data/lmm_with_ecosystem_draft_4.RData')
