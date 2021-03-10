# creation of Figure 1

# created using respiration_contribution_18O_cue.R
load('data/respiration_18O_cue_draft_4.RData')

trt_cols <- hsv(c(0,.6,.1), c(0,1,1), c(.5,.8,.8))

# get summary for B model (w/o taxonomic distinction)
resp_summary_b_se <- aggregate(cbind(co2_per_tax, g_co2_c_g_soil_week) ~ ecosystem + treatment, resp_data_b, se)
resp_summary_b <- aggregate(cbind(co2_per_tax, g_co2_c_g_soil_week) ~ ecosystem + treatment, resp_data_b, mean)

resp_summary_b <- merge(resp_summary_b, resp_summary_b_se, 
                        by=c('ecosystem', 'treatment'),
                        suffixes=c('', '_se'))

# get slope and R-squared
co2_resp_lm <- lm(qsip_resp_c_ug ~ g_co2_c_g_soil_week, resp_summary)
rsq_resp <- bquote(italic('R')^2 == .(round(summary(co2_resp_lm)$r.squared, 2)))
#
co2_resp_b_lm <- lm(co2_per_tax ~ g_co2_c_g_soil_week, resp_summary_b)
rsq_resp_b <- bquote(italic('R')^2 == .(round(summary(co2_resp_b_lm)$r.squared, 2)))

# plot respiration from total growth
ggplot(resp_summary, aes(g_co2_c_g_soil_week, qsip_resp_c_ug, fill=treatment, color=treatment, shape=ecosystem)) +
  geom_abline(slope=1, intercept=0, color=gray(.75)) +
  geom_smooth(aes(shape=NULL, group=NULL, color=NULL, fill=NULL), method='lm', formula=y ~ x, se=F, color='black', show.legend=F) +
  geom_point(stroke=0, size=2.25) +
  ylim(0, 600) + # 1100 for draft-3 model
  xlim(0, 600) + # 1100 for draft-3 model
  geom_errorbar(aes(ymin=qsip_resp_c_ug - qsip_resp_c_ug_se,
                    ymax=qsip_resp_c_ug + qsip_resp_c_ug_se), width=0) +
  geom_errorbarh(aes(xmin=g_co2_c_g_soil_week - g_co2_c_g_soil_week_se,
                     xmax=g_co2_c_g_soil_week + g_co2_c_g_soil_week_se), height=0) +
  scale_color_manual(values=trt_cols, labels=expression(''^18*O, 'C', 'C+N')) +
  scale_fill_manual(values=trt_cols, labels=expression(''^18*O, 'C', 'C+N')) +
  scale_shape_manual(values=c(24, 22, 21, 25)) +
  annotate('text', x=0, y=Inf, label=as.expression(rsq_resp), hjust=0, vjust=2, size=3) +
  annotate('text', x=575, y=575, label='1:1', color=gray(.75), vjust=1, hjust=0) +
  guides(shape=guide_legend(override.aes=list(fill='black'))) +
  labs(tag='a') +
  theme(legend.position=c(1, 0), 
        legend.justification=c(1, 0),
        legend.box='horizontal',
        legend.text=element_text(hjust=0),
        axis.title=element_blank())


# plot respiration from taxnomically-insensitive model
ggplot(resp_summary_b, aes(g_co2_c_g_soil_week, co2_per_tax, fill=treatment, color=treatment, shape=ecosystem)) +
  geom_abline(slope=1, intercept=0, color=gray(.75)) +
  geom_smooth(aes(shape=NULL, group=NULL, color=NULL, fill=NULL), method='lm', formula=y ~ x, se=F, color='black', show.legend=F) +
  geom_point(stroke=0, size=2.25) +
  ylim(0, 600) + # 1100 for draft-3 model
  xlim(0, 600) + # 1100 for draft-3 model
  geom_errorbar(aes(ymin=co2_per_tax - co2_per_tax_se,
                    ymax=co2_per_tax + co2_per_tax_se), width=0) +
  geom_errorbarh(aes(xmin=g_co2_c_g_soil_week - g_co2_c_g_soil_week_se,
                     xmax=g_co2_c_g_soil_week + g_co2_c_g_soil_week_se), height=0) +
  scale_color_manual(values=trt_cols, labels=expression(''^18*O, 'C', 'C+N')) +
  scale_fill_manual(values=trt_cols, labels=expression(''^18*O, 'C', 'C+N')) +
  scale_shape_manual(values=c(24, 22, 21, 25)) +
  annotate('text', x=0, y=Inf, label=as.expression(rsq_resp_b), hjust=0, vjust=2, size=3) +
  annotate('text', x=575, y=575, label='1:1', color=gray(.75), vjust=1, hjust=0) +
  guides(shape=guide_legend(override.aes=list(fill='black'))) +
  labs(tag='b') +
  theme(legend.position='none',
        axis.title=element_blank())