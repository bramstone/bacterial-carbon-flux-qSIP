# load data
load('data/taxonomic_data.RData')
source('figure_creation_functions.R')

tax <- as.matrix(tax)
tax[is.na(tax)] <- 'unknown'
tax <- as.data.frame(tax)

# use resp_data but combine with taxonomic info
ord_data <- merge(resp_data[,names(resp_data) %in% c('taxonID', 'ecosystem', 'treatment', 
                                                     'replicate', 'sampleID', 'activity',
                                                     'pop_c_ug', 'qsip_maint_resp_c_ug', 'net_status',
                                                     'qsip_growth_resp_c_ug', 'qsip_resp_c_ug')], 
                  tax, 
                  all.x=T)

names(ord_data)[names(ord_data) %in% 'pop_c_ug'] <- 'growth_c_ug'


# load in PCoA axes and perMANOVA tests
source('9_ordination_permanova_tests.R')

# process PCoA centroid coordinates together
pcoa <- list(bray_abund_cent$centroids[,1:2], bray_c_cent$centroids[,1:2])
pcoa <- lapply(pcoa, as.data.frame)
pcoa <- lapply(pcoa, function(x) {x$site_trt <- rownames(x); x})
pcoa <- Map(function(x, y) {x$type <- y; x}, pcoa, c('abund', 'c_use'))
pcoa <- do.call(rbind, pcoa)

#get treatment and site info
pcoa <- within(pcoa, {
  treatment <- sub('^(\\w{2})\\.(.*)', '\\2', site_trt)
  ecosystem <- sub('^(\\w{2})\\.(.*)', '\\1', site_trt)
  treatment <- factor(treatment, levels=c('18O', 'C', 'CN'))
  ecosystem <- factor(ecosystem, levels=c('MC', 'PP', 'PJ', 'GL'))
  site_trt <- NULL
})

pcoa_abund_df <- pcoa[pcoa$type=='abund',]
pcoa_c_df <- pcoa[pcoa$type=='c_use',]

# function to get ellipse data comparable to vegan's ordiellipse()
vegan_ellipse <- function(cov, center=c(0,0), conf=1, npoints=100) {
  if (conf > 1 | conf < 0) stop('confidence interval must be between 0 and 1')
  if(conf!=1) {
    conf <- sqrt(qchisq(conf, 2))
  } else {
    conf <- 1
  }
  theta <- (0:npoints) * 2 * (pi / npoints)
  Circle <- cbind(cos(theta), sin(theta))
  t(center + conf * t(Circle %*% chol(cov)))
}

# function to use in ggplot that applies vegan_ellipse to each group level
# set to standard error and 95% confidence interval
#
# CHANGE group='' ARGUMENT TO MAKE CIRCLES FORM AROUND DIFFERENT GROUPING FACTOR
#
ellipse_data <- function(x, group='ecosystem', conf=.95, kind='se') {
  kind <- match.arg(kind)
  if(!is.factor(x[,group])) x[,group] <- factor(x[,group])
  df <- data.frame()
  for(i in levels(x[,group])) {
    y <- x[x[,group]==i,]
    mat <- cov.wt(cbind(y$PCoA1, y$PCoA2), wt = rep(1/nrow(y), nrow(y)))
    if(kind=='se') mat$cov <- mat$cov * sum(mat$wt^2)
    # build data frame of ellipse coordinates using the function above
    df <- rbind(df, data.frame(vegan_ellipse(mat$cov,
                                             center=c(mean(y$PCoA1), mean(y$PCoA2)),
                                             conf=conf), 
                               group=i))
  }
  df$group <- factor(df$group, levels=levels(x[,group]))
  names(df)[1:2] <- c('PCoA1', 'PCoA2')
  names(df)[which(names(df)=='group')] <- group
  return(df)
}

# plot
pcoa_abund_plot <- ggplot(pcoa_abund_df, aes(PCoA1, PCoA2)) +
  geom_polygon(aes(x=PCoA1, y=PCoA2, group=ecosystem), color=gray(.65), fill=NA, size=1, inherit.aes=F, data=ellipse_data, show.legend=F) +
  geom_point(aes(color=treatment, fill=treatment, shape=ecosystem)) +
  xlab(paste0('Relative abundance PCoA 1 (', round(bray_abund_cent$eig[1] / sum(bray_abund_cent$eig), 2) * 100, '%)')) +
  ylab(paste0('Relative abundance PCoA 2 (', round(bray_abund_cent$eig[2] / sum(bray_abund_cent$eig), 2) * 100, '%)')) +
  scale_color_manual(values=c(trt_cols, rep(gray(.75), 4)), labels=c('Control', 'C', 'C+N')) +
  scale_fill_manual(values=trt_cols, labels=c('Control', 'C', 'C+N')) +
  scale_shape_manual(values=c(24, 22, 21, 25)) +
  labs(tag='a') +
  theme(legend.position='none')

pcoa_c_plot <- ggplot(pcoa_c_df, aes(PCoA1, PCoA2)) +
  geom_polygon(aes(x=PCoA1, y=PCoA2, group=ecosystem), color=gray(.65), fill=NA, size=1, inherit.aes=F, data=ellipse_data, show.legend=F) +
  geom_point(aes(color=treatment, fill=treatment, shape=ecosystem)) +
  xlab(paste0('Relativized carbon use PCoA 1 (', round(bray_c_cent$eig[1] / sum(bray_c_cent$eig), 2) * 100, '%)')) +
  ylab(paste0('Relativized carbon use PCoA 2 (', round(bray_c_cent$eig[2] / sum(bray_c_cent$eig), 2) * 100, '%)')) +
  scale_color_manual(values=c(trt_cols, rep(gray(.75), 4)), labels=c('Control', 'C', 'C+N')) +
  scale_fill_manual(values=trt_cols, labels=c('Control', 'C', 'C+N')) +
  scale_shape_manual(values=c(24, 22, 21, 25)) +
  guides(shape=guide_legend(override.aes=list(fill='black')), fill=F) +
  labs(tag='b') +
  theme(legend.position=c(0, 0), 
        legend.justification=c(0, 0),
        legend.spacing.x=unit(.001, units='npc'),
        legend.text=element_text(hjust=0),
        legend.box='horizontal')

grid.arrange(heights=rep(1, 2), pcoa_abund_plot, pcoa_c_plot)