# Figure 2 creation

load('data/growth_ug_c_contributions_per_asv_draft_4.RData')        # biomass growth
load('data/resp_growth_ug_c_contributions_per_asv_draft_4.RData')   # growth-respiration

# keep only ones summarized at the family level
rm(list=ls(pattern='2|4'))

# change "Abundance" column names for merging
names(growth_data5)[grep('Abund', names(growth_data5))] <- 'production'
names(resp_growth_data5)[grep('Abund', names(resp_growth_data5))] <- 'resp_growth'

# merge just growth and growth-related respiration
bar_data <- Reduce(function(x,y) merge(x,y, all.x=T), list(resp_growth_data5, growth_data5))

# make total C use from growth and resp
bar_data$total_c <- apply(bar_data[,c('production', 'resp_growth')], 1, sum, na.rm=T)
bar_data <- bar_data[bar_data$total_c > 0,]

# average per taxon contributions to each ecosystem:treatment combination
bar_data <- aggregate(total_c ~ ., 
                      data=bar_data[,!names(bar_data) %in% c('replicate', 'Sample', 'groupID', 'Species', 'Genus')], 
                      mean, na.rm=T)

# keep only select phyla
phyla_keep <- c('Actinobacteria', 'Proteobacteria', 'Acidobacteria',
                'Verrucomicrobia', 'Firmicutes', 'Chloroflexi',
                'Gemmatimonadetes', 'Bacteroidetes', 'Planctomycetes')

bar_data <- bar_data[bar_data$Phylum %in% phyla_keep,]
bar_data$Phylum <- factor(bar_data$Phylum, levels=phyla_keep)

# re-order family factor by dominant phyla, secondly by growth
bar_data <- bar_data[!is.na(bar_data$total_c),]
bar_data <- bar_data[order(as.numeric(bar_data$Phylum), bar_data$total_c, method='radix', decreasing=c(F, T)),]
family_levels <- as.character(unique(bar_data$Family))
family_levels <- family_levels[!grepl('Other', family_levels)]
family_levels <- c(family_levels, 'Other')
bar_data$Family <- factor(bar_data$Family, levels=family_levels)

# Total C use
ggplot(bar_data, aes(x=treatment, y=total_c / 3, fill=Family)) +
  geom_bar(aes(), stat="identity", position=position_stack(reverse=T), width=.5) +
  scale_fill_manual(values=rev(tax_color(bar_data$Family, 103119))[c(1:15,17)]) +
  ylab(expression('Total C use (('*mu*'g C g soil '^{'-1'}*' wk '^{'-1'}*')')) +
  xlab('') +
  scale_x_discrete(labels=NULL) +
  facet_grid(. ~ ecosystem) +
  labs(tag='a') +
  theme(legend.position='none',
        plot.margin=unit(c(5.5, 5.5, 0, 5.5), units='pt'))

# relativize total C use
bar_data_rel<- split(bar_data, interaction(bar_data$ecosystem, bar_data$treatment))
bar_data_rel <- lapply(bar_data_rel, function(x) {x$total_c <- x$total_c / sum(x$total_c); x})
bar_data_rel <- do.call(rbind, bar_data_rel)
rownames(bar_data_rel) <- NULL

# Relative C use
ggplot(bar_data_rel, aes(x=treatment, y=total_c, fill=Family)) +
  geom_bar(aes(), stat="identity", position=position_stack(reverse=T), width=.5) +
  scale_fill_manual(values=rev(tax_color(bar_data$Family, 103119))[c(1:15,17)]) +
  ylab(expression('Relative C use (g soil '^{'-1'}*' wk '^{'-1'}*')')) +
  xlab('') +
  scale_x_discrete(labels=expression(''^18*O, 'C', 'C+N')) +
  facet_grid(. ~ ecosystem) +
  labs(tag='b') +
  guides(fill=guide_legend(ncol=3)) +
  theme(strip.text=element_blank(),
        legend.spacing.y=unit(.25, 'char'), 
        legend.key.width=unit(1, 'char'), 
        legend.position='bottom',
        legend.margin=unit(0, 'cm'),
        axis.text.x=element_text(vjust=0),
        plot.margin=unit(c(0, 5.5, 5.5, 4.5), units='pt'))