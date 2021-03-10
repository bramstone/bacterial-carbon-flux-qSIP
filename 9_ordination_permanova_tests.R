# multivariate ordinations and PERMANOVA tests

install.packages('vegan')

library(vegan)

load('data/feature_table.RData')
load('data/respiration_18O_cue_draft_4.RData')
load('data/taxonomic_data.RData')
load('data/cue_18O_method.RData')
tax <- as.matrix(tax)
tax[is.na(tax)] <- 'unknown'
tax <- as.data.frame(tax)


# initial data prep---------------------------------------------------------------------

# use resp_data but combine with taxonomic info
ord_data <- merge(resp_data[,names(resp_data) %in% c('taxonID', 'ecosystem', 'treatment', 
                                                            'replicate', 'sampleID', 'activity',
                                                            'pop_c_ug', 'qsip_maint_resp_c_ug', 'net_status',
                                                            'qsip_growth_resp_c_ug', 'qsip_resp_c_ug')], 
                         tax, 
                         all.x=T)

# since we already calculated maintenance respiration based on pop levels, turn dormant growth to 0
names(ord_data)[names(ord_data) %in% 'pop_c_ug'] <- 'growth_c_ug'


# convert growth and respiration into matrices
growth_matrix <- reshape(ord_data[,c('taxonID', 'growth_c_ug', 'sampleID')],
                         idvar='taxonID',
                         timevar='sampleID',
                         direction='wide')
rownames(growth_matrix) <- growth_matrix$taxonID
growth_matrix <- as.matrix(growth_matrix[,-1])
growth_matrix[is.na(growth_matrix)] <- 0
growth_matrix <- t(growth_matrix)
rownames(growth_matrix) <- sub('growth_c_ug\\.', '', rownames(growth_matrix))

resp_matrix <- reshape(aggregate(qsip_resp_c_ug ~ sampleID + taxonID, ord_data[,c('taxonID', 'qsip_resp_c_ug', 'sampleID')], sum),
                       idvar='taxonID',
                       timevar='sampleID',
                       direction='wide')
rownames(resp_matrix) <- resp_matrix$taxonID
resp_matrix <- as.matrix(resp_matrix[,-1])
resp_matrix[is.na(resp_matrix)] <- 0
resp_matrix <- t(resp_matrix)
rownames(resp_matrix) <- sub('qsip_resp_c_ug\\.', '', rownames(resp_matrix))

# MAKE AND PLOT TOTAL C USE MATRIX, RATHER THAN GROWTH AND RESPIRATION SEPARATELY
resp_matrix <- resp_matrix[,match(colnames(growth_matrix), colnames(resp_matrix))]
resp_matrix[is.na(resp_matrix)] <- 0
c_matrix <- resp_matrix + growth_matrix

# data to color plots by (and to set contrasts)
plot_data <- unique(ord_data[,c('ecosystem', 'treatment', 'replicate', 'sampleID')])
eco_contrasts <- contrasts(plot_data$ecosystem)
eco_contrasts[1:length(eco_contrasts)] <- contr.sum(4)
contrasts(plot_data$ecosystem) <- eco_contrasts
rm(eco_contrasts)

growth_matrix <- growth_matrix[match(plot_data$sampleID, rownames(growth_matrix)),]
resp_matrix <- resp_matrix[match(plot_data$sampleID, rownames(resp_matrix)),]
c_matrix <- c_matrix[match(plot_data$sampleID, rownames(c_matrix)),]
asv <- asv[match(plot_data$sampleID, rownames(asv)),]

# perform scaling transformation (so that differences won't just be to higher C flux in MC)
# hellinger transformation divides C flux of each taxa in the sample by the sample total (relativizes) and takes the square root
# thus, differences between rare taxa are maximized, perhaps too much, use simple relative abundances instead
abund_matrix <- decostand(as(asv, 'matrix'), 'total', 1)
growth_matrix <- decostand(growth_matrix, 'total', 1)
resp_matrix <- decostand(resp_matrix, 'total', 1)
c_matrix <- decostand(c_matrix, 'total', 1)
# hellinger transformations show that differences in composition driven largely by rare species


# Bray-Curtis dissimilarities and PERMANOVA-------------------------------------

# get dissimilarity measures
bray_abund <- vegdist(abund_matrix, method='bray')
bray_growth <- vegdist(growth_matrix, method='bray')
bray_resp <- vegdist(resp_matrix, method='bray')
bray_c <- vegdist(c_matrix, method='bray')

# get replicate centroids FIRST, THEN perform perMANOVAs and plot
bray_abund_cent <- betadisper(bray_abund, interaction(plot_data$ecosystem, plot_data$treatment), type='centroid')
bray_growth_cent <- betadisper(bray_growth, interaction(plot_data$ecosystem, plot_data$treatment), type='centroid')
bray_resp_cent <- betadisper(bray_resp, interaction(plot_data$ecosystem, plot_data$treatment), type='centroid')
bray_c_cent <- betadisper(bray_c, interaction(plot_data$ecosystem, plot_data$treatment), type='centroid')

centroid_data <- data.frame(ecosystem=sub('(.*)\\.(.*)', '\\1', rownames(bray_abund_cent$centroids)),
                            treatment=sub('(.*)\\.(.*)', '\\2', rownames(bray_abund_cent$centroids)))
centroid_data$ecosystem <- factor(centroid_data$ecosystem, levels=c('MC', 'PP', 'PJ', 'GL'))

# generate perMANOVAs
(perman_abund <- adonis(bray_abund_cent$centroids ~ ecosystem + treatment, centroid_data, method='euclidean', permutations=9999))
(perman_growth <- adonis(bray_growth_cent$centroids ~ ecosystem + treatment, centroid_data, method='euclidean', permutations=9999))
(perman_resp <- adonis(bray_resp_cent$centroids ~ ecosystem + treatment, centroid_data, method='euclidean', permutations=9999))
(perman_c <- adonis(bray_c_cent$centroids ~ ecosystem + treatment, centroid_data, method='euclidean', permutations=9999))


# mantel correlation between dissimilarity measures
euc_abund <- dist(bray_abund_cent$centroids[, 1:2])
euc_growth <- dist(bray_growth_cent$centroids[, 1:2])
euc_resp <- dist(bray_resp_cent$centroids[, 1:2])
euc_c <- dist(bray_c_cent$centroids[, 1:2])

mantel(euc_abund, euc_growth, permutations=9999) # r = 0.77 ***
mantel(euc_abund, euc_resp, permutations=9999) # r = 0.77 ***
mantel(euc_abund, euc_c, permutations=9999) # r = 0.69 ***
mantel(euc_growth, euc_resp, permutations=9999) # r = 0.98 ***