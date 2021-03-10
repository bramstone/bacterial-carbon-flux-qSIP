# per-taxon contributions to carbon flow (maintenance, growth respiration, growth)
# Note: this may take a few minutes to calculate

library(phyloseq)

load('data/phyloseq_data_silva.RData')
load('data/respiration_18O_cue_draft_4.RData')

source('figure_creation_functions.R')


# Growth (ug C)---------------------------------------------------------------------

# asv table
growth_asv <- reshape_wide(resp_data, 'pop_c_ug')
growth_asv[is.na(growth_asv)] <- 0
growth_asv <- growth_asv[rowSums(growth_asv) > 0,]

# tax table
growth_tax <- tax_table(dat_silva)
growth_tax <- growth_tax[match(rownames(growth_asv), rownames(growth_tax)),]

# sample data
sam_data <- colnames(growth_asv)
sam_data <- data.frame(groupID=sam_data,
                       ecosystem=sub('(\\w{2})\\.(\\w+)\\.(\\d+)', '\\1', sam_data, perl=T),
                       treatment=sub('(\\w{2})\\.(\\w+)\\.(\\d+)', '\\2', sam_data, perl=T),
                       replicate=sub('(\\w{2})\\.(\\w+)\\.(\\d+)', '\\3', sam_data, perl=T))
sam_data$ecosystem <- factor(sam_data$ecosystem, levels=c('MC', 'PP', 'PJ', 'GL'))
rownames(sam_data) <- sam_data$groupID

# combine
growth_data <- phyloseq(otu_table(growth_asv, taxa_are_rows=T),
                        tax_table(growth_tax),
                        sample_data(sam_data))

growth_data2 <- refactor_groups(growth_data, relevel=T)
growth_data5 <- refactor_groups(growth_data, level='Family', cutoff=.9, relevel=T)
growth_data4 <- refactor_groups(growth_data, level='Order', cutoff=.95, relevel=T)



# Growth related respiration (ug C)-----------------------------------------------------

# asv table
resp_growth_asv <- reshape_wide(resp_data, 'qsip_growth_resp_c_ug')
resp_growth_asv[is.na(resp_growth_asv)] <- 0
resp_growth_asv <- resp_growth_asv[rowSums(resp_growth_asv) > 0,]

# combine
resp_growth_data <- phyloseq(otu_table(resp_growth_asv, taxa_are_rows=T),
                             tax_table(growth_tax),
                             sample_data(sam_data))

resp_growth_data2 <- refactor_groups(resp_growth_data, relevel=T)
resp_growth_data5 <- refactor_groups(resp_growth_data, level='Family', cutoff=.9, relevel=T)
resp_growth_data4 <- refactor_groups(resp_growth_data, level='Order', cutoff=.975, relevel=T)


# Maintenance related respiration (ug C)------------------------------------------------

# asv table
resp_maint_asv <- reshape_wide(resp_data, 'qsip_maint_resp_c_ug')
resp_maint_asv[is.na(resp_maint_asv)] <- 0
resp_maint_asv <- resp_maint_asv[rowSums(resp_maint_asv) > 0,]

# tax table
maint_tax <- tax_table(dat_silva)
maint_tax <- maint_tax[match(rownames(resp_maint_asv), rownames(maint_tax)),]

# combine
resp_maint_data <- phyloseq(otu_table(resp_maint_asv, taxa_are_rows=T),
                            tax_table(maint_tax),
                            sample_data(sam_data))

resp_maint_data2 <- refactor_groups(resp_maint_data, relevel=T)
resp_maint_data5 <- refactor_groups(resp_maint_data, level='Family', cutoff=.9, relevel=T)
resp_maint_data4 <- refactor_groups(resp_maint_data, level='Order', cutoff=.975, relevel=T)



# save
save(growth_data2, growth_data4, growth_data5, file='data/growth_ug_c_contributions_per_asv_draft_4.RData')
save(resp_growth_data2, resp_growth_data4, resp_growth_data5, file='data/resp_growth_ug_c_contributions_per_asv_draft_4.RData')
save(resp_maint_data2, resp_maint_data4, resp_maint_data5, file='data/resp_maint_ug_c_contributions_per_asv_draft_4.RData')
