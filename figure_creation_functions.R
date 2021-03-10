# utility functions for Figure creation analysis pipeline

# function that produces randomized color
tax_color <- function(x, y, alpha=1, include_other=T) {
  n_tax <- nlevels(x)
  if(n_tax / 2 != n_tax %/% 2) n_tax <- n_tax + 1
  set.seed(y)
  tax_colors <- hsv(seq(0,1,length.out=n_tax), rep(c(.6,1), n_tax / 2), rep(c(.6,1), n_tax / 2), alpha=alpha)
  if(include_other) {
    tax_colors <- c(gray(.5, alpha=alpha), tax_colors[sample.int(n_tax)])
  } else {
    tax_colors <- tax_colors[sample.int(n_tax)]
  }
  return(tax_colors)
}

# function that edits colors based on saturation and brightness (s & v) as well as alpha settings
color_edit <- function(x, s=1, v=1, alpha=1, mode=c('product', 'sum')) {
  if(missing(mode)) mode <- 'product' else mode <- match.arg(mode, c('product', 'sum'))
  if(is.character(x)) x <- col2rgb(x)
  x <- apply(x, 2, function(x) rgb2hsv(r=x['red'], g=x['green'], b=x['blue'], maxColorValue=255))
  rownames(x) <- c('h', 's', 'v')
  if(mode=='product') {
    return(hsv(h=x['h',], s=x['s',] * s, v=x['v',] * v, alpha=alpha))
  } else {
    return(hsv(h=x['h',], s=x['s',] + s, v=x['v',] + v, alpha=alpha))
  }
}

# function which changes scientific notation
good_digits <- function(n) {
  n <- format(n, digits=3, scientific=F)
  n <- gsub('(.*)(0{6})$', '\\1M', n)
  n <- gsub('(.*)(0{3})$', '\\1K', n)
  n <- gsub('^(.*)e\\+00', '\\1', n)
  n <- sub('^^(\\d*)\\.(0+)$' ,'\\1', n)
  if (any(grepl('^(\\d{2,})00K|^(\\d{4,})K', n))) {
    to_fix <- grepl('^(\\d{2,})00K|^(\\d{4,})K', n)
    n[to_fix] <- gsub('^(\\d?)(\\d?)(0*)K$', '\\1\\.\\2M', n[to_fix])
  }
  n
}

# function that converts matrices or data frames to long format
qsip_long <- function(x, group_name='rep_group', value='value') {
  if (is.matrix(x)) {
    y <- reshape(data.frame(x),
                 idvar='taxonID',
                 ids=rownames(x),
                 varying=list(colnames(x)),
                 times=colnames(x),
                 timevar=group_name,
                 direction='long')
    names(y)[2] <- value
  } else if (is.data.frame(x)) {
    y <- reshape(x,
                 idvar='taxonID',
                 varying=list(names(x[,-1])),
                 times=names(x[,-1]),
                 timevar=group_name,
                 direction='long')
    names(y)[length(y)] <- value
  }
  rownames(y) <- NULL
  y <- y[,c('taxonID', group_name, value)]
  y <- y[!is.na(y[value]),]
}

# use for summarizing across replicates
se <- function(x) sd(x, na.rm=T) / length(x)
ci <- function(x, ci=.95) {
  y <- t.test(x, conf.level=ci)
  y <- y$conf.int
  (y[2] - y[1])/2
}

# convert to long form
reshape_long <- function(x, value_name='value', idvar=c('ecosystem', 'treatment')) {
  x <- reshape(x,
               idvar=idvar,
               varying=list(names(x)[-(1:length(idvar))]),
               times=names(x)[-(1:length(idvar))],
               timevar='measure',
               v.names=value_name,
               direction='long')
  rownames(x) <- NULL
  return(x)
}

# binary function to get differences (replacing single NAs with 0, but double NAs should return NA)
`%+%` <- function(x, y) {
  xx <- x; yy <- y
  xx[is.na(xx)] <- 0
  yy[is.na(yy)] <- 0
  z <- xx + yy
  z[is.na(x) & is.na(y)] <- NA
  return(z)
}

# pop_flux() subtracts new pop from old pop (Nt - N0), creating negative values when death is higher than birth
# however, if reverse=TRUE, then population flux will be calculated based on population numbers after incubation
pop_flux <- function(r, N, t=7, reverse=F, exponential=T) {
  if(exponential==TRUE) {
    if(reverse==FALSE) ((N * exp(r * t)) - N)
    else (N - (N * exp(-r * t)))
  } else {
    if(reverse==FALSE) (N * r * t) - N
    else N - (N / ((r * t) + 1))
  }
}

# Function to reorder taxonomies by dominant growers, grouping all others that cumulatively account for 1 - cutoff or less
# e.g. cutoff of .9 means taxonomic groups whose contribution to total ecosystem grwoth is less than .1 will be grouped
#  together as non-dominant growers
refactor_groups <- function(x, level='Phylum', cutoff=.9, relevel=F, rename_uncultured=T) {
  y <- x
  if(is(y)!='phyloseq') stop('Must provide phyloseq object')
  if(class(level)=='numeric') level <- rank_names(y)[level]
  level_num <- which(rank_names(y)==level)
  trunc_tax <- as(y@tax_table, 'matrix')[,1:level_num]
  trunc_feature <- as(y@otu_table, 'matrix')
  if(taxa_are_rows(y)==F) trunc_feature <- t(trunc_feature)
  # remove NA occurrences
  na_position <- apply(trunc_tax, 1, function(x) {y <- which(is.na(x)); ifelse(length(y)==0, 0, y)})
  if(rename_uncultured) {
    unclt_position <- apply(trunc_tax, 1, function(x) {y <- which(grepl('uncultured', x)); ifelse(length(y)==0, 0, y)})
  }
  for(i in 1:nrow(trunc_tax)) {
    if(na_position[i]!=0) {
      trunc_tax[i,na_position[i]:ncol(trunc_tax)] <- paste0('unknown ', trunc_tax[i,(na_position[i] -1)])
      y@tax_table[i,na_position[i]:ncol(trunc_tax)] <- paste0('unknown ', trunc_tax[i,(na_position[i] -1)])
    }
    if(rename_uncultured & unclt_position[i]!=0) {
      trunc_tax[i,unclt_position[i]:ncol(trunc_tax)] <- paste0('uncultured ', trunc_tax[i,(unclt_position[i] -1)])
      y@tax_table[i,unclt_position[i]:ncol(trunc_tax)] <- paste0('uncultured ', trunc_tax[i,(unclt_position[i] -1)])
    }
  }
  # split feature rowSums by the tax rank selected
  trunc_tax <- data.frame(trunc_tax, 
                          feature=rownames(trunc_tax), 
                          mergeID=apply(trunc_tax[,1:level_num], 1, paste, collapse=':'),
                          stringsAsFactors=F)
  sort_tax <- split(rowSums(trunc_feature, na.rm=T), trunc_tax$mergeID)
  # reorder phylum as well
  if(level!='Phylum') {
    sort_phylum <- split(rowSums(trunc_feature, na.rm=T), trunc_tax[['Phylum']])
    sort_phylum <- abs(sapply(sort_phylum, sum))
    sort_phylum <- sort(sort_phylum, decreasing=T)
  }
  # calculate cumulative contributions of the taxa, keep those levels that just meet the cutoff
  sort_tax <- abs(sapply(sort_tax, sum))
  sort_tax <- sort(sort_tax, decreasing=T)
  sort_tax <- data.frame(mergeID=names(sort_tax),
                         contribution=cumsum(sort_tax) / sum(sort_tax),
                         row.names=NULL)
  sort_tax$include <- sort_tax$contribution <= cutoff
  sort_tax[sum(sort_tax$include) + 1, 'include'] <- T
  trunc_tax <- merge(trunc_tax, sort_tax, all.x=T)
  # apply cut-off renaming to tax_table
  trunc_tax <- trunc_tax[match(taxa_names(y), trunc_tax$feature),]
  y@tax_table[!trunc_tax$include, level_num] <- 'Other'
  if(relevel) {
    y <- phyloseq::psmelt(y)
    trunc_tax <- trunc_tax[!duplicated(trunc_tax$mergeID) & 
                             trunc_tax$include==T,]
    trunc_tax <- trunc_tax[order(trunc_tax$contribution),]
    tax_levels <- rev(trunc_tax[[level]])
    if(any(duplicated(tax_levels))) {
      tax_levels[duplicated(tax_levels)] <- trunc_tax[duplicated(tax_levels), 'mergeID']
    }
    y[[level]] <- factor(y[[level]], levels=c('Other', tax_levels))
    if(level!='Phylum') {
      y[['Phylum']] <- factor(y[['Phylum']], levels=names(sort_phylum))
    }
    return(y)
  } else return(y)
}

# function which converts long-format to wide-format matrix to add to phyloseq
reshape_wide <- function(x, measure=c()) {
  if(is.null(measure)) stop('must provide measure for feature table', call.=F)
  x$site_trt_rep <- factor(interaction(x$ecosystem, x$treatment, x$replicate))
  wide <- reshape(x[,c('taxonID', measure, 'site_trt_rep')],
                  v.names=measure,
                  idvar='taxonID',
                  timevar='site_trt_rep',
                  direction='wide')
  names(wide) <- sub(paste0(measure, '\\.'), '', names(wide))
  rownames(wide) <- wide$taxonID
  as.matrix(wide[,-1])
}

# function to split matrix into submatrices
split_matrix <- function(x, factor) {
  tax_names <- colnames(x)
  sam_names <- rownames(x)
  #
  x <- split(x, factor)
  sam_names <- split(sam_names, factor)
  #
  x <- lapply(x, matrix, ncol=length(tax_names), byrow=F)
  x <- Map(function(x, y) {colnames(x) <- tax_names; rownames(x) <- y; x}, x, sam_names)
  x
}

# function takes NA positions (x), their identities (y), and matches them against known taxa (z)
match_ids <- function(x, y, z, up_to_genus=F) {
  ids_found <- vector('logical', length(x))
  # create temp matrix, same number of rows as y, with each row being the sequence of 1-ncol(y)
  temp <- rep(1, nrow(y)) %o% 1:ncol(y)
  # identify where in the matrix the BEST ID position happens
  temp <- sweep(temp, 1, x, '==')
  # match with matricized y to get the IDs of the best matches
  search_id <- t(as.matrix(y))[t(temp)]
  # look at tax in the data frame that has copy no. and genome size
  if(up_to_genus) {
    copy_tax <- as.matrix(z[,c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')])
  } else {
    copy_tax <- as.matrix(z[,c('Kingdom', 'Phylum', 'Class', 'Order', 'Family')])
  }
  copy_tax <- unique(c(copy_tax))
  copy_tax <- copy_tax[!is.na(copy_tax)]
  # Are these IDs among those that we have genome and 16S copy no. data for?
  ids_found <- is.element(search_id, copy_tax)
  # Any IDs that are not need to be moved up a taxonomic level, with the upper limit the Kingdom (1)
  x[!ids_found] <- x[!ids_found] - 1
  x[x==0] <- 1
  x
}

# function that repeats the match ID search until a taxon's lineage finds a match with known data
match_repeat <- function(x, y, z, up_to_genus=F) {
  # x represents NA position, but we want to find BEST ID, which is one taxa level ABOVE the NA position
  # so the x - 1 represents a change from NA position to best ID position
  x[x > 1] <- x[x > 1] - 1
  xx <- match_ids(x, y, z, up_to_genus=up_to_genus)
  xxx <- match_ids(xx, y, z, up_to_genus=up_to_genus)
  i <- 1
  while (!isTRUE(all.equal(xx, xxx))) {
    xx <- match_ids(xx, y, z, up_to_genus=up_to_genus)
    xxx <- match_ids(xxx, y, z, up_to_genus=up_to_genus)
    i <- i + 1
    if(i > 7) {
      break
      warning('match_repeat had difficulty finding matching IDs')
    }
  }
  xxx
}

# function that transforms S4 matrix into long-form data frame
mat_to_df <- function(x) {
  # if s3 class matrix, convert
  if(grepl('matrix', class(x))) y <- Matrix::Matrix(x, sparse=TRUE)
  # convert to dgTMatrix
  y <- as(y, 'dgTMatrix')
  # convert to data.frame
  data.frame(taxonID=y@Dimnames[[2]][y@j + 1], 
             sampleID=y@Dimnames[[1]][y@i + 1],
             value=y@x,
             stringsAsFactors=T)
}