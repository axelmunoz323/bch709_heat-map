df <- read.delim('results/chr_feature_counts.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)
png('results/chr_feature_counts_plot.png', width=1800, height=1200, res=180)
par(mfrow=c(2,1), mar=c(8,5,3,1))
barplot(df$gene_per_Mb, names.arg=df$chrom, las=2, col='#2a9d8f',
        main='Gene Density by Chromosome', ylab='Genes per Mb')
mat <- t(as.matrix(df[, c('n_gene', 'n_exon_unique', 'n_tRNA', 'n_snoRNA')]))
barplot(mat, beside=FALSE, names.arg=df$chrom, las=2,
        col=c('#264653','#e9c46a','#f4a261','#e76f51'),
        main='Feature Counts by Chromosome', ylab='Count')
legend('topright', legend=c('n_gene','n_exon_unique','n_tRNA','n_snoRNA'),
       fill=c('#264653','#e9c46a','#f4a261','#e76f51'), cex=0.85)
dev.off()
cat('Saved plot to results/chr_feature_counts_plot.png\n')
