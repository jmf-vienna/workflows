
#.libPaths("/lisc/data/scratch/jmf/R/x86_64-pc-linux-gnu-library/4.5")

library(Nonpareil)
library(ggplot2)

samples <- read.table('intermediates/samples.txt', sep='\t', header=TRUE, as.is=TRUE)
attach(samples)
nps <- Nonpareil.set(File, labels=Name, plot.opts=list(plot.observed=FALSE))

png(file="intermediates/Sequencing_effort.png", height = 8, width = 16, units = "in", res = 150, pointsize = 10)
Nonpareil.set(File, labels=Name, plot.opts=list(plot.observed=FALSE))
dev.off()

detach(samples)

NPS_df <- as.data.frame(summary(nps))
write.csv(NPS_df, file="intermediates/Sequencing_coverage.csv")
