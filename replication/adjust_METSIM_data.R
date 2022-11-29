library(dplyr)
library(data.table)

#import Olink data in wide and remove Nan
targets <- read.table("GSE135134_METSIM_subcutaneousAdipose_RNAseq_log-TPM1_n434_PCA1.txt", sep='\t', header=TRUE)
targets[targets$gene.stable.ID=="", ] <- NA
targets_na <- na.omit(targets)

targets_filt <- subset(targets_na, select = -c(Gene.stable.ID, principal.component.1))
targets_scale <- data.frame(scale(targets_filt))

dat_adj <- function(df) {
	fit <- lm(df ~ principal.component.1, data = targets_na)
	return(resid(fit))
}

adj_res <- lapply(targets_filt, dat_adj)
df_adj <- data.frame(adj_res)
cols <- subset(targets_na, select = c(Gene.stable.ID))
out <- cbind(cols, df_adj)

write.table(out, "GSE135134_METSIM_subcutaneousAdipose_RNAseq_log-TPM1_n434_PCA1_adj-PC1.txt", row.names=FALSE, sep="\t", quote = FALSE)
