library(data.table)
setDTthreads(12)
df_temp <- fread("./scbs_data/scbs_filtered_kidney-atlas_matrix_var5per/methylation_fractions.csv.gz",header=TRUE,sep=",")
df_temp <- as.matrix(df_temp,rownames=1)
df_temp <- round(df_temp,4)
#logit(as.matrix(data.frame(a=c(.01,.02,.03,.1,.2,.3),b=c(0.001,0.5,0.8,0,1,0.001))),adjust=0.05)
df_temp <- logit(df_temp,adjust=0.05)
#methyl.logit.transform <- function(X) {
#	#X = 1 - X
#	delta <- 1e-3
#	U <- delta + X * (1 - 2*delta)
#	log(U/(1-U))
#}
