args = commandArgs(trailingOnly=TRUE)
filename <- args[1]

df <- read.table(filename,header = T)
df <- df[2:dim(df)[2]]
labeled_df <- data.frame(label = as.numeric(gsub(colnames(df),pattern="NZMean_",replace="")),value = t(df),row.names = NULL)
write.table(labeled_df,gsub(filename,pattern = ".txt",replacement = ".csv"),row.names = F)