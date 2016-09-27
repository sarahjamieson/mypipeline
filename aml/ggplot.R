args <- commandArgs(TRUE)
frag <- args[1]
ngs <- args[2]
itd <- args[3]

df <- data.frame(frag, ngs, itd)

library(ggplot2)

png("all_samples.png")

p <- ggplot(data = df, aes(x = frag, ngs))
p + geom_point(aes(size = itd))

dev.off()
