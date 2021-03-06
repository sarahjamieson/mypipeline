library(ggplot2)
png("/home/shjn/PycharmProjects/mypipeline/aml/static/aml/all_samples_no_outliers.png")
Frag_Analysis <- c(0.24, 0.17, 0.04, 0.21, 0.1, 0.77, 0.5, 0.01, 0.01, 0.34, 0.04, 0.16, 0.01, 0.42, 0.3, 0.04, 0.1, 0.02, 0.36, 0.4, 0.04, 1.39, 0.43)
NGS <- c(0.078, 0.085, 0.014, 0.079, 0.033, 0.373, 0.285, 0.007, 0.012, 0.284, 0.037, 0.051, 0.035, 0.228, 0.039, 0.026, 0.033, 0.010, 0.227, 0.309, 0.013, 0.609, 0.236)
ITD_size <- c(45, 57, 72, 81, 120, 72, 66, 84, 30, 207, 243, 45, 39, 57, 204, 246, 72, 27, 24, 42, 54, 201, 33)
df <- data.frame(Frag_Analysis, NGS, ITD_size)
p <- ggplot(data = df, aes(x = Frag_Analysis, NGS))
p + geom_point(aes(size=ITD_size)) + xlab("Fragment Analysis AR") + ylab("NGS AR")
dev.off()
