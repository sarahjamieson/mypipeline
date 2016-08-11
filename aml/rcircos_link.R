args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]

library(RCircos);

data(UCSC.HG19.Human.CytoBandIdeogram);
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);

out.file <- output;
png(file=out.file, height=480, width=480);

RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot();

link.data <- read.table(input, col.names=c("Chromosome", "chromStart", "chromEnd", "Chromosome.1", "chromStart.1", "chromEnd.1", "PlotColor"));
track.num <- 2;

RCircos.Link.Plot(link.data, track.num, FALSE);

dev.off()



