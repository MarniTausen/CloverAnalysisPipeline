library(ggplot2)

spectrum <- read.csv("sfs.csv")
spectrum <- spectrum[2:nrow(spectrum),1:2]

cat("Number of sites:", sum(spectrum$Count))

cat("MAF sites:", sum(spectrum$Count[1:20]))

cat("After removing MAF:", sum(spectrum$Count[-1:-20]))

spectrum$Proportion <- spectrum$Count/sum(spectrum$Count)

#spectrum <- spectrum[1:(floor(nrow(spectrum)/2)-1), 2]+spectrum[nrow(spectrum):(floor(nrow(spectrum)/2)+2), 2]

expected <- (1/1:393)/sum(1/1:393)

spectrum$Proportion[1:(ceiling(393/2)+1)] <- spectrum$Proportion[1:(ceiling(393/2)+1)] + spectrum$Proportion[394:(ceiling(393/2))]

expected <- expected[1:(ceiling(393/2)+1)]+expected[394:(ceiling(393/2))]

spectrum <- spectrum[1:(ceiling(393/2)+1), 1:3]


ggplot(spectrum, aes(x=AN, y=Proportion)) +
    geom_bar(stat="identity", fill="cyan", color="black", size=0.2) +
    labs(y="Proportion", x="Allele Frequency") + theme_classic() +
    geom_line(data=NULL, aes(x=1:198, y=expected), size=0.7, alpha=0.7)


#qplot(as.vector(rep(1:198, spectrum)), geom="histogram", binwidth=1,
#      color=I("black"), fill=I("cyan")) + theme_classic() + labs(y="Count", x="Allele Frequency")

