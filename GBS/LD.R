library(ggplot2)
library(dplyr)

LD <- read.table("LD_10k_window.geno.ld", header=TRUE, sep="\t")
#LD <- read.table("TestData.txt", header=TRUE, sep="\t")

#LD$Distance <- abs(LD$POS1-LD$POS2)

LD %>% filter(!is.nan(R.2) & R.2>0 & N_INDV>180) %>%
    mutate(Distance = abs(POS1-POS2)) -> LD

chromosomes <- sort(as.vector(unique(LD$CHR)))

for(chromosome in chromosomes[-1]){
    chrsub <- LD[LD$CHR==chromosome,]
    expmodel <- lm(log(R.2) ~ Distance, data=chrsub)
    r2 <- summary(expmodel)$r.squared
    xseq <- seq(0, 10000, by=0.01)
    reg_fit <- exp(predict(expmodel, list(Distance=xseq)))
    names(reg_fit) <- NULL
    expfit <- data.frame(Distance=xseq, R.2=reg_fit)
    cat("Length of fit:", length(reg_fit), "Length of sequence", length(xseq))
    fig <- ggplot(chrsub, aes(x=Distance, y=R.2)) + geom_point(size=0.3) +
        theme_bw() + labs(title=chromosome, x = "Distance", y = "R-squared linkage") +
        theme(text = element_text(size=8)) +
        annotate("text", label=paste0("R^2 == ", r2), x=5000, y=0.90, parse=TRUE) +
        geom_line(data=expfit, aes(x=Distance, y=R.2), colour="red") +
        scale_x_log10(limits=c(1, 10000))
    ggsave(plot=fig, file=paste0(chromosome, ".png"), width=10, height=4, unit="in")
}
