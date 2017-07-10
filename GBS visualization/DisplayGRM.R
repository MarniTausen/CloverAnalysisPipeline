library(ggplot2)
library(reshape2)
library(dplyr)

GRM <- as.matrix(read.csv("GRM.csv", header=FALSE))

Plants <- read.csv("White_Clover_Individuals_Library_Info.csv")

cat("Dim of GRM", dim(GRM))

seq_names <- c('1', '10', '100', '101', '102', '104', '105', '106', '107', '108', '109', '11', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '12', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '13', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '14', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '15', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '16', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '17', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '18', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '19', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '2', '20', '200', '201', '21', '22', '23', '25', '26', '27', '28', '29', '3', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '4', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '5', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '6', '60', '61', '62', '64', '65', '66', '67', '69', '7', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '8', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '9', '90', '91', '92', '94', '95', '96', '97', '98', '99')

cat("Length of seq_names", length(seq_names))

rownames(GRM) <- seq_names
colnames(GRM) <- seq_names

GRM <- GRM[as.character(sort(as.integer(rownames(GRM)))),
           as.character(sort(as.integer(colnames(GRM))))]

collectFastaName <- Vectorize(function(x) as.integer(unlist(strsplit(as.character(x), "\\."))[1]))

plantnames <- collectFastaName(Plants$fastq.file)

indices <- vector()
for(n in rownames(GRM)){
    indices <- c(indices, which(plantnames==n))
}

all_names <- as.character(Plants$Variety[indices])
for(i in seq_along(all_names)){
    all_names[i] <- paste(unlist(strsplit(all_names[i], "\\_"))[2:3], collapse ="_")
}

rownames(GRM) <- all_names
colnames(GRM) <- all_names

GRM <- GRM[as.character(sort(rownames(GRM))),
           as.character(sort(colnames(GRM)))]

#rownames(GRM) <- 1:nrow(GRM)
#colnames(GRM) <- 1:nrow(GRM)

GRM[GRM>1] <- 1

GRMdata <- melt(GRM)

##heatmap(GRM)

ggplot(data=GRMdata, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_classic() +
    scale_fill_gradientn(colours=c("black", "purple", "pink", "white")) +
    labs(x = "Individual", y = "Individual") + theme(axis.text.x = element_text(size=4),
                                                     axis.text.y = element_text(size=4))

## PCA

colorvector <- c("#a76d00","#4168ff","#59e069","#420097","#acd44f",
                 "#ff71f1","#00931d","#de006b","#47dada","#fe0e16",
                 "#002552","#9a9900","#70006e","#2d6800","#980016",
                 "#003003","#f6b89c","#4e0726","#ab9173","#504200")

collectGroup <- Vectorize(function(x) unlist(strsplit(as.character(x), "\\_"))[1])

PCA <- princomp(GRM, scores=TRUE)

plot(PCA)

plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC2=PCA$scores[,2],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC2),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC2, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC2, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 2")


plotdata <- data.frame(PC2=PCA$scores[,2],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,2])))

find_hull <- function(df) df[chull(df$PC2, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC2, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC2, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC2, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 2", y="Principle component 3")


plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 3")



## Fixing mismatches

## Manually fixing labels

GRM <- GRM[-which(rownames(GRM)=="Aberconcor_07"), -which(rownames(GRM)=="Aberconcor_07")]
GRM <- GRM[-which(rownames(GRM)=="Aberconcor_03"), -which(rownames(GRM)=="Aberconcor_03")]
GRM <- GRM[-which(rownames(GRM)=="Aberpearl_10"), -which(rownames(GRM)=="Aberpearl_10")]
GRM <- GRM[-which(rownames(GRM)=="Avalon_05"), -which(rownames(GRM)=="Avalon_05")]
GRM <- GRM[-which(rownames(GRM)=="Avalon_06"), -which(rownames(GRM)=="Avalon_06")]
GRM <- GRM[-which(rownames(GRM)=="Barblanca_01"), -which(rownames(GRM)=="Barblanca_01")]
GRM <- GRM[-which(rownames(GRM)=="Borek_08"), -which(rownames(GRM)=="Borek_08")]
GRM <- GRM[-which(rownames(GRM)=="Brianna_01"), -which(rownames(GRM)=="Brianna_01")]

rownames(GRM)[rownames(GRM)=="Coolfin_04"] <- "Riesling_11"
colnames(GRM)[colnames(GRM)=="Coolfin_04"] <- "Riesling_11"

GRM <- GRM[as.character(sort(rownames(GRM))),
           as.character(sort(colnames(GRM)))]

cat("Dim of GRM: ", dim(GRM))

GRMdata <- melt(GRM)

ggplot(data=GRMdata, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_classic() +
    scale_fill_gradientn(colours=c("black", "purple", "pink", "white")) +
    labs(x = "Individual", y = "Individual") + theme(axis.text.x = element_text(size=4),
                                                     axis.text.y = element_text(size=4))

collectGroup <- Vectorize(function(x) unlist(strsplit(as.character(x), "\\_"))[1])

PCA <- princomp(GRM, scores=TRUE)

plot(PCA)

plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC2=PCA$scores[,2],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC2),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC2, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC2, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 2")


plotdata <- data.frame(PC2=PCA$scores[,2],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,2])))

find_hull <- function(df) df[chull(df$PC2, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC2, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC2, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC2, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 2", y="Principle component 3")


plotdata <- data.frame(PC1=PCA$scores[,1],
                       PC3=PCA$scores[,3],
                       Groups=collectGroup(names(PCA$scores[,1])))

find_hull <- function(df) df[chull(df$PC1, df$PC3), ]

hulls <- data.frame()
for(Group in unique(plotdata$Groups)){
    Subdf <- plotdata %>% filter(Groups==Group)
    hull <- Subdf[chull(Subdf$PC1, Subdf$PC3),]
    hulls <- rbind(hulls, hull)
}

ggplot(plotdata, aes(x=PC1, y=PC3, color=Groups)) + geom_point() + theme_bw() +
    geom_polygon(data = hulls, aes(x=PC1, y=PC3, fill=Groups), alpha = 0.1) +
    scale_color_manual(values=colorvector) +
    scale_fill_manual(values=colorvector) +
    labs(x="Principle component 1", y="Principle component 3")
