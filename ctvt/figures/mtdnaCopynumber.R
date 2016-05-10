args <- commandArgs(TRUE)
allData <- read.table(args[1], header=T, sep="")
pdf(args[2],useDingbats=FALSE)
clade1="clade1"
clade2="clade2"
clade3="clade3"
clade4="clade4"
clade5="clade5"
host="host"

plot(allData$POS,allData$DEPTH, pch=16,cex=0.5, col=ifelse(allData$CLADE %in% clade1, "firebrick1",ifelse(allData$CLADE %in% clade2,"blue", ifelse(allData$CLADE %in% clade3, "green3", ifelse(allData$CLADE %in% clade4,"gold1",ifelse(allData$CLADE %in% clade5,"blueviolet","black"))))),bty='l',ylab="",yaxt="n",xlab="") #use cex to change size of points
# with labels
#dotchart(allData$DEPTH, labels=allData$NAME,cex=.4,groups=allData$CLADE, main="Estimated Relative Mitochondrial Copynumber",gcolor="black",color=allData$color, pch=15)
# without labels
#dotchart(allData$DEPTH, labels=NULL,groups=allData$CLADE, gcolor="black",color=allData$color, pch=16,bty='l',cex=0.5)
rect(-4,-100,11100,4170,col = "mistyrose",border=NA)
rect(-4,4171,11000,10560,col = "lightcyan",border=NA)
rect(-4,10561,11000,11801,col = "azure",border=NA)
rect(-4,11802,11000,12200,col = "lightgoldenrodyellow",border=NA)
rect(-4,12201,11000,12501,col = "lavender",border=NA)
rect(-4,12502,11000,21926,col = "gray97",border=NA)
par(new=TRUE) 
plot(allData$POS,allData$DEPTH, pch=16,cex=0.57, col=ifelse(allData$CLADE %in% clade1, "firebrick1",ifelse(allData$CLADE %in% clade2,"blue", ifelse(allData$CLADE %in% clade3, "green3", ifelse(allData$CLADE %in% clade4,"gold1",ifelse(allData$CLADE %in% clade5,"blueviolet","black"))))),bty='l',ylab="",yaxt="n",xlab="") #use cex to change size of points
# Exit
q()
