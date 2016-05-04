
args <- commandArgs(TRUE)
allData <- read.table(args[1], header=T, sep="")
vafData <- read.table(args[2], header=T, sep="")
pdf(args[3], useDingbats=FALSE)
cutoff <- as.numeric(args[4])
par(mar=c(5.1, 9.1, 5.1, 5.1), xpd=TRUE)

plot(allData$POS, allData$PM, xlab="", main="", ylab="", col="gray", cex=1.5, pch=19, xlim=c(0,17000), ylim=c(0,1), yaxt="n",bty='l')
par(new=TRUE)
plot(vafData$POS, vafData$PM, xlab="", main=args[5], ylab="", col="black", cex=1.5, pch=19, xlim=c(0,17000), ylim=c(0,1), yaxt="n",bty='l')
tumour_line = median(vafData$PM)
host_line = median(vafData$PM)

segments(-40,tumour_line,17000,tumour_line,col="lightgrey",lty=c(1),lwd=1,lend="square")
text(19000,cutoff,labels = 'VAF Cutoff')
axis(2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
segments(-40, cutoff, 17000, cutoff, col="black", lty=c(1), lwd=1, lend="square")
legend("topleft", inset=c(-.45,0),c("Host contamination","Tumour Variant"), pch=c(18, 18), col=c("gray", "black",bty="n") 

