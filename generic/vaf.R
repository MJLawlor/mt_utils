# mitotypus data post-processing pipeline: R script for generating VAF plots
# author:       Máire Ní Leathlobhair (ml677@cam.ac.uk)
# date:         May 2016

args <- commandArgs(TRUE)
allData <- read.table(args[1], header=T, sep="")
vafData <- read.table(args[2], header=T, sep="")
pdf(args[3], useDingbats=FALSE)
cutoff <- as.numeric(args[4])
par(mar=c(5.1, 9.1, 5.1, 5.1), xpd=TRUE)

# plot all substitutions
plot(allData$POS, allData$VAF, xlab="", main="", ylab="", col="gray", cex=1.5, pch=19, xlim=c(0,17000), ylim=c(0,1), yaxt="n",bty='l')
par(new=TRUE)
# plot only unfiltered substitutions
plot(vafData$POS, vafData$VAF, xlab="", main=args[5], ylab="", col="black", cex=1.5, pch=19, xlim=c(0,17000), ylim=c(0,1), yaxt="n",bty='l')
# label tumour and host VAF bands
tumour_line = median(vafData$VAF)
host_line = median(vafData$VAF)
if (grepl("T",args[4])){	
	text(19000,tumour_line,labels = 'Tumour line')
} else {
	text(19000,host_line,labels = 'Host line')
}
segments(-40,tumour_line,17000,tumour_line,col="lightgrey",lty=c(1),lwd=1,lend="square")
# label VAF cutoff
text(19000,cutoff,labels = 'VAF Cutoff')
segments(-40, cutoff, 17000, cutoff, col="black", lty=c(1), lwd=1, lend="square")
# set y-axis intervals
axis(2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
# add legend
legend("topleft", inset=c(-.45,0),c("Variant","Contamination"), pch=c(19, 19), col=c("black", "gray"),bty="n") 
