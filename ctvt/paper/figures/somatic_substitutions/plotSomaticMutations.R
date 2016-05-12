args <- commandArgs(TRUE)

germlineData <- read.table("graph_input_all.txt", header=T, sep="")
somaticData <- read.table("graph_input_somatic.txt", header=T, sep="")
pdf("all_mutations_per_clade.pdf",useDingbats=FALSE)

clade1="clade1"
clade2="clade2"
clade3="clade3"
clade4="clade4"
clade5="clade5"
# plot using a different colour for germline
plot(germlineData$POS, germlineData$NUM, col="grey", type="h", bty="l", xlab='',ylab='',xaxt='n',ylim=c(0,50))

par(new=TRUE)
# plot each clade using a different colour
plot(somaticData$POS, somaticData$NUM, col=ifelse(somaticData$CLADE %in% clade1, "blue", ifelse(somaticData$CLADE %in% clade2, "red1", ifelse(somaticData$CLADE %in% clade3, "green3", ifelse(somaticData$CLADE %in% clade4, "yellow1",ifelse(somaticData$CLADE %in% clade5, "blueviolet","grey"))))), type="h", bty="l", xlab='', xaxt='n',ylab="Number of Somatic Mutations",main="Somatic mutations per clade",ylim=c(0,50))
 
# Exit
