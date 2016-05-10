# We are thankful to David C. Wedge, Wellcome Trust Sanger Institute, who kindly provided the original code, which we adapted.
args <- commandArgs(TRUE)
even <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38)
data = read.table(args[1],sep="\t",header=F,row.names=NULL)
names(data) = c("chr","startpos","endpos","no.reads","bases.covered","segment.length","fractional.coverage")
data = data[order(data$chr,data$startpos),]
chr.lengths = read.table("chr_lengths.txt",sep="\t",header=T)
cum.lengths = cumsum(chr.lengths$length)
data$genome.pos = data$startpos

for(c in 2:nrow(chr.lengths)){
	data$genome.pos[data$chr==chr.lengths$chr[c]] = data$genome.pos[data$chr==chr.lengths$chr[c]] + cum.lengths[c-1]
}
png(args[2],width=2000,height=400)
#plot(data$genome.pos,data$no.reads, main=args[3],pch=20,col="red",ylim=c(0,1000),cex=0.5,xlab="",ylab="depth/10kb",xaxt="n")
plot(data$genome.pos,data$no.reads, main=args[3],pch=20,col=ifelse(data$chr %in% even, "dodgerblue", "darkorchid"),ylim=c(0,400),cex=0.5,xlab="",ylab="depth/10kb",xaxt="n")
abline(v=0,lty=1,col="lightgrey",lwd=2)
for(c in 1:nrow(chr.lengths)){
	abline(v=cum.lengths[c],lty=1,col="lightgrey",lwd=2)
	if(c==1){
		text(cum.lengths[c]/2,100,chr.lengths$chr[c],pos=1,cex=2)
	}else{
		text((cum.lengths[c-1]+cum.lengths[c])/2,100,chr.lengths$chr[c],pos=1,cex=2)
	}
}
#lines(x0=1000,125610000,y=415,415,col="gray44",lwd=20,lty=1)
#segments(x0=1000,y0=415,x1=125610000,y1=415,col="gray44",lwd=20,lty=1)
#segments(x0=12561100,y0=415,x1=214020000,y1=415,col="gray72",lwd=20,lty=1)
dev.off()
