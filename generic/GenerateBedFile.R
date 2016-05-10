# We are thankful to David C. Wedge, Wellcome Trust Sanger Institute, who kindly provided the 
# original code, which we adapted for our own purposes.

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
snps = read.table(args[2],sep="\t",header=F,stringsAsFactors=F,row.names=NULL)
names(snps)[2:3] = c("chr","pos")
chrs = unique(snps$chr)
maxpos = vector(mode="integer",length=length(chrs))
bed.table = array(NA,c(0,3))
windowsize = 10000
for(c in 1:length(chrs)){	
	chr=chrs[c]
	maxpos[c] = windowsize*ceiling(max(snps$pos[snps$chr==chr])/windowsize)
	bed.table = rbind(bed.table,cbind(chr,format(seq(0,maxpos[c]-windowsize,windowsize),scientific=F,trim=T),format(seq(windowsize,maxpos[c],windowsize),scientific=F,trim=
T)))
}
write.table(cbind(chrs,maxpos),args[3],col.names = c("chr","length"),row.names=F,quote=F,sep="\t")
write.table(bed.table,paste("windows_",windowsize/1000,"kb.bed",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
