args <- commandArgs(TRUE)
p_angle1 <- function(a1) {
#########pi for labels to run counterclockwise;-pi for labels to run clockwise
p_a1=a1/16727*2*-pi+pi/2
return(p_a1)
}

p_angle <- function(a1,a2) {
p_a1=a1/16727*2*-pi+pi/2
p_a2=a2/16727*2*-pi+pi/2
return(c(seq(p_a1,p_a2,by=1/120*2*-pi),p_a2))
}
###################
#
# Text file containing specifications for main mitochondrial figure
#
##################
genelist<-as.matrix(read.table("dogGene_chrMTs.txt",header=FALSE,sep="\t"))
len_genes=length(genelist[,1])

r1=9
r2=10
plot(c(-15,15),c(-15,15),axes=FALSE,xlab="",ylab="",main=NULL)

###############
#
# Draw lines on central mitochondrial ring marking intergenic regions
#
###############
pos1=5730
pos2=5760
polar_pos=p_angle(pos1,pos2)
x1<-r1*cos(polar_pos)
x2<-r2*cos(rev(polar_pos))
y1<-r1*sin(polar_pos)
y2<-r2*sin(rev(polar_pos))
polygon(c(x1,x2),c(y1,y2),col="#66B2FF")
pos1=8270
pos2=8293
polar_pos=p_angle(pos1,pos2)
x1<-r1*cos(polar_pos)
x2<-r2*cos(rev(polar_pos))
y1<-r1*sin(polar_pos)
y2<-r2*sin(rev(polar_pos))
polygon(c(x1,x2),c(y1,y2),col="#66B2FF")
pos1=5892
pos2=5902
polar_pos=p_angle(pos1,pos2)
x1<-r1*cos(polar_pos)
x2<-r2*cos(rev(polar_pos))
y1<-r1*sin(polar_pos)
y2<-r2*sin(rev(polar_pos))
polygon(c(x1,x2),c(y1,y2),col="#66B2FF")

###################
#
# Creates central map of mitochondria
#
###################

for (n in 1:len_genes) {
        start=as.integer(genelist[n,5])
        stop=as.integer(genelist[n,6])
        polar_pos=p_angle(start,stop)
        x1<-r1*cos(polar_pos)
        x2<-r2*cos(rev(polar_pos))
        y1<-r1*sin(polar_pos)
        y2<-r2*sin(rev(polar_pos))
        polygon(c(x1,x2),c(y1,y2),col=genelist[n,17])
        if ((stop-start) > 500){
                fontangle=(min(polar_pos)+max(polar_pos))/2*180/pi-90
                if (fontangle > 90 & fontangle < 270){
                        fontangle=fontangle+180
                }
                par(cex=0.45,font=2, srt=fontangle)
                text((r1+r2)/2*cos((min(polar_pos)+max(polar_pos))/2),(r1+r2)/2*sin((min(polar_pos)+max(polar_pos))/2),genelist[n,13],adj=0.5)
        }
        if ((stop-start) <=500 ) {
                fontangle=(min(polar_pos)+max(polar_pos))/2*180/pi

                fontangle1=0
                if (nchar(genelist[n,13]) > 1) fontangle1=fontangle
                this_adj=1
                if (fontangle1 > 90 & fontangle1 < 270) {
                        fontangle1=fontangle1+180
                        this_adj=0
                }
                this_offset_y=0
#               if (n==30) this_offset_y=-0.3
#               if (n==32) this_offset_y= 0.3
                ##########labels all of the small tRNA segments
                par(cex=0.2, srt=fontangle1)
#               if (genelist[n,4]=="+") text((r1-0.5)*cos(fontangle*pi/180),(r1-0.5)*sin(fontangle*pi/180)+this_offset_y,genelist[n,13],adj=this_adj)
#               if (genelist[n,4]=="-") text((r2+0.5)*cos(fontangle*pi/180),(r2+0.5)*sin(fontangle*pi/180)+this_offset_y,genelist[n,13],adj=(1-this_adj))
                if (genelist[n,4]=="+") text((r1+r2)/2*cos(fontangle*pi/180),(r1+r2)/2*sin(fontangle*pi/180)+this_offset_y,genelist[n,13],adj=0.5)
                if (genelist[n,4]=="-") text((r1+r2)/2*cos(fontangle*pi/180),(r1+r2)/2*sin(fontangle*pi/180)+this_offset_y,genelist[n,13],adj=0.5)
        }
}


tickpos<-c(1,seq(1000,16000,by=1000))
p_tickpos=p_angle1(tickpos)
par(cex=0.4, srt=0)
text((r2+0.7)*cos(p_tickpos),(r2+0.5)*sin(p_tickpos),tickpos)
###switching cos and sin rotates position of 1 clockwise by 90 degrees

for (n in 1:17) {
        lines(c((r2+0.05)*cos(p_tickpos[n]),(r2+0.15)*cos(p_tickpos[n])),c((r2+0.05)*sin(p_tickpos[n]),(r2+0.15)*sin(p_tickpos[n])))
}


#######Specify name of file containing list of mutations
varlist<-as.matrix(read.table("muts_inds.dat",header=TRUE,sep="\t"))
mutationData <- read.table("muts_inds.dat", header=T, sep="")
#######################################################
C2A="C>A"
C2G="C>G"
C2T="C>T"
T2A="T>A"
T2C="T>C"
T2G="T>G"
A2G="A>G"
A2C="A>C"
A2T="A>T"
G2A="G>A"
G2C="G>C"
G2T="G>T"
IN="C>CTT"
varlist[,2][varlist[,2]=="C>A"]<-5
varlist[,2][varlist[,2]=="C>G"]<-1
varlist[,2][varlist[,2]=="C>T"]<-2
varlist[,2][varlist[,2]=="T>A"]<-6
varlist[,2][varlist[,2]=="T>C"]<-4
varlist[,2][varlist[,2]=="T>G"]<-3
varlist[,2][varlist[,2]=="C>CTT"]<-7

varlist[,5][varlist[,5]=="3"]<-1 #syn
varlist[,5][varlist[,5]=="8"]<-0 #
varlist[,5][varlist[,5]=="4"]<-6 #noncoding
varlist[,5][varlist[,5]=="2"]<-4 #nonsyn


angle_forCircle<-p_angle(1,16727)

polygon((r1-4)*cos(angle_forCircle),(r1-4)*sin(angle_forCircle),border=8)
polygon((r1-3)*cos(angle_forCircle),(r1-3)*sin(angle_forCircle),border=8)
polygon((r1-2)*cos(angle_forCircle),(r1-2)*sin(angle_forCircle),border=8)
polygon((r1-1)*cos(angle_forCircle),(r1-1)*sin(angle_forCircle),border=8)
polygon((r2+1)*cos(angle_forCircle),(r2+1)*sin(angle_forCircle),border=8)
polygon((r2+2)*cos(angle_forCircle),(r2+2)*sin(angle_forCircle),border=8)
polygon((r2+3)*cos(angle_forCircle),(r2+3)*sin(angle_forCircle),border=8)
polygon((r2+4)*cos(angle_forCircle),(r2+4)*sin(angle_forCircle),border=8)
lines(c(0,0),c(r2+1,r2+4))

############Draw little lines marking the clade rings
#lines(c(0,0.2),c(r2+2,r2+2))
#lines(c(0,0.2),c(r2+3,r2+3))
############Text describing which ring corresponds to which clade
text(0.2,r2+1.2,"Clade 4",adj=0,cex=1.1,srt=0)
text(0.2,r2+2.2,"Clade 3",adj=0,cex=1.1,srt=0)
text(0.2,r2+3.2,"Clade 2",adj=0,cex=1.1,srt=0)
text(0.2,r2+4.2,"Clade 1",adj=0,cex=1.1,srt=0)
text(0,r2+5.0,"H strand",cex=1.6,srt=0)

lines(c(0,0),c(r1-1,r1-4))


#lines(c(0,0.2),c(r1-2,r1-2))
#lines(c(0,0.2),c(r1-1,r1-1))
text(0.2,r1-3.8,"Clade 4",adj=0,cex=1.1,srt=0)
text(0.2,r1-2.8,"Clade 3",adj=0,cex=1.1,srt=0)
text(0.2,r1-1.8,"Clade 2",adj=0,cex=1.1,srt=0)
text(0.2,r1-0.8,"Clade 1",adj=0,cex=1.1,srt=0)
text(0,r1-5.0,"L strand",cex=1.6,srt=0)

varlist_p<-varlist[varlist[,3]=="+",]
varlist_n<-varlist[varlist[,3]=="-",]
#################Plot variant points
points((r1-3+2*as.double(varlist_p[,4]))*cos(p_angle1(as.integer(varlist_p[,1]))),(r1-3+2*as.double(varlist_p[,4]))*sin(p_angle1(as.integer(varlist_p[,1]))),col=ifelse(mutationData$MUT %in% C2A, "cyan", ifelse(mutationData$MUT %in% C2T, "red", ifelse(mutationData$MUT %in% C2G, "black", ifelse(mutationData$MUT %in% T2A,"magenta", ifelse(mutationData$MUT %in% T2G,"green",ifelse(mutationData$MUT %in% T2C, "blue",ifelse(mutationData$MUT %in% G2T, "cyan", ifelse(mutationData$MUT %in% G2A, "red", ifelse(mutationData$MUT %in% G2C, "black", ifelse(mutationData$MUT %in% A2T,"magenta", ifelse(mutationData$MUT %in% A2C,"green",ifelse(mutationData$MUT %in% A2G, "blue",ifelse(mutationData$MUT %in% IN,"goldenrod1","goldenrod1"))))))))))))),pch=as.integer(varlist_p[,5]),cex=1)
#original points
#points((r1-3+2*as.double(varlist_p[,4]))*cos(p_angle1(as.integer(varlist_p[,1]))),(r1-3+2*as.double(varlist_p[,4]))*sin(p_angle1(as.integer(varlist_p[,1]))),col=as.integer(varlist_p[,2]),pch=as.integer(varlist_p[,5]),cex=1)
######points((r2+1+2*as.double(varlist_n[,4]))*cos(p_angle1(as.integer(varlist_n[,1]))),(r2+1+2*as.double(varlist_n[,4]))*sin(p_angle1(as.integer(varlist_n[,1]))),col=as.integer(varlist_n[,2]),pch=as.integer(varlist_n[,5]),cex=1)
#points((r2+1+2*as.double(varlist_n[4]))*cos(p_angle1(as.integer(varlist_n[,1]))),(r2+1+2*as.double(varlist_n[4]))*sin(p_angle1(as.integer(varlist_n[,1]))),col=as.integer(varlist_n[,2]),pch=as.integer(varlist_n[,5]),cex=1)
points((r2+1+2*as.double(varlist_n[4]))*cos(p_angle1(as.integer(varlist_n[,1]))),(r2+1+2*as.double(varlist_n[4]))*sin(p_angle1(as.integer(varlist_n[,1]))),col=ifelse(mutationData$MUT %in% C2A, "cyan", ifelse(mutationData$MUT %in% C2T, "red", ifelse(mutationData$MUT %in% C2G, "black", ifelse(mutationData$MUT %in% T2A,"magenta", ifelse(mutationData$MUT %in% T2G,"green",ifelse(mutationData$MUT %in% T2C, "blue",ifelse(mutationData$MUT %in% G2T, "cyan", ifelse(mutationData$MUT %in% G2A, "red", ifelse(mutationData$MUT %in% G2C, "black", ifelse(mutationData$MUT %in% A2T,"magenta", ifelse(mutationData$MUT %in% A2C,"green",ifelse(mutationData$MUT %in% A2G, "blue",ifelse(mutationData$MUT %in% IN,"goldenrod1","goldenrod1"))))))))))))),pch=as.integer(varlist_n[,5]),cex=1)

##############
#
#Different legends classifying mutation type
#
##############
#####need to reset par here or everything will be off
par(cex=0.45,font=2,srt=0)
#legend("topleft",pch=15,col=c("#66B2FF","#009900","#FF6666","#FFFF33","#88FF00","#FFA166"),legend=c("intergenic","CDS(+)","tRNA(+)","rRNA(+)","CDS(-)","tRNA(-)"),cex=2)
#legend("topleft",pch=16,col=c(5,1,2),legend=c("Germline","Somatic","Potential Somatic"),cex=1.5,bty="n")
legend("topleft",pch=16,col=c(5,1,2,6,4,3,7),legend=c("A>C T>G","G>C C>G","T>C A>G","A>T T>A","C>T G>A","G>T C>A"),cex=1.7,bty="n")
#legend("bottomright",pch=c(4,1,6,0),col=1,legend=c("non-syn","syn","non-coding","intergenic"),cex=1.5,bty="n")
legend("bottomright",pch=c(4,6,0,1),col=1,legend=c("Somatic","Potential Somatic","Germline","Indel"),cex=1.5,bty="n")

text(0,1.4,"Dog mtDNA",cex=2.5)
text(0,0,"16,727 bp",cex=2)
#text(0,-1.4,"Somatic Substitutions (N=X)",cex=1.5)
par(bg="white")
