#!/usr/bin/r

library("ggplot2")
library("gridExtra")
library("grid")
library("cowplot")
library("dplyr")
library("RColorBrewer")
library("patchwork")
'%ni%' <- Negate('%in%')

## ape virome spreadsheet
gav<-read.table("virome_samples.tsv",sep="\t",header=T,fill=T,as.is=T, quote="")
## simplify species names
gasp<-ifelse(grepl("orilla",gav[,3]),"Gorilla",gav[,3])
gasp<-ifelse(grepl("ongo",gasp),"Pongo",gasp)
gasp<-ifelse(grepl("panisc",gasp),"Bonobo",gasp)
gasp<-ifelse(grepl("Pan",gasp),"Chimpanzee",gasp)

## define things needed later on
specs=c("Chimpanzee","Bonobo","Gorilla","Pongo")
spcol=colorRampPalette(c("grey35", "grey75"))(length(specs));names(spcol)<-specs
fosi=15

##### number of specimens per species
tonu<-as.data.frame(table(gasp))
tonu<-data.frame(Species=tonu$gasp, Value=tonu$Freq)
tonu$Species = factor(tonu$Species, levels=specs)
tonu<-tonu[order(tonu$Species),]

a0<- tonu %>% ggplot(aes(x=Species, y=Value, fill=Species)) +  geom_bar(stat="identity",position="dodge") +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face="bold",hjust=0.5,size=15),axis.text.x = element_text(angle = 45, vjust = 1.3, hjust=1,size=fosi),axis.text.y = element_text(size=fosi),legend.position="none") + ggtitle(label="Libraries per species") +ylab("Number of libraries") +xlab("") + scale_x_discrete(drop=F) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=85)) +  scale_y_continuous(breaks=c(0,20,40,60,80)) + scale_fill_manual(values=spcol)

pdf(paste("plots/1_s-number.pdf",sep=""),6,6)
a0
dev.off()


## number of reads kept
kread<-read.table("code/Violin-Plot_Capture.txt",sep="\t",header=T,fill=T,as.is=T, quote="")
bb<-data.frame(Type=kread[,2],Value=kread[,1])
bb$Type<-gsub(1,"Raw reads",bb$Type)
bb$Type<-gsub(2,"Reads after Trimming",bb$Type)
bb$Type<-gsub(3,"Unique reads",bb$Type)
bb$Type = factor(bb$Type, levels=c("Raw reads","Reads after Trimming","Unique reads"))

c1<-
  ggplot(bb) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=fosi*0.6)) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=45076183))  + ggtitle(label="Number of reads per library") +xlab("") + ylab("Number of reads")

c2<-
  ggplot(bb) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=fosi*0.6),) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=4000000))  + ggtitle(label="Number of reads per library (without tail)") +xlab("") + ylab("Number of reads") +  coord_cartesian(ylim = c(0,3819402*1.05))

pdf(paste("plots/2_filtered-reads.pdf",sep=""),8,5)
c1 + c2 + plot_annotation(tag_levels = "A") 
dev.off()


###### number of reads assigned to domains
reads<-read.table("code/dnareads.tsv",sep="\t",header=T,fill=T,as.is=T, quote="")
aa<-data.frame(Type=c(rep("Bacteria",214),rep("Homo",214),rep("Viruses",214)),Value=c(reads[,2],reads[,3],reads[,4]))
aa$Type = factor(aa$Type, levels=c("Homo","Bacteria","Viruses"))
ab<-data.frame(Type=c(rep("Bacteria",214),rep("Homo",214),rep("Viruses",214)),Value=c(reads[,5],reads[,6],reads[,7]))
ab$Type = factor(ab$Type, levels=c("Homo","Bacteria","Viruses"))
ac<-ab
ac$Value<-log(ac$Value)

# percentage
b1<-  ggplot(aa) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=12)) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=100))  + ggtitle(label="Percentage of assigned reads") +xlab("") + ylab("Abundance per library (%)")

# number
b2<-  ggplot(ab) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=12)) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=100))  + ggtitle(label="Number of assigned reads") +xlab("") + ylab("Abundance per library")

# log number
b3<-  ggplot(ac) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(), plot.title = element_text(face="bold",hjust=0.5,size=12))+  geom_segment(aes(x=0.0,y=0,xend=0.0,yend=17))  + ggtitle(label="Number of assigned reads") +xlab("") + ylab("Abundance per library (log-scale)")

### tbd: add A-C
pdf(paste("plots/3_assigned-reads.pdf",sep=""),12,5)
b2 + b3 + b1 +  plot_annotation(tag_levels = "A") 
dev.off()
