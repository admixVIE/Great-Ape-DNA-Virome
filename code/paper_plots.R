#!/usr/bin/r

library("ggplot2")
library("gridExtra")
library("grid")
library("cowplot")
library("dplyr")
library("RColorBrewer")
library("patchwork")
'%ni%' <- Negate('%in%')
options(scipen=100)

## ape virome spreadsheet
gav<-read.table("virome_samples.tsv",sep="\t",header=T,fill=T,as.is=T, quote="")
## simplify species names
gasp<-ifelse(grepl("orilla",gav[,3]),"Gorilla",gav[,3])
gasp<-ifelse(grepl("ongo",gasp),"Orangutan",gasp)
gasp<-ifelse(grepl("panisc",gasp),"Bonobo",gasp)
gasp<-ifelse(grepl("Pan",gasp),"Chimpanzee",gasp)

## define things needed later on
specs=c("Chimpanzee","Bonobo","Gorilla","Orangutan")
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
fosi=12
kread<-read.table("code/Violin-Plot_Capture.txt",sep="\t",header=T,fill=T,as.is=T, quote="")
bb<-data.frame(Type=kread[,2],Value=kread[,1])
bb$Type<-gsub(1,"Raw reads",bb$Type)
bb$Type<-gsub(2,"Reads after trimming",bb$Type)
bb$Type<-gsub(3,"Unique reads",bb$Type)
bb$Type = factor(bb$Type, levels=c("Raw reads","Reads after trimming","Unique reads"))

c1<-  ggplot(bb) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(),axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1,size=fosi),axis.title.x=element_text(size=fosi), plot.title = element_text(face="bold",hjust=0.5,size=12)) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=45076183))  + ggtitle(label="Number of reads per library") +xlab("") + ylab("Number of reads")

c2<-  ggplot(bb) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1,size=fosi),axis.title.x=element_text(size=fosi), plot.title = element_text(face="bold",hjust=0.5,size=12)) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=4000000))  + ggtitle(label="Number of reads per library (without tail)") +xlab("") + ylab("Number of reads") +  coord_cartesian(ylim = c(0,3819402*1.05))


###### number of reads assigned to domains
reads<-read.table("code/dnareads.tsv",sep="\t",header=T,fill=T,as.is=T, quote="")
aa<-data.frame(Type=c(rep("Bacteria",214),rep("Homo",214),rep("Viruses",214)),Value=c(reads[,2],reads[,3],reads[,4]))
aa$Type = factor(aa$Type, levels=c("Homo","Bacteria","Viruses"))
ab<-data.frame(Type=c(rep("Bacteria",214),rep("Homo",214),rep("Viruses",214)),Value=c(reads[,5],reads[,6],reads[,7]))
ab$Type = factor(ab$Type, levels=c("Homo","Bacteria","Viruses"))
ac<-ab
ac$Value<-log(ac$Value)

# percentage
b1<-  ggplot(aa) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=fosi), axis.title.x=element_text(size=fosi), plot.title = element_text(face="bold",hjust=0.5,size=12)) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=100))  + ggtitle(label="Percentage of assigned reads") +xlab("") + ylab("Abundance per library (%)")

# number
b2<-  ggplot(ab) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=fosi),axis.title.x=element_text(size=fosi), plot.title = element_text(face="bold",hjust=0.5,size=12)) + geom_segment(aes(x=0.0,y=0,xend=0.0,yend=8000000))  + ggtitle(label="Number of assigned reads") +xlab("") + ylab("Abundance per library")

# log number
b3<-  ggplot(ac) + theme_bw() + geom_violin(mapping=aes(x=Type,y=Value, fill=Type),adjust=1.0,draw_quantiles = c(0.5), scale="width", na.rm=T  ) + geom_jitter(mapping=aes(x=Type,y=Value), height = 0, width = 0.2,na.rm=T,inherit.aes=F) + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "white"), legend.position="none",axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=fosi),axis.title.x=element_text(size=fosi), plot.title = element_text(face="bold",hjust=0.5,size=12))+  geom_segment(aes(x=0.0,y=0,xend=0.0,yend=17))  + ggtitle(label="Number of assigned reads") +xlab("") + ylab("Abundance per library (log-scale)")

#plot
pdf(paste("plots/2_read_numbers.pdf",sep=""),9,9)
(c1 | c2 )  / ( b2 | b3 | b1 ) + plot_annotation(tag_levels = "A") 
dev.off()



######## overview of viruses from kraken2 output
library("ggplot2")
library("gridExtra")
library("grid")
library("cowplot")
library("dplyr")
library("RColorBrewer")
'%ni%' <- Negate('%in%')

krk<-read.table("TableS4.txt",sep="\t",header=T,fill=T,as.is=T, quote="")

# A) stacked barplot of numbers of reads assigend to specific types
kk1<-colSums(krk[,-1])
relms<-c("Monodnaviria","Duplodnaviria","Varidnaviria","Riboviria")
uncl<-c("unclassified.Viruses", "Anelloviridae")
relms<-cbind(rep("Realm",length(relms)+1),c(names(kk1[which(names(kk1)%in%relms)]),"Other/uncertain"),c(kk1[which(names(kk1)%in%relms)],sum(kk1[which(names(kk1)%in%uncl)])))
relms<-cbind(relms,as.numeric(relms[,3])/sum(as.numeric(relms[,3])))
fnam<-c("Herpesviridae", "Adenoviridae", "Circoviridae", "Hepadnaviridae", "Genomoviridae","Papillomaviridae", "Parvoviridae", "Poxviridae","Polyomaviridae","Retroviridae")
madna<-cbind(rep("Family",length(fnam)),names(kk1[which(names(kk1)%in%fnam)]),kk1[which(names(kk1)%in%fnam)])
madna<-cbind(madna,as.numeric(madna[,3])/sum(as.numeric(madna[,3])))
kk2<-rbind(relms,madna)

av<-data.frame(grp=kk2[,1],typ=kk2[,2],val=as.numeric(kk2[,4]))
av$typ = factor(av$typ, levels=c("Monodnaviria","Duplodnaviria","Varidnaviria","Riboviria","Other/uncertain",fnam))
av$grp = factor(av$grp, levels=c("Realm","Family"))
av_colors <-  setNames(c(brewer.pal(7,'BrBG')[-c(4:5)],brewer.pal(10,'Set3')), levels(av$typ))

p3f<-
  ggplot(data=av, aes(x=grp, y=val,fill=typ)) +   geom_bar(stat="identity",position="stack") + theme_minimal(base_size=18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face="bold",hjust=0.5,size=15)) + ggtitle(label="Virus types in data") +ylab("Proportion of viral reads in category") +xlab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust=1)) +  scale_fill_manual(values=av_colors) + theme(legend.position="none") 

p3a<-  av %>%  filter(grp %in% c("Realm")) %>% ggplot(aes(x=grp, fill=typ)) +   geom_bar() + scale_fill_manual(values=av_colors[which(av$grp %in% c("Realm"))],name="Realm") + theme(legend.title = element_text(size=16),legend.text = element_text(size=14))
leg1<-get_legend(p3a)
p3a<-p3a + theme(legend.position="none")
p3c<-  av %>%  filter(grp %in% c("Family")) %>% ggplot(aes(x=grp, fill=typ)) +   geom_bar() + scale_fill_manual(values=av_colors[which(av$grp %in% c("Family"))],name="Family (Mammalian)") + theme(legend.title = element_text(size=16),legend.text = element_text(size=14))
leg2<-get_legend(p3c)
p3c<-p3c + theme(legend.position="none")
p3<-plot_grid(p3f, plot_grid( leg1, leg2, nrow = 2), ncol = 2, rel_widths = c(7,3))
f1 <- arrangeGrob(p3, top = textGrob("A", x = unit(0.02, "npc"), y   = unit(0.8, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=18)))

# B) & C) number of libraries with reads in a given category
kk3<-ifelse(krk[,-c(1:8,15,20)]>0 ,1,0)
kk4<-ifelse(krk[,-c(1:8,15,20)]>24,1,0)
av1<-data.frame(nams=colnames(kk3),val=colSums(kk3));av1<-av1[which(av1$val>0),];av1$nams=factor(av1$nams,levels=av1$nams[order(av1$val,decreasing=T)])
av2<-data.frame(nams=colnames(kk4),val=colSums(kk4));av2<-av2[which(av2$val>0),];av2$nams=factor(av2$nams,levels=av2$nams[order(av2$val,decreasing=T)])
p1<-
  ggplot(data=av1, aes(x=nams, y=val,fill=nams)) +   geom_bar(stat="identity",position="stack") + theme_minimal(base_size = 18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face="bold",hjust=0.5,size=15)) + ggtitle(label="Libraries with any virus-assigned read") +ylab("Number of libraries") +xlab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)) +  scale_fill_manual(values=rep("orange",nrow(av1))) + theme(legend.position = "none")

p2<-
  ggplot(data=av2, aes(x=nams, y=val,fill=nams)) +   geom_bar(stat="identity",position="stack") + theme_minimal(base_size = 18) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face="bold",hjust=0.5,size=15)) + ggtitle(label="Libraries with at least 25 virus-assigned read") +ylab("Number of libraries") +xlab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)) +  scale_fill_manual(values=rep("red",nrow(av2))) + theme(legend.position = "none")

f3 <- arrangeGrob(p1, top = textGrob("B", x = unit(0, "npc"), y   = unit(0.8, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=18)))
f4 <- arrangeGrob(p2, top = textGrob("C", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=18)))


pdf(paste("plots/3-abundance.pdf",sep=""),12,8)
grid.arrange(f1,arrangeGrob(f3,f4, nrow = 2),
 ncol = 2,widths=c(1,0.9))         
dev.off()
