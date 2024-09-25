rm(list=ls())

library(readxl)
library(dplyr)
library(stringr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(colorspace)
library(plotly)
library(tidyr)
library(magrittr)
library(audiometry)
library(trackViewer)

####Bar_plot####
fulldata <- read_excel("D:\\프로젝트\\MPZL2\\Database_Filtering Condition1_0920.xlsx",sheet=3)


familycount <- fulldata %>% 
  group_by(Gene) %>% 
  summarize(v1=n())
familycount <- familycount[order(familycount$v1,decreasing = T),]

names <- familycount$Gene
number <- familycount[!duplicated(familycount$v1),]

pal <- c('21'='#86D8A4','12'='#86D8A4','7'='#FCAEBB','6'='#86D8A4','4'='#86D8A4','3'='#86D8A4','2'='#86D8A4','1'='#86D8A4')

ggplot(data=familycount,mapping=aes(x=reorder(Gene,`v1`),y=v1))+
  scale_x_discrete(limits=rev(names))+
  scale_y_continuous(limits=c(0,22.5),expand=c(0,0))+
  geom_bar(mapping=aes(fill=as.factor(v1)),stat="identity")+
  coord_flip()+
  geom_label(aes(label=v1),nudge_y=0.7,nudge_x=0,size=4,data=number)+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

sum(familycount$v1)

rm(names,pal,familycount,number,fulldata)

####variant type bar plot####
variant <- c("c.220C>T/c.220C>T","c.220C>T/c.463del","c.220C>T/c.68del","c.463del/c.463del")
n <- c(19,3,1,1)
data <- data.frame(variant,n)

pal <- c('c.220C>T/c.220C>T'='#d2ffd2','c.220C>T/c.463del'='#28e7ff','c.220C>T/c.68del'='#ffe900','c.463del/c.463del'='#ffdcff')

ggplot(data=data,mapping=aes(x=variant,y=n))+
  scale_x_discrete(limits=variant)+
  scale_y_continuous(limits=c(0,20),expand=c(0,0))+
  geom_bar(mapping=aes(fill=as.factor(variant)),stat="identity")+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

rm(variant,pal,n,data)


####Pie_Simple ver####
#14명, 28variant
data_vr2=data.frame(
  varianttype2=c("Frameshift","Nonsense"),
  v1=c(1,27)
)

plot_ly(data_vr2,labels=~varianttype2,values=~v1,textposition='none',
        marker=list(colors=c('#FFB3B5','#86D8A4')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)


data_vr2=data.frame(
  varianttype2=c("short indel","SNV"),
  v1=c(1,27)
)

plot_ly(data_vr2,labels=~varianttype2,values=~v1,textposition='none',
        marker=list(colors=c('#FFF780','#6CD6E4')),sort=FALSE) %>% 
  add_pie(hole=0.5) %>% 
  layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         showlegend=FALSE)
rm(data_vr2)

####Audio gram####
audio <- read_excel("D:\\프로젝트\\MPZL2\\MPZL2_Audio.xlsx",sheet=1)

#Audio 1 Age
x <- audio[59,9:14]
x <- t(x)
x <- unname(x[,1])

y <- audio[59,16:21]
y <- t(y)
y <- unname(y[,1])

z <- append(x,y)

fre <- c(250,500,1000,2000,4000,8000)
sp <- c(1,1,1,1,1,1,4,4,4,4,4,4)

audio_set <- data.frame(time=gl(2,6),
                        f=rep(fre,2),
                        sp=rep(sp,2),
                        t=z)

gg_pta(audio_set,xlab="Frequency (Hz)",ylab="Hearing level (dB)") +
  geom_point(mapping=aes(x = f, y = t, color = time, shape=factor(sp)), size = 5, stroke=1.5) +
  geom_line(aes(x = f, y = t, color = time), lwd=1)+
  scale_color_manual(values=c("red","blue"))+
  scale_shape_manual(values=c(1,4))+
  theme(legend.position="none")

rm(y,x,z,audio_set)

#Audio 2 Age
x <- audio[48,9:14]
x <- t(x)
x <- unname(x[,1])

x2 <- audio[49,9:14]
x2 <- t(x2)
x2 <- unname(x2[,1])

y <- audio[48,16:21]
y <- t(y)
y <- unname(y[,1])

y2 <- audio[49,16:21]
y2 <- t(y2)
y2 <- unname(y2[,1])

z <- append(x,y)
z <- append(z,x2)
z <- append(z,y2)

fre <- c(250,500,1000,2000,4000,8000)
sp <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4)

audio_set <- data.frame(time=rep(sp,1),
                        f=rep(fre,4),
                        t=z)

gg_pta(audio_set,xlab="Frequency (Hz)",ylab="Hearing level (dB)") +
  geom_point(mapping=aes(x = f, y = t, color = factor(time), shape=factor(time)), size = 5, stroke=1.5) +
  geom_line(aes(x = f, y = t, color = factor(time),linetype=factor(time)), lwd=1)+
  scale_color_manual(values=c("red","blue","red","blue"))+
  scale_shape_manual(values=c(1,4,1,4))+
  scale_linetype_manual(values = c(1,1,2,2))+
  theme(legend.position="none")

rm(y,x,z,x2,y2,audio_set)

####domain map####
#Domain & AminoAcid
mut <- c("p.Met1Ile","p.Leu18Phe","p.Pro23LeufsTer2","p.Ile24MetfsTer22",
         "p.Gln74Ter","p.Arg94Trp","p.Ala155LeufsTer10")
loc <- c(1,18,23,24,74,94,155)
col=c("#6CD6E4","#6CD6E4","#D1ACF9","#6CD6E4","#D1ACF9","#6CD6E4","#D1ACF9")

features <- GRanges("chr1",IRanges(c(27),width=c(113),names=c("Ig-like V-type")),
                    fill = c("#ccfccf"),
                    height=0.05)

MPZL2.gr <- GRanges("chr1",IRanges(loc,width=1,names=mut),
                    color=col,label.parameter.rot=c(60))
MPZL2.gr$SNPsideID <- c("bottom","bottom","bottom","bottom","bottom","bottom","bottom")
MPZL2.gr$label.parameter.gp <- gpar(fontsize=11)


lolliplot(MPZL2.gr,features,ranges=GRanges("chr1",IRanges(1,215)),xaxis=c(1,215),ylab=FALSE)

rm(mut,loc,features,MPZL2.gr)

#MPZL2 Exon-Intron
mut1 <- c("c.3G>T","c.52C>T","c.68del","c.72del","c.220C>T","c.280C>T","c.463del")

loc1 <- c(147,196,1210,1214,1362,1704,4123)

col=c("#6CD6E4","#6CD6E4","#D1ACF9","#6CD6E4","#D1ACF9","#6CD6E4","#D1ACF9")


features <- GRanges("chr1",
                    IRanges(c(1,145,1201,1650,4097,6985,7049,9065),
                            width=c(144,58,167,211,148,64,12,1818)),
                    fill = c("gray","white","white","white","white","white","gray","gray"),
                    height=0.05)

mpzl2.gr <- GRanges("chr1",IRanges(loc1,width=1,names=mut1),
                    color=col,
                    label.parameter.rot=c(60))

mpzl2.gr$label.parameter.gp <- gpar(fontsize=11)

lolliplot(mpzl2.gr,features,ranges=GRanges("chr1",IRanges(1,10882)),xaxis=c(1,10882),ylab=FALSE)

#####decade audiogram_right&left#####
#unknown age
audio <- audio[-54,]

audio0 <- subset(audio,`나이`<10)
a <- audio0[,9:14]
b <- audio0[,16:21]
names(a) <- c("250","500","1000","2000","4000","8000")
names(b) <- c("250","500","1000","2000","4000","8000")
audio0 <- rbind(a,b)
x <- data.frame(apply(audio0,2,mean,na.rm=TRUE))
x <- x[,1]

audio10 <- subset(audio,9<`나이`& 20>`나이`)
a <- audio10[,9:14]
b <- audio10[,16:21]
names(a) <- c("250","500","1000","2000","4000","8000")
names(b) <- c("250","500","1000","2000","4000","8000")
audio10 <- rbind(a,b)
y <- data.frame(apply(audio10,2,mean,na.rm=TRUE))
y <- y[,1]

audio20 <- subset(audio,19<`나이`& 30>`나이`)
a <- audio20[,9:14]
b <- audio20[,16:21]
names(a) <- c("250","500","1000","2000","4000","8000")
names(b) <- c("250","500","1000","2000","4000","8000")
audio20 <- rbind(a,b)
z <- data.frame(apply(audio20,2,mean,na.rm=TRUE))
z <- z[,1]

audio30 <- subset(audio,29<`나이`& 40>`나이`)
a <- audio30[,9:14]
b <- audio30[,16:21]
names(a) <- c("250","500","1000","2000","4000","8000")
names(b) <- c("250","500","1000","2000","4000","8000")
audio30 <- rbind(a,b)
k <- data.frame(apply(audio30,2,mean,na.rm=TRUE))
k <- k[,1]

audio70 <- subset(audio,69<`나이`& 80>`나이`)
a <- audio70[,9:14]
b <- audio70[,16:21]
names(a) <- c("250","500","1000","2000","4000","8000")
names(b) <- c("250","500","1000","2000","4000","8000")
audio70 <- rbind(a,b)
l <- data.frame(apply(audio70,2,mean,na.rm=TRUE))
l <- l[,1]


a <- append(x,y)
a <- append(a,z)
a <- append(a,k)
a <- append(a,l)

fre <- c(250,500,1000,2000,4000,8000)
sp <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5)
audio_set <- data.frame(time=rep(sp,1),
                        f=rep(fre,5),
                        t=a)

gg_pta(audio_set,xlab="Frequency (Hz)",ylab="Hearing level (dB)") +
  geom_point(mapping=aes(x = f, y = t, color = factor(time), shape=factor(time)), size = 5, stroke=1.5) +
  geom_line(aes(x = f, y = t, color = factor(time),linetype=factor(time)), lwd=1)+
  scale_color_manual(values=c("red","blue","green","orange","purple"))+
  scale_shape_manual(values=c(1,1,1,1,1))+
  scale_linetype_manual(values = c(1,1,1,1,1))+
  theme(legend.position="none")
ggsave("Right_ver1.png", dpi=1000)

gg_pta(audio_set,xlab="Frequency (Hz)",ylab="Hearing level (dB)") +
  geom_line(aes(x = f, y = t, color = factor(time),linetype=factor(time)), lwd=1)+
  scale_color_manual(values=c("red","blue","green","orange","purple"))+
  scale_linetype_manual(values = c(1,1,1,1,1))+
  theme(legend.position="none")
ggsave("Right_ver2.png", dpi=1000)

####Hearing loss & Age####
audio <- read_excel("D:\\프로젝트\\MPZL2\\MPZL2_Audio.xlsx",sheet=1)
audio <- audio[-33,]
audio <- audio[-54,]
#Low/Mid/High
a <- audio[,c(4,13)]
names(a) <- c("age","dB")
b <- audio[,c(4,14)]
names(b) <- c("age","dB")
c <- audio[,c(4,20)]
names(c) <- c("age","dB")
d <- audio[,c(4,21)]
names(d) <- c("age","dB")
data <- rbind(a,b,c,d)
rm(a,b,c,d)
data <- na.omit(data)

a <- lm(dB~age,data=data)
summary(a)
sqrt(mean(a$residuals^2))

ggplot(data,aes(x = age , y = dB)) +
  geom_point() +
  coord_cartesian(ylim=c(100,10))+
  scale_y_continuous(breaks=seq(100,10,-10))+
  geom_smooth(colour="red",fill="#69B3E7",method = 'lm') +
  theme_classic()

rm(a,data)


#L,M,H age&HL
ca <- c("Low",NA,"Middle",NA,"High")
for(i in c(1,3,5)){
  a <- audio[,c(4,i+8)]
  names(a) <- c("age","dB")
  b <- audio[,c(4,i+15)]
  names(b) <- c("age","dB")
  data1 <- rbind(a,b)
  rm(a,b)
  data1 <- na.omit(data1)
  
  a <- audio[,c(4,i+9)]
  names(a) <- c("age","dB")
  b <- audio[,c(4,i+16)]
  names(b) <- c("age","dB")
  data2 <- rbind(a,b)
  rm(a,b)
  data2 <- na.omit(data2)
  
  data3 <- rbind(data1,data2)
  data3$class <- ca[i]
  
  assign(paste0("x",ca[i]),data3)
  rm(data1,data2,data3)
}

data_fin <- rbind(xHigh,xLow,xMiddle)

COLS = c("skyblue","green","orange")
names(COLS) = c("Low","Middle","High")

ggplot(data_fin,aes(x = age , y = dB)) +
  coord_cartesian(ylim=c(100,10))+
  scale_y_continuous(breaks=seq(100,10,-10))+
  geom_smooth(aes(colour=class,fill=class),method = 'lm') +
  theme_classic()+
  scale_fill_manual(name="class",values=COLS)+
  scale_color_manual(name="class",values=COLS)+
  theme(legend.position = "none")

rm(ca,i,COLS,data,data_fin,xHigh,xLow,xMiddle)

####Variant_country####
data <- read_excel("C:\\Users\\user\\Desktop\\기타\\MPZL2\\MPZL2_variant.xlsx",sheet=2)

variant <- c('c.3G>T','c.52C>T','c.68del','c.72del','c.220C>T','c.280C>T','c.463del')
country_rank <- c('Korea','Chinese','Turkey','Iran','Dutch','Georgian Ukrainian','Moroccan','White','Eastern European Jewish')

data$variant <- factor(data$variant,levels=variant)

color1=c('#FFB3B5','#86D8A4','#d2ffd2','#6CD6f1','#FDCA96','#ffdcff','#ffee66')

ggplot(data,aes(x=country,y=count,fill=factor(`variant`)))+
  scale_x_discrete(limits=country_rank)+
  scale_y_continuous(limits=c(0,50),expand=c(0,0))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color1)+
  theme_classic()+
  theme(legend.position=c(0.8,0.7),legend.title=element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank())



#variant_country2
#ver1
data <- read_excel("C:\\Users\\user\\Desktop\\기타\\MPZL2\\MPZL2_variant.xlsx",sheet=4)

variant <- c('c.3G>T','c.52C>T','c.68del','c.72del','c.220C>T','c.280C>T','c.463del')
country_rank <- c('EA','Other')

data$variant <- factor(data$variant,levels=variant)

color1=c('#FFB3B5','#86D8A4','#d2ffd2','#6CD6f1','#FDCA96','#ffdcff','#ffee66')

ggplot(data,aes(x=country,y=count,fill=factor(`variant`)))+
  scale_x_discrete(limits=country_rank)+
  scale_y_continuous(limits=c(0,85),expand=c(0,0))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color1)+
  theme_classic()+
  theme(legend.position=c(0.8,0.8),legend.title=element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank())

#ver2
x <- data.frame(variant = c('c','c.220c>T','c','c.220c>T'),
                country= c('EA','EA','Other','Other'),
                count = c(10,72,43,3))

color1=c('#aaaaaa','#FDCA96')

ggplot(x,aes(x=country,y=count,fill=factor(`variant`)))+
  scale_x_discrete(limits=country_rank)+
  scale_y_continuous(limits=c(0,85),expand=c(0,0))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color1)+
  theme_classic()+
  theme(legend.position=c(0.8,0.8),legend.title=element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank())

#ver3
x <- data.frame(variant = c('c','c.220c>T','c','c.220c>T'),
                country= c('EA','EA','Other','Other'),
                count = c(12.2,87.8,86.96,13.04))

color1=c('#aaaaaa','#FDCA96')

ggplot(x,aes(x=country,y=count,fill=factor(`variant`)))+
  scale_x_discrete(limits=country_rank)+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color1)+
  theme_classic()+
  theme(legend.position="none",legend.title=element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_blank())

####circos 2####
### You need several libraries
library(chorddiag)

m <- matrix(c(0,0,0,0,0,1,1,3,0,87,0,10,
              0,0,0,0,0,0,0,0,23,3,0,0,
              0,0,0,0,0,0,0,0,4,0,0,0,
              0,0,0,0,0,0,0,0,10,0,2,0,
              0,0,0,0,0,0,0,0,3,0,1,0,
              1,0,0,0,0,0,0,0,0,0,0,0,
              1,0,0,0,0,0,0,0,0,0,0,0,
              3,0,0,0,0,0,0,0,0,0,0,0,
              0,23,4,10,3,0,0,0,0,0,0,0,
              87,3,0,0,0,0,0,0,0,0,0,0,
              0,0,0,2,1,0,0,0,0,0,0,0,
              10,0,0,0,0,0,0,0,0,0,0,0
),
byrow = TRUE,
nrow = 12, ncol = 12)
haircolors <- c("East Asian","Middle Eastern","African","European","White",
                "c.3G>T","c.52C>T","c.68del","c.72del","c.220C>T","c.280C>T","c.463del")
dimnames(m) <- list(have = haircolors,
                    prefer = haircolors)

groupColors <- c("#c0ffff", "#FFDD89", "#957244", "#F26223","#d2ff88",
                 "#00bfff","#9063CD","#FDFDaf","#FF9DFF","#D1F0DC","#568256","#ffe900")
chorddiag(m,groupColors = groupColors, groupnamePadding = 10,groupnameFontsize = 10,showTicks=FALSE)
#groupThickness = 0,
m2 <- matrix(c(1,1,3,0,87,0,10,
               0,0,0,23,3,0,0,
               0,0,0,4,0,0,0,
               0,0,0,10,0,2,0,
               0,0,0,3,0,1,0),
             byrow=TRUE,
             nrow=5,ncol=7)
dimnames(m2) <- list(have = c("East Asian","Middle Eastern","African","European","White"),
                     prefer = c("c.3G>T","c.52C>T","c.68del","c.72del","c.220C>T","c.280C>T","c.463del"))
groupColors <- c("#c0ffff", "#FFDD89", "#957244", "#F26223","#d2ff88")

chorddiag(m2,type="bipartite", groupColors = groupColors, groupnamePadding = 10,groupnameFontsize = 10,showTicks=FALSE,
          categoryNames = FALSE)


####Gene Count_240912####
Genecount <- read_excel("D:\\프로젝트\\MPZL2\\MPZL2_Gene_SNUH_SHU.xlsx",sheet=1)
Genecount <- Genecount[order(Genecount$Gene),]

familycount <- Genecount[,c(1,4)]
familycount <- familycount[order(familycount$Total,decreasing = T),]

names <- familycount$Gene
number <- familycount[!duplicated(familycount$Total),]

pal <- c('57'='#86D8A4','28'='#86D8A4','14'='#FCAEBB','9'='#86D8A4','6'='#86D8A4',
         '4'='#86D8A4','3'='#86D8A4','2'='#86D8A4','1'='#86D8A4')

ggplot(data=familycount,mapping=aes(x=reorder(Gene,`Total`),y=Total))+
  scale_x_discrete(limits=rev(names))+
  scale_y_continuous(limits=c(0,60),expand=c(0,0))+
  geom_bar(mapping=aes(fill=as.factor(Total)),stat="identity")+
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())

#count over 2 
familycount2 <- familycount[familycount$Total>=2,]
names <- familycount2$Gene
number <- familycount2[!duplicated(familycount2$Total),]

pal <- c('57'='#86D8A4','28'='#86D8A4','14'='#FCAEBB','9'='#86D8A4','6'='#86D8A4',
         '4'='#86D8A4','3'='#86D8A4','2'='#86D8A4','1'='#86D8A4')

ggplot(data=familycount2,mapping=aes(x=reorder(Gene,`Total`),y=Total))+
  scale_x_discrete(limits=rev(names))+
  scale_y_continuous(limits=c(0,60),expand=c(0,0))+
  geom_bar(mapping=aes(fill=as.factor(Total)),stat="identity")+
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values=pal)+
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())


rm(names,pal,familycount,familycount2,number,Genecount)


sum(Genecount$Total)


