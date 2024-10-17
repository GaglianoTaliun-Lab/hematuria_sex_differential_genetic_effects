# hematuria prevalence by sex and age 

# Load Packages -----------------------------

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(rstatix)

# Arguments ---------------------------------

project_dir = "/Users/frida/Documents/research-projects/col4a2_hematuria"

# Main --------------------------------------

hem <- read.table(here(project_dir, "case.control-wb.withcovariates-4regenie.txt"), sep = "\t", header = T)
age <- read.table(here(project_dir, "f21022.txt"), sep = "\t", header = T) %>%
  dplyr::rename(FID = f.eid, age = f.21022.0.0)

hem_all <- inner_join(hem, age, by = "FID")

subset=hem_all
age_col=16
pheno_col=3
AGE_col=age_col
prev_col=pheno_col
df=subset
qtile=10 #number of quantiles


df_female<-subset(df, df$f.22001.0.0==0)
df_male<-subset(df, df$f.22001.0.0==1)

if (!sum(unique(df[[prev_col]])==c(0,1))==2) {
  print("Column for calculating prevalence of trait must be a binary variable. Expects 0 (controls) and 1 (cases).")
}
if (sum(qtile)<2*length(qtile)){ #check qtile
  print("q-quantiles should be number of divisions for data set and must be greater than 1")
}

##initialize data structures
p<-(100/qtile)/100
index<-c(seq(from=0,to=1,by=p)*100)
prevalences<-rep(NA,qtile+1) #initialize prevalence vector
ns<-rep(NA,qtile+1) #initialize count vector
ses<-rep(NA,qtile+1)#initialize se vector
tiles<-quantile(df[[AGE_col]],seq(from=0,to=1,by=p)) #quantile values
for (i in 1:length(index)-1) {
  print(i)
}
i=1
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalence
i=2
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen

i=3
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen

i=4
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=5
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=6
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=7
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=8
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=9
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=10
prev_list<-subset(df_female, df_female[AGE_col] > tiles[i] & df_female[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen

##create object
pq<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=tiles)
class(pq)<-"prev_quantile_obj"

#1-based indices for your grs and pheno column

pqdf1<-data.frame(prev=pq$prev,se=pq$se,i=pq$i,n=pq$n,tiles=pq$tiles)
pqdf1$Group="female"

i=1
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalence
i=2
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen

i=3
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen

i=4
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=5
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=6
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=7
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=8
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=9
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen
i=10
prev_list<-subset(df_male, df_male[AGE_col] > tiles[i] & df_male[AGE_col] <= tiles[i+1])
prevalences[i]<-sum(prev_list[prev_col])/nrow(prev_list) #how many affected in given quantile
ns[i]<-nrow(prev_list)
ses[i]<-sqrt((prevalences[i]*(1-prevalences[i]))/nrow(prev_list)) #what is SE for this prevalen

##create object
pq<-list(prev=prevalences,se=ses,i=index,n=ns,tiles=tiles)
class(pq)<-"prev_quantile_obj"

#1-based indices for your grs and pheno column

pqdf2<-data.frame(prev=pq$prev,se=pq$se,i=pq$i,n=pq$n,tiles=pq$tiles)
pqdf2$Group<-"male"

pqdfboth<-rbind(pqdf1, pqdf2)

pqdf<-pqdfboth


#pqdf$frac=pqdf$i/100
pqdf$frac=seq(1,i+1,1)
#pqdf<-pqdf[pqdf$frac!=1.00,]
pqdf<-pqdf[pqdf$frac!=i+1,]
pqdf$ub<-pqdf$prev+(1.96*pqdf$se)
pqdf$lb<-pqdf$prev-(1.96*pqdf$se)
#ymax<-max(pqdf$prev) 
ymax<-max(pqdf$ub)
ymin<-min(pqdf$lb)

main<-"Prevalence in Age deciles"
xlab<-"Age decile"
ylab<-"Hematuria prevalence"


#pdf(file="noref_prev_plot.pdf",height=5,width=5,useDingbats=FALSE)
#ggplot(pqdf,aes(x=frac,y=prev,color=as.factor(1))) + geom_point() +
#           scale_color_manual(values=c("grey")) +
#           geom_errorbar(aes(ymin=pqdf$lb,ymax=pqdf$ub),color="grey")  +
#              theme_bw() + labs(title=main) + xlab(xlab) + ylab(ylab)  + 
#                     coord_cartesian(ylim=c(0,ymax)) +
#                     scale_x_discrete(limits=c(pqdf$frac)) +
#                     theme(legend.text=element_text(color = "white"), legend.title = element_text(color = "white"), legend.key = element_rect(fill = "white"))
# dev.off()


pdf(file="female_male_prev_plot.pdf",height=5,width=5,useDingbats=FALSE)
ggplot(pqdf,aes(x=frac,y=prev, group=Group, color=Group)) + geom_point() +
  #scale_color_manual(values=c("grey")) +
  geom_errorbar(aes(ymin=pqdf$lb,ymax=pqdf$ub))  +
  theme_bw() + labs(title=main) + xlab(xlab) + ylab(ylab)  + 
  coord_cartesian(ylim=c(ymin,ymax)) +
  scale_x_discrete(limits=c(pqdf$frac)) +
  theme(legend.text=element_text(color = "black"), legend.title = element_text(color = "black"), legend.key = element_rect(fill = "white"))
dev.off()

## statistics

# t-test to compare means across groups (sexes):
sink(here(project_dir, "age_sex_analyses", "t.test_prev_between_sexes.txt"))
t.test(pqdf$prev[1:10], pqdf$prev[11:20], alternative = "two.sided")
sink()
