library(tidyverse)
library(data.table)
pre <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/softclipping_pre.txt") 
post <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/softclipping_post.txt") 
post_split <- fread("/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/ONT_preprocessing/scripts/softclipping_post_split.txt") 

boxplot(c("pre"=pre$V1, "post"=post$V1, "split"=post_split$V1), outline = F)

summary(pre)
summary(post)
summary(post_split)

df <- data.frame(softclippedbases=c(pre$V1,post$V1, post_split$V1), Group=rep(c("pre","post", "post_split"), times=c(length(pre$V1), length(post$V1), length(post_split$V1))))


theme <- theme_minimal() + theme(axis.ticks = element_line(linewidth = 0.2, color = "black"), axis.text = element_text(size = 11, color="black"),
                                 axis.title = element_text(size=12, vjust = -0.5, color = "black"),
                                 legend.text = element_text(size=12), legend.title = element_text(size = 12, face = "bold"),
                                 panel.border = element_rect(linewidth = 0.4, fill = FALSE), panel.background = element_blank(),   panel.grid = element_line(linewidth =0.2),
                                 strip.text = element_text(size=11),   strip.background = element_blank(),
                                 legend.margin = margin(r = 10, l = 5, t = 5, b = 2),
                                 legend.key.size = unit(15, "pt"))

ggplot(df,aes(y=softclippedbases, fill=Group, x=Group))+
  geom_boxplot(outliers=F)+
  scale_x_discrete(limits=c("pre", "post", "post_split"))+
  theme
