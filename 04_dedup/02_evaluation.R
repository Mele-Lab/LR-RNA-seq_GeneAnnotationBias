library(tidyverse)

data2ed <- data.frame("Repeated_UMI"=paste0("x",1:4), "Count"=c(4654791, 28800, 3, 8))
data3ed <- data.frame("Repeated_UMI"=paste0("x",1:4), "Count"=c(4654203, 29388, 3, 8))


ggplot(data3ed)+
  geom_col(aes(x=Repeated_UMI, y=Count, fill=Repeated_UMI))+
  theme_bw()+
  scale_fill_manual(values=c("#329614", "#dbaf0f", "#db720f", "#db200f" ))
