library(tidyverse)
library(data.table)
library(scales)

data <- data.frame("UMI-length"=12:16,
                   "Number of different UMI"=4^c(12:16),
                   "Edit Distance in UMI-16"=4:0,
                   "1 False Positive Duplicate in"= 20000*4^c(12:16))
data2 <- melt(data, measure.vars= c("Number.of.different.UMI", "X1.False.Positive.Duplicate.in"), variable.name="metric"  )
ggplot(data2)+
  geom_col(aes(x=UMI.length, y=value, fill=metric), position="dodge")+
  scale_fill_manual(values=c("#0c49ab", "#4287f5"), labels=c("Possible different UMI", "Possible UMI-gene pairs"))+
  geom_rect(aes(xmin=11.5,
                xmax=16.5, ymin=3000000, 
                ymax=15000000), fill="grey", alpha=0.1)+
  scale_y_continuous(trans="log10",expand = c(0, 0), 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides="l")+  
  theme_bw()+
  xlab("UMI length - Edit distance")+
  ylab("Combinatorics: Possible UMI*Genes")+
  theme(legend.title= element_blank(),
        legend.position = "top",
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,4,1,1), "lines"))+
  annotate("text", x=17, y=10000000, label="ONT\nthroughput")+
  coord_cartesian(clip = 'off', xlim = c(11.5, 16.5))
  

