## 0---------------------------------HEADER------------------------------------0
##
## AUTHOR:  Pau Clavell-Revelles
## EMAIL:   pauclavellrevelles@gmail.com
## CREATION DATE: 2024-06-18
## NOTES:
## SETTINGS:  
##    write path relative to 
##    /gpfs/projects/bsc83/ or /home/pclavell/mounts/mn5/
relative_path <- "Projects/pantranscriptome/pclavell"
##
##    load setup sources
mn5source <- "/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/Z_resources/myutils.R"
localsource <-"/home/pclavell/mounts/mn5/Projects/pantranscriptome/pclavell/Z_resources/myutils.R"
ifelse(file.exists(mn5source), source(mn5source), source(localsource))
rm(mn5source, localsource)
##    sets wd, load tidyverse, data.table, sets DT threads
setup_script(relative_path, 3, 48)
##
##    catch arguments from standard input 
##      catch_args(number_of_args, "obj_name1", "obj_name2")
catch_args(0)
##
## 0----------------------------END OF HEADER----------------------------------0

data <- fread("04_transcriptome_assembly/04_evaluation/05_mastertable/data/29102024_PODER_mastertable.tsv")
matched_cols <- grep("^[A-Z]{3}[0-9]", colnames(data), value = TRUE)
long <-melt(data[, .SD, .SDcols = c("isoform", matched_cols)], id.vars = "isoform", value.name = "detected", variable.name = "sample")
long <- long[detected==1]

samples <- unique(long$sample)
results <- c()
for(i in 1:43){
  mysample <- samples[1:i]
  results <- c(results, long[sample%in%mysample, uniqueN(isoform)])
}

mydf <- data.frame("uniqueiso"=results, "n_samples"=1:43)

# Step 2: Plot
ggplot(mydf, aes(x = n_samples, y = uniqueiso)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Samples", y = "Unique Isoforms Discovered") +
  mythemen

# Compute differences
increase <- diff(results)


n_fun <- function(x, y){
  return(data.frame(y = y, label = paste0("n = ",length(x))))
}
median_fun <- function(x,y) {
  return(data.frame(y = y, label =  paste0("median =",round(median(x), 0))))
}

# Convert to dataframe correctly
df <- data.frame("increase"=as.data.frame(increase), "n_samples"=1:42)
ggplot(df, aes(x = n_samples, y = increase)) +
  geom_point() +
  labs(x = "Number of Samples", y = "Unique Isoforms Discovered") +
  mythemen+
  scale_y_continuous(limits =c(1,17000))+
  geom_hline(yintercept=min(df$increase), color="darkgrey", linetype="dashed")+
  annotate(geom="text", x=30, y=7500, label="Each sample adds at least\n329 transcripts")

# fit a line
fit_poly <- lm(df$increase ~ poly(df$n_samples, 2, raw=TRUE))
summary(fit_poly)

ggplot(df, aes(x = "", y = increase)) +
  geom_boxplot()


ggplot(df, 
       aes(x="", y=increase))+
  geom_boxplot()+
  stat_summary(fun.data = n_fun, geom = "text", fun.args = list(y = 0, size=5*0.35)) +
  stat_summary(fun.data = median_fun, geom = "text", fun.args = list(y = -1000), size=6*0.35)+
  mythemen+
  labs(x="", y="Increase in Number of\nDiscovered Transcripts per Sample")



















a_start <- max(mydf$uniqueiso, na.rm = TRUE)  # Max value as asymptote
b_start <- 1 / mean(mydf$n_samples, na.rm = TRUE)  # Rough initial guess

exp_model <- minpack.lm::nlsLM(
  uniqueiso ~ a * (1 - exp(-b * n_samples)), 
  start = list(a = a_start, b = b_start), 
  data = mydf
)


data$fit_exp <- predict(exp_model)

ggplot(mydf, aes(n_samples, uniqueiso)) +
  geom_point() +
  geom_smooth(method = "nls", formula = y ~ a * (1 - exp(-b * x)), 
              method.args = list(start = list(a = a_start, b = b_start)), 
              se = FALSE, color = "red") +
  theme_minimal()
