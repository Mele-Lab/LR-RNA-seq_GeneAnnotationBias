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

data <- fread("../../../Data/GWAScatalog/gwas_catalog_v1.0-associations_e113_r2024-10-21.tsv")


ancestries <- unique(data$`INITIAL SAMPLE SIZE`)

# try to parse as many labels as possible
ancestries[!grepl("Africa|Black|Maya|Asia|Chin|Japan|Turk|Vietn|Korea|Bangladesh|Pakistan|Iranian|Tanzania|Indian|Taiwan|Mexican|American|Saudi|Cameroonian|Qatari|Ugandan|Ghanaian|Aborig|Arab|Malagasy|Ethiopia|Sri Lankan|Peruvian|Thai|Monoglian|Filipino|Punjabi", ancestries) & 
             !grepl("Eur|Iceland|British|Jewish|Finnish|White|Sardinian|Caucasian|Dutch|Amish|Basque|Italian|Spanish|Iberian|German|Russian|Danish|Greek|Hutter|Swedish|Austria|Norweg|French|Irish|Polish", ancestries)]

# with these filters i can tag 95% of all studies ( I won't be tagging Hispanic/Latino because it is very difficult to associate this label with genetics)
noneuropeantags <- ancestries[grepl("Africa|Black|Maya|Asia|Chin|Japan|Turk|Vietn|Korea|Bangladesh|Pakistan|Iranian|Tanzania|Indian|Taiwan|Mexican|American|Saudi|Cameroonian|Qatari|Ugandan|Ghanaian|Aborig|Arab|Malagasy|Ethiopia|Sri Lankan|Peruvian|Thai|Monoglian|Filipino|Punjabi", ancestries)]
europeantags <- ancestries[grepl("Eur|Iceland|British|Jewish|Finnish|White|Sardinian|Caucasian|Dutch|Amish|Basque|Italian|Spanish|Iberian|German|Russian|Danish|Greek|Hutter|Swedish|Austria|Norweg|French|Irish|Polish", ancestries)]

# tag data
data[`INITIAL SAMPLE SIZE`%in%noneuropeantags, ancestry:="Diverse"]
data[`INITIAL SAMPLE SIZE`%in%europeantags, ancestry:="European"]
data$ancestry <- factor(data$ancestry)
