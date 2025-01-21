# This file contains resources useful for other scripts
# prepare setup

setup_script <-function(relative_path, DIGITS=3, MNTHREADS=48){
  ## set working directory for local and mn5
  cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
  
  machine <- ifelse(Sys.info()[7]=="pclavell", "local",
                    ifelse(grepl("bsc", Sys.info()[7]), "mn5",
                           "ERROR"))
  if(machine=="ERROR"){
    stop("ERROR: User not recognized, this is not your laptop nor your MN5 login")
  }
  cat(paste0("YOU ARE WORKING IN ", machine, "... \n\n"))
  
  if(machine=="local"){
    mn5projects <- "/home/pclavell/mounts/mn5/"
  }else if(machine=="mn5"){
    mn5projects <- "/gpfs/projects/bsc83/"
    .libPaths( c( .libPaths(), "/gpfs/projects/bsc83/utils/Rpackages") )

  }
  
  wd <- paste0(mn5projects, relative_path, "/")
  if(!dir.exists(wd)){dir.create(wd, recursive = T)}
  setwd(wd)
  cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, "... \n\n",sep = "")
  # set non-scientific notation with chosen # of DIGITS
  #options(scipen = 6, digits = DIGITS) # non-scientific notation
  
  ## load up the packages that we always use
  cat("LOADING LIBRARIES", sep = "")
  library(tidyverse)
  library(data.table)
  source(paste0(mn5projects, "Projects/pantranscriptome/pclavell/Z_resources/mytheme.R"))
  # Set Data.table threads depending on machine
  if (machine == "mn5") {
    setDTthreads(threads = MNTHREADS)
  } else if (machine == "local") {
    setDTthreads(threads = 4)
  }
  cat("\n\n")}


# LOAD ARGUMENTS FUNCTION
catch_args <- function(num, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) {
  if (num > 0) {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (num >= 1 && length(args) >= 1) {
      first <- args[1]
      assign(arg1, first, envir = .GlobalEnv)
    }
    if (num >= 2 && length(args) >= 2) {
      second <- args[2]
      assign(arg2, second, envir = .GlobalEnv)
    }
    if (num >= 3 && length(args) >= 3) {
      third <- args[3]
      assign(arg3, third, envir = .GlobalEnv)
    }
    if (num >= 4 && length(args) >= 4) {
      fourth <- args[4]
      assign(arg4, fourth, envir = .GlobalEnv)
    }
    if (num >= 5 && length(args) >= 5) {
      fifth <- args[5]
      assign(arg5, fifth, envir = .GlobalEnv)
    }
    if (num >= 6 && length(args) >= 6) {
      sixth <- args[6]
      assign(arg6, sixth, envir = .GlobalEnv)
    }
    if (num >= 7 && length(args) >= 7) {
      seventh <- args[7]
      assign(arg7, seventh, envir = .GlobalEnv)
    }
    if (num >= 8 && length(args) >= 8) {
      eighth <- args[8]
      assign(arg8, eighth, envir = .GlobalEnv)
    }
    if (num >= 9 && length(args) >= 9) {
      ninth <- args[9]
      assign(arg9, ninth, envir = .GlobalEnv)
    }
    if (num >= 10 && length(args) >= 10) {
      tenth <- args[10]
      assign(arg10, tenth, envir = .GlobalEnv)
    }
  }
}