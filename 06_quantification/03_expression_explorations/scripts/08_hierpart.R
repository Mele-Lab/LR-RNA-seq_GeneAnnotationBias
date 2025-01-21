# Tens un parell de plots també!
# 11:34
# merge_re son els valor d'expressió de cada variant (event) cada row es una variant i cada columna un donor
# 11:35
# I en el cas que hagis d'executar Tissue com a random :
# Aquesta es la funció modificada de hier.part per afegir Tissue as random
# hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
#   pcan <- dim(xcan)[2]
#   gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
#   hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
#   names(xcan)[names(xcan)=="Tissue"] <- "(1 | Tissue)"
#   params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")),
#                  family = family, link = link, gof = gof)
#   if(sum(hp$IJ$I<0)>0){
#     NA
#   }else{
#     list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
#   }
# }


individual_traits <- c("Tissue", "Age", "Sex")

n_cores <- 2

# Perform hier.part modified
current.model.mod <- function (y, current.comb, xcan, SStot=0,family = c("gaussian","quasibinomial"), 
                               link = c("logit"), gof = c("Rsqu","RMSPE"), ...)  {
  comb.data <- data.frame(xcan[, current.comb])
  colnames(comb.data) <- colnames(xcan)[current.comb]
  data <- data.frame(y, comb.data)
  depv <- names(data)[1]
  n.comb <- dim(comb.data)[2]
  xs <- vector("character", n.comb)
  for (i in 1:(n.comb - 1)) xs[i] <- paste(names(comb.data)[i], 
                                           "+", sep = "")
  xs[n.comb] <- names(comb.data)[n.comb]
  xss <- paste(xs, collapse = " ", sep = "")
  formu <- stats::formula(paste(depv, "~", xss, sep = ""))
  if (gof == "RMSPE") gf <- sqrt(sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu") {
    if (family == "quasibinomial") 
      gf <- (SStot-sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))/SStot
    if (family == "gaussian") 
      gf <- summary(stats::lm(formu, data = data))$r.squared
  }
  gf
}
all.regs.mod <- function (y, xcan, family = c("gaussian", "quasibinomial"), link = c("logit"), gof = c("Rsqu","RMSPE"),...) { 
  if (sum(is.na(xcan)) > 0) {
    missing <- is.na(apply(xcan, 1, FUN = sum))
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in xcan\n"), 
            call. = FALSE)
  }
  if (sum(is.na(y)) > 0) {
    missing <- is.na(y)
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in y\n"), 
            call. = FALSE)
  }
  pcan <- dim(xcan)[2]
  n <- (2^pcan) - 1
  combs <- combos1(pcan)$ragged
  SStot <- sum((y-mean(y))^2)
  
  if (gof == "RMSPE")  gfs <- sqrt(sum((stats::glm(y ~ 1, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu")   gfs <- 0
  
  for (i in 1:n) {
    if (i%%500 == 0) 
      cat(i, "regressions calculated:", n - i, "to go...\n")
    current.comb <- as.vector(combs[i, ][combs[i, ] > 0]) 
    combn <- paste(names(data.frame(xcan)[current.comb]), "", collapse = "")
    if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    gfs <- c(gfs, new.line)
  }
  gfs
}

hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
  pcan <- dim(xcan)[2]
  gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
  hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
  
  params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                 family = family, link = link, gof = gof)
  if(sum(hp$IJ$I<0)>0){
    NA
  }else{
    list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
  }
}

metadata_info <- metadata_info %>%
  filter(ID != "ENCSR301QXH_2")

hier.part.results <- mclapply(rownames(merged_re), function(event)
  hier.part.mod(y=as.numeric(merged_re[event,]), x=metadata_info[,individual_traits],
                fam = "quasibinomial", link = "logit", gof = "Rsqu",control = list(maxit = 100)), mc.cores = n_cores)
names(hier.part.results) <- rownames(merged_re)


### JUST IN CASE YOU HAVE ERRORS RUNNING THE COMMAND ABOVE ###

# Initialize a list to store results
hier.part.results <- list()

# Iterate over each row in merged_re
for (event in rownames(merged_re)) {
  cat("Processing event:", event, "\n")  # Print progress
  
  # Use tryCatch to handle potential errors
  hier.part.results[[event]] <- tryCatch({
    # Run the function for the current row
    hier.part.mod(
      y = as.numeric(merged_re[event, ]),
      x = metadata_info[, individual_traits],
      fam = "quasibinomial",
      link = "logit",
      gof = "Rsqu",
      control = list(maxit = 100)
    )
  }, error = function(e) {
    # Log the error and continue
    cat("Error in event:", event, "\n")
    print(e)
    NULL  # Store NULL for events that cause errors
  })
}

# Check which events returned NA
failed_events <- names(hier.part.results)[sapply(hier.part.results, function(x) any(is.na(x)))]

cat("Failed events (returned NA):", failed_events, "\n")

hier.part.results <- hier.part.results[!sapply(hier.part.results, function(x) any(is.na(x)))]

common_names <- intersect(rownames(merged_re), names(hier.part.results))

#### UNTIL  HERE  ###

# If you were able to run it properly with the first command then substitute 
# common_names by rownames() and the dataframe where you store the expression values

rsq <- sapply(common_names, function(event) 
  sum(hier.part.results[[event]]$IJ[,1])) 
names(rsq) <- common_names
rel_perc <- do.call(rbind.data.frame,
                    lapply(names(hier.part.results), function(event) 
                      as.numeric(unlist(hier.part.results[[event]]$I.perc))))
rownames(rel_perc) <- common_names
colnames(rel_perc) <- individual_traits
abs_perc <-  do.call(rbind.data.frame,
                     lapply(names(hier.part.results), function(event) 
                       hier.part.results[[event]]$IJ[,1])
)
rownames(abs_perc) <- common_names
colnames(abs_perc) <- individual_traits
hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
colnames(hier_data) <- c("R2",paste0(individual_traits,"_rel"),paste0(individual_traits,"_abs"))


# Barplot 
# Compute means for `_rel` columns
rel_means <- hier_data %>%
  select(ends_with("_rel")) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

# Reshape data from wide to long format
rel_means_long <- rel_means %>%
  pivot_longer(cols = everything(), names_to = "Trait", values_to = "Mean")

rel_means_long$Trait <- c("Tissue", "Age", "Sex")

# Create the stacked barplot
ggplot(rel_means_long, aes(x = 1, y = Mean, fill = Trait)) +
  geom_bar(stat = "identity") +  # Stacked barplot
  scale_x_continuous(breaks = NULL) +  # Remove x-axis labels
  labs(x = NULL, y = "Mean", title = "Variance explained by each Trait") +
  theme_bw() +
  theme(legend.title = element_blank())


# Violin and Boxplot
abs_values <- hier_data %>%
  select(ends_with("_abs")) 

colnames(abs_values) <- c("Tissue", "Age", "Sex")

# Reshape data from wide to long format
abs_values_long <- abs_values %>%
  pivot_longer(cols = everything(), names_to = "Trait", values_to = "Absolute")

ggplot(abs_values_long, aes(x = Trait, y = Absolute)) +
  geom_boxplot() +  # Boxplot layer
  theme_bw() +  # Clean theme
  ylab("Absolute Value") +  # y-axis label
  xlab("Trait")  # x-axis label