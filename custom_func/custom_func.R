#'---------------------------------------------------------------#'
#' Functions used in the methodological paper 
#' Quantitative synthesis: a practical guide on meta-analysis, meta-regression, and publication bias tests for environmental sciences
#' Oct 2022
#'---------------------------------------------------------------#'


############################################################################
#' @title cv_avg
#' @description Calculates the weighted average CV^2 within a study and the weighted average CV^2 across a study
cv_avg <- function(x, sd, n, group, data, label = NULL, sub_b = TRUE, cv2=FALSE){

  # Check if the name is specified or not. If not, then assign it the name of the mean, x, variable input in the function. https://stackoverflow.com/questions/60644445/converting-tidyeval-arguments-to-string
  if(is.null(label)){
    label <- purrr::map_chr(enquos(x), rlang::as_label)
  }

  # Calculate between study CV. Take weighted mean CV within study, and then take a weighted mean across studies of the within study CV. Weighted based on sample size and pooled sample size.
  b_grp_cv_data <- data                                             %>%
    dplyr::group_by({{group}})                            %>%
    dplyr::mutate(   w_CV2 = weighted_CV({{sd}}, {{x}}, {{n}}, cv2=cv2),
                     n_mean = mean({{n}}, na.rm = TRUE))   %>%
    dplyr::ungroup(.)                                     %>%
    dplyr::mutate(b_CV2 = weighted.mean(w_CV2, n_mean, na.rm = TRUE), .keep = "used")

  # Make sure that label of the calculated columns is distinct from any other columns
  names(b_grp_cv_data) <- paste0(names(b_grp_cv_data), "_", label)

  # Append these calculated columns back to the original data and return the full dataset.
  if(sub_b){
    b_grp_cv_data <- b_grp_cv_data %>% dplyr::select(grep("b_", names(b_grp_cv_data)))
    dat_new <- cbind(data, b_grp_cv_data)
  } else {
    dat_new <- cbind(data, b_grp_cv_data)
  }

  return(data.frame(dat_new))
}

############################################################################

#' @title weighted_CV
#' @description Calculates the weighted average CV^2 or CV followed by squaring within a study and the weighted averages CV^2 across a studies
weighted_CV <- function(sd, x, n, cv2=FALSE){
  if(cv2){
    weighted.mean(na_if((sd / x)^2, Inf), n, na.rm = TRUE)
  }else{
    weighted.mean(na_if((sd / x), Inf), n, na.rm = TRUE)^2
  }
}

############################################################################

#' @title get_est
#' @description Extracts estimates from rma.mv and rma model objects
get_est <- function(model){
  est <- coef(model)
  ci.lb <- model$ci.lb
  ci.ub <- model$ci.ub
  se <- model$se
  return(data.frame(Est. = est, SE=se, "95% LCI" = ci.lb, "95% UCI" = ci.ub, check.names = FALSE))
}

############################################################################

#' @title lnrr_laj
#' @description Calculates log response ratio based on Taylor expansion from Jajeunesse 2011
lnrr_laj <- function(m1, m2, cv1_2, cv2_2, n1, n2){
  log(m1 / m2) + 0.5*((cv1_2 / n1) - (cv2_2 / n2))
}

############################################################################

#' @title v_lnrr_laj
#' @description Calculates the sampling variance for log response ratio based on second order Taylor expansion proposed by Lajeunesse 2011
v_lnrr_laj <- function(cv1_2, cv2_2, n1, n2){
  ((cv1_2) / n1) + ((cv2_2) / n2) +
    ((cv1_2)^2 / (2*n1)^2) + ((cv2_2)^2 / (2*n2)^2)
}

############################################################################

#' @title gen.miss
#' @description Generates random missing data in two columns of a meta-analytic dataset
gen.miss <- function(data, missVar, missCol2, n_miss){
  data[sample(rownames(data), n_miss), missVar] <- NA
  data[is.na(data[,missVar]), missCol2] <- NA
  return(data)
}

############################################################################
#' @title estimates.model
#' @description Extract mean and CI from each metafor model
estimates.model <- function(model){
  db.mf <- data.frame(round(model$b, 3),row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,round(model$ci.lb, 3),round(model$ci.ub,3),row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}

############################################################################

mod_results <- function(model, mod = "1", group, data, N = NULL,  weights = "prop", by = NULL, at = NULL, subset = FALSE, upper = TRUE, ...){
  
  if(any(grepl("-1|0", as.character(model$formula.mods)))){
    warning("It is recommended that you fit the model with an intercept. Unanticipated errors can occur otherwise.")
  }
  
  if(missing(model)){
    stop("Please specify the 'model' argument by providing rma.mv or rma model object. See ?mod_results")
  }
  
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv, rma, or robust.rma")}
  
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  
  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  
  if(is.null(stats::formula(model))){
    #model <- stats::update(model, "~1")
    model$formula.mods <- ~ 1
    dat_tmp <- data$`1` <- "Intrcpt"
    model$data <- dat_tmp
  } else{
    model$data <- data
  }
  
  if(model$test == "t"){
    df_mod = as.numeric(model$ddf[[1]])
  } else{
    df_mod = 1.0e6 # almost identical to z value
  }
  
  if(is.character(data[[mod]]) | is.factor(data[[mod]]) | is.null(data[[mod]])) {
    grid <- emmeans::qdrg(formula = stats::formula(model), at = at, data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1) ## NOTE: Added data argument emmeans >vers 1.7.4. Object is unstable so feeding in the relevant arguments from model object directly. Note, we should think about df!
    mm <- emmeans::emmeans(grid, specs = mod, df = df_mod, by = by, weights = weights, ...)
    
    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    
    if(is.null(by)){
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
      
    } else{
      mod_table <- data.frame(name = firstup(as.character(mm_pi[,1]), upper = upper),
                              condition = mm_pi[,2], estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }
    
    # Extract data
    data2 <- get_data_raw(model, mod, group, N, data, at = at, subset)
    
    mod_table$name <- factor(mod_table$name,
                             levels = mod_table$name,
                             labels = mod_table$name)
    
  } else{
    at2 <- list(mod = seq(min(data[,mod], na.rm = TRUE), max(data[,mod], na.rm = TRUE), length.out = 100))
    names(at2) <- mod
    grid <- emmeans::qdrg(formula =  stats::formula(model), data = model$data, coef = model$b,
                          vcov = stats::vcov(model), df = model$k-1, at = c(at2, at))  # getting 100 points. Fixing this to make it more general
    mm <- emmeans::emmeans(grid, specs = mod, by = c(mod, by), weights = weights, df = df_mod)
    
    # getting prediction intervals
    mm_pi <- pred_interval_esmeans(model, mm, mod = mod)
    
    if(is.null(by)){
      mod_table <- data.frame(moderator = mm_pi[,1],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    } else{
      mod_table <- data.frame(moderator = mm_pi[,1],
                              condition = mm_pi[,2],
                              estimate = mm_pi[,"emmean"],
                              lowerCL = mm_pi[,"lower.CL"],
                              upperCL = mm_pi[,"upper.CL"],
                              lowerPR = mm_pi[,"lower.PI"],
                              upperPR = mm_pi[,"upper.PI"])
    }
    
    # extract data
    data2 <- get_data_raw_cont(model, mod, group, N, data, by = by)
    
  }
  
  
  output <- list(mod_table = mod_table,
                 data = data2)
  
  class(output) <- c("orchard", "data.frame")
  
  return(output)
}




############# Key Sub-functions #############

#' @title pred_interval_esmeans
#' @description Function to get prediction intervals (credibility intervals) from \code{esmeans} objects (\pkg{metafor}).
#' @param model \code{rma.mv} object.
#' @param mm result from \code{emmeans::emmeans} object.
#' @param mod Moderator of interest.
#' @param ... other arguments passed to function.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export


pred_interval_esmeans <- function(model, mm, mod, ...){
  
  tmp <- summary(mm)
  tmp <- tmp[ , ]
  test.stat <- stats::qt(0.975, tmp$df[[1]])
  
  if(length(model$tau2) <= 1){ # including gamma2
    sigmas <- sum(model$sigma2)
    PI <- test.stat * base::sqrt(tmp$SE^2 + sigmas)
  } else {
    sigmas <- sum(model$sigma2)
    taus   <- model$tau2
    gammas <- model$gamma2
    w_tau <- model$g.levels.k
    w_gamma <- model$g.levels.k
    
    if(mod == "1"){
      tau <- weighted_var(taus, weights = w_tau)
      gamma <- weighted_var(gamma, weights = w_gamma)
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + tau + gamma)
      
    } else {
      PI <- test.stat * sqrt(tmp$SE^2 + sigmas + taus + gammas)
    }
  }
  
  tmp$lower.PI <- tmp$emmean - PI
  tmp$upper.PI <- tmp$emmean + PI
  
  # renaming "overall" to ""
  if(tmp[1,1] == "overall"){tmp[,1] <- "intrcpt"}
  
  return(tmp)
}

#' @title get_data_raw
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor}.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or whatever other grouping variable one wishes to present sample sizes.
#' @param data The data frame used to fit the \code{rma.mv} model object.
#' @param N The name of the column in the data specifying the sample size, N. Defaults to \code{NULL}, so precision is plotted instead of sample size.
#' @param at List of moderators. If \code{at} is equal to \code{mod} then levels specified within \code{at} will be used to subset levels when \code{subset = TRUE}. Otherwise, it will marginalise over the moderators at the specified levels.
#' @param subset Whether or not to subset levels within the \code{mod} argument. Defaults to \code{FALSE}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#' @examples \dontrun{
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#'  test <- get_data_raw(model, mod = "trait.type", group = "group_ID", data = warm_dat, at = list(trait.type = c("physiology", "morphology")))
#'  test2 <- get_data_raw(model, mod = "1", group = "group_ID", data = warm_dat)
#'
#'  data(english)
#'  # We need to calculate the effect sizes, in this case d
#'  english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"), data = english)
#'  model <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#'  test3 <-  get_data_raw(model, mod = "1", group = "StudyNo", data = english)}

get_data_raw <- function(model, mod, group, N = NULL, data, at = NULL, subset = TRUE){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  # Extract data
  # Check first if missing data exists
  if(length(attr(model$X, "dimnames")[[1]]) > 0){
    # full model delete missing values so need to adjust
    position <- attr(model$X, "dimnames")[[1]]
    data <- data[position, ] }
  if(!is.null(at) & subset){
    # Find the at slot in list that pertains to the moderator and extract levels
    at_mod <- at[[mod]]
    position2 <- which(data[,mod] %in% at_mod)
    # Subset the data to only the levels in the moderator
    data <- data[position2,]
    yi <- model$yi[position2]
    vi <- model$vi[position2]
    type <- attr(model$yi, "measure")
  } else {
    # Extract effect sizes
    yi <- model$yi
    vi <- model$vi
    type <- attr(model$yi, "measure")
  }
  if(mod == "1"){
    moderator <- "Intrcpt"
  }else{
    # Get moderator
    moderator <- as.character(data[[mod]]) # Could default to base instead of tidy
    moderator <- firstup(moderator)
  }
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, stdy, type)
  #names(data_reorg)[4] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)
  
  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }
  
  return(data_reorg)
}

#' @title get_data_raw_cont
#' @description Collects and builds the data used to fit the \code{rma.mv} or \code{rma} model in \pkg{metafor} when a continuous variable is fit within a model object.
#' @param model \code{rma.mv} object.
#' @param mod the moderator variable.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes.
#' @param N  The name of the column in the data specifying the sample size, N. Defaults to \code{NULL} so that precision is plotted instead of sample size.
#' @param data The data frame used to fit the \code{rma.mv} model object.
#' @param by Character name(s) of the 'condition' variables to use for grouping into separate tables.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export

#TODO what if there is no "by"

get_data_raw_cont <- function(model, mod, group, N = NULL, data, by){
  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?mod_results")
  }
  if(missing(data)){
    stop("Please specify the 'data' argument by providing the data used to fit the model. See ?mod_results")
  }
  # Extract data
  # Check first if missing data exists
  if(length(attr(model$X, "dimnames")[[1]]) > 0){
    # full model delete missing values so need to adjust
    position <- attr(model$X, "dimnames")[[1]]
    data <- data[position, ] }
  # Extract effect sizes
  yi <- model$yi
  vi <- model$vi
  type <- attr(model$yi, "measure")
  # Get moderator
  moderator <- data[[mod]] # Could default to base instead of tidy
  #names(moderator) <  "moderator"
  if(is.null(by)){
    condition <- data[ , by]
  }else{
    condition <- data[[by]]
  }
  #names(condition) <  "condition"
  # Extract study grouping variable to calculate the
  stdy <- data[[group]] # Could default to base instead of tidy
  data_reorg <- data.frame(yi, vi, moderator, condition, stdy, type)
  # if(!is.na(names(data_reorg)[names(data_reorg) == by]) == TRUE) {  ## FAILING HERE
  #   names(data_reorg)[names(data_reorg) == by] <- "condition"
  # }
  #names(data_reorg)[5] <- "stdy" # sometimes stdy gets replaced by group's names
  row.names(data_reorg) <- 1:nrow(data_reorg)
  
  if(is.null(N) == FALSE){
    data_reorg$N <- data[ ,N]
  }
  
  return(data_reorg)
}

############# Helper-functions #############

#' @title firstup
#' @description Uppercase moderator names
#' @param x a character string
#' @param upper logical indicating if the first letter of the character string should be capitalized. Defaults to TRUE.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a character string with all combinations of the moderator level names with upper case first letters
#' @export

firstup <- function(x, upper = TRUE) {
  if(upper){
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  } else{ x }
}


#' @title print.orchard
#' @description Print method for class 'orchard'
#' @param x an R object of class orchard
#' @param ... Other arguments passed to print
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#'

print.orchard <- function(x, ...){
  return(print(x$mod_table))
}

#' @title weighted_var
#' @description Calculate weighted variance
#' @param x A vector of tau2s to be averaged
#' @param weights Weights, or sample sizes, used to average the variance
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a vector with a single weighted variance
#' @export
#'

weighted_var <- function(x, weights){
  weight_var <- sum(x * weights) / sum(weights)
  return(weight_var)
}


#' @title num_studies
#' @description Computes how many studies are in each level of categorical moderators of a \code{rma.mv} model object.
#' @param mod Character string describing the moderator of interest.
#' @param data Raw data from object of class "orchard"
#' @param group A character string specifying the column name of the study ID grouping variable.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a table with the number of studies in each level of all parameters within a \code{rma.mv} or \code{rma} object.
#' @export
#' @examples
#' \dontrun{data(fish)
#'warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list( ~1 | es_ID,~1 | group_ID), mods = ~experimental_design-1, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' num_studies(model$data, experimental_design, group_ID)
#' }

num_studies <- function(data, mod, group){
  
  # Summarize the number of studies within each level of moderator
  table <- data        %>%
    dplyr::group_by({{mod}}) %>%
    dplyr::summarise(stdy = length(unique({{group}})))
  
  table <- table[!is.na(table$moderator),]
  # Rename, and return
  colnames(table) <- c("Parameter", "Num_Studies")
  return(data.frame(table))
  
}

############################################################################
#' @title orchard_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, it creates an orchard plot from mean effect size estimates for all levels of a given categorical moderator, and their corresponding confidence and prediction intervals.
#' @param object model object of class \code{rma.mv}, \code{rma}, or \code{orchard} table of model results.
#' @param mod the name of a moderator. Defaults to \code{"1"} for an intercept-only model. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param data The data frame used to fit the \code{rma.mv} model object.  Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding varaibles in 'by'. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param xlab The effect size measure label.
#' @param N The name of the column in the data specifying the sample size, N. Defaults to \code{NULL}, so that precision is plotted instead of sample size.
#' @param alpha The level of transparency for effect sizes represented in the orchard plot.
#' @param angle The angle of y labels. The default is 90 degrees.
#' @param cb If \code{TRUE}, it uses 20 colour blind friendly colors.
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD). Defaults to \code{"none"}.
#' @param condition.lab Label for the condition being marginalized over.
#' @param trunk.size Size of the mean, or central point.
#' @param branch.size Size of the confidence intervals.
#' @param twig.size Size of the prediction intervals.
#' @param legend.pos Where to place the legend. To remove the legend, use \code{legend.pos = "none"}.
#' @param k.pos Where to put k (number of effect sizes) on the plot.
#' @param colour Colour of effect size shapes. By default, effect sizes are colored according to the \code{mod} argument. If \code{TRUE}, they are colored according to the grouping variable
#' @param fill If \code{TRUE}, effect sizes will be filled with colours. If \code{FALSE}, they will not be filled with colours.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param upper Logical, defaults to \code{TRUE}, indicating that the first letter of the character string for the moderator variable should be capitalized.
#' @return Orchard plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au

orchard_plot <- function(object, mod = "1", group, data, xlab, N = NULL,
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE, g = TRUE,
                         trunk.size = 3, branch.size = 1.2, twig.size = 0.5,
                         transfm = c("none", "tanh"), condition.lab = "Condition",
                         legend.pos = c("bottom.right", "bottom.left",
                                        "top.right", "top.left",
                                        "top.out", "bottom.out",
                                        "none"), # "none" - no legends
                         k.pos = c("right", "left", "none"),
                         colour = FALSE,
                         fill = TRUE,
                         weights = "prop", by = NULL, at = NULL, upper = TRUE)
{
  ## evaluate choices, if not specified it takes the first choice
     transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
       k.pos <- match.arg(NULL, choices = k.pos)

	if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){

	    if(mod != "1"){
	    #results <-  orchaRd::mod_results(object, mod, group, data, N, by = by, at = at, weights = weights, upper = upper)
	    # remove orchaRd:: in front of mod_results as it introduce extra trouble!!!!
	    results <-  mod_results(object, mod, group, data, N,
	                                        by = by, at = at, weights = weights, upper = upper)
	  } else {
	    # results <-  orchaRd::mod_results(object, mod = "1", group, data, N,by = by, at = at, weights = weights, upper = upper)
	    results <-  mod_results(object, mod = "1", group, data, N,
	                                        by = by, at = at, weights = weights, upper = upper)
	    }
	  }

	if(any(class(object) %in% c("orchard"))) {
			results <- object
	}

	mod_table <- results$mod_table

  data_trim <- results$data
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)

  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
	legend <- "Precision (1/SE)"

	if(any(N != "none")){
	  data_trim$scale <- data_trim$N
		  legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
		  #latex2exp::TeX()
	}

	if(transfm == "tanh"){
		                   cols <- sapply(mod_table, is.numeric)
		mod_table[,cols] <- Zr_to_r(mod_table[,cols])
		data_trim$yi <- Zr_to_r(data_trim$yi)
		                  label <- xlab
	}else{
		label <- xlab
	}

	# Add in total effect sizes for each level
	 mod_table$K <- as.vector(by(data_trim, data_trim[,"moderator"], function(x) length(x[,"yi"])))

	# Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
	 mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])

	 # the number of groups in a moderator & data points
	 group_no <- length(unique(mod_table[, "name"]))

	 #data_no <- nrow(data)

	# colour blind friendly colours with grey
	 cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

	 # setting fruit colour
	 if(colour == TRUE){
	   color <- as.factor(data_trim$stdy)
	   color2 <- NULL
	 }else{
	   color <- data_trim$mod
	   color2 <- mod_table$name
	 }

	 # whether we fill fruit or not
	 if(fill == TRUE){
	   fill <- color
	 }else{
	     fill <- NULL
	   }

	 # whether marginal
	 if(names(mod_table)[2] == "condition"){

	   # the number of levels in the condition
	   condition_no <- length(unique(mod_table[, "condition"]))

	   plot <- ggplot2::ggplot() +
	     # pieces of fruit (bee-swarm and bubbles)
	     ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +

	     ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
	     # creating CI
	     ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL),
	                             size = branch.size, position = ggplot2::position_dodge2(width = 0.3)) +
	     # drowning point estimate and PI
	     ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name, ymin = lowerPR, ymax = upperPR,  shape = as.factor(condition), fill = color2), size = twig.size, position = ggplot2::position_dodge2(width = 0.3), fatten = trunk.size) +
	     # this will only work for up to 5 different conditions
	     # flipping things around (I guess we could do use the same geoms but the below is the original so we should not change)
	     ggplot2::scale_shape_manual(values =  20 + (1:condition_no)) + ggplot2::coord_flip() +
	     ggplot2::theme_bw() +
	     ggplot2::guides(fill = "none", colour = "none") +
	     ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
	     ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
	     ggplot2::theme(legend.direction="horizontal") +
	     ggplot2::theme(legend.background = ggplot2::element_blank()) +
	     ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
	     ggplot2::labs(shape = condition.lab) +
	     ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
	                                                        hjust = 0.5,
	                                                        angle = angle))

	 } else {

	  plot <- ggplot2::ggplot() +
	    # pieces of fruit (bee-swarm and bubbles)
	    ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = color, fill = fill), alpha=alpha, shape = 21) +

	    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = alpha) +
	    # creating CI
	    ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size) +
	    # drowning point estimate and PI
	    ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name,  ymin = lowerPR, ymax = upperPR, fill = color2), size = twig.size, fatten = trunk.size, shape = 21) +
	    ggplot2::coord_flip() +
	    ggplot2::theme_bw() +
	    ggplot2::guides(fill = "none", colour = "none") +
	    ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
	    ggplot2::theme(legend.direction="horizontal") +
	    ggplot2::theme(legend.background = ggplot2::element_blank()) +
	    ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
	    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
	                                                       hjust = 0.5,
	                                                       angle = angle)) #+
	    #ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))

	 }

	   # adding legend
	 if(legend.pos == "bottom.right"){
	   plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
	 } else if ( legend.pos == "bottom.left") {
	   plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
	 } else if ( legend.pos == "top.right") {
	   plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
	 } else if (legend.pos == "top.left") {
	   plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
	 } else if (legend.pos == "top.out") {
	   plot <- plot + ggplot2::theme(legend.position="top")
	 } else if (legend.pos == "bottom.out") {
	   plot <- plot + ggplot2::theme(legend.position="bottom")
	 } else if (legend.pos == "none") {
	   plot <- plot + ggplot2::theme(legend.position="none")
	 }

	  # putting colors in
	  if(cb == TRUE){
	    plot <- plot +
	      ggplot2::scale_fill_manual(values = cbpl) +
	      ggplot2::scale_colour_manual(values = cbpl)
	  }

	  # putting k and g in
	  if(k == TRUE && g == FALSE && k.pos == "right"){
	    plot <- plot +
	      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                        label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = 3.5)
	  } else if(k == TRUE && g == FALSE && k.pos == "left") {
	    plot <- plot +  ggplot2::annotate('text', y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                                      label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "left", size = 3.5)
	  } else if (k == TRUE && g == TRUE && k.pos == "right"){
	    # get group numbers for moderator
	    plot <- plot + ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                        label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
	                        parse = TRUE, hjust = "right", size = 3.5)
	  } else if (k == TRUE && g == TRUE && k.pos == "left"){
	    # get group numbers for moderator
	    plot <- plot + ggplot2::annotate('text',  y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
	                        label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
	                        parse = TRUE, hjust = "left", size = 3.5)
	  }


	  return(plot)
}


#' @title Zr_to_r
#' @description Converts Zr back to r (Pearson's correlation coefficient)
#' @param df data frame of results of class 'orchard'
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals with estimates converted back to r
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

Zr_to_r <- function(df){
	return(sapply(df, tanh))
}
############################################################################
#' @title bubble_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, the \code{bubble_plot} function creates a bubble plot from slope estimates. In cases when a model includes interaction terms, this function creates panels of bubble plots.
#' @param object Model object of class \code{rma}, \code{rma.mv}, or \code{orchard} table of model results
#' @param mod The name of a continuous moderator, to be plotted on the x-axis of the bubble plot.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param data The data frame used to fit the \code{rma.mv} model object. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding variables in \code{by}. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights How to marginalize categorical variables; used when one wants marginalised means. The default is \code{weights = "prop"}, which weights means for moderator levels based on their proportional representation in the data. For example, if \code{"sex"} is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal when, for example, males and females are typically roughly equally prevalent in a population. In cases such as these, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param xlab Moderator label.
#' @param ylab Effect size measure label.
#' @param k.pos The position of effect size number, k.
#' @param N The vector of sample size which an effect size is based on. Defaults to precision (the inverse of sampling standard error).
#' @param alpha The level of transparency for pieces of fruit (effect size).
#' @param cb If \code{TRUE}, it uses a colourblind-friendly palette of 20 colours (do not make this \code{TRUE}, when colour = \code{TRUE}).
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param est.lwd Size of the point estimate.
#' @param ci.lwd Size of the confidence interval.
#' @param pi.lwd Size of the prediction interval.
#' @param est.col Colour of the point estimate.
#' @param ci.col Colour of the confidence interval.
#' @param pi.col Colour of the prediction interval.
#' @param condition.nrow Number of rows to plot condition variable.
#' @param legend.pos Where to place the legend, or not to include a legend ("none").
#'
#' @return Orchard plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(lim)
#' lim[, "year"] <- as.numeric(lim$year)
#' lim$vi<- 1/(lim$N - 3)
#' model<-metafor::rma.mv(yi=yi, V=vi, mods= ~Environment*year,
#' random=list(~1|Article,~1|Datapoint), data=na.omit(lim))
#' test <- orchaRd::mod_results(model, mod = "year", group = "Article", data = lim, weights = "prop", by = "Environment")
#' orchaRd::bubble_plot(test, mod = "year", data = lim, group = "Article",legend.pos = "top.left")
#' # Or just using model directly
#' orchaRd::bubble_plot(model, mod = "year", legend.pos = "top.left", data = lim, group = "Article", weights = "prop", by = "Environment")
#'
#' }
#' @export


bubble_plot <- function(object, mod, group = NULL, data,
                        xlab = "Moderator", ylab = "Effect size", N = "none",
                        alpha = 0.5, cb = TRUE, k = TRUE, g = FALSE,
                        est.lwd = 1, ci.lwd = 0.5, pi.lwd = 0.5,
                        est.col = "black", ci.col = "black", pi.col = "black",
                        legend.pos = c("top.left", "top.right",
                                       "bottom.right", "bottom.left",
                                       "top.out", "bottom.out",
                                       "none"),
                        k.pos = c("top.right", "top.left",
                                  "bottom.right", "bottom.left",
                                  "none"),
                        condition.nrow = 2,
                         #condition.lab = "Condition",
                        weights = "prop", by = NULL, at = NULL)
  {
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  #facet <- match.arg(NULL, choices = facet)

  if(missing(data)){
         stop("Please specify the 'data' argument by providing the data used to fit the model. See ?bubble_plot")
  }

  if(any(grepl(mod, colnames(data))) == FALSE){
    error("The moderator specified is not found in your data. Did you transform the variable in the model and forget to add it to your dataframe?")
  }

  if(missing(group)){
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?bubble_plot")
  }

  if(is.numeric(by)){
   k = FALSE
   g = FALSE
  }


  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){

    if(mod != "1"){
      # change orchaRd::mod_results to mod_results to avide trouble of R version
      results <-  mod_results(object, mod, group, data,
                                       by = by, at = at, weights = weights)
    } else {
      results <-  mod_results(object, mod = "1", group, data,
                                       by = by, at = at, weights = weights)
    }
  }

  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }

  mod_table <- results$mod_table

  data_trim <- results$data
  #data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)

  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision (1/SE)"

  if(any(N != "none")){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
  }

  label <- xlab
  # if(transfm == "tanh"){
  #   cols <- sapply(mod_table, is.numeric)
  #   mod_table[,cols] <- Zr_to_r(mod_table[,cols])
  #   data_trim$yi <- Zr_to_r(data_trim$yi)
  #   label <- xlab
  # }else{
  #   label <- xlab
  # }

  if(is.null(data_trim$condition) == TRUE){
  # the number of effect sizes
  effect_num <- nrow(data_trim)

  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  group_num <- length(unique(data_trim$stdy))

  dat_text <- data.frame(K = effect_num, G = group_num)

  }else{
  effect_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(x[,"yi"])))

  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  #group_num <- c(2,4)
  group_num <- as.vector(by(data_trim, data_trim[,"condition"], function(x) base::length(base::unique(x[,"stdy"]))))

  dat_text <- data.frame(K = effect_num, G = group_num, condition = as.vector(base::unique(data_trim$condition)))
  }
  # the number of groups in a moderator & data points
  #group_no <- length(unique(mod_table[, "name"]))

  #data_no <- nrow(data)

  # # colour blind friendly colours with grey
  # cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

  if(is.null(data_trim$condition) == TRUE){
   plot <-ggplot2::ggplot() +
    # putting bubbles
     ggplot2::geom_point(data = data_trim, ggplot2::aes(x = moderator, y = yi, size = scale), shape = 21, alpha = alpha, fill = "grey90" ) +
    # prediction interval
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x, se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
     # confidence interval
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x, se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x, se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
     # main line
     ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
    #facet_grid(rows = vars(condition)) +
     ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
     ggplot2::guides(fill = "none", colour = "none") +
    # themes
     ggplot2::theme_bw() +
    #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
     ggplot2::theme(legend.direction="horizontal") +
    #theme(legend.background = element_rect(fill = "white", colour = "black")) +
     ggplot2::theme(legend.background = ggplot2::element_blank()) +
     ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))
  } else if(is.character(data_trim$condition) == TRUE || is.factor(data_trim$condition) == TRUE){
    plot <-ggplot2::ggplot() +
      # putting bubbles
      ggplot2::geom_point(data = data_trim, ggplot2::aes(x = moderator, y = yi, size = scale, fill = condition), shape = 21, alpha = alpha) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x,se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x,se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x,se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themses
      ggplot2::theme_bw() +
      #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
      ggplot2::theme(legend.direction="horizontal") +
      #theme(legend.background = element_rect(fill = "white", colour = "black")) +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))

    # if(facet == "rows"){
    #   plot <- plot + facet_grid(rows = vars(condition))
    # } else{
    #   plot <- plot + facet_grid(cols = vars(condition))
    # }


  } else{
    plot <-ggplot2::ggplot() +
      # putting bubbles
      #geom_point(data = data, aes(x = moderator, y = yi, size = scale), shape = 21, alpha = alpha, fill = "grey90" ) +
      # prediction interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerPR), method =  "loess", formula = y~x, se = FALSE, lty =  "dotted", lwd = pi.lwd, colour = pi.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperPR), method =  "loess", formula = y~x,se = FALSE, lty = "dotted", lwd = pi.lwd, colour = pi.col) +
      # confidence interval
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = lowerCL), method =  "loess", formula = y~x,se = FALSE,lty = "dashed", lwd = ci.lwd, colour = ci.col) +
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = upperCL), method =  "loess", formula = y~x,se = FALSE, lty ="dashed", lwd = ci.lwd, colour = ci.col) +
      # main line
      ggplot2::geom_smooth(data = mod_table, ggplot2::aes(x = moderator, y = estimate), method =  "loess", formula = y~x, se = FALSE, lwd = est.lwd, colour = est.col) +
      ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) +
      ggplot2::labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
      ggplot2::guides(fill = "none", colour = "none") +
      # themses
      ggplot2::theme_bw() # +
      #theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
      # theme(legend.direction="horizontal") +
      # #theme(legend.background = element_rect(fill = "white", colour = "black")) +
      # theme(legend.background = element_blank()) +
      # theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))
  }

  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }

  # putting k and g in
  # c("top.right", "top.left", "bottom.right", "bottom.left","none")
  if(k == TRUE && g == FALSE && k.pos == "top.right"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                        mapping = ggplot2::aes(x = Inf, y = Inf),
                        label =  paste("italic(k)==", rev(dat_text$K)),
                        parse = TRUE,
                        hjust   = 2,
                        vjust   = 2.5
                        )

  } else if(k == TRUE && g == FALSE && k.pos == "top.left") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2.5
      )
  } else if(k == TRUE && g == FALSE && k.pos == "bottom.right") {
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = 2,
                         vjust   = -1.5
      )
  } else if (k == TRUE && g == FALSE && k.pos == "bottom.left"){
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==", rev(dat_text$K)),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -1.5
      )
    # below get g ----

  } else if (k == TRUE && g == TRUE && k.pos == "top.right"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                                   mapping = ggplot2::aes(x = Inf, y = Inf),
                                   label =  paste("italic(k)==",
                                                  rev(dat_text$K),
                                                         "~","(", rev(dat_text$G), ")"),
                                   parse = TRUE,
                                   hjust   = 1.5,
                                   vjust   = 2)

  } else if (k == TRUE && g == TRUE && k.pos == "top.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = Inf),
                         label =  paste("italic(k)==",
                                        rev(dat_text$K),
                                        "~","(", rev(dat_text$G), ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = 2)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.right"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = Inf, y = -Inf),
                         label =  paste("italic(k)==",
                                        rev(dat_text$K),
                                        "~","(", rev(dat_text$G), ")"),
                         parse = TRUE,
                         hjust   = 1.5,
                         vjust   = -0.5)
  } else if (k == TRUE && g == TRUE && k.pos == "bottom.left"){
    # get group numbers for moderator
    plot <- plot +
      ggplot2::geom_text(data = dat_text,
                         mapping = ggplot2::aes(x = -Inf, y = -Inf),
                         label =  paste("italic(k)==",
                                        rev(dat_text$K),
                                        "~","(", rev(dat_text$G), ")"),
                         parse = TRUE,
                         hjust   = -0.5,
                         vjust   = -0.5)
  }

  # # putting colors in
  # if(cb == TRUE){
  #   plot <- plot +
  #     ggplot2::scale_fill_manual(values=cbpl) +
  #     ggplot2::scale_colour_manual(values=cbpl)
  # }

  return(plot)
}
############################################################################
#' @title i2_ml
#' @description I2 (I-squared) for mulilevel meta-analytic models, based on Nakagawa & Santos (2012). Under multilevel models, we can have multiple I2 (see also Senior et al. 2016). Alternatively, the method proposed by Wolfgang Viechtbauer can also be used.
#' @param model Model object of class \code{rma.mv} or \code{rma}. Currently only model objects using the \code{mods} argument work (e.g., \code{mod = ~1}).
#' @param method Method used to calculate I2. Two options exist: a ratio-based calculation proposed by Nakagawa & Santos (\code{"ratio"}), or Wolfgang Viechtbauer's matrix method (\code{"matrix"}).
#' @param data Data frame used to fit the model.
#' @param boot Number of simulations to run to produce 95 percent confidence intervals for I2. Default is \code{NULL}, where only the point estimate is provided.
#' @return A data frame containing all the model results including mean effect size estimate, confidence, and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' # IMPORTANT NOTE ** boot = 10 is set LOW deliberately to make the models run fast. You should always run for at least boot = 1000
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E,
#' m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng_1 <- i2_ml(english_MA, data = english, boot = 10)
#' I2_eng_2 <- i2_ml(english_MA, data = english, method = "ratio")
#' I2_eng_3 <- i2_ml(english_MA, data = english, method = "matrix")
#'
#' ## Fish example
#' data(fish)
#' warm_dat <- fish
#' model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
#' random = list(~1 | group_ID, ~1 | es_ID),
#' mods = ~ experimental_design + trait.type + deg_dif + treat_end_days,
#' method = "REML", test = "t", data = warm_dat,
#' control=list(optimizer="optim", optmethod="Nelder-Mead"))
#' I2_fish_1 <- i2_ml(model, data = warm_dat, boot = 10)
#' I2_fish_2 <- i2_ml(model, method = c("matrix"),data = warm_dat)
#' I2_fish_2 <- i2_ml(model, method = c("ratio"),data = warm_dat)
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence.
#' The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' I2_lim_1 <- i2_ml(lim_MR, data=lim, boot = 10)
#' I2_lim_2 <- i2_ml(lim_MR, data=lim)
#' }
#' @references Senior, A. M., Grueber, C. E., Kamiya, T., Lagisz, M., O’Dwyer, K., Santos, E. S. A. & Nakagawa S. 2016. Heterogeneity in ecological and evolutionary meta-analyses: its magnitudes and implications. Ecology 97(12): 3293-3299.
#'  Nakagawa, S, and Santos, E.S.A. 2012. Methodological issues and advances in biological meta-analysis.Evolutionary Ecology 26(5): 1253-1274.
#' @export

i2_ml <- function(model, method = c("ratio", "matrix"), data, boot = NULL) {

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  if(any(model$tau2 > 0)) { stop("Sorry. At the moment i2_ml cannot take models with heterogeneous variance.")}

  ## evaluate choices
  method <- match.arg(method)

  if (method == "matrix") {
    # Wolfgang Viechtbauer's method
    I2s <- matrix_i2(model)
  } else {
    # Nakagawa & Santos (2012)
    I2s <- ratio_i2(model)
  }

  if(!is.null(boot)){
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot)

    # Get formula from model object.
    random_formula <- model$random
      mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
                vi <- model$vi

    # Parametric bootstrap
                pb <- progress::progress_bar$new(total = boot,
                                                 format = "Bootstrapping [:bar] :percent ETA: :eta",
                                                 show_after = 0)

     I2_each <- sapply(sim, function(ysim) {

              # The model
             tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                             mods = mods_formula,
                           random = random_formula,
                             data = data))
             pb$tick()
             Sys.sleep(1 / boot)

            if(method == "matrix"){
              I2 <- matrix_i2(tmp)
            } else {
              I2 <- ratio_i2(tmp)
            }

             return(I2) })

      # Summarise the bootstrapped distribution.
       I2s_each_95 <- data.frame(t(apply(I2_each, 1, stats::quantile, probs=c(0.5, .025, .975))))
               I2s <-  round(I2s_each_95, digits = 3)
      colnames(I2s) = c("Est.", "2.5%", "97.5%")
  }

  return(I2s)
}

#' @title matrix_i2
#' @description Calculated I2 (I-squared) for mulilevel meta-analytic models, based on a matrix method proposed by Wolfgang Viechtbauer.
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng <- i2_ml(english_MA, data = english, method = "matrix")
#'
#' # Lim example
#' data(lim)
#' # Add in the sampling variance
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' # Lets fit a meta-regression - I will do Article non-independence. The phylogenetic model found phylogenetic effects, however, instead we could fit Phylum as a fixed effect and explore them with an Orchard Plot
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' I2_lim <- i2_ml(lim_MR, data=lim, method = "matrix")
#' }
#' @export
matrix_i2 <- function(model){

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

    W <- solve(model$V)
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- 100* (sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
    I2_each <- 100* (model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P))))
    names(I2_each) <- paste0("I2_", model$s.names)
    names(I2_total) <- "I2_Total"
    I2s <- c(I2_total, I2_each)
    return(I2s)
}


#' @title ratio_i2
#' @description I2 (I-squared) for mulilevel meta-analytic models based on Nakagawa & Santos (2012). Under multilevel models, we can have a multiple I2 (see also Senior et al. 2016).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @examples
#' \dontrun{
#' # English example
#' data(english)
#' english <- escalc(measure = "SMD", n1i = NStartControl,
#' sd1i = SD_C, m1i = MeanC, n2i = NStartExpt,
#' sd2i = SD_E, m2i = MeanE, var.names=c("SMD","vSMD"),data = english)
#' english_MA <- rma.mv(yi = SMD, V = vSMD,
#' random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
#' I2_eng_1 <- i2_ml(english_MA, data = english, boot = 1000)
#' I2_eng_2 <- i2_ml(english_MA, data = english, method = "ratio")
#' }
#' @export
ratio_i2 <- function(model){

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / model$vi) * (model$k - 1) /
              (sum(1 / model$vi)^2 - sum((1 / model$vi)^2))

  # s^2_t = total variance
  I2_total <- 100 * (sum(model$sigma2) / (sum(model$sigma2) + sigma2_v))
   I2_each <- 100 * (model$sigma2 / (sum(model$sigma2) + sigma2_v))
  names(I2_each) <- paste0("I2_", model$s.names)
  names(I2_total) <- "I2_Total"

  I2s <- c(I2_total, I2_each)
  return(I2s)
}

############################################################################
#' @title r2_ml
#' @description R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @param data Data frame used to fit the \code{rma.mv} or \code{rma} model object
#' @param boot The number of parametric bootstrap iterations, if desired. Defaults to \code{NULL}. A setting of 1000 is recommended as a minimum number of iterations.
#' @return A data frame containing all model results, including: mean effect size estimate, confidence and prediction intervals, with estimates converted back to r.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @references Nakagawa, S, and Schielzeth, H. 2013. A general and simple method for obtaining R2 from generalized linear mixed‐effects models. *Methods in Ecology and Evolution* 4(2): 133-142.
#' @examples
#' \dontrun{
#' data(lim)
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
#' R2 <- r2_ml(lim_MR,data=lim, boot = 10)
#' }
#' @export
r2_ml <- function(model, data, boot = NULL) {

  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}

  if(any(model$tau2 > 0)) { stop("Sorry. At the moment r2_ml cannot take models with heterogeneous variance.")}

  R2 <- R2_calc(model)

  if(!is.null(boot)){

    if(any(class(model) %in% c("robust.rma")) == TRUE){stop("Sorry, bootstrapping currently doesn't work for robust.rma objects. Please use rma.mv instead.")}
    # Simulate the vector of effect sizes
    sim <- metafor::simulate.rma(model, nsim=boot)

    # Get formula from model object.
    random_formula <- model$random
    mods_formula <- metafor::formula.rma(model, type = "mods") #in case moderators
    vi <- model$vi

    pb <- progress::progress_bar$new(total = boot,
                                     format = "Bootstrapping [:bar] :percent ETA: :eta",
                                     show_after = 0)
    # Parametric bootstrap
    R2 <- sapply(sim, function(ysim) {
      # The model
      tmp <- tryCatch(metafor::rma.mv( ysim, vi,
                     mods = mods_formula,
                     random = random_formula,
                     data = data))
      R2s <- R2_calc(tmp)
      pb$tick()
      Sys.sleep(1 / boot)
      return(R2s)
    })

    # Summarise the bootstrapped distribution.
    R2 <- data.frame(t(apply(R2, 1, stats::quantile, probs=c(0.5, .025, .975))))
    R2 <-  round(R2, digits = 3)
    colnames(R2) = c("Est.", "2.5%", "97.5%")
}

return(R2)

}

#' @title R2_calc
#' @description Calculated R2 (R-squared) for mixed (mulitlevel) models, based on Nakagawa & Schielzeth (2013).
#' @param model Model object of class \code{rma.mv} or \code{rma}.
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(lim)
#' lim$vi<-(1/sqrt(lim$N - 3))^2
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
#' random=list(~1|Article, ~1|Datapoint), data=lim)
#' R2 <- R2_calc(lim_MR)
#' }
#' @export

R2_calc <- function(model){
  if(all(class(model) %in% c("robust.rma", "rma.mv", "rma", "rma.uni")) == FALSE) {stop("Sorry, you need to fit a metafor model of class robust.rma, rma.mv, rma, rma.uni")}
  # fixed effect variance
  fix <- stats::var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))

  # marginal
  R2m <- fix / (fix + sum(model$sigma2))

  # conditional. Need to remove 'residual' variance; assume this is the sigma level with the largest k. Really the only way we can get that to work.
  R2c <- (fix + sum(model$sigma2) - model$sigma2[which(model$s.nlevels.f == max(model$s.nlevels.f))]) /
    (fix + sum(model$sigma2))

  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  return(R2s)
}
############################################################################