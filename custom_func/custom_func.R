#'---------------------------------------------------------------#'
#' Functions used in the methodological paper 
#' Quantitative synthesis: a practical guide on meta-analysis, meta-regression, and publication bias tests for environmental sciences
#' Oct 2022
#'---------------------------------------------------------------#'


############################################################################
#' @title cv_avg
#' @description Calculates the weighted average CV^2 within a study and the weighted average CV^2 across a study
#' @param x Mean of an experimental group
#' @param sd Standard deviation of an experimental group
#' @param n The sample size of an experimental group
#' @param group Study, grouping or cluster variable one wishes to calculate the within and between weighted CV^2. In meta-analysis this will most likely be 'study'.
#' @param data The dataframe containing the mean, sd, n and grouping variables
#' @param label A character string specifying the label one wishes to attach to columns to identify the treatment. Otherwise, if not specified it will default to the variable name for x
#' @param sub_b A logical indicating whether the between study CV^2 (b_CV2) should be appended to the data only ('TRUE') or whether both within study CV^2 (w_CV2), mean sample size (n_mean) and between study CV^2 (b_CV2) should all be appended to the data only ('FALSE')
#' @param cv2 A logical indicating whether one should take the weighted average of CV2 or the weighted average of CV followed by squaring this average. Default to FALSE.
#' @example \dontrun{
#' # test data for cv_avg function
#' library(tidyverse)
#' set.seed(76)
#' x1 = rnorm(16, 6, 1)
#' x2 = rnorm(16, 6, 1)
#' test_dat <- data.frame(stdy = rep(c(1,2,3,4), each = 4),x1 = x1,sd1 = exp(log(x1)*1.5 + rnorm(16, 0,
#' sd = 0.10)),n1 = rpois(16, 15),x2 = x2,sd2 = exp(log(x2)*1.5 + rnorm(16, 0, sd = 0.10)),n2 = rpois(16, 15))
#' rm(list = c("x1", "x2"))
#' # # Now generate some missing data
#' t2 <- gen.miss(test_dat, "sd1", "sd2", 6)
#' t2_cv2 <- cv_avg(x = x1, sd = sd1, n = n1, stdy, data =  t2, sub_b = FALSE, cv2 = TRUE)
#' t2_cv2 <- cv_avg(x2, sd2, n2, stdy, label = "2", data =  t2_cv, sub_b = FALSE)
#' # Check calculations are correct. All match what is expected
#' test <- t2_cv %>%  filter(stdy == "1")
#' # Within
#' mean(test$n1) # Matches 14.75
#' mean(test$n2) # Matches 14
#' # CV^2
#' weighted.mean((test$sd1 / test$x1)^2, test$n1, na.rm = T)
#' weighted.mean((test$sd2 / test$x2)^2, test$n2, na.rm = T)
#' # mean(CV)^2
#' t2_cv <- cv_avg(x = x1, sd = sd1, n = n1, stdy, data =  t2, sub_b = FALSE, cv2 = FALSE)
#' weighted.mean((test$sd1 / test$x1), test$n1, na.rm = T)^2
#' # Between
#' wCV1 = unique(t2_cv2$w_CV2_x1)
#' w_nt1 = c(59,58,58,50)
#' weighted.mean(wCV1, w_nt1)
#' wCV2 = unique(t2_cv2$w_CV2_2)
#' w_nt2 = c(56, 56, 72, 63)
#' weighted.mean(wCV2, w_nt2)
#' }

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
#' @param sd Standard deviation of an experimental group
#' @param x Mean of an experimental group
#' @param n The sample size of an experimental group
#' @param cv2 Logical indicating whether the weighted average of CV^2 or CV should be taken (followed by squaring weighted average CV). Defaults to weighted average of CV.

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
#' @param model The rma.mv model object

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
#' @param m1 Mean of treatment group 1
#' @param m2 Mean of treatment group 2
#' @param cv1_2 Coefficient of variation squared (CV^2) for treatment group 1
#' @param cv2_2 Coefficient of variation squared (CV^2) for treatment group 2
#' @param n1 Sample size for treatment group 1
#' @param n2 Sample size for treatment group 2
#'
lnrr_laj <- function(m1, m2, cv1_2, cv2_2, n1, n2){
  log(m1 / m2) + 0.5*((cv1_2 / n1) - (cv2_2 / n2))
}

############################################################################

#' @title v_lnrr_laj
#' @description Calculates the sampling variance for log response ratio based on second order Taylor expansion proposed by Lajeunesse 2011
#' @param cv1_2 Coefficient of variation squared (CV^2) for treatment group 1
#' @param cv2_2 Coefficient of variation squared (CV^2) for treatment group 2
#' @param n1 Sample size for treatment group 1
#' @param n2 Sample size for treatment group 2
v_lnrr_laj <- function(cv1_2, cv2_2, n1, n2){
  ((cv1_2) / n1) + ((cv2_2) / n2) +
    ((cv1_2)^2 / (2*n1)^2) + ((cv2_2)^2 / (2*n2)^2)
}

############################################################################

#' @title gen.miss
#' @description Generates random missing data in two columns of a meta-analytic dataset
#' @param data The dataframe containing the mean, sd, n and grouping variables
#' @param missVar A character string specifying the first column one wishes to have missing data. Example, "sd1" column.
#' @param missCol2 A character string specifying the second column one wishes to have missing data. Example, "sd2" column.
#' @param n_miss The number of missing data (NA) one wishes to introduce to the missVar and missCol2
gen.miss <- function(data, missVar, missCol2, n_miss){
  data[sample(rownames(data), n_miss), missVar] <- NA
  data[is.na(data[,missVar]), missCol2] <- NA
  return(data)
}


############################################################################
# Effect size (lnRR and SMD) for 2 main effects and interaction effect


interactive_es <- function(CC_n, CC_mean, CC_SD,
                           EC_n, EC_mean, EC_SD,
                           CS_n, CS_mean, CS_SD,
                           ES_n, ES_mean, ES_SD)
  {
    # lnRR
    # main effect for environmental group
    lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
      log(0.5*(CS_mean+ CC_mean))
    
    lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) + 
      (1/(CS_mean + CC_mean))^2*(CS_SD^2 / CS_n + CC_SD^2 / CC_n)
    
    # main effect for stress group
    lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
      log(0.5*(EC_mean+ CC_mean))
    
    lnRRV_S <- lnRRV_E
    
    # interaction between experimental and stress groups
    
    lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
      (log(EC_mean) - log(CC_mean))
    
    
    lnRRV_ES <- 
      (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
         ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
         ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
         ((CC_SD)^2 / ((CC_mean)^2*CC_n)))

    effect <- tibble(
      # lnRR
      lnRR_E = lnRR_E,
      lnRRV_E = lnRRV_E, 
      lnRR_S = lnRR_S, 
      lnRRV_S = lnRRV_S,
      lnRR_ES =lnRR_ES, 
      lnRRV_ES = lnRRV_ES
    )
    effect
}


############################################################################

# meta-analysis of magnitude
## folded effect size
folded_es <-function(mean, variance){ # the sampling variance of magnitude   
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_mu
}
## folded variance
folded_var <- function(mean, variance){ # the sampling variance of magnitude   
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_se <- sqrt(mu^2 + sigma^2 - fold_mu^2)
  # adding se to make bigger mean
  fold_v <- fold_se^2
  fold_v
}

############################################################################

# custom function for extracting mean and CI from each metafor model
estimates.model <- function(model){
  db.mf <- data.frame(round(model$b, 3),row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,round(model$ci.lb, 3),round(model$ci.ub,3),row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}

############################################################################
