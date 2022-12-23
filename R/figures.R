# creating figures

# packages

library(tidyverse)
library(orchaRd)
library(metafor)
library(here)
library(patchwork)
library(grDevices)
library(cowplot)

# datasets

dat_Midolo_2019 <- read.csv(here("data","Midolo_2019_Global Change Biology.csv"))
dat_Midolo_2019 <- dat_Midolo_2019[,which(colnames(dat_Midolo_2019) %in% c("Study_ID", "study_name", "species", "trait", "treatment", "control", "sd_treatment", "sd_control", "n_treatment", "n_control", "elevation", "elevation_log"))]

# adding lnRR

lnRR <- escalc(measure = "ROM",  # "ROM" means ratio of means; lnRR is specified to be calculated;
               m1i = treatment, # mean of group 1 (e.g., environmental stressor) 
               m2i = control, # mean of group 2 (e.g., control)
               sd1i = sd_treatment, # standard deviation of mean of group 1 (e.g., environmental stressor)
               sd2i = sd_control, # standard deviation of group 2 (e.g., control) 
               n1i = n_treatment, # sample size of group 1 (e.g., environmental stressor) 
               n2i = n_control, # sample size of group 2 (e.g., control) 
               data = dat_Midolo_2019, # dataset containing the above information (here is the dataset of our working example)
               digits = 3,
               append = FALSE)

### bind the four sets of effect sizes into one dataframe
metrics_set <- data.frame(lnRR = lnRR$yi, lnRRV = lnRR$vi)
dat <- cbind(dat_Midolo_2019, metrics_set) # bind_cols()


# getting year
dat <- dat %>% mutate(Year = as.integer(str_extract(study_name,"\\d+")))
# centered
dat$Year.c <- dat$Year - mean(dat$Year) 

# effect sample size
# write a help function to calculate effective size based sampling variance - tilde n
ess.var_cal <- function(dat){1/dat$n_control + 1/dat$n_treatment} 
# calculate tilde N
dat$ess.var <- ess.var_cal(dat)
# calculate effective size based sampling error - tilde square root n
dat$ess.se <- sqrt(dat$ess.var)


# ES_ID

dat$ES_ID <-1:nrow(dat)

# looking at data

#View(dat)


###########
# Figure 4
###########

# data for forest plot

# unique_values <- unique(dat$Study_ID)
# unique_positions <- sapply(unique_values, function(x) which(dat$Study_ID == x)[1])

dat_short <- dat %>% group_by(Study_ID) %>% slice(1) %>% 
  ungroup(Study_ID) %>% slice(1:20)

# random effects model
ma_short <- rma(yi = lnRR, # observed effect sizes / estimates of SMD; the outputs of escalc() function; 
                vi = lnRRV, # the estimates of sampling variance of SMD; 
                method = "REML", # setting restricted maximum likelihood estimator as the estimators for the amount of heterogeneity;
                data = dat_short # the dataset
)


# 100 effect sizes

dat_mid <- dat  %>% slice(1:100)

# 3 level model

ma <- rma.mv(yi = lnRR, 
             V = lnRRV, 
             random = list(~1 | Study_ID, # allows true effect sizes to vary among different primary studies - account for the between-study effect and quantify between-study heterogeneity;
                           ~1 | ES_ID), # allows true effect sizes to vary within primary studies - account for the with-study effect and quantify with-study heterogeneity;
             method = "REML", # REML is assigned as the estimator for variance components as suggested;
             test = "t", # t-distribution is specified for the tests of model coefficients and CIs;
             dfs = "contain", # the methods to adjust the (denominator) degrees of freedom;
             data = dat_mid, # our dataset
             sparse = TRUE
)



mr1 <- rma.mv(yi = lnRR, 
              V = lnRRV, 
              mod = ~ trait,
              random = list(~1 | Study_ID, # allows true effect sizes to vary among different primary studies - account for the between-study effect and quantify between-study heterogeneity;
                            ~1 | ES_ID), # allows true effect sizes to vary within primary studies - account for the with-study effect and quantify with-study heterogeneity;
              method = "REML", # REML is assigned as the estimator for variance components as suggested;
              test = "t", # t-distribution is specified for the tests of model coefficients and CIs;
              dfs = "contain", # the methods to adjust the (denominator) degrees of freedom;
              data = dat, # our dataset
              sparse = TRUE
)




mr2 <- rma.mv(yi = lnRR, 
              V = lnRRV, 
              mod = ~ elevation_log,
              random = list(~1 | Study_ID, # allows true effect sizes to vary among different primary studies - account for the between-study effect and quantify between-study heterogeneity;
                            ~1 | ES_ID), # allows true effect sizes to vary within primary studies - account for the with-study effect and quantify with-study heterogeneity;
              method = "REML", # REML is assigned as the estimator for variance components as suggested;
              test = "t", # t-distribution is specified for the tests of model coefficients and CIs;
              dfs = "contain", # the methods to adjust the (denominator) degrees of freedom;
              data = dat, # our dataset
              sparse = TRUE
)



# A

pdf(NULL)
dev.control(displaylist="enable")
forest(ma_short, 
       annotate = FALSE,
       xlab = "lnRR (log response ratio)",
       slab = dat_short$study_name,
       cex = 0.8)
a <- recordPlot()
invisible(dev.off())


# B

b <- caterpillars(ma, mod = "1", 
                  xlab = "lnRR (log response ratio)", 
                  group = "Study_ID", 
                  data = dat_mid)

# C

c <- orchard_plot(mr1, mod = "trait", 
                  group = "Study_ID", 
                  data = dat, 
                  xlab = "lnRR (log response ratio)")  

# D

d <- bubble_plot(mr2, mod = "elevation_log", 
                group = "Study_ID", 
                data = dat, 
                ylab = "lnRR (log response ratio)", 
                xlab = "log(elevation)", 
                g = TRUE,
                legend.pos = "bottom.left")  


(ggdraw(a) + plot_grid(b))/(plot_grid(c) + plot_grid(d)) + plot_annotation(tag_levels = 'A')


###########
# Figure 5
###########

dat_short2 <- dat %>% group_by(Study_ID) %>% slice(1) %>% 
  ungroup(Study_ID) %>% slice(1:30)

# random effects model
ma_short2 <- rma(yi = lnRR, # observed effect sizes / estimates of SMD; the outputs of escalc() function; 
                vi = lnRRV, # the estimates of sampling variance of SMD; 
                method = "REML", # setting restricted maximum likelihood estimator as the estimators for the amount of heterogeneity;
                data = dat_short2 # the dataset
)




mr3 <- rma.mv(yi = lnRR, 
              V = lnRRV, 
              mod = ~ ess.se,
              random = list(~1 | Study_ID, # allows true effect sizes to vary among different primary studies - account for the between-study effect and quantify between-study heterogeneity;
                            ~1 | ES_ID), # allows true effect sizes to vary within primary studies - account for the with-study effect and quantify with-study heterogeneity;
              method = "REML", # REML is assigned as the estimator for variance components as suggested;
              test = "t", # t-distribution is specified for the tests of model coefficients and CIs;
              dfs = "contain", # the methods to adjust the (denominator) degrees of freedom;
              data = dat, # our dataset
              sparse = TRUE
)


mr4 <- rma.mv(yi = lnRR, 
              V = lnRRV, 
              mod = ~ Year.c,
              random = list(~1 | Study_ID, # allows true effect sizes to vary among different primary studies - account for the between-study effect and quantify between-study heterogeneity;
                            ~1 | ES_ID), # allows true effect sizes to vary within primary studies - account for the with-study effect and quantify with-study heterogeneity;
              method = "REML", # REML is assigned as the estimator for variance components as suggested;
              test = "t", # t-distribution is specified for the tests of model coefficients and CIs;
              dfs = "contain", # the methods to adjust the (denominator) degrees of freedom;
              data = dat, # our dataset
              sparse = TRUE
)



# A
pdf(NULL)
dev.control(displaylist="enable")
funnel(ma_short2, yaxis = "seinv",
       #level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, #legend=TRUE,
      xlab = "Model residuals of lnRR (log response ratio)")
A <- recordPlot()
invisible(dev.off())

# B

B <- bubble_plot(mr3, mod = "ess.se", 
                 group = "Study_ID", 
                 data = dat, 
                 ylab = "lnRR (log response ratio)", 
                 xlab = "Square root of 'effective sample size'", 
                 g = TRUE,
                 legend.pos = "bottom.left")  


# C

C <- bubble_plot(mr4, mod = "Year.c", 
                 group = "Study_ID", 
                 data = dat, 
                 ylab = "lnRR (log response ratio)", 
                 xlab = "Publication Year (centered)", 
                 g = TRUE,
                 legend.pos = "bottom.left")  

ggdraw(A) + plot_grid(B) + plot_grid(C)  + plot_annotation(tag_levels = 'A') + plot_layout(nrow = 3, byrow = FALSE)
