#' Topics in Microeconometrics Final Project (FGV-EESP, 2025)
#' Students: Gabriela Setti, Loreta Guerra, Marcelo Alonso and Vinícius Nery
#' Professor: Marinho Bertanha
#' 
#' *Project Objective*: propose a bootstrap-based approach not to construct
#' confidence intervals, but rather to choose the bandwidth of linear local/
#' quadratic regressions. See ~line 330 for more details.
#' 
#' _Goal of This Script_: implement a routine that, given bandwidth h, runs a
#' wild bootstrap in the spirit of Bartalotti, Calhoun and He (2017).
#' See ~line 330, 352 for more details.

# Setup ========================================================================
## Libraries -------------------------------------------------------------------
# Data
library(magrittr)
library(glue)
library(purrr)
library(janitor)
library(tictoc)

# Tables + Graphs
library(tinytable)
library(tidyverse)

# Estimation
library(rdrobust)

## Environment -----------------------------------------------------------------
# Useful commands for cleaning the working environment
rm(list = ls())
cat("\014")

# Directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
str_base_directory = dirname(rstudioapi::getActiveDocumentContext()$path)
str_figure_directory = paste(str_base_directory, "figures", sep = "/")
str_table_directory = paste(str_base_directory, "tables", sep = "/")

# Seed
set.seed(121019)

# Theme
my_theme = function(){
   theme_classic() +
      theme(legend.position = 'top', 
            legend.text = element_text(size = 12),
            axis.title.x = element_text(size = 14, colour = "#1d1d1d"),
            axis.title.y = element_text(size = 14, colour = "#1d1d1d"),
            axis.text.x = element_text(colour = "#1d1d1d", size = 12),
            axis.text.y = element_text(colour = "#1d1d1d", size = 12),
            strip.text = element_text(colour = "#1d1d1d", size = 14),
            legend.title = element_blank())
}

# Simulations ==================================================================
# Number of simulations
M = 1

# Sample sizes N: same as Long and Rooklyn (2024)
## 50,000 might be infeasible: too much memory when storing the simulation data
## Can do 500, 10_000, for example and call it a day.
sample_sizes = c(500, 5000, 50000)
num_sample_sizes = length(sample_sizes)

## DGPs ------------------------------------------------------------------------
#' Will use the DGPs used by IK (2012), CCT (2014), PLCW (2022) and LR (2024).
#' It is faster (although more memory intensive) to draw everything in advance,
#' as we need only one function call to draw all the random numbers.

num_dgps = 2

generate_dgp1 = function(N, M, error=TRUE){
   # Skeleton: each array is 2 (number of variables: Y, x) x N x M
   data = array(dim = c(2, N, M), dimnames = list(c("Y", "x"), c(1:N), c(1:M)))
   
   # For each simulation,
   for (m in 1:M){
      # Draw the covariates \in [-1, 1]
      data["x", , m] = 2 * rbeta(N, shape1 = 4, shape2 = 2) - 1
      
      # Indices for x < 0 and x >= 0
      idx_left  <- (data["x", , m] < 0)
      idx_right <- (data["x", , m] >= 0)
      
      # Generate Y
      # since x = 0 is the cut-off, the treatment effect is .52 - .48 = +.04
      ## Left of cut-off
      data["Y", idx_left, m] = 
         0.48 +
         1.27 * data["x", idx_left, m] +
         7.18 * data["x", idx_left, m]^2 +
         20.21 * data["x", idx_left, m]^3 +
         21.54 * data["x", idx_left, m]^4 +
         7.33 * data["x", idx_left, m]^5
      
      ## Right of cut-off
      data["Y", idx_right, m] = 
         0.52 +
         0.84 * data["x", idx_right, m] +
         -3.00 * data["x", idx_right, m]^2 +
         7.99 * data["x", idx_right, m]^3 +
         -9.01 * data["x", idx_right, m]^4 +
         3.56 * data["x", idx_right, m]^5
      
      ## Error
      if (error){
         # Draw the errors, which are independent from x
         e = rnorm(N, mean = 0, sd = .1295)
         
         data["Y", , m] = data["Y", , m] + e
      }
   }
   
   return(data)
}

generate_dgp2 = function(N, M, error=TRUE){
   # Skeleton: each array is 2 (number of variables: Y, x) x N x M
   data = array(dim = c(2, N, M), dimnames = list(c("Y", "x"), c(1:N), c(1:M)))
   
   # For each simulation,
   for (m in 1:M){
      # Draw the covariates \in [-1, 1]
      data["x", , m] = 2 * rbeta(N, shape1 = 4, shape2 = 2) - 1
      
      # Indices for x < 0 and x >= 0
      idx_left  <- (data["x", , m] < 0)
      idx_right <- (data["x", , m] >= 0)
      
      # Generate Y
      # since x = 0 is the cut-off, the treatment effect is .26 - 3.71 = -3.45
      ## Left of cut-off
      data["Y", idx_left, m] = 
         3.71 +
         2.30 * data["x", idx_left, m] +
         3.28 * data["x", idx_left, m]^2 +
         1.45 * data["x", idx_left, m]^3 +
         0.23 * data["x", idx_left, m]^4 +
         0.03 * data["x", idx_left, m]^5
      
      ## Right of cut-off
      data["Y", idx_right, m] = 
         0.26 +
         18.49 * data["x", idx_right, m] +
         -54.81 * data["x", idx_right, m]^2 +
         74.30 * data["x", idx_right, m]^3 +
         -45.02 * data["x", idx_right, m]^4 +
         9.83 * data["x", idx_right, m]^5
      
      ## Error
      if (error){
         # Draw the errors, which are independent from x
         e = rnorm(N, mean = 0, sd = .1295)
         
         data["Y", , m] = data["Y", , m] + e
      }
   }
   
   return(data)
}

# Matrix of DGPs
## Skeleton
dgp_data = matrix(list(), nrow = num_dgps, ncol = num_sample_sizes,
                 dimnames = list(c("DGP1", "DGP2"), sample_sizes))

## Populating
for (i in seq_along(sample_sizes)) {
   dgp_data["DGP1", i] <- list(generate_dgp1(N = sample_sizes[i], M = M))
   dgp_data["DGP2", i] <- list(generate_dgp2(N = sample_sizes[i], M = M))
}

# Example of indexing: DGP1, smaller sample size, Y, 3rd obs., 1st simulation
dgp_data[["DGP1", 1]]["Y", 3, 1]


## Plots -----------------------------------------------------------------------
#' Let's plot the true DPGs (without errors) to see them in action

## Functions
fun_dgp1 = function(x){
   y = ifelse(
      x < 0,
      0.48 + 1.27 * x + 7.18 * x^2 + 20.21 * x^3 + 21.54 * x^4 + 7.33 * x^5,
      0.52 + 0.84 *x - 3.00 * x^2 + 7.99 * x^3 - 9.01 * x^4 + 3.56 * x^5
   )
   
   y[x == 0] = NA
   
   return(y)
}

fun_dgp2 = function(x){
   y = ifelse(
      x < 0,
      3.71 + 2.30 * x + 3.28 * x^2 + 1.45 * x^3 + 0.23 * x^4 + 0.03 * x^5,
      0.26 + 18.49 *x - 54.81 * x^2 + 74.30 * x^3 - 45.02 * x^4 + 9.83 * x^5
   )
   
   y[x == 0] = NA
   
   return(y)
}

## Plots
ggplot() +
   xlim(-1, 1) + 
   geom_function(fun = fun_dgp1, linewidth = 1) +
   geom_vline(xintercept = 0, linewidth = .3, linetype = 'dashed') +
   labs(x = expression(x), y = paste0(expression(y = f(x)), " (DGP 1)")) + 
   my_theme()
ggsave("figures/true_dgp1.pdf", height = 6, width = 9, dpi = 600)

ggplot() +
   xlim(-1, 1) + 
   geom_function(fun = fun_dgp2, linewidth = 1) +
   geom_vline(xintercept = 0, linewidth = .3, linetype = 'dashed') +
   labs(x = expression(x), y = paste0(expression(y = f(x)), " (DGP 2)")) + 
   my_theme()
ggsave("figures/true_dgp2.pdf", height = 6, width = 9, dpi = 600)

# Plots of example data
## DGP 1
ggplot(data.frame(Y = dgp_data[["DGP1", 1]]["Y", , 1], 
                  x = dgp_data[["DGP1", 1]]["x", , 1]), 
       aes(x = x, y = Y)) + 
   geom_point() +
   geom_vline(xintercept = 0, linewidth = .3, linetype = 'dashed') +
   labs(x = expression(x), y = paste0(expression(y = f(x)), " (DGP 1)")) + 
   my_theme()

## DGP 2
ggplot(data.frame(Y = dgp_data[["DGP2", 1]]["Y", , 1], 
                  x = dgp_data[["DGP2", 1]]["x", , 1]), 
       aes(x = x, y = Y)) + 
   geom_point() +
   geom_vline(xintercept = 0, linewidth = .3, linetype = 'dashed') +
   labs(x = expression(x), y = paste0(expression(y = f(x)), " (DGP 2)")) + 
   my_theme()

## RDDs ========================================================================
#' Will only do one time here; have to put this inside an 1:M loop afterwards;
#' in each loop, need to find h via optimization; taken it as given here.
#' In each optimization iteration, will need to perform the following bootstrap.

#' To make things faster, I think we can set the initial guess of the bandwidth
#' in the optimization to h(CCT), since we are going to estimate it anyway to
#' compare the methods. Have to see % of the time the two coincide.

### CCT (2014) -------------------------------------------------------------
#' p is the order of the polynomial used to construct point-estimator,
#' q is the order of the polynomial used to construct bias correction for 
#' valid inference; following CCT Remark 7, q = p + 1 when rho = 1. 
#' In their simulations, p = 1, q = 2.
#' Using heteroskedasticity-robust nearest neighbor variance estimator 
#' (nn with J = 3, default) and bandwidth is selected to minimize common MSE
#' of the treatment effect estimator (mserd, default).

# Array to store results: 2 DGPs x 3 Sample Sizes x M simulations
array_cct = array(list(), dim = c(num_dgps, num_sample_sizes, M), 
                  dimnames = list(c("DGP1", "DGP2"), sample_sizes, c(1:M)))

# Running CCTs
tic()
for (m in 1:M){
   for (dgp in 1:num_dgps){
      for (size in 1:num_sample_sizes){
         
         # Procedure with p = 1, q = 2, as default and in CCT simulations
         # Imposing same bandwidth (rho = 1) to be comparable with our procedure
         # Note that this only changes the bandwidth used for bias-correctness,
         # not the 'usual' bandwidth. By default, h is the same on both sides.
         array_cct[dgp, size, m] = list(rdrobust(
            y = dgp_data[[dgp, size]]["Y", , m],
            x = dgp_data[[dgp, size]]["x", , m], c = 0,
            p = 1, q = 2, kernel = "triangular", rho = 1
         ))
         
         # # Special case where h = b
         # ## Compute optimal bandwidth with p = 1
         # h_p1 = rdbwselect(
         #    y = dgp_data[[dgp, size]]["Y", , m],
         #    x = dgp_data[[dgp, size]]["x", , m], c = 0,
         #    p = 1, kernel = "triangular"
         # )$bws[1, 1]
         # 
         # ## Estimate a local quadratic regression with this h and 
         # ## imposing h = b, which is the same as rho = 1
         # array_cct[dgp, size, m] = rdrobust(
         #    y = dgp_data[[dgp, size]]["Y", , m],
         #    x = dgp_data[[dgp, size]]["x", , m], c = 0,
         #    p = 2,, kernel = "triangular", h = h_p1, rho = 1
         # )
         # 
         #'`array_cct[dgp, size, m]$coef["Conventional", 1]` will be the same
         #' #' as the `$coef["Robust", 1]` or `$coef["Bias-Corrected", 1]` of 
         # rdrobust(
         #    y = dgp_data[[dgp, size]]["Y", , m],
         #    x = dgp_data[[dgp, size]]["x", , m], c = 0,
         #    p = 1, q = 2, kernel = "triangular", rho = 1
         # )
      }
   }
}
toc() 
#' 2.5-3.5s for all sample sizes and 1 simulation with only 83kb of memory
#' 1000 simulations: 3500s, 1 hour, potentially 83Mb of memory

# To retrieve an object of DGP1, second smallest sample size and 1st simulation
summary(array_cct[["DGP1", 2, 1]])

# Bandwidth of the point estimator ("b": bias-correction)
array_cct[["DGP1", 2, 1]]$bws["h", 1]

# Estimate ("Bias-Corrected", "Conventional", "Robust" = BC)
array_cct[["DGP1", 2, 1]]$coef["Robust", 1]

# Confidence intervals ("Bias-Corrected", "Conventional", "Robust")
array_cct[["DGP1", 2, 1]]$ci["Robust", ]

# Conventional estimates are the differences in intercepts
array_cct[["DGP1", 2, 1]]$coef["Conventional", 1] ==
   array_cct[["DGP1", 2, 1]]$beta_Y_p_r[1] - 
   array_cct[["DGP1", 2, 1]]$beta_Y_p_l[1]

# Sizes to left or right of cutoff
array_cct[["DGP1", 2, 1]]$N_h

### Wild Bootstrap --------------------------------------------------------
#' Will follow the spirit of the wild bootstrap procedure of
#' Bartalotti, Calhoun and He (2017). However, the goal here is different: we
#' want to select an optimal bandwidth, while their goal is to conduct valid
#' inference using bootstraps. In their simulations, their starting point is 
#' CCT's AMS-optimal bandwidth, which is fixed given the simulation.

#' Thus, GIVEN h (which will be optimized over), will first run a L*L*R in the 
#' whole sample and do B wild bootstraps, each with L*L*R as well (different
#' from Bartalotti et al. (2017), which use LQR in the whole sample and LLR 
#' in the bootstraps). Under the distribution induced by the bootstrap, the 
#' estimated effect of the whole sample (based on LLR) is the true treatment 
#' effect, so can use the *Bootstrap Principle* to justify our method.
#' This procedure will be optimized over h to find the bandwidth that
#' minimizes the MSE between the estimated value in the whole sample and
#' the bootstrap values.

#' Now, with the optimal _h*_ in hand, we can run a L*Q*R in the whole sample
#' that will produce bias-corrected inference and point estimates under 
#' the assumptions of CCT's Remark 7, since h* was calculated using p = 1.
#' BI (2020), for instance, also work under the assumptions of this remark.
#' We can also use Bartalotti et al. (2017)'s bootstrap, again assuming h = b.
#' Final output would be to compare the MSE of the bias-corrected estimate 
#' of the L*Q*R using our bandwidth with CCT's and IK's bandwidth (both 
#' generated in the same way: choose h using LLR, final estimate using LQR
#' and same bandwidth).
#' Comparing coverage rates (% of the time the bias-corrected robust 
#' confidence interval contains the true parameter) would also be interesting


#' Will fix h and the data just to write the code: in the final version, 
#' this bootstrap needs to be inside a function which will optimize over h.
#' The function itself should be inside a loop of M simulations, which
#' should contain various sample size or have a loop inside/outside to go over
#' the various sample sizes.

h = array_cct[["DGP1", 2, 1]]$bws["h", 1]
data_h = as.data.frame(t(dgp_data[["DGP1", 2]][, , 1]))

#### Fitted Values ----------------------------------------------------------
#' A nice thing is that I am not worried about confidence intervals, 
#' only caring about the fitted values and the residuals.
#' Thus, can run a simple kernel-weighted OLS with observations in the region
#' defined by the bandwidth. See Bertanha and Imbens (2018, p. 13).

#' Kernels: heavy inspiration from the `rdrobust` source code. 
#' Returns a vector of weights for each observation
kernel_weights = function(x, h, c = 0, kernel = "triangular"){
   # Define relative distance
   u = (x - c) / h
   
   # Kernel weights
   if (kernel == 'triangular'){
      w = ((1 - abs(u)) * (abs(u) <= 1)) / h
   } else if (kernel == "uniform"){
      w = (.5 * (abs(u) <= 1)) / h
   } else { # epanechnikov
      w = (0.75 * (1 - u^2)*(abs(u) <= 1)) / h
   }
   
   return(w)
}

# Calculating weights
data_h$w = kernel_weights(x = data_h$x, h = h)

# Estimating LLR on the left side of the cutoff
(lm(formula = Y ~ x, data = data_h,
    subset = (x > -h & x < 0), weights = data_h$w))
array_cct[["DGP1", 2, 1]]$beta_Y_p_l

# Estimating LLR on the right side of the cutoff
(lm(formula = Y ~ x, data = data_h,
    subset = (x >= 0 & x < h), weights = data_h$w))
array_cct[["DGP1", 2, 1]]$beta_Y_p_r

#' The above was a test: we get the same values as `rdrobust`. 
#' Now, make it a function!

do_local_linear_reg_cutoff = function(data_frame, h, boot = FALSE,
                                      c = 0, kernel = "triangular"){
   # Calculating kernel weights
   data_frame$w = kernel_weights(x = data_frame$x, h = h, c = c, kernel = kernel)
   
   # Estimating LLR at the left and right of cut-off
   if (boot){
      llr_left = lm(formula = Y_wild_boot ~ x, data = data_frame, 
                    subset = (x > -h & x < 0), weights = data_frame$w)
      
      llr_right = lm(formula = Y_wild_boot ~ x, data = data_frame, 
                     subset = (x >= 0 & x < h), weights = data_frame$w)
   } else {
      llr_left = lm(formula = Y ~ x, data = data_frame, 
                    subset = (x > -h & x < 0), weights = data_frame$w)
      
      llr_right = lm(formula = Y ~ x, data = data_frame, 
                     subset = (x >= 0 & x < h), weights = data_frame$w)
   }
   
   # Effect estimate
   effect_estimate = llr_right$coefficients[[1]] - llr_left$coefficients[[1]]
   
   # Fitted values
   left_fit = as.vector(llr_left$fitted.values)
   right_fit = as.vector(llr_right$fitted.values)
   
   # Residuals
   left_residuals = as.vector(llr_left$residuals)
   right_residuals = as.vector(llr_right$residuals)
   
   # List
   out = list(llr_left = llr_left,
              llr_right = llr_right,
              effect_estimate = effect_estimate,
              left_Nh = length(left_fit),
              right_Nh = length(right_fit),
              left_fit = left_fit,
              right_fit = right_fit,
              left_residuals = left_residuals,
              right_residuals = right_residuals)
}

#### Bootstrap -----------------------------------------------------------------
# Number of bootstrap replications
num_bootstrap = 500 # same as B1 in Bartalotti et al. (2016)

# List to store coefficients
tau_wild_bootstrap <- rep(NA, num_bootstrap)

# Initializing auxiliary columns 
data_h$Y_wild_boot = 0

# Treatment effect in the whole sample
llr_results = do_local_linear_reg_cutoff(data_h, h)

# Looping over the bootstrap samples
tic()
for (boot_sample in 1:num_bootstrap){
   # Data
   data_h_left = data_h[data_h$x < 0 & data_h$x > -h, ]
   data_h_right = data_h[data_h$x >= 0 & data_h$x < h, ]
   
   # Sample sizes
   left_Nh = llr_results$left_Nh
   right_Nh = llr_results$right_Nh
   
   # Generating vector of sign changes
   left_sign = ifelse(runif(left_Nh, min = -1, max = 1) < 0, -1, 1)
   right_sign = ifelse(runif(right_Nh, min = -1, max = 1) < 0, -1, 1)
   
   # Sampling indexes of residuals
   left_indexes = sample(left_Nh, left_Nh, replace = TRUE)
   right_indexes = sample(right_Nh, right_Nh, replace = TRUE)
   
   # Generating the bootstrap Y
   data_h_left$Y_wild_boot = predict(llr_results$llr_left, data_h_left) +
      left_sign * llr_results$left_residuals[left_indexes]
   
   data_h_right$Y_wild_boot = predict(llr_results$llr_right, data_h_right) +
      right_sign * llr_results$right_residuals[right_indexes]
   
   # Estimating LLR on the bootstrap
   llr_results_boot = do_local_linear_reg_cutoff(
      data_frame = rbind(data_h_left, data_h_right),
      h = h, boot = TRUE
   )
   
   # Store the coefficient
   tau_wild_bootstrap[boot_sample] = llr_results_boot$effect_estimate
}
toc()  # 3.3-4.6s
#' This needs to be optimized maybe: in each simulation, this bootstrap will
#' run for different sample sizes and over many iterations of h...
#' If we do two sample sizes, 1000 simulations and in each it takes 
#' 20 iterations to find the optimal h, quick math with 4s gives 44h...
#' Don't think will be that many, since we can use the optimal h of previous
#' simulations to speed things up, but still
#' Maybe can reduce by sampling indexes and sign changes in advance (at each
#' simulation m), but at the trade-off of memory allocation (would need a 
#' B x N for each, since samples need to be independent across bootstraps).
#' Another alternative is to use low-level interface o `lm`, `lm.fit`. See
#' https://chatgpt.com/share/67d7fe1f-6bd4-8013-9b2d-fcfd01fd15ef

# Calculating rootMSE
rmse = sqrt(1 / num_bootstrap * 
   sum((tau_wild_bootstrap - llr_results$effect_estimate)^2))
