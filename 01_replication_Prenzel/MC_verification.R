##### MC verification of variance estimates for Replication_Prenzel ############
# Martin Möhler, 2026
############################################################################## #

library(ptable)
library(data.table)

## helper functions

# apply cell key perturbation to a single count value
get_perturbation_value <- function(x_ck, pt) {
  
  if(is.na(x_ck[1])) { 
    NA 
  } else {
    max_x <- max(pt$i)
    pt$v[pt$i == min(max_x, x_ck[1]) & data.table::between(x_ck[2], pt$p_int_lb, pt$p_int_ub)]
  }
}

# apply random perturbation to a vector of counts
noise_x <- function(x) {
  
  ck <- runif(length(x))
  v <- apply(cbind(x, ck), 1, get_perturbation_value, pt = pt@pTable)
  x + v
}


## index & variance functions

get_hhi <- function(x, N) { sum((x/N)^2) }
get_ent <- function(x, N) { -sum((x/N) * ifelse(x == 0, 0, log(x/N))) }

get_div1 <- function(x, N) { (1 - get_hhi(x, N)) / (1 - 1/length(x)) * 100}
get_div2 <- function(x, N) { get_ent(x, N) / log(length(x)) * 100}

hhi_getV <- function(x, N) { (sum(x^2) * 4) / N^4 }
ent_getV <- function(x, N) { sum((log(N) - ifelse(x == 0, 0, log(x)) - 1)^2) / N^2 }

div1_getV <- function(x, N) { hhi_getV(x, N) * (1 / (1 - 1/length(x)))^2 * 10^4 }
div2_getV <- function(x, N) { ent_getV(x, N) / log(length(x))^2 * 10^4 }


## Simulation parameters

pt <- create_cnt_ptable(D = 5, V = 1) # noise distribution
R <- 1e+4                             # no. of Monte-Carlo repetitions


# ----- Simulation -----

# elementary regions
N_even <- rep(100, 10)
N_skew <- c(991, rep(1, 9))
N_unev <- c(480, 230, 100, 70, 40, 20, 20, 20, 10, 10) 
# composite region
N_tota <- N_even + N_skew + N_unev


## results from unperturbed data

get_hhi(N_even, sum(N_even))  # HHI in region with even distribution
get_hhi(N_skew, sum(N_skew))  # HHI in region with highly skewed distribution
get_hhi(N_unev, sum(N_unev))  # HHI in region with somewhat uneven distribution
get_hhi(N_tota, sum(N_tota))  # HHI in composite region

get_div1(N_even, sum(N_even)) # DIV1 in region with even distribution
get_div1(N_skew, sum(N_skew)) # DIV1 in region with highly skewed distribution
get_div1(N_unev, sum(N_unev)) # DIV1 in region with somewhat uneven distribution
get_div1(N_tota, sum(N_tota)) # DIV1 in composite region

get_ent(N_even, sum(N_even))  # E in region with even distribution
get_ent(N_skew, sum(N_skew))  # E in region with highly skewed distribution
get_ent(N_unev, sum(N_unev))  # E in region with somewhat uneven distribution
get_ent(N_tota, sum(N_tota))  # E in composite region

get_div2(N_even, sum(N_even)) # DIV2 in region with even distribution
get_div2(N_skew, sum(N_skew)) # DIV2 in region with highly skewed distribution
get_div2(N_unev, sum(N_unev)) # DIV2 in region with somewhat uneven distribution
get_div2(N_tota, sum(N_tota)) # DIV2 in composite region


## variance estimates based on Monte-Carlo simulation

resV <- data.frame(r = 1:R,
                   hhi_even  = NA, hhi_skew  = NA, hhi_unev  = NA,
                   ent_even  = NA, ent_skew  = NA, ent_unev  = NA,
                   div1_even = NA, div1_skew = NA, div1_unev = NA,
                   div2_even = NA, div2_skew = NA, div2_unev = NA,
                   var_hhi_even  = NA, var_hhi_skew  = NA, var_hhi_unev  = NA, var_hhi_tota  = NA,
                   var_ent_even  = NA, var_ent_skew  = NA, var_ent_unev  = NA, var_ent_tota  = NA,
                   var_div1_even = NA, var_div1_skew = NA, var_div1_unev = NA, var_div1_tota = NA,
                   var_div2_even = NA, var_div2_skew = NA, var_div2_unev = NA, var_div2_tota = NA)
set.seed(42)
for(r in 1:R) {
  
  ckm_even <- noise_x(N_even)
  ckm_skew <- noise_x(N_skew)
  ckm_unev <- noise_x(N_unev)
  ckm_tota <- ckm_even + ckm_skew + ckm_unev
  
  resV$hhi_even[r] <- get_hhi(ckm_even, sum(N_even))
  resV$hhi_skew[r] <- get_hhi(ckm_skew, sum(N_skew))
  resV$hhi_unev[r] <- get_hhi(ckm_unev, sum(N_unev))
  resV$hhi_tota[r] <- get_hhi(ckm_tota, sum(N_tota))
  resV$ent_even[r] <- get_ent(ckm_even, sum(N_even))
  resV$ent_skew[r] <- get_ent(ckm_skew, sum(N_skew))
  resV$ent_unev[r] <- get_ent(ckm_unev, sum(N_unev))
  resV$ent_tota[r] <- get_ent(ckm_tota, sum(N_tota))
  
  resV$div1_even[r] <- get_div1(ckm_even, sum(N_even))
  resV$div1_skew[r] <- get_div1(ckm_skew, sum(N_skew))
  resV$div1_unev[r] <- get_div1(ckm_unev, sum(N_unev))
  resV$div1_tota[r] <- get_div1(ckm_tota, sum(N_tota))
  resV$div2_even[r] <- get_div2(ckm_even, sum(N_even))
  resV$div2_skew[r] <- get_div2(ckm_skew, sum(N_skew))
  resV$div2_unev[r] <- get_div2(ckm_unev, sum(N_unev))
  resV$div2_tota[r] <- get_div2(ckm_tota, sum(N_tota))
  
  resV$var_hhi_even[r] <- var(resV$hhi_even[1:r])
  resV$var_hhi_skew[r] <- var(resV$hhi_skew[1:r])
  resV$var_hhi_unev[r] <- var(resV$hhi_unev[1:r])
  resV$var_hhi_tota[r] <- var(resV$hhi_tota[1:r])
  resV$var_ent_even[r] <- var(resV$ent_even[1:r])
  resV$var_ent_skew[r] <- var(resV$ent_skew[1:r])
  resV$var_ent_unev[r] <- var(resV$ent_unev[1:r])
  resV$var_ent_tota[r] <- var(resV$ent_tota[1:r])
  
  resV$var_div1_even[r] <- var(resV$div1_even[1:r])
  resV$var_div1_skew[r] <- var(resV$div1_skew[1:r])
  resV$var_div1_unev[r] <- var(resV$div1_unev[1:r])
  resV$var_div1_tota[r] <- var(resV$div1_tota[1:r])
  resV$var_div2_even[r] <- var(resV$div2_even[1:r])
  resV$var_div2_skew[r] <- var(resV$div2_skew[1:r])
  resV$var_div2_unev[r] <- var(resV$div2_unev[1:r])
  resV$var_div2_tota[r] <- var(resV$div2_tota[1:r])
}


## variance estimates from formula

Vtheo_hhi_even <- hhi_getV(N_even, sum(N_even)) * pt@pParams@V
Vtheo_hhi_skew <- hhi_getV(N_skew, sum(N_skew)) * pt@pParams@V
Vtheo_hhi_unev <- hhi_getV(N_unev, sum(N_unev)) * pt@pParams@V
Vtheo_hhi_tota <- hhi_getV(N_tota, sum(N_tota)) * pt@pParams@V * 3

Vtheo_ent_even <- ent_getV(N_even, sum(N_even)) * pt@pParams@V
Vtheo_ent_skew <- ent_getV(N_skew, sum(N_skew)) * pt@pParams@V
Vtheo_ent_unev <- ent_getV(N_unev, sum(N_unev)) * pt@pParams@V
Vtheo_ent_tota <- ent_getV(N_tota, sum(N_tota)) * pt@pParams@V * 3

Vtheo_div1_even <- div1_getV(N_even, sum(N_even)) * pt@pParams@V
Vtheo_div1_skew <- div1_getV(N_skew, sum(N_skew)) * pt@pParams@V
Vtheo_div1_unev <- div1_getV(N_unev, sum(N_unev)) * pt@pParams@V
Vtheo_div1_tota <- div1_getV(N_tota, sum(N_tota)) * pt@pParams@V * 3

Vtheo_div2_even <- div2_getV(N_even, sum(N_even)) * pt@pParams@V
Vtheo_div2_skew <- div2_getV(N_skew, sum(N_skew)) * pt@pParams@V
Vtheo_div2_unev <- div2_getV(N_unev, sum(N_unev)) * pt@pParams@V
Vtheo_div2_tota <- div2_getV(N_tota, sum(N_tota)) * pt@pParams@V * 3


# ----- Coherence of formula & MC results -----

## check elementary regions

# variance estimate for HHI
base::plot(resV$r, resV$var_hhi_even, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_hhi_even, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_hhi_skew, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_hhi_skew, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_hhi_unev, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_hhi_unev, lty = "dashed", col = "blue")

# variance estimate for Entropy
base::plot(resV$r, resV$var_ent_even, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_ent_even, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_ent_skew, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_ent_skew, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_ent_unev, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_ent_unev, lty = "dashed", col = "blue")

# variance estimate for DIV1
base::plot(resV$r, resV$var_div1_even, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div1_even, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_div1_skew, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div1_skew, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_div1_unev, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div1_unev, lty = "dashed", col = "blue")

# variance estimate for DIV2
base::plot(resV$r, resV$var_div2_even, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div2_even, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_div2_skew, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div2_skew, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_div2_unev, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div2_unev, lty = "dashed", col = "blue")


## check composite regions

base::plot(resV$r, resV$var_hhi_tota, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_hhi_tota, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_ent_tota, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_ent_tota, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_div1_tota, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div1_tota, lty = "dashed", col = "blue")
base::plot(resV$r, resV$var_div2_tota, xlab = "MC iteration", ylab = "Var")
abline(h = Vtheo_div2_tota, lty = "dashed", col = "blue")


# ----- Bias of plugin-in estimation -----

K <- 5  # no. of elementary regions
M <- 10 # no. of categories per elementary region

R <- 5e+4 # no. of MC repetitions

set.seed(42)
X <- matrix(sample(10:100, K*M, replace = TRUE), nrow = M, ncol = K)
N_j <- rowSums(X)
N   <- sum(N_j)

# variance estimates based on unperturbed data
targ <- c(hhi  = hhi_getV(rowSums(X), N),
          ent  = ent_getV(rowSums(X), N),
          div1 = div1_getV(rowSums(X), N),
          div2 = div2_getV(rowSums(X), N))

# variance estimates based on perturbed data
test <- bias <- rbias <- matrix(nrow = R, ncol = 4)
for(r in 1:R) {
  
  Y <- matrix(noise_x(c(X)), nrow = M, ncol = K)
  
  test[r, 1] <- hhi_getV(rowSums(Y), N)
  test[r, 2] <- ent_getV(rowSums(Y), N)
  test[r, 3] <- div1_getV(rowSums(Y), N)
  test[r, 4] <- div2_getV(rowSums(Y), N)
}

for(i in 1:4) { 
  bias[, i]  <- cumsum(test[, i] - targ[i]) / (1:R) 
  rbias[, i] <- bias[, i] / targ[i]
}


## check for coherence

# HHI
base::plot(bias[, 1], ylim = c(-1e-10, 1e-10))
abline(h = 0)
abline(h = (4 * K * M * pt@pParams@V) / N^4, col = "blue", lty = "dashed")

# DIV1
base::plot(bias[, 3], ylim = c(-1e-6, 1e-6))
abline(h = 0)
abline(h = (4 * K * 10^4 * M * pt@pParams@V) / ((1 - 1/M)^2 * N^4), col = "blue", lty = "dashed")

