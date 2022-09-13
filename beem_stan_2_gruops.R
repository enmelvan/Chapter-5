library(cmdstanr)
library(beemStatic)
data(beemDemo)
library(rstan)

stan_beem <- 
'
functions {
 
  // dirichlet-multinomial probability mass function
  real dm_lpmf(int[] a_o, vector log_a_p, real alpha0, int S) {
    vector[S] a_summand = to_vector(a_o);
    vector[S] alpha_summand = softmax(log_a_p) * alpha0;
    real sum1 = sum(a_summand);
    real sum2 = sum(log(a_summand));
    real sum3 = sum(lbeta(alpha_summand, a_summand));
    return log(sum1) + lbeta(alpha0, sum1) - sum2 - sum3;
  }

// alternative expression for dirichlet-multinomial probability mass function
  real dm1_lpmf(int[] a_o, vector log_a_p, real alpha0, int S) {
    vector[S] a_summand = to_vector(a_o);
    vector[S] alpha_summand = softmax(log_a_p) * alpha0;
    real sum1 = sum(a_summand);
    real sum2 = lgamma(sum1 + 1) - sum(lgamma(a_summand + 1));
    real sum3 = sum(lgamma(alpha_summand + a_summand));
    real sum4 = sum(lgamma(alpha_summand));
    real sum5 = lgamma(alpha0 + sum1) - lgamma(alpha0);
    return sum2 + sum3 - sum4 - sum5;
  }

  // predict abundances calculated from b matrix
  // values <=0 are reassigning values <=1e-10 to penalise likelihood
  vector pred_a(matrix b, int S, int[] txn) {
    matrix[S, S] b_sub;
    vector[S] negs;
    vector[S] p;
    vector[S] log_p;
    real my_pred;
    for(j in 1:S) {
        for(i in 1:S) {
          b_sub[i, j] = b[txn[i], txn[j]];
        }
    }
    negs = rep_vector(-1.0, S);
    p = b_sub \\ negs;
    for(i in 1:S) {
      my_pred = p[i];
      if(my_pred <= 1e-10) {
        log_p[i] = (my_pred > -23) ? -23 + my_pred : -46;
      } else {
        log_p[i] = log(my_pred);
      }
    }
    return log_p;
  }

}

data {

  // abundance data
  int M; // # of samples
  int N; // # of taxa
  int G; // # of groups
  int E; // # of non-zero abundance values
  int P_dia; // # of posthoc tests (diagonal elements of b)
  int P_off; // # of posthoc tests (b off-diagonal elements of b)
  int S[M]; // # of taxa per sample
  int grp[M]; // group assignment for each sample
  real inv_med_a_o_tot; // inverse of median total abundance per sample
  int r[M, 2]; // index to first/last entry of txn/a_o
  int txn[E]; // list of taxa present in each sample
  int a_o[E]; // abundance data
  int a_o_raw[M, N]; // raw abundance data (not used)
  int b_dia_n[G]; // total number of taxa present in each group
  int b_dia_x[G, N]; // b_vec index for taxa present/absent
  int b_dia_t[G, N]; // list of taxa present in each group
  int b_off_n[G]; // # of taxon pairs present in each group
  int b_off_x[G, N*(N-1)]; // b_vec index to pairs present/absent
  int post_dia[P_dia]; // index for post hoc tests of carrying capacities
  int post_off[P_off,3]; // index for post hoc tests of interactions

}

parameters {

  vector<lower=1e-6>[b_dia_n[1]] b1_dia; // negative of b diagonal elements
  vector<lower=1e-6>[b_dia_n[2]] b2_dia; // negative of b diagonal elements
  vector[b_off_n[1]] b1_offa; // off-diagonal auxilliary var. 1
  vector[b_off_n[2]] b2_offa; // off-diagonal auxilliary var. 1
  vector<lower=1e-6>[b_off_n[1]] b1_offb; // off-diagonal auxilliary var. 2
  vector<lower=1e-6>[b_off_n[2]] b2_offb; // off-diagonal auxilliary var. 2
  real<lower=1> alpha0; // dm concentration parameter
  real<lower=1e-3> s; // scale parameter for cauchy distribution

}

transformed parameters {

  vector[N*N] b_vec[G]; // vector of unnormalised b elements
  vector[N] b_diag[G]; // b diagonal elements
  matrix[N, N] b[G]; // normalised interaction matrix b
  vector[E] log_a_p; // predicted abundances

  b_vec[1,b_dia_x[1,1:b_dia_n[1]]] = rep_vector(-1.0, b_dia_n[1]);
  b_diag[1,b_dia_t[1,1:b_dia_n[1]]] = b1_dia[1:b_dia_n[1]];
  if(b_dia_n[1] < N) {
    b_vec[1,b_dia_x[1,(b_dia_n[1] + 1):N]] = rep_vector(0, N - b_dia_n[1]);
    b_diag[1,b_dia_t[1,(b_dia_n[1] + 1):N]] = rep_vector(0, N - b_dia_n[1]);
  }
  b_vec[1,b_off_x[1,1:b_off_n[1]]] = b1_offa[1:b_off_n[1]] .*
                                       sqrt(b1_offb[1:b_off_n[1]]); // cauchy
  if(b_off_n[1] < N*(N-1)) {
      b_vec[1,b_off_x[1,(b_off_n[1] + 1):(N*(N-1))]] = 
        rep_vector(0, N*(N-1) - b_off_n[1]);
  }
  
  b_vec[2,b_dia_x[2,1:b_dia_n[2]]] = rep_vector(-1.0, b_dia_n[2]);
  b_diag[2,b_dia_t[2,1:b_dia_n[2]]] = b2_dia[1:b_dia_n[2]];
  if(b_dia_n[2] < N) {
    b_vec[2,b_dia_x[2,(b_dia_n[2] + 1):N]] = rep_vector(-999, N - b_dia_n[2]);
    b_diag[2,b_dia_t[2,(b_dia_n[2] + 1):N]] = rep_vector(-999, N - b_dia_n[2]);
  }
  b_vec[2,b_off_x[2,1:b_off_n[2]]] = b2_offa[1:b_off_n[2]] .*
                                       sqrt(b2_offb[1:b_off_n[2]]); // cauchy
  if(b_off_n[2] < N*(N-1)) {
      b_vec[2,b_off_x[2,(b_off_n[2] + 1):(N*(N-1))]] = 
        rep_vector(0, N*(N-1) - b_off_n[2]);
  }

  for(i in 1:G) {  
    for(k in 1:N) {
      for(j in 1:N) {
        b[i, j, k] = b_vec[i, (k - 1)*N + j] * b_diag[i, j];
      }
    }
  }

  // calculate abundance estimates for each taxon in each sample
  for(i in 1:M) {
      log_a_p[r[i,1]:r[i,2]] = pred_a(b[grp[i]], S[i], txn[r[i,1]:r[i,2]]);
  }

}

model {

  // priors
  alpha0 ~ exponential(inv_med_a_o_tot);
  b1_dia ~ inv_gamma(1.0, 1.0);
  b2_dia ~ inv_gamma(1.0, 1.0);
  s ~ exponential(1.0);
  b1_offa ~ std_normal();
  b2_offa ~ std_normal();
  b1_offb ~ inv_gamma(0.5,0.5*s*s);
  b2_offb ~ inv_gamma(0.5,0.5*s*s);

  // dirichlet multinomial likelihood
  for(i in 1:M) {
    a_o[r[i,1]:r[i,2]] ~ dm(log_a_p[r[i,1]:r[i,2]], alpha0, S[i]);
  }

}

generated quantities {

  vector<lower=0>[M] a_p_tot;
  vector[M] deviance;
  vector[P_dia] posthoc_dia;
  vector[P_off] posthoc_off;

  for(i in 1:M) {
    a_p_tot[i] = sum(exp(log_a_p[r[i,1]:r[i,2]]));
    deviance[i] = 2.0*(dm_lpmf(a_o[r[i,1]:r[i,2]] | 
                               log(to_vector(a_o[r[i,1]:r[i,2]])), alpha0, S[i]) -
                       dm_lpmf(a_o[r[i,1]:r[i,2]] | 
                               log_a_p[r[i,1]:r[i,2]], alpha0, S[i]));
  }
  
  for(i in 1:P_dia) {
    posthoc_dia[i] = 1/b_diag[1,post_dia[i]] - 1/b_diag[2,post_dia[i]];
  }
  
  for(i in 1:P_off) {
    posthoc_off[i] = b_vec[1,post_off[i,1]] - b_vec[2,post_off[i,1]];
  }

}
'


prep_func <- function(a_o) {

  # remove taxa/samples with no counts
  grp <- a_o[,1]
  a_o <- a_o[,-1]
  a_o <- a_o[,apply(a_o>0,2,sum)>0]
  xx <- apply(a_o>0,1,sum)>0
  a_o <- a_o[xx,]
  grp <- as.numeric(as.factor(grp[xx]))
  
  # change abundance to long format
  a_o_mlt <- data.frame(grp = rep(grp,times=ncol(a_o)),
                        smp = rep(1:nrow(a_o),times=ncol(a_o)),
                        txn = rep(1:ncol(a_o),each=nrow(a_o)),
                        a_o = as.vector(a_o))
  a_o_mlt <- a_o_mlt[a_o_mlt$a_o > 0,]
  a_o_mlt <- a_o_mlt[order(a_o_mlt$smp, a_o_mlt$txn),]
  
  b_dia <- as.vector(diag(TRUE, nrow = ncol(a_o))) == TRUE
  b_dia_n <- numeric()
  b_dia_x <- matrix(NA, nrow = ncol(a_o), ncol = 2)
  b_dia_t <- matrix(NA, nrow = ncol(a_o), ncol = 2)
  b_off_n <- numeric()
  b_off_x <- matrix(NA, nrow = ncol(a_o)*(ncol(a_o) - 1), ncol = 2)
  b_smp <- list()
  b_off <- as.vector(diag(TRUE, nrow = ncol(a_o))) == FALSE
  
  for(i in 1:length(unique(grp))) {
    tmp <- a_o[grp==i,]
    # index to off-diagonal entries of b; taxa and taxon pairs that aren't 
    # observed in a group are excluded from model fitting
    b_smp[[i]] <- matrix(NA, nrow = ncol(a_o), ncol = ncol(a_o))
    for(k in 1:ncol(a_o)) {
      for(j in k:ncol(a_o)) {
        b_smp[[i]][j,k] <- b_smp[[i]][k,j] <- sum(tmp[,j]>0 & tmp[,k]>0)>0
      }
    }
    b_dia_n[i] <- length(which(b_dia & as.vector(b_smp[[i]])))
    b_dia_t[1:b_dia_n[i],i] <- which(diag(b_smp[[i]]))
    if(b_dia_n[i] < ncol(a_o)) {
      b_dia_t[(b_dia_n[i] + 1):ncol(a_o),i] <- which(!diag(b_smp[[i]]))
    }
    tmp2 <- as.vector(b_smp[[i]])
    b_dia_x[1:b_dia_n[i],i] <- which(b_dia & tmp2)
    if(b_dia_n[i] < ncol(a_o)) {
      b_dia_x[(b_dia_n[i] + 1):ncol(a_o),i] <- which(b_dia & !tmp2)
    }
    b_off_n[i] <- length(which(b_off & tmp2))
    b_off_x[1:b_off_n[i],i] <- which(b_off & tmp2)
    if(sum(b_off) > b_off_n[i]) {
      b_off_x[(b_off_n[i] + 1):sum(b_off),i] <- which(b_off & !tmp2)
    }
  }

  # post hoc stuff
  post_dia <- which(diag(b_smp[[1]]) & diag(b_smp[[2]]));
  post_off_idx <- which(as.vector(b_smp[[1]]) & as.vector(b_smp[[2]]) & b_off)
  post_off_row <- as.vector(matrix(rep(1:ncol(a_o),times=ncol(a_o)),
                                   ncol=ncol(a_o)))[post_off_idx]
  post_off_col <- as.vector(matrix(rep(1:ncol(a_o),each=ncol(a_o)),
                                   ncol=ncol(a_o)))[post_off_idx]
  post_off <- cbind(post_off_idx,post_off_row,post_off_col)
 
  return(list(M = nrow(a_o),
              N = ncol(a_o),
              G = max(grp),
              E = nrow(a_o_mlt),
	            P_dia = length(post_dia),
              P_off = nrow(post_off),
              S = apply(a_o>0,1,sum),
              grp = grp,
              inv_med_a_o_tot = 1/median(apply(a_o,1,sum)),
              r = cbind(tapply(1:nrow(a_o_mlt),a_o_mlt$smp,min),
                        tapply(1:nrow(a_o_mlt),a_o_mlt$smp,max)),
              txn = a_o_mlt$txn,
              a_o = a_o_mlt$a_o,
              a_o_raw = a_o,
              b_dia_n = b_dia_n,
              b_dia_x = t(b_dia_x),
              b_dia_t = t(b_dia_t),
              b_off_n = b_off_n,
              b_off_x = t(b_off_x),
              post_dia = post_dia,
	            post_off = post_off))

}


init_func <- function(stan_data,chains=1) {
  
  df <- t(apply(stan_data$a_o_raw[stan_data$grp==1,], 1, function(x) x/sum(x)))
  bdg1 <- apply(df,2,function(x) quantile(x[x>0],0.9))
  df <- t(apply(stan_data$a_o_raw[stan_data$grp==2,], 1, function(x) x/sum(x)))
  bdg2 <- apply(df,2,function(x) quantile(x[x>0],0.9))
  
  ans <- list()
  for(i in 1:chains) {
    ans[[i]] <- list(b1_dia = c(sd(na.omit(bdg1))/na.omit(bdg1)),
                     b2_dia = c(sd(na.omit(bdg2))/na.omit(bdg2)),
                     b1_offa = rep(0,sum(stan_data$b_off_n[1])),
                     b2_offa = rep(0,sum(stan_data$b_off_n[2])),
                     b1_offb = rep(1,sum(stan_data$b_off_n[1])),
                     b2_offb = rep(1,sum(stan_data$b_off_n[2])),
                     alpha0 = 1/stan_data$inv_med_a_o_tot,
                     s=1)
  }
  return(ans)
  
}


# create database with 2 groups (numbered 1 and 2) in the first column
#a_o <- cbind(rep(c(1,2),times=250),t(beemDemo$dat.w.noise))

a_o <- as.matrix(read.csv('/Users/enamelvan/Desktop/PhD/zadnji_kod_rezultati/HMP_merged_datasets_two_groups.csv')[,-1])

# analyse full set of noisy beemStatic data (groups split in half)
stan_data <- prep_func(a_o)
stan_init <- init_func(stan_data)
f <- write_stan_file(stan_beem,getwd())
mod <- cmdstan_model(f)

# fit model using penalised maximum likelihood
fit_optim <- mod$optimize(data = stan_data,
                          seed = 123,
                          refresh = 1000,
                          init = stan_init,
                          iter = 1e6)

stan_init[[1]]$b1_dia <- data.frame(fit_optim$summary('b1_dia'))[,2]
stan_init[[1]]$b2_dia <- data.frame(fit_optim$summary('b2_dia'))[,2]
stan_init[[1]]$b1_offa <- data.frame(fit_optim$summary('b1_offa'))[,2]
stan_init[[1]]$b2_offa <- data.frame(fit_optim$summary('b2_offa'))[,2]
stan_init[[1]]$b1_offb <- data.frame(fit_optim$summary('b1_offb'))[,2]
stan_init[[1]]$b2_offb <- data.frame(fit_optim$summary('b2_offb'))[,2]
stan_init[[1]]$alpha0 <- data.frame(fit_optim$summary('alpha0'))[,2]
stan_init[[1]]$s <- data.frame(fit_optim$summary('s'))[,2]

# fit model using variational bayes
fit_varia <- mod$variational(data = stan_data,
                             seed=123,
                             refresh = 100,
                             init = stan_init)

fit_sampl <- mod$sample(data = stan_data,
                             seed=123,
                             refresh = 100,
                             init = stan_init,
                            chains = 1)


# posthoc results
d <- data.frame(fit_varia$summary('posthoc_dia'))[c(1,2,6,7)]
mat_b <- cbind(stan_data$post_off,data.frame(fit_varia$summary('posthoc_off'))[c(2,6,7)])

# good predictions on biomass (up to multiplicative constant)
# because parallel to 1:1 line on log-log scale
#plot(beemDemo$biomass.true,
#     data.frame(fit_optim$summary('a_p_tot'))[,2],
#     xlab='actual total biomass',ylab='predicted total biomass',log='xy')
#abline(0,1)
#cor.test(beemDemo$biomass.true,
#         data.frame(fit_varia$summary('a_p_tot'))[,2])
