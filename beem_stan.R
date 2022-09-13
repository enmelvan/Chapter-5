library(cmdstanr)
library(beemStatic)
data(beemDemo)
library(rstan)

stan_beem <- 
'
functions {
 
  // probability mass function for dirichlet-multinomial
  real dm_lpmf(int[] a_o, vector a_p, real alpha0, 
               real a_o_tot, int S, int[] idx) {
    vector[S] a_summand = to_vector(a_o[idx]);
    vector[S] alpha = softmax(log(a_p[idx]));
    vector[S] alpha_summand = alpha * alpha0;
    real sum1 = sum(lgamma(a_summand + alpha_summand));
    real sum2 = sum(lgamma(alpha_summand));
    real sum3 = sum(lgamma(a_summand + 1.0));
    return lgamma(alpha0) + lgamma(a_o_tot + 1.0) + sum1
           - lgamma(a_o_tot + alpha0) - sum2 - sum3;
  }

  // predict abundance per species using b matrix; values <=0 are
  // reassigning values <=1e-10 to penalised likelihood
  vector pred_a(matrix b, int S, int[] idx, int N) {
    matrix[S, S] b_sub;
    vector[S] negs;
    vector[S] p_sub;
    vector[N] p;
    real my_pred;
    for(j in 1:S) {
        for(i in 1:S) {
          b_sub[i, j] = b[idx[i], idx[j]];
        }
    }
    negs = rep_vector(-1.0, S);
    p_sub = b_sub \\ negs;
    for(i in 1:S) {
      if(p_sub[i] <= 0) {
        my_pred = p_sub[i];
        p_sub[i] = 10^(-10 + my_pred);
      }
    }
    p = rep_vector(0.0, N);
    p[idx] = p_sub;
    return p;
  }

}

data {

  // abundance data
  int M; // number of samples
  int N; // number of taxa
  int S[M]; // number of taxa per sample
  vector<lower=1>[M] a_o_tot; // observed total abundance per sample
  real inv_med_a_o_tot; // inverse of median total abundance per sample
  int taxa_N; // number of non-zero taxa entries in abundance data
  int taxa_list[taxa_N]; // list of taxa present in each sample
  int taxa_idx[M]; // index to first entry in taxa_list for each sample
  int a_o[M,N]; // observed abundances
  int b_diag_idx[N]; // index to diagonal B elements
  int b_off1_N; // # of fitted off-diagonal B elements
  int b_off1_idx[b_off1_N]; // index to fitted off-diagonal B elements
  int b_off0_N; // # of excluded off-diagonal B elements
  int b_off0_idx[b_off0_N]; // index to excluded off-diagonal B elements
  vector[b_off0_N] b_off0; // vector of 0s

}

transformed data {

  vector[N] neg_ones = rep_vector(-1.0, N);

}

parameters {

  vector<lower=1e-10>[N] b_diag; // negative of diagonal b elements
  vector[b_off1_N] b_off1_a; // off-diagonal b elements
  vector<lower=1e-10>[b_off1_N] b_off1_b; // off-diagonal b elements
  real<lower=1e-10> alpha0; // dm concentration parameter
  real<lower=1e-10> s; // square of scale parameter for cauchy distribution

}

transformed parameters {

  matrix[N, N] b; // normalised interaction matrix b
  vector[N*N] b_vec; // vector of b elements prior to normalisation
  vector[N] a_p[M]; // predicted abundances

  // b elements prior to normalisation by focal species carrying capacities
  b_vec[b_diag_idx] = neg_ones;
  b_vec[b_off1_idx] = b_off1_a .* sqrt(b_off1_b); // cauchy conversion
  b_vec[b_off0_idx] = b_off0;

  // normalised b matrix 
  for(j in 1:N) {
    for(i in 1:N) {
      b[i, j] = b_vec[(j - 1)*N + i] * b_diag[i];
    }
  }
  
  // calculate abundance estimates for each taxon in each sample
  for(i in 1:M) {
      a_p[i] = pred_a(b, S[i],
                      segment(taxa_list, taxa_idx[i], S[i]), N);
  }

}

model {

  // priors
  alpha0 ~ exponential(inv_med_a_o_tot);
  b_diag ~ inv_gamma(1.0, 1.0);
  s ~ exponential(1.0);
  b_off1_a ~ std_normal();
  b_off1_b ~ inv_gamma(0.5,0.5*s);

  // dirichlet multinomial likelihood
  for(i in 1:M) {
    a_o[i] ~ dm(a_p[i], alpha0, a_o_tot[i], S[i],
                segment(taxa_list, taxa_idx[i], S[i]));
  }

}

generated quantities {

  vector<lower=0>[M] a_p_tot;
  vector[M] deviance;
  
  for(i in 1:M) {
    a_p_tot[i] = sum(a_p[i]);
    deviance[i] = 2.0*(dm_lpmf(a_o[i] | to_vector(a_o[i]), alpha0,
                     a_o_tot[i], S[i],
                     segment(taxa_list, taxa_idx[i], S[i])) -
                     dm_lpmf(a_o[i] | a_p[i], alpha0,
                     a_o_tot[i], S[i],
                     segment(taxa_list, taxa_idx[i], S[i])));
  }

}
'


prep_func <- function(a_o) {

  # remove taxa/samples with no counts
  a_o <- a_o[,apply(a_o>0,2,sum)>0]
  a_o <- a_o[apply(a_o>0,1,sum)>0,]
  rownames(a_o) <- paste0('smp_',formatC(1:nrow(a_o),width=3,flag='0'),
                          '_',rownames(a_o))
  colnames(a_o) <- paste0('txn_',formatC(1:ncol(a_o),width=3,flag='0'),
                          '_',colnames(a_o))
  
  # create list of taxa present in each sample
  a_o_melt <- data.frame(smp = rep(rownames(a_o),times=ncol(a_o)),
                         txn = rep(colnames(a_o),each=nrow(a_o)),
                         cnt = as.vector(a_o),stringsAsFactors = TRUE)
  a_o_melt <- a_o_melt[order(a_o_melt[,1],a_o_melt[,2]),]
  a_o_melt[,2] <- as.numeric(a_o_melt[,2])
  a_o_melt <- a_o_melt[a_o_melt[,3]>0,]
  taxa_N <- nrow(a_o_melt)
  taxa_list <- a_o_melt[,2]
  taxa_idx <- tapply(1:nrow(a_o_melt),a_o_melt[,1],min)
  
  # index to diagonal entries of interaction matrix b
  b_diag_idx <- which(as.vector(diag(TRUE,nrow=ncol(a_o)))==TRUE)
  
  # index to off-diagonal entries of b; taxa pairs that don't 
  # co-occur are excluded from model fitting
  b_off_idx <- matrix(NA,ncol(a_o),ncol(a_o))
  for(i in 1:ncol(a_o)) {
    for(j in i:ncol(a_o)) {
      b_off_idx[i,j] <- b_off_idx[j,i] <- sum(a_o[,i]>0 & a_o[,j]>0)>0
    }
  }
  b_off1_idx <- which(as.vector(b_off_idx & 
                                  diag(TRUE,nrow=ncol(a_o))==FALSE)==TRUE)
  b_off0_idx <- which(as.vector(!b_off_idx &
                                  diag(TRUE,nrow=ncol(a_o))==FALSE)==TRUE)
  
  return(list(M = nrow(a_o),
              N = ncol(a_o),
              S = apply(a_o>0,1,sum),
              a_o_tot = apply(a_o,1,sum),
              inv_med_a_o_tot = 1/median(apply(a_o,1,sum)),
              taxa_N = taxa_N,
              taxa_list = taxa_list,
              taxa_idx  = taxa_idx,
              a_o = a_o,
              b_diag_idx = b_diag_idx,
              b_off1_N = length(b_off1_idx),
              b_off1_idx = b_off1_idx,
              b_off0_N = length(b_off0_idx),
              b_off0_idx = b_off0_idx,
              b_off0 = rep(0,length = length(b_off0_idx))))

}


init_func <- function(stan_data,chains=1) {
  
  df <- t(apply(stan_data$a_o, 1, function(x) x/sum(x)))
  bdg <- 1/apply(df,2,function(x) quantile(x[x>0],0.9))
  ans <- list()
  for(i in 1:chains) {
    ans[[i]] <- list(b_diag = bdg/sd(bdg),
                     b_off1_a = rep(0,stan_data$b_off1_N),
                     b_off1_b = rep(1,stan_data$b_off1_N),
                     alpha0 = 1/stan_data$inv_med_a_o_tot,
                     s=0.1)
  }
  return(ans)

}


# analyse full set of noisy beemStatic data
stan_data <- prep_func(t(beemDemo$dat.w.noise))

#ADD MY DATA - HMP FIRST ONE
my_data <- readRDS("/Users/enamelvan/Desktop/PhD/data/HMP_final_final.RDS")
(my_phylum <- tax_glom(my_data, taxrank="Rank1"))
(my_phylum <- tax_glom(my_data, taxrank="Phylum"))
(my_order <- tax_glom(my_data, taxrank="Rank3"))
(my_family <- tax_glom(my_data, taxrank = "Rank4"))
(my_genus <- tax_glom(my_data, taxrank = "Rank5"))
(my_species <- tax_glom(my_data, taxrank = "Rank6"))
(my_strain <- tax_glom(my_data, taxrank = "Rank7"))

my_prep_data <- taxon_sub

lean_names    <- rownames(my_data@sam_data)[my_data@sam_data$Disease == "OB"]
mySubset      <- my_prep_data@otu_table@.Data[,colnames(my_prep_data@otu_table@.Data) %in% lean_names]
t <- cbind(1, t(mySubset))
lean_names    <- rownames(my_data@sam_data)[my_data@sam_data$Disease == "OB"]
mySubset      <- my_prep_data@otu_table@.Data[,colnames(my_prep_data@otu_table@.Data) %in% lean_names]
q <- cbind(2, t(mySubset))
m  <- rbind(t,q)
m <- t(m)

stan_data <- prep_func(t(m))




stan_init <- init_func(stan_data)
f <- write_stan_file(stan_beem,getwd())
mod <- cmdstan_model(f)


# fit model using penalised maximum likelihood
# runs fast (3 minutes) because we only find modes of posteriors
# new code models cauchy as scale mixture of normal distributions
fit_optim <- mod$optimize(data = stan_data,
                          seed=123,
                          refresh = 1000,
                          init=stan_init,
                          iter=1e6)

stan_init[[1]]$bdiag <- data.frame(fit_optim$summary('b_diag'))[,2]
stan_init[[1]]$b_off1_a <- data.frame(fit_optim$summary('b_off1_a'))[,2]
stan_init[[1]]$b_off1_b <- data.frame(fit_optim$summary('b_off1_b'))[,2]
stan_init[[1]]$alpha0 <- data.frame(fit_optim$summary('alpha0'))[,2]
stan_init[[1]]$s <- data.frame(fit_optim$summary('s'))[,2]

#Do this instead of optim!!!!
fit_varia <- mod$variational(data = stan_data,
                             seed=123,
                             refresh = 1000,
                             init=stan_init,
                             iter=1e6)

#Do this instead of optim!!!!
#fit_varia <- mod$variational(data = stan_data)

#LOAD FIT OPTIM HMP
fit_optim <- readRDS("/Users/enamelvan/Desktop/PhD/zadnji_kod_rezultati/hmp_genus_fit_optim.RDS")
# good predictions on biomass (up to multiplicative constant)
# because parallel to 1:1 line on log-log scale
plot(beemDemo$biomass.true,
     data.frame(fit_optim$summary('a_p_tot'))[,2],
     xlab='actual total biomass',ylab='predicted total biomass',log='xy')
abline(0,1)
cor.test(beemDemo$biomass.true,
         data.frame(fit_optim$summary('a_p_tot'))[,2])

# good predictions on biomass (up to multiplicative constant)
# because parallel to 1:1 line on log-log scale
plot(my_prep_data@otu_table@.Data,
     data.frame(fit_optim$summary('a_p_tot'))[,2],
     xlab='actual total biomass',ylab='predicted total biomass',log='xy')
abline(0,1)
cor.test(my_data@otu_table@.Data,
         data.frame(fit_optim$summary('a_p_tot'))[,2])

# good fit on diagonal of b matrix
# because parallel to 1:1 line on log-log scale
plot(beemDemo$scaled.params$a.truth,
     1/data.frame(fit_optim$summary('b_diag'))[,2],
     xlab='actual carrying capacity',ylab='predicted carrying capacity',log='xy')
abline(0,1)
cor.test(beemDemo$scaled.params$a.truth,
         1/data.frame(fit_optim$summary('b_diag'))[,2])

# bayesian fit
fit <- mod$sample(data = stan_data,
                  seed = 123,
                  refresh = 200,
                  init = stan_init,  #change init to output from optim, change the name
                  chains = 1,
                  max_treedepth = 15)


plot(beemDemo$biomass.true[idx.r],
     data.frame(fit_optim$summary('a_p_tot'))[,2],
     xlab='observed (z score)',ylab='predicted (z score)',log='xy')
abline(0,1)


plot(beemDemo$scaled.params$a.truth,
     1/data.frame(fit$summary('b_diag'))[,2],
     xlab='observed (z score)',ylab='predicted (z score)',log='xy')
abline(0,1)
