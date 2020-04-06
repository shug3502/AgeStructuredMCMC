library(rstan)
library(dplyr)
library(ggplot2)
library(patchwork)
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write=TRUE)
set.seed(42)

#beta = get_beta_from_prem
#the following are from the table in keeling write up dated 2020-03-22
#TODO: reconcile the 17 age classes in these tables with the 21 categories described for Prem contact matrices
# prop_hospitalised <- 0.01*c(3.8,3.8,2.6,2.6,2.8,2.8,2.7,2.7,
#                             5.4,5.4,12.6,12.6,19.7,19.7,28.7,
#                             28.7,27.3,27.3,27.3,27.3,27.3)
# prop_die <- 0.01*c(3.8,3.8,3.8,3.8,3.8,3.8,3.8,4.0,4.5,5.6,7.8,
#                    11.3,16.9,23.2,29.1,34.8,53.5,53.5,53.5,53.5,53.5)

prop_hospitalised <- 0.01*c(3.8,19.7,28.7)
prop_die <- 0.01*c(3.8,16.9,53.5)
T <- 30
Tlockdown <- 7
Tfinal <- 50
nAges <- 3
nCompartments <- 16
beta <- purrr::map(1:4,function(x) runif(nAges^2) %>% matrix(.,ncol=nAges))

theta <-    c(alpha=1,
              tau=0.5,
              sigma=0.5,
              gamma=1/3,
              H=0.75,
              phi=0.7,
              det_rate=c(rep(0.2,2),rep(0.5,1)),
              delta=1) #alpha, tau, beta_scaling, gamma, H, phi, det_rate[a], delta
  N <- rep(10^6,nAges)
  y0 <- rep(0,nAges)
  for (a in 1:nAges){
    y0[a] = N[a] - 100;
    y0[nAges+a] = 100;
    for (i in 3:nCompartments){
      y0[(i-1)*nAges+a] = 0;
    }
  }
stan_sim_data <- list(T = T,
                      Tlockdown = Tlockdown,
                      ts1 = seq_len(Tlockdown),
                      ts2 = seq(from=Tlockdown+1,to=T),
                      nAges=nAges,
                      nCompartments=nCompartments,
                      beta = beta,
                      N = N,
                      prop_hospitalised=prop_hospitalised,
                      prop_die=prop_die,
                      theta = theta,
                      obs_noise = 1,
                      y0=y0)

#simulate from version of model
samples <- stan(file = 'src/stan_files/generate_age_structured_data.stan',
                data = stan_sim_data,
                seed = 42,
                chains = 1,
                iter = 100,
                algorithm = "Fixed_param"
)

#extract and summarise results in sensible format
z_df <- tidybayes::spread_draws(samples,z[t,ind])
z_median_df <- z_df %>% summarise(z = median(z)) %>%
  mutate(a = ((ind-1) %% nAges) + 1,
         compartment = floor((ind-1) / nAges) + 1)

#make nice plot
compartment_names <- c("Susceptible","Exposed (F)","Exposed (SD)",
                       "Exposed (SU)","Detected (F)",
                       "Detected (SD)", "Detected (SU)",
                       "Undetected (F)", "Undetected (S)",
                       "Exposed (Q)","Detected (QF)","Detected (QS)",
                       "Undetected (Q)","Hospitalised","Recovered",
                       "Deaths")
compartment_names_factor <- factor(compartment_names,levels = compartment_names)
p1 <- ggplot(data=z_median_df,aes(t,z,color=factor(a))) +
    geom_line() +
    geom_vline(data=tibble(Tlockdown=rep(Tlockdown,nCompartments),compartment=seq_len(nCompartments)),
               aes(xintercept=Tlockdown),linetype="dashed",color="grey") +
    facet_wrap(.~compartment_names_factor[compartment],scales="free_y") +
    theme_bw() +
    theme(legend.position = "None") +
    labs(x="Time (days)",
         y=NULL)
p1
ggsave(here::here(paste0("plots/simulated_data",Sys.time(),".eps")),device=cairo_ps)

sim_df <- tidybayes::spread_draws(samples,hospitalised[t,a],deaths[t,a])
sim_median_df <- sim_df %>% summarise(hospitalised=median(hospitalised),
                     deaths=median(deaths)) %>% tidyr::gather(compartment,num,-t,-a)
 p2 <- ggplot(sim_median_df,aes(t,num,color=factor(a))) +
    geom_line() +
    facet_wrap(.~compartment) +
   theme_bw() +
   theme(legend.position="None") +
   labs(x="Time (days)", y=NULL)
 p2

hospitalised <- (sim_median_df %>% filter(compartment=="hospitalised") %>%
  pull(num) %>% matrix(.,nrow=nAges)) %>% round() %>% t() #NB rounding is required because taking median produces fractional values
deaths <-  (sim_median_df %>% filter(compartment=="deaths") %>%
              pull(num) %>% matrix(.,nrow=nAges)) %>% round() %>% t()

############################
#inference using simulated data
stan_data <- list(T = T,
                  Tlockdown = Tlockdown,
                  Tfinal = Tfinal,
                  ts1 = seq_len(Tlockdown),
                  ts2 = seq(from=Tlockdown+1,to=T),
                  ts3 = seq(fromT+1,to=Tfinal),
                  nAges=nAges,
                  nCompartments=nCompartments,
                  beta = beta,
                  N = N,
                  prop_hospitalised=prop_hospitalised,
                  prop_die=prop_die,
                  y0=y0,
                  hospitalised = hospitalised,
                  deaths=deaths)
initF <- function() list(alpha=1,
              tau=0.5,
              sigma=1,
              gamma=1/3,
              H=0.75,
              phi=0.7,
              det_rate=c(rep(0.2,2),rep(0.5,1)),
              delta=1,
              obs_noise=1)

fit <- stan(file = 'src/stan_files/age_structured_seir.stan',
                data = stan_data,
                seed = 42,
                chains = 4,
                iter = 50,
                init=initF
)
saveRDS(fit,file=here::here('fits/age_structured_fit_synthetic_data.rds'))

draws <- as.data.frame(fit)
g1 <- draws %>% dplyr::select(alpha,tau,sigma,gamma,H,phi,delta) %>%
  bayesplot::mcmc_recover_intervals(., theta[c(1:6,10)]) +
  scale_y_continuous(limits=c(0,3)) +
  labs(title="Estimates based on data (Posterior)")
draws %>% dplyr::select(alpha,tau,sigma,gamma,H,phi,delta,obs_noise) %>%
  bayesplot::mcmc_pairs()
#bayesplot::mcmc_trace(draws,pars=names(theta)[1:6])

prior_samples <- stan(file = 'src/stan_files/priors.stan',
                            data=stan_data,
                             seed = 42,
                             chains = 4,
                             iter = 50
)
g2 <- prior_samples %>% as.data.frame() %>% dplyr::select(alpha,tau,sigma,gamma,H,phi,delta) %>%
  bayesplot::mcmc_recover_intervals(., theta[c(1:6,10)]) +
  scale_y_continuous(limits=c(0,3)) +
  labs(title="Estimates before using data (Prior)")

g1+g2+plot_layout(ncol=1)
ggsave(here::here(paste0("plots/estimates_from_prior_and_posterior",Sys.time(),".eps")),device=cairo_ps)
