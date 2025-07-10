pkgs = c("ggplot2", "tidyr", "dplyr", "stringr" ,"rstanarm", "rstan", "bayesrules", "tidybayes", "bayesplot", "gridExtra", "extraDistr", "latex2exp")
lapply(pkgs, library, character.only = TRUE)

# load models and data from phase2

load("../scripts/shedding_df.Rdata")
load("../scripts/model_phase_2.Rdata")
model_inital_pool = model_sim
load("../scripts/df_counts.Rdata")

#phase3 data

phase_3A_df = read.csv("../data/phase_3A.csv")
names(phase_3A_df)[8] = "blue_synthetic_A"
phase_3C_df = read.csv("../data/phase_3C.csv")

set.seed(0)

# clean phase 3 data
phase_3_longer = function(df){
  longer = df %>% 
    pivot_longer(cols = c(8:10), names_to = "fibre", values_to = "count")
  return(longer)
}

phase_3_counts = function(df_long){
  counts = df_long %>% 
    group_by(Participant, replicate, role, fibre) %>% 
    summarize(sum = sum(count))
  return(counts)
}

# join both a and C tables in long form
phase_3_long = rbind(phase_3_longer(phase_3A_df),phase_3_longer(phase_3C_df))

# Calculate sums
phase_3_counts_seperate = phase_3_counts(phase_3_long)

# self join to put inter and receiver in same row total sums and proportions
phase_3_counts = phase_3_counts_seperate %>% 
  merge(phase_3_counts_seperate,
        by = c("Participant", "replicate", "fibre")) %>% 
  filter (role.x == "intermediate", role.y =="receiver") %>% 
  mutate(sum = sum.x + sum.y, prop_inter = round(sum.x/sum,3), prop_rec = round(sum.y/sum,3), phase = "after", 
         intermediate = case_when(replicate == 1|replicate == 2|replicate == 3 ~"smooth",
                                  .default = "rough"))

phase_3_counts$fibre = factor(phase_3_counts$fibre, levels = unique(phase_3_counts$fibre))

# colour selection
colour_fibres = c("blue", "purple", "red", "black", "blue", "red")

phase3_total_boxplot = phase_3_counts %>% 
  ggplot(aes(x = fibre, y = log10(sum), fill = fibre))+
  geom_boxplot()+
  scale_fill_manual(values = colour_fibres)

phase3_total_boxplot

phase_2_counts = df_counts %>% 
  filter(replicate != "01") %>% 
  filter((participant == "A"& grepl("_A", fibre))| (participant =="C" & grepl("_C", fibre))) %>% 
  mutate(phase = "before")

combine = rbind(phase_2_counts[,-c(1,2)], phase_3_counts[,c(3,8,11)])

complare_total_boxplot_pub = 
  combine %>% 
  mutate(phase = factor (phase, levels = c("before", "after"))) %>% 
  ggplot(aes(y = log10(sum), x = phase, fill = fibre))+
  geom_boxplot(size = 2)+
  scale_fill_manual(values = colour_fibres)+
  facet_wrap(vars(fibre))+
  xlab("intermediate (before transfer) VS intermediate + recipient (after transfer)")+
  ylab(TeX("$log_{10}(count)$"))+
  theme(text = element_text(size = 30),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(l = 55))

#compare_totals_boxplot width = 1500
complare_total_boxplot_pub

combine %>% 
  group_by(fibre, phase) %>% 
  summarize(mean = mean(sum), sd = sd(sum), median = median(sum), min = min(sum), max = max(sum))


prop_boxplot = phase_3_counts %>% 
  ggplot(aes(x = fibre, y = prop_rec, fill = fibre))+
  geom_boxplot()+
  scale_fill_manual(values = colour_fibres)

prop_boxplot

phase_3_counts %>% 
  summarise(mean = mean(prop_rec), sd = sd(prop_rec), median = median(prop_rec), min = min(prop_rec), max = max(prop_rec))

prop_roughness_boxplot = phase_3_counts %>% 
  ggplot(aes(x = intermediate, y = prop_rec, fill = fibre))+
  geom_boxplot()+
  scale_fill_manual(values = colour_fibres)

prop_roughness_boxplot

prop_roughness_boxplot_pub = phase_3_counts %>% 
  ggplot(aes(x = fibre, y = prop_rec, fill = fibre))+
  geom_boxplot(size = 2)+
  facet_grid(col = vars(intermediate))+
  scale_fill_manual(values = colour_fibres)+
  ylab("Proportion transferred onto recipient")+
  xlab("Target fibre")+
  theme(text = element_text(size = 30),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(l = 55))

#prop_roughness_boxplot 1500 x1400
prop_roughness_boxplot_pub

phase_3_counts %>% 
  group_by(intermediate) %>% 
  summarize(mean = mean(prop_rec), sd = sd(prop_rec), median = median(prop_rec), min = min(prop_rec), max = max(prop_rec))


phase_3_counts %>% 
  summarize(mean = mean(sum.y), sd = sd(sum.y), median = median(sum.y), min = min(sum.y), max = max(sum.y))

phase_3_counts %>% 
  group_by(intermediate) %>% 
  summarize(mean = mean(sum.y), sd = sd(sum.y), median = median(sum.y), min = min(sum.y), max = max(sum.y))


# Model: scenario 1
# Count based

mean_prior = 5
sd_prior = 5

p = mean_prior /sd_prior^2
r = mean_prior^2/(sd_prior^2-mean_prior)

alpha_nbinom = 3
beta_nbinom = 10

plot(x = seq(0,1, by = 0.01), y = dbeta(x = seq(0,1, by = 0.01), 
                                        shape1 = alpha_nbinom, 
                                        shape2 = beta_nbinom))
# transfer data
mean = 29
sd = 17
n = 36

# post

alpha_star_nbinom = alpha_nbinom + r*n
beta_star_nbinom = beta_nbinom + sum(as.numeric(phase_3_counts$sum.y))


p_star_nbinom = qbeta(p = 0.5, shape1 = alpha_star_nbinom, 
                      shape2 = beta_star_nbinom)
y_bar_star_nbinom = qnbinom(p = 0.5, prob = p_star_nbinom, size = r)

# post pred
total_nb_beta_sim = rbnbinom(n=1000, size = 1, 
                             alpha = y_bar_star_nbinom, 
                             beta = beta_star_nbinom)


# Prop based
# phase 3 total parameters
prop_ybar_total = 0.35
prop_sigma_total = 0.08

# priors p around 0.3
alpha_binom = 5
beta_binom = 10


alpha_star_binom = alpha_binom + sum(as.numeric(phase_3_counts$sum.y))
beta_star_binom = beta_binom + sum(as.numeric(phase_3_counts$sum.x))


p_draws_binom = rbeta(1000, shape1 = alpha_star_binom, shape2 = beta_star_binom)


draw_phase_2 = rnbinom(n = 1000, mu = 99, size = 2)

total_binom = mapply(rbinom, n=1, p = p_draws_binom, size = draw_phase_2)

total_sim_df = data.frame(model = c(rep("data", nrow(phase_3_counts)),
                                    rep("proportion-based", 1000),
                                    rep ("count-based", 1000)),
                          sum = as.numeric(
                            c(phase_3_counts$sum.y, 
                              total_binom,
                              total_nb_beta_sim)))


density_total_model_prop_pub = total_sim_df %>% 
  ggplot(aes(x = sum, fill = model))+
  geom_density(alpha = 0.5, size = 1)+
  ylab(TeX("Density"))+
  xlab("Number of fibres transferred")+
  labs(fill = "Model:")+
  scale_fill_discrete()+
  theme(text = element_text(size = 30),
        legend.position = "top",
        legend.title = element_text(size = 25))

plot_total_compare_model_data_plot = 
  total_sim_df %>% 
  ggplot(aes(x = model, y = log10(sum)))+
  geom_boxplot(fill = "skyblue", size = 1)+
  ylab(TeX("$log_{10}$(count)"))+
  xlab("Model")+
  theme(text = element_text(size = 30))


# positional distribution
positional_df = phase_3_long %>% 
  mutate(position = paste0(left.right, top.bottom)) %>% 
  group_by(Participant, replicate, fibre, role, position) %>% 
  summarise (sum = sum(count))

positional_merge_df = positional_df %>% 
  merge(positional_df,
        by = c("Participant", "replicate", "fibre", "position"))%>%
  filter(role.x == "receiver", role.y == "intermediate") %>% 
  mutate(prop = sum.x/(sum.y + sum.x))

positional_plot = positional_merge_df %>% 
  ggplot(aes(x = position, y = prop))+
  geom_boxplot()

positional_corr_plot = positional_merge_df %>% 
  ggplot(aes(x = sum.x, y = sum.y))+
  geom_point(aes(colour = position))+
  geom_smooth(method = "lm", formula =  y~x)


positional_individual_mean_df = positional_merge_df %>% 
  ungroup() %>% 
  group_by(Participant, replicate, fibre) %>% 
  summarise(mean = mean(sum.x))

positional_individual_means_tbl = cbind(positional_individual_mean_df$mean,
                                        positional_individual_mean_df$mean,
                                        positional_individual_mean_df$mean)


positional_individual_df = positional_merge_df %>% 
  ungroup() %>% 
  group_by(Participant, replicate, fibre) %>%
  mutate(difference = sum.x - positional_individual_mean_df$mean)

positional_corr_plot_pub = positional_merge_df %>% 
  mutate(surface = ifelse(replicate == "1" |replicate == "2" |replicate == "3", "smooth", "rough")) %>% 
  ggplot(aes(x = sum.x, y = sum.y))+
  geom_point(size = 3, aes(colour = position))+
  geom_smooth(method = "lm", size = 2)+
  # facet_grid(rows = vars(surface))+
  xlab("Number of recovered fibres on recipient")+
  ylab("Number of recovered fibres on intermediate")+
  labs(colour = "Zone")+
  scale_color_manual(labels = c ("1", "2", "3", "4"), 
                     values = c("#000000", "#E69F00", "#56B4E9", "#009E73"))+
  theme(text = element_text(size = 20))

# position 1000
positional_corr_plot_pub


# model 1, width 1750
grid.arrange(density_total_model_prop_pub, plot_total_compare_model_data_plot, nrow = 1)


# Model: scenario 2

# sum by smooth/rough
sum2 = phase_3_counts %>% 
  group_by(intermediate) %>% 
  summarise(sum.y = sum(sum.y), sum.x= sum(sum.x), sum = sum(sum))

# posteriors
alpha_star_nbinom2 = alpha_nbinom + r*(n/2)
beta_star_nbinom2_smooth = beta_nbinom + sum2$sum.y[2]
beta_star_nbinom2_rough = beta_nbinom + sum2$sum.y[1]

p_star_nbinom2_smooth = qbeta(p = 0.5, shape1 = alpha_star_nbinom2, 
                              shape2 = beta_star_nbinom2_smooth)
y_bar_star_nbinom2_smooth = qnbinom(p = 0.5, prob = p_star_nbinom2_smooth,
                                    size = r)

p_star_nbinom2_rough = qbeta(p = 0.5, shape1 = alpha_star_nbinom2, 
                             shape2 = beta_star_nbinom2_rough)
y_bar_star_nbinom2_rough = qnbinom(p = 0.5, prob = p_star_nbinom2_rough,
                                   size = r)

#pred 
smooth_nb_beta_sim = rbnbinom(n=1000, size = 1, 
                              alpha = y_bar_star_nbinom2_smooth, 
                              beta = beta_star_nbinom2_smooth)
rough_nb_beta_sim = rbnbinom(n=1000, size = 1, 
                             alpha = y_bar_star_nbinom2_rough, 
                             beta = beta_star_nbinom2_rough)

alpha_star_binom2_smooth = 
  alpha_binom + sum2$sum.y[2]
beta_star_binom2_smooth  = 
  beta_binom + sum2$sum.x[2]

alpha_star_binom2_rough = 
  alpha_binom + sum2$sum.y[1]
beta_star_binom2_rough  = 
  beta_binom + sum2$sum.x[1]

p_draws_binom2_smooth = rbeta(1000, 
                              shape1 = alpha_star_binom2_smooth, 
                              shape2 = beta_star_binom2_smooth)

p_draws_binom2_rough = rbeta(1000, 
                             shape1 = alpha_star_binom2_rough, 
                             shape2 = beta_star_binom2_rough)

total_binom2_smooth = mapply(rbinom, n=1, p = p_draws_binom2_smooth, size = draw_phase_2)
total_binom2_rough = mapply(rbinom, n=1, p = p_draws_binom2_rough, size = draw_phase_2)

phase_3_counts_texture = phase_3_counts %>% 
  group_by(intermediate) %>% 
  arrange(intermediate)

text_sim_df = data.frame(model = c(rep("data", nrow(phase_3_counts)),
                                   rep("proportion-based", 2000),
                                   rep ("count-based", 2000)),
                         texture = c(rep("smooth", nrow(phase_3_counts)/2),
                                     rep("rough", nrow(phase_3_counts)/2),
                                     rep("smooth", 1000),
                                     rep("rough", 1000),
                                     rep("smooth", 1000),
                                     rep("rough", 1000)),
                         sum = as.numeric(
                           c(
                             phase_3_counts_texture$sum.y[(nrow(phase_3_counts)/2+1): nrow(phase_3_counts)],
                             phase_3_counts_texture$sum.y[1:(nrow(phase_3_counts)/2)],
                             total_binom2_smooth,
                             total_binom2_rough,
                             smooth_nb_beta_sim,
                             rough_nb_beta_sim
                           )))


density_text_model_prop_pub = text_sim_df %>% 
  ggplot(aes(x = sum, fill = model))+
  geom_density(alpha = 0.5, size = 1)+
  facet_grid(cols = vars(texture))+
  ylab(TeX("Density"))+
  xlab("Number of fibres transferred")+
  labs(fill = "Model:")+
  scale_fill_discrete()+
  theme(text = element_text(size = 30),
        legend.position = "top",
        legend.title = element_text(size = 25))

plot_total_text_compare_model_data_pub = 
  text_sim_df %>% 
  ggplot(aes(x = model, y = log10(sum)))+
  geom_boxplot(fill = "skyblue", size = 1)+
  facet_grid(cols = vars(texture))+
  ylab(TeX("$log_{10}$(count)"))+
  xlab("Model")+
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1))

# model 2, height 1250
grid.arrange(density_text_model_prop_pub, plot_total_text_compare_model_data_pub, nrow = 2)

## Model: scenario 3

#increase size of p_draws
p_draws_binom = rbeta(20000, shape1 = alpha_star_binom, shape2 = beta_star_binom)
filtered_shedding_df = shedding_df %>% 
  filter(grepl("_A", target)|grepl("_C", target)) # filter for a and b

# predict draws for size from model 1
shedding_newdata = data.frame(average = filtered_shedding_df$average)
size_draws = posterior_predict(model_sim,newdata = shedding_newdata)



draw_binom_by_shedding = function(size_draws, p_draws){
  mapply(rbinom, size = size_draws, 
         prob = p_draws, 
         MoreArgs = list(n = 1))
}


draws_binom_3_df = data.frame(apply(size_draws, 2, draw_binom_by_shedding, p_draws = p_draws_binom))



# manula reorder
draws_binom_3_df =draws_binom_3_df[,c(1,3,2,4,5)]

names(draws_binom_3_df) = levels(phase_3_counts$fibre)[-4]

draws_binom_3_df_long = pivot_longer(
  draws_binom_3_df, cols = c(1:5), names_to = "fibre", values_to = "sum.y"
) %>% 
  mutate(model = "model")


scenario3_df = phase_3_counts %>% 
  ungroup() %>% 
  select(c(fibre, sum.y)) %>% 
  filter(fibre != "black_animal_C") %>%
  mutate(model = "data") %>% 
  rbind(draws_binom_3_df_long)


scenario3_df_dens_plot_pub = scenario3_df %>% 
  ggplot(aes(x = log10(sum.y), fill = model))+
  geom_density(alpha = 0.5)+
  facet_wrap(~fibre)+
  ylab(TeX("Density"))+
  xlab(TeX("$log_{10}$(number of fibres transferred)"))+
  scale_fill_discrete()+
  theme(text = element_text(size = 30),
        legend.position = "top",
        legend.title = element_blank())

scenario3_df_plot_pub = scenario3_df %>% 
  ggplot(aes(x = fibre, y = log10(sum.y), fill = model))+
  geom_boxplot(size = 1, alpha = 0.75)+
  ylab(TeX("$log_{10}(count)$"))+
  xlab("Fibre")+
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# model 3 width = 1750
grid.arrange(scenario3_df_dens_plot_pub, scenario3_df_plot_pub, nrow = 1)


# Model: scenario 4

p_draws_binom2_smooth = rbeta(20000, 
                              shape1 = alpha_star_binom2_smooth, 
                              shape2 = beta_star_binom2_smooth)

p_draws_binom2_rough = rbeta(20000, 
                             shape1 = alpha_star_binom2_rough, 
                             shape2 = beta_star_binom2_rough)
draws_binom_4_df_smooth = data.frame(apply(size_draws, 2, draw_binom_by_shedding, p_draws = p_draws_binom2_smooth))
draws_binom_4_df_rough = data.frame(apply(size_draws,2, draw_binom_by_shedding, p_draws = p_draws_binom2_rough))

draws_binom_4_df_smooth = draws_binom_4_df_smooth[,c(1,3,2,4,5)]
draws_binom_4_df_rough = draws_binom_4_df_rough[,c(1,3,2,4,5)]

names(draws_binom_4_df_smooth) = levels(phase_3_counts$fibre)[-4]
names(draws_binom_4_df_rough) = levels(phase_3_counts$fibre)[-4]

draws_binom_4_df_long_smooth = pivot_longer(
  draws_binom_4_df_smooth, cols = c(1:5), names_to = "fibre", values_to = "sum.y"
) %>% 
  mutate(intermediate = "smooth", model = "model")

draws_binom_4_df_long_rough = pivot_longer(
  draws_binom_4_df_rough, cols = c(1:5), names_to = "fibre", values_to = "sum.y"
) %>% 
  mutate(intermediate = "rough", model = "model")


scenario4_df_dens_plot_pub = scenario4_df %>% 
  ggplot(aes(x = log10(sum.y), fill = model))+
  geom_density(alpha = 0.5)+
  facet_grid(intermediate~fibre)+
  ylab(TeX("Density"))+
  xlab(TeX("$log_{10}$(number of fibres transferred)"))+
  scale_fill_discrete()+
  theme(text = element_text(size = 30),
        legend.position = "top",
        legend.title = element_blank())


scenario4_df_plot_pub = scenario4_df %>% 
  ggplot(aes(x = fibre, y = log10(sum.y), fill = model))+
  geom_boxplot(size = 1, alpha = 0.75)+
  facet_grid(~intermediate)+
  ylab(TeX("$log_{10}(count)$"))+
  xlab("Fibre")+
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# model 4 height = 1250
grid.arrange(scenario4_df_dens_plot_pub, scenario4_df_plot_pub,
             nrow = 2)

# model evaluations

total_sim_df %>% 
  filter(model != "1data") %>% 
  group_by(model) %>% 
  summarize (mean = mean(sum), sd = sd(sum), error = abs(mean((sum-30))), 
             error_scaled = error/sd(sum),
             ci_upper = quantile(sum, 0.975),
             ci_lower = quantile(sum, 0.025))

text_sim_df %>% 
  filter(model != "1data") %>% 
  group_by(model, texture) %>% 
  summarize (mean = mean(sum), sd = sd(sum), error = 
               case_when(texture == "smooth" ~ abs(mean(sum-35)),
                         texture == "rough" ~ abs(mean-23)), 
             error_scaled = error/sd(sum),
             ci_upper = quantile(sum, 0.975),
             ci_lower = quantile(sum, 0.025)) %>% 
  unique()

fibre_means = scenario3_df %>% 
  filter(model!= "model") %>% 
  group_by(fibre) %>% 
  summarise(mean(sum.y))

scenario3_df %>% 
  filter(model != "data") %>% 
  group_by(fibre) %>% 
  summarize (mean = mean(sum.y), sd = sd(sum.y), 
             ci_upper = quantile(sum.y, 0.975),
             ci_lower = quantile(sum.y, 0.025)) %>% 
  mutate(error = abs(mean-fibre_means$`mean(sum.y)`), error_scaled = error/sd) %>% 
  unique()


fibre_means_text = scenario4_df %>% 
  filter(model!= "model") %>% 
  group_by(fibre,intermediate) %>% 
  summarise(mean(sum.y))

scenario4_df %>% 
  filter(model != "data") %>% 
  group_by(fibre, intermediate) %>% 
  summarize (mean = mean(sum.y), sd = sd(sum.y), 
             ci_upper = quantile(sum.y, 0.975),
             ci_lower = quantile(sum.y, 0.025)) %>% 
  ungroup() %>% 
  mutate(error = abs(mean-fibre_means_text$`mean(sum.y)`), error_scaled = error/sd) %>% 
  unique()

# compare to first pub
dir = c(13,12,4,1)
interval = c("none","small","medium","large")
draws_old = rmultinom(1000, sum(dir),prob =dir)

draws_binom = mapply(rbinom, n=sum(dir), p = p_draws_binom, size = draw_phase_2)


none = function(x){
  return(sum(x == 0))
}
small = function(x){
  return(sum(x > 1&x<=5))
}
medium = function(x){
  return(sum(x > 5&x<=50))
}

large = function(x){
  return(sum(x>50))
}

draws_new_none =  apply(draws_binom*0.1, 2, none)
draws_new_small = apply(draws_binom*0.1, 2, small)
draws_new_medium = apply(draws_binom*0.1, 2, medium)
draws_new_large = apply(draws_binom*0.1, 2, large)


draws_new =cbind(draws_new_none, draws_new_small,draws_new_medium, draws_new_large, rep("new", 1000))


# p values

list = vector(mode = "list", length = 10000)

n = nrow(phase_3_counts)

observed = c(median(phase_3_counts$sum.y), min(phase_3_counts$sum.y), max(phase_3_counts$sum.y))

count_sim_generator = function(x){
  sim = rbnbinom(n=n, size = 1, 
                 alpha = y_bar_star_nbinom, 
                 beta = beta_star_nbinom)
  
  return(c(median(sim), min(sim), max(sim)))
}

prop_sim_generator = function (X){
  p_draws_binom = rbeta(n, shape1 = alpha_star_binom, shape2 = beta_star_binom)
  n_draws = rnbinom(n = n, mu = 99, size = 2)
  
  sim = mapply(rbinom, n=1, p = p_draws_binom, size = n_draws)
  
  return(c(median(sim), min(sim), max(sim)))
}

count_based_1_sims = lapply(list, count_sim_generator) %>% 
  unlist() %>% 
  matrix(ncol = 10000, nrow = 3) %>% 
  t() %>% 
  data.frame()
prop_based_1_sims = lapply(list, prop_sim_generator)%>% 
  unlist() %>% 
  matrix(ncol = 10000, nrow = 3) %>% 
  t() %>% 
  data.frame()


# p values 

pppv_count_1 = count_based_1_sims %>%
  mutate(mean = X1 - observed[1], min = X2 - observed[2],
         max = X3-observed[3]) %>% 
  summarise (median = sum(mean > 0)/10000, min = sum(min>0)/10000,
             max = sum(max > 0)/10000) %>% 
  print()



pppv_prop_1 = prop_based_1_sims %>%
  mutate(mean = X1 - observed[1], min = X2 - observed[2],
         max = X3-observed[3]) %>% 
  summarise (median = sum(mean > 0)/10000, min = sum(min>0)/10000,
             max = sum(max > 0)/10000) %>% 
  print()


n = 18

observed_smooth = phase_3_counts %>% 
  filter (intermediate == "smooth") %>% 
  summarise( median = median(sum.y), min = min(sum.y), max = max(sum.y))

observed_rough = phase_3_counts %>% 
  filter (intermediate == "rough") %>% 
  summarise( median = median(sum.y), min = min(sum.y), max = max(sum.y))

count_sim_generator_smooth = function(x){
  sim = rbnbinom(n=n, size = 1, 
                 alpha = y_bar_star_nbinom2_smooth, 
                 beta = beta_star_nbinom2_smooth)
  return(c(median(sim), min(sim), max(sim)))
}

count_sim_generator_rough = function(x){
  sim = rbnbinom(n=n, size = 1, 
                 alpha = y_bar_star_nbinom2_rough, 
                 beta = beta_star_nbinom2_rough)
  return(c(median(sim), min(sim), max(sim)))
}

prop_sim_generator_smooth = function (X){
  p_draws_binom = rbeta(n, 
                        shape1 = alpha_star_binom2_smooth, 
                        shape2 = beta_star_binom2_smooth)
  
  n_draws = rnbinom(n = n, mu = 99, size = 2)
  
  sim = mapply(rbinom, n=1, p = p_draws_binom, size = n_draws)
  
  return(c(median(sim), min(sim), max(sim)))
}
prop_sim_generator_rough = function (X){
  p_draws_binom = rbeta(n, 
                        shape1 = alpha_star_binom2_rough, 
                        shape2 = beta_star_binom2_rough)
  
  n_draws = rnbinom(n = n, mu = 99, size = 2)
  
  sim = mapply(rbinom, n=1, p = p_draws_binom, size = n_draws)
  
  return(c(median(sim), min(sim), max(sim)))
}

count_based_smooth_sims = lapply(list, count_sim_generator_smooth) %>% 
  unlist() %>% 
  matrix(ncol = 10000, nrow = 3) %>% 
  t() %>% 
  data.frame()

count_based_rough_sims = lapply(list, count_sim_generator_rough) %>% 
  unlist() %>% 
  matrix(ncol = 10000, nrow = 3) %>% 
  t() %>% 
  data.frame()

prop_based_smooth_sims = lapply(list, prop_sim_generator_smooth)%>% 
  unlist() %>% 
  matrix(ncol = 10000, nrow = 3) %>% 
  t() %>% 
  data.frame()

prop_based_rough_sims = lapply(list, prop_sim_generator_rough)%>% 
  unlist() %>% 
  matrix(ncol = 10000, nrow = 3) %>% 
  t() %>% 
  data.frame()



# p values 

pppv_count_smooth = count_based_smooth_sims %>%
  mutate(mean = X1 - as.numeric(observed_smooth[1]), 
         min = as.numeric(observed_smooth[2]),
         max = X3-as.numeric(observed_smooth[3])) %>% 
  summarise (median = sum(mean > 0)/10000, min = sum(min>0)/10000,
             max = sum(max > 0)/10000) %>% 
  print()

pppv_count_rough = count_based_rough_sims %>%
  mutate(mean = X1 - as.numeric(observed_rough[1]), 
         min = X2 - as.numeric(observed_rough[2]),
         max = X3 - as.numeric(observed_rough[3])) %>% 
  summarise (median = sum(mean > 0)/10000, min = sum(min>0)/10000,
             max = sum(max > 0)/10000) %>% 
  print()


pppv_prop_smooth = prop_based_smooth_sims %>%
  mutate(mean = X1 - as.numeric(observed_smooth[1]), 
         min = X2 - as.numeric(observed_smooth[2]),
         max = X3 - as.numeric(observed_smooth[3])) %>% 
  summarise (median = sum(mean > 0)/10000, min = sum(min>0)/10000,
             max = sum(max > 0)/10000) %>% 
  print()

pppv_prop_rough = prop_based_rough_sims %>%
  mutate(mean = X1 - as.numeric(observed_rough[1]),
         min = X2 - as.numeric(observed_rough[2]),
         max = X3 - as.numeric(observed_rough[3])) %>% 
  summarise (median = sum(mean > 0)/10000, min = sum(min>0)/10000,
             max = sum(max > 0)/10000) %>% 
  print()


n = 6
list = vector(mode = "list", length = 1000)

observed_shedding = phase_3_counts %>% 
  filter((Participant == "A"|Participant == "C"), fibre != "black_animal_C") %>%   group_by(Participant, fibre) %>% 
  summarise(median = median(sum.y), min = min(sum.y), max = max(sum.y))

observed_shedding_vector = unlist(as.numeric(t(observed_shedding[,c(3:5)])))

prop_based_generator_shedding = function(x){
  p_draws_binom = rbeta(n, shape1 = alpha_star_binom, 
                        shape2 = beta_star_binom)
  shedding_newdata = data.frame(average = filtered_shedding_df$average)
  
  size_draws = posterior_predict(model_sim,newdata = shedding_newdata, draws = n)
  
  sim = data.frame(apply(size_draws, 2, draw_binom_by_shedding, p_draws = p_draws_binom))
  sim = sim[,c(1,3,2,4,5)]
  names(sim) = levels(phase_3_counts$fibre)[-4]
  
  sim %>% 
    summarise(median_bsa = median(blue_synthetic_A),
              min_bsa = min(blue_synthetic_A),
              max_bsa = max(blue_synthetic_A),
              median_pca = median(purple_cotton_A),
              min_pca = min(purple_cotton_A),
              max_pca = max(purple_cotton_A),
              median_rsa = median(red_synthetic_A),
              min_rsa = min(red_synthetic_A),
              max_rsa = max(red_synthetic_A),
              median_bac = median(blue_animal_C),
              min_bac = min(blue_animal_C),
              max_bac = max(blue_animal_C),
              median_rac = median(red_animal_C),
              min_rac = min(red_animal_C),
              max_rac = max(red_animal_C)) %>% 
    return()
  
}

shedding_sims = lapply(list, prop_based_generator_shedding) %>% 
  unlist() %>% 
  matrix(ncol = 1000, nrow = 15) %>% 
  t() %>% 
  data.frame()


observed_shedding_matrix = t(matrix(c(rep(observed_shedding_vector, 1000)),
                                    nrow = 15, ncol = 1000))

pppv_shedding = as.data.frame(shedding_sims - observed_shedding_matrix > 0)%>%   colSums()/1000 %>% 
  print()

n = 3
list = vector(mode = "list", length = 1000)

observed_shedding_intermediate = phase_3_counts %>% 
  filter((Participant == "A"|Participant == "C"), fibre != "black_animal_C") %>%   group_by(intermediate, Participant, fibre) %>% 
  summarise(median = median(sum.y), min = min(sum.y), max = max(sum.y))

observed_shedding_smooth_vector = unlist(as.numeric(t(observed_shedding_intermediate[c(6:10),c(4:6)])))
observed_shedding_rough_vector = unlist(as.numeric(t(observed_shedding_intermediate[c(1:5),c(4:6)])))

prop_based_generator_shedding_smooth = function(x){
  p_draws_binom = rbeta(n, shape1 = alpha_star_binom2_smooth, 
                        shape2 = beta_star_binom2_smooth)
  shedding_newdata = data.frame(average = filtered_shedding_df$average)
  
  size_draws = posterior_predict(model_sim,newdata = shedding_newdata, draws = n)
  
  sim = data.frame(apply(size_draws, 2, draw_binom_by_shedding, p_draws = p_draws_binom))
  sim = sim[,c(1,3,2,4,5)]
  names(sim) = levels(phase_3_counts$fibre)[-4]
  
  sim %>% 
    summarise(median_bsa = median(blue_synthetic_A),
              min_bsa = min(blue_synthetic_A),
              max_bsa = max(blue_synthetic_A),
              median_pca = median(purple_cotton_A),
              min_pca = min(purple_cotton_A),
              max_pca = max(purple_cotton_A),
              median_rsa = median(red_synthetic_A),
              min_rsa = min(red_synthetic_A),
              max_rsa = max(red_synthetic_A),
              median_bac = median(blue_animal_C),
              min_bac = min(blue_animal_C),
              max_bac = max(blue_animal_C),
              median_rac = median(red_animal_C),
              min_rac = min(red_animal_C),
              max_rac = max(red_animal_C)) %>% 
    return()
  
}

prop_based_generator_shedding_rough = function(x){
  p_draws_binom = rbeta(n, shape1 = alpha_star_binom2_rough, 
                        shape2 = beta_star_binom2_rough)
  shedding_newdata = data.frame(average = filtered_shedding_df$average)
  
  size_draws = posterior_predict(model_sim,newdata = shedding_newdata, draws = n)
  
  sim = data.frame(apply(size_draws, 2, draw_binom_by_shedding, p_draws = p_draws_binom))
  sim = sim[,c(1,3,2,4,5)]
  names(sim) = levels(phase_3_counts$fibre)[-4]
  
  sim %>% 
    summarise(median_bsa = median(blue_synthetic_A),
              min_bsa = min(blue_synthetic_A),
              max_bsa = max(blue_synthetic_A),
              median_pca = median(purple_cotton_A),
              min_pca = min(purple_cotton_A),
              max_pca = max(purple_cotton_A),
              median_rsa = median(red_synthetic_A),
              min_rsa = min(red_synthetic_A),
              max_rsa = max(red_synthetic_A),
              median_bac = median(blue_animal_C),
              min_bac = min(blue_animal_C),
              max_bac = max(blue_animal_C),
              median_rac = median(red_animal_C),
              min_rac = min(red_animal_C),
              max_rac = max(red_animal_C)) %>% 
    return()
  
}

shedding_sims_smooth = lapply(list, prop_based_generator_shedding_smooth) %>% 
  unlist() %>% 
  matrix(ncol = 1000, nrow = 15) %>% 
  t() %>% 
  data.frame()
shedding_sims_rough = lapply(list, prop_based_generator_shedding_rough) %>% 
  unlist() %>% 
  matrix(ncol = 1000, nrow = 15) %>% 
  t() %>% 
  data.frame()


observed_shedding_matrix_smooth = t(matrix(c(rep(observed_shedding_smooth_vector, 1000)),
                                           nrow = 15, ncol = 1000))
observed_shedding_matrix_rough = t(matrix(c(rep(observed_shedding_rough_vector, 1000)),
                                          nrow = 15, ncol = 1000))

pppv_shedding_smooth = as.data.frame(shedding_sims_smooth - observed_shedding_matrix_smooth > 0)%>%  
  colSums()/1000 %>% 
  print()

pppv_shedding_rough = as.data.frame(shedding_sims_rough - observed_shedding_matrix_rough > 0)%>%  
  colSums()/1000 %>% 
  print()
