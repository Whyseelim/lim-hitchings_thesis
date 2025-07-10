pkgs = c("ggplot2", "tidyr", "dplyr", "stringr" ,"rstanarm", "rstan", "bayesrules", "tidybayes", "bayesplot", "gridExtra")
lapply(pkgs, library, character.only = TRUE)

# function to separate phase, participant

seperate_participant = function(x){
  phase_1 = x[which(startsWith(x,"01") == TRUE)]
  phase_2 = x[which(startsWith(x,"01") == FALSE)]
  participant_1 = str_extract(phase_1, "(?<=-).*|(?<=_).*")
  participant_2 = str_extract(phase_2, "(?<=-).*(?=-)|(?<=_).*(?=_)")
  return(c(participant_1,participant_2))
  
}

seperate_replicate = function(x){
  phase_1 = x[which(startsWith(x,"01") == TRUE)]
  phase_2 = x[which(startsWith(x,"01") == FALSE)]
  rep_1 = str_extract(phase_1, ".*(?=-)|.*(?=_)")
  rep_2 = str_extract(phase_2, "(?<=-..-).*|(?<=_.._).*")
  return(c(rep_1,rep_2))
  
}

# import
df_raw = read.csv("../data/phase_1_2_data.csv")

# correctly labeling replicate 2
df_raw$item = str_replace(df_raw$item, "(?<=02_..)_01", "_02")

df_clean = 
  df_raw %>% 
  # delete all NA
  na.omit() %>% 
  # drop phase number
  mutate(participant = seperate_participant(item), replicate = seperate_replicate(item)) %>% 
  # get rid of front_back and item column since useless
  select(-c(front_back,item))

# Convert to alphabet

convert_to_alphabet = function(x){
  last_2 =  substr(x, nchar(x)-1, nchar(x))
  replace = case_when(last_2 == "01" ~ "A",
                      last_2 == "04" ~ "B",
                      last_2 == "05" ~ "C",
                      last_2 == "14" ~ "D",
                      last_2 == "15" ~ "E",
                      .default = last_2)
  new = str_replace(x, "(?<=_)\\d\\d", replace)
  return(new)
}


# long table
df_clean_long = pivot_longer(df_clean, !c(left_right, top_bottom, grid, participant, replicate), values_to = "count", names_to = "fibre") %>% 
  mutate(participant = case_when(participant == "01" ~ "A",
                                 participant == "04" ~ "B",
                                 participant == "05" ~ "C",
                                 participant == "14" ~ "D",
                                 participant == "15" ~ "E",
                                 .default = "12"),
         fibre = convert_to_alphabet(fibre))



# counts of each fibre type on each participant/replicate
df_counts = df_clean_long %>%
  group_by(participant, replicate, fibre) %>% 
  summarize (sum = sum(count)) 

#as factors

df_counts$fibre = factor(df_counts$fibre, levels = unique(df_clean_long$fibre))
df_counts$participant = factor(df_counts$participant)
df_counts$replicate = factor(df_counts$replicate)

# select colors
unique_fibres = levels(df_counts$fibre)

colours = ifelse(grepl("blue", unique_fibres),"#3467e0",
                 ifelse(grepl("purple",unique_fibres),"#613485",
                        ifelse(grepl("red",unique_fibres), "#c23030",
                               ifelse(grepl("white",unique_fibres),
                                      "grey",ifelse(grepl("black",unique_fibres), "#1f1f1f",ifelse(grepl("orange", unique_fibres), "#f59842",ifelse(grepl("yellow", unique_fibres), "#f5dd42",ifelse(grepl("green", unique_fibres),"#6d996b","#30c2bd"))))))))


colours_rmv = colours[-c(10,11,12,13)]

plot_bar_rep1 = df_counts %>% 
  filter(replicate == "01", participant !=12, !grepl("_12", fibre)) %>% 
  ggplot(mapping = aes(x = fibre, y = log10(sum), fill = fibre))+
  geom_bar(stat = "identity", show.legend = FALSE)+
  scale_fill_manual(values = colours_rmv)+
  facet_grid(rows = vars(participant))+
  geom_vline(xintercept = c(3.5,6.5,9.5,13.5), size = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.spacing.y = unit(5,"mm"),
        text = element_text(size = 40), plot.margin = margin(l = 65))+
  ylab("log(count)") +
  xlab("Target fibre")

plot_bar_rep1_removed = df_counts %>% 
  filter(replicate == "01", participant !=12, sum < 500, !grepl("_12", fibre)) %>% 
  ggplot(mapping = aes(x = fibre, y = sum, fill = fibre))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colours_rmv)+
  facet_grid(rows = vars(participant))+
  geom_vline(xintercept = c(3.5,6.5,9.5,13.5,17.5))+
  theme(axis.text.x = element_text(angle = 90))

plot_native_vs_non = df_counts %>% 
  mutate(native = str_sub(fibre, -1,-1) == participant) %>% 
  filter(participant !=12,  !grepl("_12", fibre)) %>% 
  ggplot(mapping = aes(x = native, y = log10(sum)))+
  geom_boxplot(fill ="skyblue", size = 2)+
  theme(text = element_text (size = 40))+
  scale_x_discrete(name = "Target fibre", labels = c("Non-Native", "Native"))+
  scale_y_continuous(name = "log(count)")

#vector of participants remove 4
participants = unique(df_counts$participant)[-1]

# function to separate each participant and their native fibres and plot
plot_rep = function(value){
  # set colour
  colour_participants = colours[which(str_detect(levels(df_counts$fibre), paste0("_", value, "$")))]
  #plot
  df_counts %>%
    filter(participant == value, str_detect(fibre, paste0("_", value, "$"))) %>% 
    ggplot(df_counts, mapping = aes(x = replicate, y = log10(sum), fill = fibre))+
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE)+
    scale_y_continuous(name = "log(count)", limits = c(0,3))+
    xlab(paste0("Replicates (", value, ")"))+
    scale_fill_manual(values = colour_participants)+
    theme(text = element_text(size = 20))
  
}

# get list of plots, 1 for each participant 
plots_allrep = lapply(X = as.character(participants), FUN=plot_rep)

grid.arrange(plots_allrep[[1]], plots_allrep[[2]], plots_allrep[[3]], plots_allrep[[4]], plots_allrep[[5]])


# classify garment
shirt_type = function(x){
  if(x == "01" | x == "02"){
    return("own")
  }else{
    if(x == "03" | x == "04"){
      return("smooth")
    }
    return("rough")
  }
}

# function to figure out differences and prop change from first collection
count_diff = function(value){
  # create vector of first rep results that lines up with sum column
  first_rep = rep(df_counts$sum[which(df_counts$replicate =="01" & df_counts$participant==value & str_detect(df_counts$fibre, paste0("_", value, "$"))) ],6)
  x = df_counts %>%
    ungroup() %>% 
    filter(participant == value, str_detect(fibre, paste0("_", value, "$"))) %>% 
    mutate(difference = sum - first_rep, type = sapply(replicate, shirt_type)) %>% 
    mutate(perc_change = (difference/first_rep) *100)
  return(x)
}

df_diff = bind_rows(lapply(participants, count_diff))

plot_diff_from_rep1 = df_diff %>% 
  filter(replicate !="01") %>% 
  ggplot(aes(x = type, y = difference))+
  geom_boxplot(fill = "skyblue")


plot_sum_group_participants_fibres = ggplot(df_diff, aes(y = log10(sum), x = participant, fill = fibre))+
  geom_boxplot(, size = 2, outlier.size = 4)+
  geom_point(size = -1, aes(fill = fibre))+
  theme(text = element_text (size = 40), legend.text = element_text(size = 20), legend.title = element_blank())+
  xlab("Participant")+
  scale_y_continuous(name = "log(count)", limits = c(0,3.05))+
  scale_fill_manual(values = colours_rmv)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), size = 1)

plot_dist_part = df_diff %>% 
  filter(replicate != "01") %>% 
  ggplot(aes(x = sum, fill = participant ))+
  geom_density(alpha = 0.5)

plot_dist_fibre = df_diff %>% 
  filter(replicate != "01") %>% 
  ggplot(aes(x = sum, fill = fibre))+
  scale_fill_manual(values  = colours_rmv)+
  geom_density(alpha = 0.5)+
  facet_grid(rows = vars(participant))

plot_hist_total = df_diff %>% 
  filter(fibre != "white_cat_D", replicate != "01") %>% 
  ggplot(aes(x = sum))+
  geom_histogram(fill = "skyblue", bins = 18)

mean = round(mean(df_diff_outlier_rmv$sum),0)
size = 2
x = c(0:500)
dens = dnbinom(x, mu = mean, size = size)

model = data.frame(x = x, dens = dens)

plot_dens_total_vs_model = df_diff %>% 
  filter(fibre != "white_cat_D", replicate != "01") %>% 
  ggplot(aes(x = sum))+
  geom_density(fill = "skyblue", alpha = 0.5, size = 2)+
  geom_line(data = model, aes(x = x, y = dens), linetype = "dashed", size = 2)+
  annotate("text", x = 400, y = 0.001, size = 10, label = paste0("~NegBin(", mean, ",", size,")"))+
  xlab("Count")+
  ylab("Density")+
  theme(text = element_text (size = 30))

shedding_df = read.csv("../data/shedding.csv") %>% 
  rowwise() %>% 
  mutate(average = mean(c(shedding.1, shedding.2, shedding.3, shedding.4, shedding.5)), target = convert_to_alphabet(target)) 

save(shedding_df, file = "../scripts/shedding_df.Rdata")

# add to df

df_counts_shedding = df_counts %>% 
  mutate(native = str_sub(fibre, -1,-1) == participant) %>% 
  filter(replicate!="01", native == TRUE , fibre %in% shedding_df$target) %>% 
  inner_join(y = shedding_df, by = join_by(fibre == target))

plot_shedding = df_counts_shedding %>%
  # remove outlier of high shedder
  #filter(average < 100) %>% 
  ggplot(aes(x = average, y = log10(sum)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x)

plot_shedding_use = df_counts_shedding %>%
  # remove outlier of high shedder
  filter(average < 500) %>% 
  ggplot(mapping = aes(x = average, y = log10(sum)))+
  geom_point(mapping = aes(colour = use), size = 3)+
  geom_smooth(method = "lm", formula = y~x, mapping = aes(x = average, y = log10(sum), colour = use), se = FALSE, linetype = "dashed", size = 3)+
  geom_smooth(method = "lm", formula = y~x, colour = "black", size = 3, alpha = 0.6)+
  scale_x_continuous(name = "Shedding", limits = c(8,90))+
  ylab("log(count)")+
  labs(colour = "Use")+
  theme(text = element_text (size = 30), legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 25), axis.title.x = element_blank())+
  scale_color_manual(values = c("#009E73", "#D55E00"))

plot_shedding_store = df_counts_shedding %>%
  # remove outlier of high shedder and not storage together as we only have 2
  filter(average < 500, storage != "together") %>% 
  ggplot(aes(x = average, y = log10(sum), colour = storage))+
  geom_point(mapping = aes(colour = storage), size = 3)+
  geom_smooth(method = "lm", formula = y~x, mapping = aes(x = average, y = log10(sum), colour = storage), se = FALSE, linetype = "dashed", size = 3)+
  geom_smooth(method = "lm", formula = y~x, colour = "black",size = 3, alpha = 0.6)+
  scale_x_continuous(name = "Shedding", limits = c(8,90))+
  ylab("log(count)")+
  labs(colour = "Storage")+
  theme(text = element_text (size = 30), legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 25))+
  scale_color_manual(values = c("#56B4E9", "#CC79A7"))

grid.arrange(plot_shedding_use, plot_shedding_store)

plot_shedding_use_store = df_counts_shedding %>%
  # remove outlier of high shedder
  filter(average < 100, storage != "together") %>% 
  mutate(use_storage = paste0(use, "_", storage)) %>% 
  ggplot(aes(x = average, y = log10(sum), colour = use_storage))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x)

plot_shedding_use_store 

seed = 484468

df_counts_shedding_notogether = df_counts_shedding %>% 
  filter(storage != "together")

model_sim <- stan_glm(
  sum ~ average, 
  data = df_counts_shedding_notogether, family = neg_binomial_2,
  prior_intercept = normal(1.5, 2, autoscale = FALSE),
  prior = normal(0, 1, autoscale = FALSE), 
  prior_aux = exponential(1, autoscale = FALSE),
  chains = 4, iter = 5000*2, seed = seed)

ppcheck1= pp_check(model_sim)+
  theme(text = element_text(size = 30), legend.text = element_text(size = 20))
ppcheck1_graph = 
  ppcheck1+
  annotate("text", x = 700, y = 0.01, label = "Model 1", size = 8)+
  xlim(0,800)+
  ylim (0,0.014)

model_sim_2 <- stan_glm(
  sum ~ average + use + storage, 
  data = df_counts_shedding_notogether, family = neg_binomial_2,
  prior_intercept = normal(10, 2, autoscale = TRUE),
  prior = normal(0, 1, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = seed)

ppcheck2 = pp_check(model_sim_2)+
  theme(text = element_text(size = 30), legend.text = element_text(size = 20))

ppcheck2_graph =
  ppcheck2+
  annotate("text", x = 700, y = 0.01, label = "Model 2", size = 8)+
  xlim(0,800)+
  ylim (0,0.014)

model_sim_3 <- stan_glm(
  sum ~ average + use, 
  data = df_counts_shedding_notogether, family = neg_binomial_2,
  prior_intercept = normal(10, 2, autoscale = TRUE),
  prior = normal(0, 1, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = seed)

ppcheck3 = pp_check(model_sim_3)+
  theme(text = element_text(size = 30), legend.text = element_text(size = 20))+
  xlab("Count")

ppcheck3_graph = 
  ppcheck3+
  annotate("text", x = 700, y = 0.01, label = "Model 3", size = 8)+
  xlim(0,800)+
  ylim (0,0.014)

model_sim_4 <- stan_glm(
  sum ~ average + storage, 
  data = df_counts_shedding_notogether, family = neg_binomial_2,
  prior_intercept = normal(10, 2, autoscale = TRUE),
  prior = normal(0, 1, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = seed)

ppcheck4=pp_check(model_sim_4)+
  theme(text = element_text(size = 30), legend.text = element_text(size = 20))+
  xlab("Count")

ppcheck4_graph = 
  ppcheck4 +
  annotate("text", x = 700, y = 0.01, label = "Model 4", size = 8)+
  xlim(0,800)+
  ylim (0,0.014)

grid.arrange(ppcheck1_graph, ppcheck2_graph, ppcheck3_graph, ppcheck4_graph)

predictions <- posterior_predict(model_sim_4, newdata = df_counts_shedding_notogether)
dim(predictions)

ppc_intervals(df_counts_shedding_notogether$sum, yrep = predictions, x = df_counts_shedding_notogether$average, 
              prob = 0.5, prob_outer = 0.95)

model_accuracy = lapply(list(model_sim, model_sim_2, model_sim_3, model_sim_4),prediction_summary, data = df_counts_shedding_notogether[,c(4,12,14,15)], prob_inner = 0.68, prob_outer = 0.95, stable = TRUE)

model_accuracy_table = signif(do.call(rbind.data.frame, model_accuracy),2)
model_accuracy_table

waic = lapply(lapply(list(model_sim, model_sim_2, model_sim_3, model_sim_4),FUN = log_lik ), FUN = waic)

print(waic)