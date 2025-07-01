load("../inter-framework-functions.RData")

# Background class

# FFG data

vector_f = c(0.0158, 0.0317, 0.0635, 0.127, 0.254, 0.508) # priors

vector_f_resample = c(rep("0", 1), 
                      rep("1-10", 2), 
                      rep("11-20", 2), 
                      rep("21-30", 4), 
                      rep("31-40", 8), 
                      rep(">40", 16))

resample_vector = rep(list(vector_f_resample), 1000)

vector_f_resampled = sapply(resample_vector, resample_f)

# Size resample
vector_s = c(0.76, 0.20, 0.04)

vector_s_resample = c(rep("small", 38), rep("medium", 10), rep("large", 2))

vector_s_resampled = sapply(resample_vector, resample_s)


# Occurrence assignment
gamma_common = 0.14
gamma_rare = 0.00171

# Transfer class

# occurrence assignment
gamma_prime_common = 0.00566
gamma_prime_rare = 0.000855

# Primary transfer assignment
vector_t_primary = c (0.001, 0.009, 0.05, 0.94)

# Secondary transfer

vector_t = c(0.43, 0.4, 0.13, 0.04)

vector_t_resample = c(rep("none", 13), 
                      rep("small", 12), 
                      rep("medium", 4), 
                      rep("large", 1))

resample_vector = rep(list(vector_t_resample), 1000)

vector_t_resampled = sapply(resample_vector, resample_t)

# Condition 1
n_condition_1 = 4

c_condition_1_4R = c(rep(FALSE, 4))
gamma_condition_1_4R = assign_gamma(c_condition_1_4R)
gamma_prime_condition_1_4R = assign_gamma(c_condition_1_4R, 
                                          common = common, 
                                          rare = rare)

c_condition_1_4C = c(rep(TRUE, 4))
gamma_condition_1_4C = assign_gamma(c_condition_1_4C)
gamma_prime_condition_1_4C = assign_gamma(c_condition_1_4C, 
                                          common = common, 
                                          rare = rare)

c_condition_1_1R3C = c(FALSE, rep(TRUE, 3))
gamma_condition_1_1R3C = assign_gamma(c_condition_1_1R3C)
gamma_prime_condition_1_1R3C = assign_gamma(c_condition_1_1R3C, 
                                            common = common, 
                                            rare = rare)

c_condition_1_2R2C = c(rep(FALSE, 2), rep(TRUE, 2))
gamma_condition_1_2R2C = assign_gamma(c_condition_1_2R2C)
gamma_prime_condition_1_2R2C = assign_gamma(c_condition_1_2R2C, 
                                            common = common, 
                                            rare = rare)

c_condition_1_3R1C = c(rep(FALSE, 3), TRUE, 2)
gamma_condition_1_3R1C = assign_gamma(c_condition_1_3R1C)
gamma_prime_condition_1_3R1C = assign_gamma(c_condition_1_3R1C, 
                                            common = common, 
                                            rare = rare)


vector_t_list = rep(list(vector_t), 4)
cpts_t_condition_1 = lapply(vector_t_list, write_cpt_t)



nocontact = matrix(rep(c(1,0,0,0),1000),nrow = 1000, ncol = 4, byrow = TRUE)
vector_ts = cbind(resample_2t, nocontact)
cpts_t_condition_1 = apply(vector_ts,1, into_t_list)

## create oobns
bn_condition_1_4R_resamples = mapply(build_bn, cpts_t = cpts_t_condition_1, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_1, gamma = gamma_condition_1_4R, gamma_prime = gamma_prime_condition_1_4R), SIMPLIFY = FALSE)
bn_condition_1_4C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_1, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_1, gamma = gamma_condition_1_4C, gamma_prime = gamma_prime_condition_1_4C), SIMPLIFY = FALSE)
bn_condition_1_1R3C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_1, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_1, gamma = gamma_condition_1_1R3C, gamma_prime = gamma_prime_condition_1_1R3C), SIMPLIFY = FALSE)
bn_condition_1_2R2C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_1, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_1, gamma = gamma_condition_1_2R2C, gamma_prime = gamma_prime_condition_1_2R2C), SIMPLIFY = FALSE)
bn_condition_1_3R1C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_1, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_1, gamma = gamma_condition_1_3R1C, gamma_prime = gamma_prime_condition_1_3R1C), SIMPLIFY = FALSE)

## Evaluate LR medium sized
evidence_1 = list(E1 = c(0, 0, 1, 0))
evidence_2 = list(E1 = c(0, 0, 1, 0), 
                  E2 = c(0, 0, 1, 0))
evidence_3 = list(E1 = c(0, 0, 1, 0), 
                  E2 = c(0, 0, 1, 0), 
                  E3 = c(0, 0, 1, 0))
evidence_4 = list(E1 = c(0, 0, 1, 0), 
                  E2 = c(0, 0, 1, 0), 
                  E3 = c(0, 0, 1, 0), 
                  E4 = c(0, 0, 1, 0))


LR_condition_1_1R_resamples = sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_1), get_LR)

LR_condition_1_1C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_1), get_LR)

LR_condition_1_2R_resamples = sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_1_2C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_1_1R1C_resamples = 
  sapply(lapply(bn_condition_1_1R3C_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_1_2R1C_resamples = 
  sapply(lapply(bn_condition_1_2R2C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_1_3R_resamples = 
  sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_3), get_LR) 

LR_condition_1_3C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_1_1R2C_resamples = 
  sapply(lapply(bn_condition_1_1R3C_resamples, query_node, evidence = evidence_3), get_LR)


LR_condition_1_4R_resamples = 
  sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_4C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_1R3C_resamples = 
  sapply(lapply(bn_condition_1_1R3C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_2R2C_resamples = 
  sapply(lapply(bn_condition_1_2R2C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_3R1C_resamples = 
  sapply(lapply(bn_condition_1_3R1C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_medium_resamples = 
  signif(c(LR_condition_1_1C_resamples, 
           LR_condition_1_2C_resamples, 
           LR_condition_1_3C_resamples, 
           LR_condition_1_4C_resamples, 
           LR_condition_1_1R_resamples, 
           LR_condition_1_2R_resamples, 
           LR_condition_1_3R_resamples, 
           LR_condition_1_4R_resamples,
           LR_condition_1_1R1C_resamples, 
           LR_condition_1_1R2C_resamples,
           LR_condition_1_1R3C_resamples,
           LR_condition_1_2R1C_resamples, 
           LR_condition_1_2R2C_resamples, 
           LR_condition_1_3R1C_resamples), 2)


## Evaluate LR small sized
evidence_1 = list(E1 = c(0, 1, 0, 0))
evidence_2 = list(E1 = c(0, 1, 0, 0), 
                  E2 = c(0, 1, 0, 0))
evidence_3 = list(E1 = c(0, 1, 0, 0), 
                  E2 = c(0, 1, 0, 0), 
                  E3 = c(0, 1, 0, 0))
evidence_4 = list(E1 = c(0, 1, 0, 0), 
                  E2 = c(0, 1, 0, 0), 
                  E3 = c(0, 1, 0, 0), 
                  E4 = c(0, 1, 0, 0))

LR_condition_1_1R_resamples = sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_1), get_LR)

LR_condition_1_1C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_1), get_LR)

LR_condition_1_2R_resamples = sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_1_2C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_1_1R1C_resamples = 
  sapply(lapply(bn_condition_1_1R3C_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_1_2R1C_resamples = 
  sapply(lapply(bn_condition_1_2R2C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_1_3R_resamples = 
  sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_3), get_LR) 

LR_condition_1_3C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_1_1R2C_resamples = 
  sapply(lapply(bn_condition_1_1R3C_resamples, query_node, evidence = evidence_3), get_LR)


LR_condition_1_4R_resamples = 
  sapply(lapply(bn_condition_1_4R_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_4C_resamples = 
  sapply(lapply(bn_condition_1_4C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_1R3C_resamples = 
  sapply(lapply(bn_condition_1_1R3C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_2R2C_resamples = 
  sapply(lapply(bn_condition_1_2R2C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_3R1C_resamples = 
  sapply(lapply(bn_condition_1_3R1C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_1_small_resamples = 
  signif(c(LR_condition_1_1C_resamples, 
           LR_condition_1_2C_resamples, 
           LR_condition_1_3C_resamples, 
           LR_condition_1_4C_resamples, 
           LR_condition_1_1R_resamples, 
           LR_condition_1_2R_resamples, 
           LR_condition_1_3R_resamples, 
           LR_condition_1_4R_resamples, 
           LR_condition_1_1R1C_resamples, 
           LR_condition_1_1R2C_resamples,
           LR_condition_1_1R3C_resamples,
           LR_condition_1_2R1C_resamples, 
           LR_condition_1_2R2C_resamples, 
           LR_condition_1_3R1C_resamples), 2)


# plots
names = c("1C", "2C", "3C", "4C", "1R", "2R", "3R", "4R", "1R1C", "1R2C", "1R3C", "2R1C", "2R2C", "3R1C")

names_vector = c(sapply(names, rep, times = 1000))

LR_condition_1_medium_resamples_df = cbind(signif(log10(LR_condition_1_medium_resamples[,-1]),3), names_vector, rep("medium", 14000))
LR_condition_1_small_resamples_df = cbind(signif(log10(LR_condition_1_small_resamples[,-1]),3), names_vector, rep("small", 14000))
LR_condition_1_df = as.data.frame(rbind(LR_condition_1_medium_resamples_df, LR_condition_1_small_resamples_df))
names(LR_condition_1_df) = c("log10_LR", "Findings", "Size_recovered")
LR_condition_1_df$Findings = factor(LR_condition_1_df$Findings, levels = names)
LR_condition_1_df$Size_recovered = factor(LR_condition_1_df$Size_recovered, levels = c("small", "medium"))
LR_condition_1_df$log10_LR = as.numeric(LR_condition_1_df$log10_LR)

ggplot(LR_condition_1_df, aes(x = Findings, y = log10_LR, fill = Size_recovered))+
  geom_boxplot(alpha=0.8)+
  ylab("Log(LR)")+
  geom_hline(yintercept = 0, color = "darkgreen")+
  geom_hline(yintercept = 1, color = "darkgreen")+
  annotate("text", x = 17, y= 0.8, label = "Weak", size = 10)+
  geom_hline(yintercept = 2, color = "darkgreen")+
  annotate("text", x = 17, y= 1.8, label = "Moderate", size = 10)+
  geom_hline(yintercept = 3, color = "darkgreen")+
  annotate("text", x = 17, y= 2.8, label = "Strong", size = 10)+
  geom_hline(yintercept = 4, color = "darkgreen")+
  annotate("text", x = 17, y= 3.8, label = "Very strong", size = 10)+
  annotate("text", x = 17.1, y= 4.8, label = "Extremely strong", size = 10)+
  annotate("text", x = 20, y= 4.8, label = "")+
  scale_fill_manual(values = c("#56B4E9", "#D55E00"))+
  theme(axis.title = element_text(size = 30), axis.text.x = element_text(size = 24, angle = 45), axis.text.y = element_text(size = 24), legend.title = element_text(size = 26), legend.text = element_text( size = 24))


ggsave("../figures/plot_condition_1_resamples.png", height =8, width = 15)


# Condition 2

n_condition_2 = 5

c_condition_2_4R = c(TRUE, rep(FALSE, 4))
gamma_condition_2_4R = assign_gamma(c_condition_2_4R)
gamma_prime_condition_2_4R = assign_gamma(c_condition_2_4R,
                                          common = common, 
                                          rare = rare)

c_condition_2_4C = c(TRUE, rep(TRUE, 4))
gamma_condition_2_4C = assign_gamma(c_condition_2_4C)
gamma_prime_condition_2_4C = assign_gamma(c_condition_2_4C, 
                                          common = common, 
                                          rare = rare)

c_condition_2_1R3C = c(TRUE, FALSE, rep(TRUE, 3))
gamma_condition_2_1R3C = assign_gamma(c_condition_2_1R3C)
gamma_prime_condition_2_1R3C = assign_gamma(c_condition_2_1R3C, 
                                            common = common, 
                                            rare = rare)

c_condition_2_2R2C = c(TRUE, rep(FALSE, 2), rep(TRUE, 2))
gamma_condition_2_2R2C = assign_gamma(c_condition_2_2R2C)
gamma_prime_condition_2_2R2C = assign_gamma(c_condition_2_2R2C, 
                                            common = common, 
                                            rare = rare)

c_condition_2_3R1C = c(TRUE, rep(FALSE, 3), TRUE, 2)
gamma_condition_2_3R1C = assign_gamma(c_condition_2_3R1C)
gamma_prime_condition_2_3R1C = assign_gamma(c_condition_2_3R1C, 
                                            common = common, 
                                            rare = rare)



## set transfer cpt
vector_t_primary = c (0.001, 0.009, 0.05, 0.94)
cpts_t_condition_2 = lapply(cpts_t_condition_1, append, values = list(c(vector_t_primary, 1,0,0,0)), after = 0)

## create oobns
bn_condition_2_4R_resamples = mapply(build_bn, cpts_t = cpts_t_condition_2, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_2, gamma = gamma_condition_2_4R, gamma_prime = gamma_prime_condition_2_4R), SIMPLIFY = FALSE)
bn_condition_2_4C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_2, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_2, gamma = gamma_condition_2_4C, gamma_prime = gamma_prime_condition_2_4C), SIMPLIFY = FALSE)
bn_condition_2_1R3C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_2, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_2, gamma = gamma_condition_2_1R3C, gamma_prime = gamma_prime_condition_2_1R3C), SIMPLIFY = FALSE)
bn_condition_2_2R2C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_2, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_2, gamma = gamma_condition_2_2R2C, gamma_prime = gamma_prime_condition_2_2R2C), SIMPLIFY = FALSE)
bn_condition_2_3R1C_resamples = mapply(build_bn, cpts_t = cpts_t_condition_2, vector_s = as.list(data.frame(t(resample_size))), vector_f = as.list(data.frame(t(resample_ffg))), MoreArgs = list(n = n_condition_2, gamma = gamma_condition_2_3R1C, gamma_prime = gamma_prime_condition_2_3R1C), SIMPLIFY = FALSE)

## Evaluate LR medium sized
evidence_1 = list(E1 = c(0, 0, 0, 1))
evidence_2 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 0, 1, 0))
evidence_3 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 0, 1, 0), 
                  E3 = c(0, 0, 1, 0))
evidence_4 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 0, 1, 0), 
                  E3 = c(0, 0, 1, 0), 
                  E4 = c(0, 0, 1, 0))
evidence_5 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 0, 1, 0), 
                  E3 = c(0, 0, 1, 0), 
                  E4 = c(0, 0, 1, 0), 
                  E5 = c(0, 0, 1, 0))


LR_condition_2_base_resamples = sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_1), get_LR)

LR_condition_2_1R_resamples = sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_2_1C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_2_2R_resamples = sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_2_2C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_2_1R1C_resamples = 
  sapply(lapply(bn_condition_2_1R3C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_2_2R1C_resamples = 
  sapply(lapply(bn_condition_2_2R2C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_2_3R_resamples = 
  sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_4), get_LR) 

LR_condition_2_3C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_2_1R2C_resamples = 
  sapply(lapply(bn_condition_2_1R3C_resamples, query_node, evidence = evidence_4), get_LR)


LR_condition_2_4R_resamples = 
  sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_4C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_1R3C_resamples = 
  sapply(lapply(bn_condition_2_1R3C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_2R2C_resamples = 
  sapply(lapply(bn_condition_2_2R2C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_3R1C_resamples = 
  sapply(lapply(bn_condition_2_3R1C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_medium_resamples = 
  signif(c(LR_condition_2_1C_resamples, 
           LR_condition_2_2C_resamples, 
           LR_condition_2_3C_resamples, 
           LR_condition_2_4C_resamples, 
           LR_condition_2_1R_resamples, 
           LR_condition_2_2R_resamples, 
           LR_condition_2_3R_resamples, 
           LR_condition_2_4R_resamples,
           LR_condition_2_1R1C_resamples, 
           LR_condition_2_1R2C_resamples,
           LR_condition_2_1R3C_resamples,
           LR_condition_2_2R1C_resamples, 
           LR_condition_2_2R2C_resamples, 
           LR_condition_2_3R1C_resamples), 2)


## Evaluate LR small sized
evidence_1 = list(E1 = c(0, 1, 0, 0))
evidence_2 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 1, 0, 0))
evidence_3 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 1, 0, 0), 
                  E3 = c(0, 1, 0, 0))
evidence_4 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 1, 0, 0), 
                  E3 = c(0, 1, 0, 0), 
                  E4 = c(0, 1, 0, 0))
evidence_5 = list(E1 = c(0, 0, 0, 1), 
                  E2 = c(0, 1, 0, 0), 
                  E3 = c(0, 1, 0, 0), 
                  E4 = c(0, 1, 0, 0), 
                  E5 = c(0, 1, 0, 0))

LR_condition_2_1R_resamples = sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_2_1C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_2), get_LR)

LR_condition_2_2R_resamples = sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_2_2C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_2_1R1C_resamples = 
  sapply(lapply(bn_condition_2_1R3C_resamples, query_node, evidence = evidence_3), get_LR)

LR_condition_2_2R1C_resamples = 
  sapply(lapply(bn_condition_2_2R2C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_2_3R_resamples = 
  sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_4), get_LR) 

LR_condition_2_3C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_4), get_LR)

LR_condition_2_1R2C_resamples = 
  sapply(lapply(bn_condition_2_1R3C_resamples, query_node, evidence = evidence_4), get_LR)


LR_condition_2_4R_resamples = 
  sapply(lapply(bn_condition_2_4R_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_4C_resamples = 
  sapply(lapply(bn_condition_2_4C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_1R3C_resamples = 
  sapply(lapply(bn_condition_2_1R3C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_2R2C_resamples = 
  sapply(lapply(bn_condition_2_2R2C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_3R1C_resamples = 
  sapply(lapply(bn_condition_2_3R1C_resamples, query_node, evidence = evidence_5), get_LR)

LR_condition_2_small_resamples = 
  signif(c(LR_condition_2_1C_resamples, 
           LR_condition_2_2C_resamples, 
           LR_condition_2_3C_resamples, 
           LR_condition_2_4C_resamples, 
           LR_condition_2_1R_resamples, 
           LR_condition_2_2R_resamples, 
           LR_condition_2_3R_resamples, 
           LR_condition_2_4R_resamples, 
           LR_condition_2_1R1C_resamples, 
           LR_condition_2_1R2C_resamples,
           LR_condition_2_1R3C_resamples,
           LR_condition_2_2R1C_resamples, 
           LR_condition_2_2R2C_resamples, 
           LR_condition_2_3R1C_resamples), 2)



# plots
names = c("1C", "2C", "3C", "4C", "1R", "2R", "3R", "4R", "1R1C", "1R2C", "1R3C", "2R1C", "2R2C", "3R1C")

names_vector = c(sapply(names, rep, times = 1000))

LR_condition_2_base_resamples_df = cbind(signif(log10(LR_condition_2_base_resamples[,-1]),3), rep("base", 1000), rep("base", 1000))
LR_condition_2_medium_resamples_df = cbind(signif(log10(LR_condition_2_medium_resamples[,-1]),3), names_vector, rep("medium", 14000))
LR_condition_2_small_resamples_df = cbind(signif(log10(LR_condition_2_small_resamples[,-1]),3), names_vector, rep("small", 14000))
LR_condition_2_df = as.data.frame(rbind(LR_condition_2_base_resamples_df, LR_condition_2_medium_resamples_df, LR_condition_2_small_resamples_df))
names(LR_condition_2_df) = c("log10_LR", "Findings", "Size_recovered")
LR_condition_2_df$Findings = factor(LR_condition_2_df$Findings, levels = c("base",names))
LR_condition_2_df$Size_recovered = factor(LR_condition_2_df$Size_recovered, levels = c("small", "base", "medium"))
LR_condition_2_df$log10_LR = as.numeric(LR_condition_2_df$log10_LR)

ggplot(LR_condition_2_df, aes(x = Findings, y = log10_LR, fill = Size_recovered))+
  geom_boxplot(alpha=0.8)+
  ylab("Log(LR)")+
  geom_hline(yintercept = 0, color = "darkgreen")+
  geom_hline(yintercept = 1, color = "darkgreen")+
  annotate("text", x = 18, y= 0.8, label = "Weak", size = 10)+
  geom_hline(yintercept = 2, color = "darkgreen")+
  annotate("text", x = 18, y= 1.8, label = "Moderate", size = 10)+
  geom_hline(yintercept = 3, color = "darkgreen")+
  annotate("text", x = 18, y= 2.8, label = "Strong", size = 10)+
  geom_hline(yintercept = 4, color = "darkgreen")+
  annotate("text", x = 18, y= 3.8, label = "Very strong", size = 10)+
  annotate("text", x = 18.1, y= 4.8, label = "Extremely strong", size = 10)+
  annotate("text", x = 21, y= 4.8, label = "")+
  scale_fill_manual(values = c("#56B4E9", "#009E73" ,"#D55E00"))+
  theme(axis.title = element_text(size = 30), axis.text.x = element_text(size = 24, angle = 45), axis.text.y = element_text(size = 24), legend.title = element_text(size = 26), legend.text = element_text( size = 24))


ggsave("../images/plot_condition_2_resamples.png", height =8, width = 15)

# Sensitivity analysis
df_resampling = as.data.frame(rbind(cbind(log10(LR_resampleT_4R[,1]), rep("Secondary transfer", 1000)), cbind(log10(LR_resampleS_4R[,1]), rep("Size of FFGs", 1000)), cbind(log10(LR_resampleF_4R[,1]), rep("Number of FFGs", 1000))))

colnames(df_resampling) = c("log10_LR", "Data")
df_resampling$log10_LR = as.numeric(df_resampling$log10_LR)
df_resampling$Data = factor(df_resampling$Data, levels = c("Secondary transfer", "Size of FFGs", "Number of FFGs"))



plot_resampling = ggplot(data = df_resampling, aes(x = Data, y = log10_LR))+
  ylab("Log(LR)")+
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 45))+
  theme(legend.position = "none")+
  geom_boxplot(fill = "skyblue")
plot_resampling


ggsave("../inter-framework/images/resample.png")