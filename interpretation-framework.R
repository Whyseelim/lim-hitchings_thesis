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


bn_condition_1_4R = build_bn(n = n_condition_1, 
                             gamma = gamma_condition_1_4R, 
                             gamma_prime =  gamma_prime_condition_1_4R, 
                             vector_s = vector_s, cpts_t = cpts_t_condition_1)


bn_condition_1_4C = build_bn(n = n_condition_1, 
                             gamma = gamma_condition_1_4C,
                             gamma_prime =  gamma_prime_condition_1_4C, 
                             vector_s = vector_s, cpts_t = cpts_t_condition_1)

bn_condition_1_1R3C = build_bn(n = n_condition_1, 
                               gamma = gamma_condition_1_1R3C, 
                               gamma_prime =  gamma_prime_condition_1_1R3C, 
                               vector_s = vector_s, 
                               cpts_t = cpts_t_condition_1)


bn_condition_1_2R2C = build_bn(n = n_condition_1, 
                               gamma = gamma_condition_1_2R2C, 
                               gamma_prime =  gamma_prime_condition_1_2R2C, 
                               vector_s = vector_s, 
                               cpts_t = cpts_t_condition_1)

bn_condition_1_3R1C = build_bn(n = n_condition_1, 
                               gamma = gamma_condition_1_3R1C, 
                               gamma_prime =  gamma_prime_condition_1_3R1C, 
                               vector_s = vector_s, 
                               cpts_t = cpts_t_condition_1)