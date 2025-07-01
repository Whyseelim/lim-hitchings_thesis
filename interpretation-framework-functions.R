pkgs=c("bnlearn","dplyr","stringr","Rgraphviz", "DiagrammeR","ggplot2","gRain", "tidyr", "kableExtra")
lapply(pkgs, library, character.only = TRUE)

# initialize an empty cpt

create_empty_cpts=function(nodes){
  cpts=vector(mode="list", length = length(nodes))
  names(cpts)=nodes
  return(cpts)
}

# Building a cpt

build_cpt = function (node, list, arcs){
  index=which(names(list)==substr(node,1,1))
  
  dimensions = as.vector(sapply(list[[index]]$var_names, length))
  
  if (nchar(node)>1){
    cpt_index = as.numeric(substr(node,nchar(node),nchar(node)))
  } else {
    cpt_index = 1
  }
  
  array = array(dimnames = list[[index]]$var_names, dim = dimensions, data = list[[index]]$cpt[[cpt_index]])
  
  dependencies = c(list[[index]]$dep)
  
  names = vector(mode = "character", length = length(dependencies))
  names[which(nchar(dependencies)>1)] = paste0(sub("_.*", "", dependencies[which(nchar(dependencies)>1)]), substr(node, 2, 2))
  names[which(nchar(dependencies)==1)] = dependencies[which(nchar(dependencies)==1)]
  
  dimnames=c(node, names)
  names(dimnames(array))=dimnames
  return(array)
}

# assigning gamma

assign_gamma = function(c, common = 0.014, rare = 0.00171){
  gamma = vector(mode = "numeric", length = length(c))
  gamma[which(c == TRUE)] = common
  gamma[which(c == FALSE)] = rare
  return(gamma)
}

# Background: resample f

resample_f = function (vector) {
  resample = sample(vector, size = length(vector), replace = TRUE)
  probability_vector = c((length(resample[which(resample == "0")]))/length(resample), (length(resample[which(resample == "1-10")]))/length(resample), (length(resample[which(resample == "11-20")]))/length(resample), (length(resample[which(resample == "21-30")]))/length(resample), (length(resample[which(resample == "31-40")]))/length(resample), (length(resample[which(resample == ">40")]))/length(resample))
  return(probability_vector)
}

# Background: resample s

resample_s = function (vector) {
  resample = sample(vector, size = length(vector), replace = TRUE)
  probability_vector = c((length(resample[which(resample == "small")]))/length(resample), (length(resample[which(resample == "medium")]))/length(resample), (length(resample[which(resample == "large")]))/length(resample))
  return(probability_vector)
}


# Build background class

build_B = function(n, gamma, vector_f = c(0.0158, 0.0317, 0.0635, 0.127, 0.254, 0.508), vector_s = c(0.767, 0.207, 0.026)){
  #Create node strings for each target fiber and dependencies
  nodes_s = paste0("[S", c(1:n), "]", collapse="")
  nodes_m = paste0("[M",c(1:n),"|F]", collapse = "")
  nodes_b = paste0("[B", c(1:n), "|M", c(1:n), ":S", c(1:n), "]", collapse="")
  
  #Build empty network
  
  bn_b = model2network(paste0("[F]",nodes_s, nodes_m, nodes_b))
  graphviz.plot(bn_b)
  
  # Dependency list
  names_s = c("small", "medium", "large")
  names_f = c("0", "1-10", "11-20", "21-30", "31-40", ">40")
  list_b = list(F = list(), S = list(), M = list(), B = list())
  
  list_b$F = list(dep = NULL, var_names = list(F = names_f), cpt = list(vector_f))
  
  cpts_s = rep(list(vector_s), n)
  
  list_b$S = list(dep = NULL, var_names = list(S = names_s), cpt = cpts_s)
  
  get_match = function (g, ffg){
    if(ffg == 0){
      mu = 0
    }else{
      k = c(1:ffg)
      mu = g*sum((1-g)^(k-1))
    }
    return(mu)
  }
  
  get_match_vector = function(f_vector, g){
    mu_avg = mean(mapply(get_match, g=g, ffg=f_vector))
    return(mu_avg)
  }
  
  get_match_list = function(g, f_list){
    mu_avgs = sapply(f_list, get_match_vector, g = g)
    return (mu_avgs)
  }
  
  
  f_dep_list = list(f_0 = 0, f_1 = 1:10, f_2 = 11:20, f_3 = 21:30, f_4 = 31:40, f_5 = 40:50)
  match = sapply(gamma, get_match_list, f_list = f_dep_list)
  no_match = 1-match
  match_to_cpt = function (index, match, no_match){
    matrix = as.vector(rbind(match[,index], no_match[,index]))
    return(matrix)
  }
  index=c(1:n)
  cpt_match = lapply(index, match_to_cpt, match, no_match)
  
  #simplified this to just match
  list_b$M = list(dep = c("F"), var_names = list (M = c("match", "no match"), F = names_f), cpt = cpt_match)
  
  names_x = c("none", "small", "medium", "large")
  
  list_b$B = list(dep = c("S_index", "M_index"), var_names = list (B = names_x, S = names_s, M = c("match", "no match")), cpt = rep(list(c(as.vector(rbind(c(rep(0,3)), diag(1, 3, 3))), rep (c(1, 0, 0, 0), 3))), n))
  
  # Define CPTs
  
  cpts_B = lapply(nodes(bn_b), build_cpt, list = list_b, arcs = arcs(bn_b))
  names(cpts_B) = nodes(bn_b)
  
  # fit
  bn_b_fit = custom.fit(bn_b, cpts_B)
  graphviz.plot(bn_b_fit)
  
  return(bn_b_fit)
}


# Transfer

# Transfer resample
resample_t = function (vector) {
  resample = sample(vector, size = length(vector), replace = TRUE)
  probability_vector = c((length(resample[which(resample == "none")]))/length(resample), (length(resample[which(resample == "small")]))/length(resample), (length(resample[which(resample == "medium")]))/length(resample), (length(resample[which(resample == "large")]))/length(resample))
  return(probability_vector)
}

# Write transfer vectors into CPTs

write_cpt_t = function (vector_t){
  cpt_vector = c(vector_t, 1, 0, 0, 0)
  return(cpt_vector)
}

# Build T
build_T = function(n, gamma_prime, cpts){
  #Create node strings for each target fiber and dependencies
  nodes_m_prime = paste0("[Mprime", c(1:n), "|H]", collapse = "")
  nodes_t = paste0("[T",c(1:n),"|H:Mprime", c(1:n), "]", collapse = "")
  
  #Build empty network
  
  bn_t = model2network(paste0("[H]", nodes_m_prime, nodes_t))
  
  # Dependency list
  names_x = c("none","small", "medium", "large")
  list_t = list(H = list(), M_prime = list(), T = list())
  names(list_t) = c ("H", "M", "T")
  
  list_t$H = list(dep = NULL, var_names = list(H = c("H1", "H2")), cpt = list(c(0.5, 0.5)))
  
  match = gamma_prime
  no_match = 1-match
  match_to_cpt_prime = function (index, match, no_match){
    matrix = as.vector(c(1, 0, match[index], no_match[index]))
    return(matrix)
  }
  index=c(1:n)
  cpt_match_prime = lapply(index, match_to_cpt_prime, match, no_match)
  
  list_t$M = list(dep = c("H"), var_names = list(M_prime = c("match", "no match"), H = c("H1", "H2")), cpt = cpt_match_prime)
  
  list_t$T = list(dep = c("Mprime_index", "H"), var_names = list (T = names_x, M_prime = c("match", "no match"), H = c("H1", "H2")), cpt = cpts)
  
  # Define CPTs
  
  cpts_T = lapply(nodes(bn_t), build_cpt, list = list_t, arcs = arcs(bn_t))
  names(cpts_T) = nodes(bn_t)
  
  # fit
  bn_t_fit = custom.fit(bn_t, cpts_T)
  graphviz.plot(bn_t_fit)
  
  return(bn_t_fit)
}

# Evidence CPT
#Given no background

B0 = diag(1, nrow = 4, ncol = 4)

# Small background

B1T0 = c(0, 1, 0, 0)


# Create function to sum up combinations

sum_comb = function (x, range){
  sum = x + range
  return (sum)
}

B_small = c(1:5)
sum_b1t1 = sapply(B_small, sum_comb, range = c(1:5))
sum_b1t2 = sapply(B_small, sum_comb, range = c(6:50))

B1T1 = c(0,(1-(length(which(sum_b1t1 > 5))/25)), (length(which(sum_b1t1 > 5))/25), 0)
B1T2 = c(0, 0, 1 - (length(which(sum_b1t2 > 50))/225), (length(which(sum_b1t2 > 50))/225))
B1T3 = c (0, 0, 0, 1)

B2T0 = c(0, 0, 1, 0)
B_medium = c(6:50)
sum_b2t1 = sapply(B_medium, sum_comb, range = c(1:5))
sum_b2t2 = sapply(B_medium, sum_comb, range = c(6:50))

B2T1 = c(0, 0, 1 - (length(which(sum_b2t1 > 50))/225), (length(which(sum_b2t1 > 50))/225))
B2T2 = c(0, 0, 1 - (length(which(sum_b2t2 > 50))/2025), (length(which(sum_b2t2 > 50))/2025))
B2T3 = c(0, 0, 0, 1)

# Big background

B3 = c(rep(c(0,0,0,1),4))

vector_e = c(as.vector(B0), B1T0, B1T1, B1T2 ,B1T3, B2T0, B2T1, B2T2, B2T3, B3)

# Build E class
build_E = function(n){
  #Create node strings for each target fiber and dependencies
  nodes_b = paste0("[B",c(1:n), "]", collapse = "")
  nodes_t = paste0("[T",c(1:n), "]", collapse = "")
  nodes_e = paste0("[E",c(1:n),"|T", c(1:n), ":", "B", c(1:n), "]", collapse = "")
  
  #Build empty network
  
  bn_e = model2network(paste0(nodes_b, nodes_t, nodes_e))
  
  
  # Dependency list
  names_x = c("none", "small", "medium", "large")
  list_e = list(B = list(), T = list(), E = list())
  
  list_e$B = list(dep = NULL, var_names = list(B = names_x), cpt = rep(list(c(rep(1/4,4))), n))
  
  
  list_e$T = list(dep = NULL, var_names = list(T = names_x), cpt = rep(list(c(rep(1/4,4))), n))
  
  list_e$E = list(dep = c("T_index", "B_index"), var_names = list (E= names_x, T = names_x, B = names_x), cpt = rep(list(vector_e),n))
  
  # Define CPTs
  
  cpts_E = lapply(nodes(bn_e), build_cpt, list = list_e, arcs = arcs(bn_e))
  names(cpts_E) = nodes(bn_e)
  
  # fit
  bn_e_fit = custom.fit(bn_e, cpts_E)
  graphviz.plot(bn_e_fit)
  
  return(bn_e_fit)
}

# Join network
join_bn = function (B, T, E){
  
  nodes = unique(c(nodes(B), nodes(T), nodes(E)))
  arcs = rbind(arcs(B), arcs(T), arcs(E))
  
  bn = empty.graph(nodes)
  arcs(bn) = arcs
  
  get_cpt = function(node, graph){
    index = which(nodes(graph) == node)
    return(graph[[index]]$prob)
  }
  
  cpt_b = lapply(nodes(B), get_cpt, graph = B)
  names(cpt_b) = nodes(B)
  
  cpt_t = lapply(nodes(T), get_cpt, graph = T)
  names(cpt_t) = nodes(T)
  
  E_nodes_index = which(substr(nodes(E), 1, 1) == "E")
  cpt_e = lapply(nodes(E)[E_nodes_index], get_cpt, graph = E)
  names(cpt_e) = nodes(E)[E_nodes_index]  
  
  cpts = c(cpt_b, cpt_t, cpt_e)
  
  bn_fit = custom.fit(bn, cpts)
  graphviz.plot(bn_fit)
  return(bn_fit)
  
}

# overall build BN
build_bn = function (n, gamma, 
                     vector_f = c(0.0158, 0.0317, 0.0635, 0.127, 0.254, 0.508),
                     vector_s = c(0.767, 0.207, 0.026), gamma_prime, cpts_t){
  E = build_E(n)
  B = build_B(n, gamma = gamma, vector_f = vector_f, vector_s = vector_s)
  T = build_T(n, gamma_prime = gamma_prime, cpts = cpts_t)
  
  bn = join_bn(B = B, T = T, E = E)
  return(bn)
  
}

# Evalauting the BN for LR
query_node = function (bn, nodes = "H", evidence){
  convert = compile(as.grain(bn))
  results = querygrain(convert, nodes = nodes, evidence = evidence)
  return(results)
}

get_LR = function (results){
  return(as.numeric(results[[1]][1]/results[[1]][2]))
}

get_LR_direct = function (bn, evidence){
  convert = compile(as.grain(bn))
  results = querygrain(convert, nodes = "H", evidence = evidence)
  return(as.numeric(results[[1]][1]/results[[1]][2]))
}


