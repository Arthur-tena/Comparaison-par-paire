library(survival)

library(doParallel)
library(parallel)
library(foreach)

# générer des variables de type time-to-event
# nécessite le package "survival"
# rentre en argument mu le paramètre de la censure ; lambda le paramètre des valeurs observées 
# donne un dataframe comportant les données tte
gener_tte=function(mu, lambda=1){
  X=round(rexp(n, rate=lambda),3)
  C=round(rexp(n, rate = mu),3)
  TT=pmin(X,C)
  delta=as.numeric(X==TT)
  delta=abs(1-delta)
  S=Surv(time = TT, event = delta)
  S=as.data.frame(S)
  return(S)
}

# générer des variables de type continue sous la forme de valeur absolue de loi normale
# rentre en argument mean la moyenne et sd l'écart-type de la loi normale
# donne un dataframe comportant les données continue
gener_continue=function(mean=1, sd=0){
  X=round(abs(rnorm(n, mean=mean, sd=sd)),3)
  return(data.frame(X))
}

# générer des variables de type binaire
# rentre en argument prob, la probabilité de tiré 1 pour une loi binomiale
# donne un dataframe comportant les données binaire
gener_binom=function(prob){
  X=rbinom(n, size = 1, prob = prob)
  X=factor(X)
  return(data.frame(X))
}

# rentre en argument treatmentdata le traitement que l'on veut étudier et L le nombre d'outcome (colonne)
# donne un vecteur de taille L comportant le type de variable de treatmentdata pour chaque outcome 
type_variable = function(treatmentdata, L) {
  type = rep(0, L)
  
  type = sapply(treatmentdata, function(x) {
    classe = class(x)
    
    if ("Surv" %in% classe) {
      return("tte")
    } else if ("numeric" %in% classe) {
      return("continue")
    } else {
      return("binaire")
    }
  })
  
  return(unname(type))
}

# extract_tte ne fonctionne que pour les variables de type tte
# rentre en argument treatmentdata le traitement que l'on veut étudier et l une outcome
# donne un sortie 2 vecteurs, le vecteur de temps observé "t_obs" et le vecteur de censure "censure"
extract_tte=function(treatmentdata,l){
  t_obs=treatmentdata[[l]][,1]
  censure=treatmentdata[[l]][,2]
  return(cbind(t_obs,censure))
}

# affect_crit se base sur les tableaux donné dans le document sur la GPC de Marc Buyse (2010)
# rentre en argument treatmentdata le nouveau traitement, controldata le traitement de contrôle et threshold le seuil
# donne en sortie une matrice de taille n1*n2xL de paire suivante si la paire est favorable, défavorablle, neutre ou non-informative

affect_crit=function(treatmentdata,controldata,threshold=0){
  n1 = nrow(treatmentdata)
  n2 = nrow(controldata)
  L = ncol(treatmentdata)
  
  if (ncol(treatmentdata) != ncol(controldata)) {
    stop("il n'y a pas le même nombre d'outcomes")
  }
  
  type1 = type_variable(treatmentdata, L)
  
  paire = matrix("", nrow = n1 * n2, ncol = L)
  
  p = 1 
  
  groupe = as.factor(rep(c("T", "C"), c(n1, n2)))
  comp = rbind(treatmentdata, controldata)
  comp = data.frame(groupe = groupe, outcome = comp)
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      for (l in 1:L) {
        if (type1[l] == "tte") {
          
          t_obs1 = extract_tte(comp[groupe=="T",], l+1)[,1]
          censure1 = extract_tte(comp[groupe=="T",], l+1)[,2]
          t_obs2 = extract_tte(comp[groupe=="C",], l+1)[,1]
          censure2 = extract_tte(comp[groupe=="C",], l+1)[,2]
          
          diff_tte = t_obs1[i] - t_obs2[j]
          
          if (censure1[i] == 0 && censure2[j] == 0) {
            paire[p, l] = ifelse(diff_tte > threshold, "favorable",
                                 ifelse(diff_tte < -threshold, "défavorable", "neutre"))
          } else if (censure1[i] == 1 && censure2[j] == 0) {
            paire[p, l] = ifelse(diff_tte > threshold, "favorable", "non-informative")
          } else if (censure1[i] == 0 && censure2[j] == 1) {
            paire[p, l] = ifelse(diff_tte > threshold, "non-informative", "défavorable")
          } else {
            paire[p, l] = "non-informative"
          }
        } 
        
        if (type1[l] == "continue") {
          diff = comp[groupe=="T", l+1][i] - comp[groupe=="C", l+1][j]
          paire[p, l] = ifelse(diff > threshold, "favorable",
                               ifelse(diff < -threshold, "défavorable", "neutre"))
        }
        if (type1[l] == "binaire") {
          paire[p, l] = ifelse(comp[groupe=="T", l+1][i] == 1 & comp[groupe=="C", l+1][j] == 0, "favorable",
                               ifelse(comp[groupe=="T", l+1][i] == 0 & comp[groupe=="C", l+1][j] == 1, "défavorable", "neutre"))
        }
      }
      p = p + 1 
    }
  }
  return(paire)
}

# rentre en argument une matrice de paire comportant des valeurs de type charactère donné par la fonction affect_crit
# donne en sortie le nombre de win, de lose et de tie
calcul_stat=function(paire){
  n1 = sqrt(nrow(paire))
  n2 = sqrt(nrow(paire))
  L = ncol(paire)
  
  N_w = 0
  N_l = 0
  N_t = 0
  
  for (p in 1:(n1 * n2)) {
    for (l in 1:L) {
      if (paire[p, l] != "non-informative") {
        if (paire[p, l] == "favorable") {
          N_w = N_w + 1
        } else if (paire[p, l] == "défavorable") {
          N_l = N_l + 1
        } else if (paire[p, l] == "neutre") {
          N_t = N_t + 1
        }
        break
      }
    }
  }
  return(c(N_w,N_l,N_t))
}

# GPC_WO_WR necesiite les package doParallele, parallele et foreach
# rentre en argument treatmentdata le nouveau traitrement, controldata le traitement de contrôle, threshold le seuil, 
#   p.val le test unilatéral ou bilatéral et n_perm le nombre de permutation
# donne en sortie une liste de 3 dataframe avec les résultats de la GPC, des WR et des WO et leur p-valeur,
#    l'intervalle de confiance pour ces 3 valeurs et le nombre de win,lose et tie 
GPC_WO_WR = function(treatmentdata, controldata, threshold = 0, p.val = c("one.sided", "two.sided"), n_perm = 1000) {
  
  n_cores = detectCores()-3  
  cl = makeCluster(n_cores)
  registerDoParallel(cl)
  
  n1 = nrow(treatmentdata)
  n2 = nrow(controldata)
  L = ncol(treatmentdata)
   
   groupe = as.factor(rep(c("T", "C"), c(n1, n2)))
   comp = rbind(treatmentdata, controldata)
   comp = data.frame(groupe = groupe, outcome = comp)

  paire = affect_crit(treatmentdata, controldata, threshold)
  stat_init = calcul_stat(paire)

  N_w = stat_init[1]
  N_l = stat_init[2]
  N_t = stat_init[3]

  Delta_obs = round((N_w - N_l) / (n1 * n2), 3)
  WR_obs = round(N_w/N_l,3)
  WO_obs = round((N_w+0.5*N_t)/(N_l+0.5*N_t),3)
  
  
  N_w_perm=rep(0, n_perm)
  N_l_perm=rep(0, n_perm)
  Delta_perm=rep(0, n_perm)
  WR_perm=rep(0, n_perm)
  WO_perm=rep(0, n_perm)
  
  Delta_perm_res = foreach(s = 1:n_perm, .combine = rbind, .packages = c("dplyr", "survival"), 
                           .export = c("affect_crit", "calcul_stat", "type_variable", "extract_tte")) %dopar% {
  
                             
                             comp_perm = comp
                             comp_perm$groupe = sample(comp$groupe)
                             
                             compT = subset(comp_perm, groupe == "T")[,-1]
                             compC = subset(comp_perm, groupe == "C")[,-1]
                             
                             paire_perm = affect_crit(compT, compC, threshold) 
                             stat_perm = calcul_stat(paire_perm)  
                             
                             N_w_perm = stat_perm[1] 
                             N_l_perm = stat_perm[2]  
                             N_t_perm = stat_perm[3]
                             
                             
                             Delta_perm = (N_w_perm - N_l_perm) / (n1 * n2)
                             WR_perm = N_w_perm/N_l_perm
                             WO_perm = (N_w_perm+0.5*N_t_perm)/(N_l_perm+0.5*N_t_perm)
                             
                             
                             return(unname(c(N_w_perm,N_l_perm,Delta_perm, WR_perm,WO_perm)))
                           }
  
  N_w_perm = Delta_perm_res[, 1]
  N_l_perm = Delta_perm_res[, 2]
  Delta_perm = Delta_perm_res[, 3]
  WR_perm = Delta_perm_res[, 4]
  WO_perm = Delta_perm_res[, 5]
  
  CI_GPC = quantile(Delta_perm, c(0.025, 0.975), na.rm = TRUE)
  #CI_GPC = c(Delta_obs - CI_GPC[1], Delta_obs + CI_GPC[2])
  
  CI_WR = quantile(WR_perm, c(0.025, 0.975), na.rm = TRUE)
  #CI_WR = c(WR_obs - CI_WR[1], WR_obs + CI_WR[2])
  
  CI_WO = quantile(WO_perm, c(0.025, 0.975), na.rm = TRUE)
  #CI_WO = c(WO_obs - CI_WO[1], WO_obs + CI_WO[2])
  
  hist(Delta_perm, breaks = 30, main = "Distribution de Δ sous H0 (permutation), GPC", 
       xlab = "Δ permuté", col = "lightblue", border = "black", xlim=c(-1, 1))
  abline(v = Delta_obs, col = "red", lwd = 2, lty = 2)
  abline(v=CI_GPC[1], col = "green", lwd = 2, lty = 2)
  abline(v=CI_GPC[2], col = "green", lwd = 2, lty = 2)
  
  hist(WR_perm, breaks = 30, main = "Distribution de Δ sous H0 (permutation), WR", 
       xlab = "Δ permuté", col = "lightblue", border = "black", xlim=c(0, 5))
  abline(v = WR_obs, col = "red", lwd = 2, lty = 2)
  abline(v=CI_WR[1], col = "green", lwd = 2, lty = 2)
  abline(v=CI_WR[2], col = "green", lwd = 2, lty = 2)
  
  hist(WO_perm, breaks = 30, main = "Distribution de Δ sous H0 (permutation), WO", 
       xlab = "Δ permuté", col = "lightblue", border = "black", xlim=c(0, 5))
  abline(v = WO_obs, col = "red", lwd = 2, lty = 2)
  abline(v=CI_WO[1], col = "green", lwd = 2, lty = 2)
  abline(v=CI_WO[2], col = "green", lwd = 2, lty = 2)
  
  # boxplot(WO_perm, main = "Distribution de WO sous H0", col = "lightblue")
  # points(1, WO_obs, col = "red", pch = 19, cex = 1.5)
  
  sigma_GPC = sd(Delta_perm)
  z_GPC = (Delta_obs) / sigma_GPC
  
  sigma_WR = sd(WR_perm)
  z_WR = (WR_obs - mean(WR_perm)) / sigma_WR
  
  sigma_WO = sd(WO_perm)
  z_WO = (WO_obs - mean(WO_perm)) / sigma_WO
  
  p_value_GPC = ifelse(p.val == "one.sided", mean(Delta_perm >= Delta_obs), 2*mean(abs(Delta_perm) >= abs(Delta_obs)))
  p_value_WR = ifelse(p.val == "one.sided", mean(WR_perm >= WR_obs), 2*mean(abs(WR_perm) >= abs(WR_obs)))
  p_value_WO = ifelse(p.val == "one.sided", mean(WO_perm >= WO_obs),2* mean(abs(WO_perm) >= abs(WO_obs)))
  
  
  signif_GPC = dplyr::case_when(
    p_value_GPC < 0.001 ~ "***",
    p_value_GPC < 0.01 ~ "**",
    p_value_GPC < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  signif_WR = dplyr::case_when(
    p_value_WR < 0.001 ~ "***",
    p_value_WR < 0.01 ~ "**",
    p_value_WR < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  signif_WO = dplyr::case_when(
    p_value_WO < 0.001 ~ "***",
    p_value_WO < 0.01 ~ "**",
    p_value_WO < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  data1 <- data.frame(
    Method = c("GPC", "Win Ratio (WR)", "Win Odds (WO)"),
    Estimate = c(Delta_obs, WR_obs, WO_obs),
    Z_score = c(z_GPC, z_WR, z_WO),
    P_value = c(p_value_GPC, p_value_WR, p_value_WO),
    Signif. = c(signif_GPC, signif_WR, signif_WO)
  )
  
  data2 = data.frame(
    Method = c("GPC", "WR", "WO"),
    CI_lower = c(CI_GPC[1], CI_WR[1], CI_WO[1]),
    CI_upper = c(CI_GPC[2], CI_WR[2], CI_WO[2])
  )
  
  #data3 = data.frame(Nb_win = N_w, Nb_lose = N_l, Nb_tie = N_t, row.names = "")
  
  stopCluster(cl)
  
  return(list(results = data1, confidence_intervals = data2))
}
