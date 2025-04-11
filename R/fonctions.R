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
    
    if ("Surv" %in% classe | "matrix" %in% classe) {
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
eval_diff = function(i, j, l, type1, comp, threshold) {
  indices_T = which(comp$groupe == "T")
  indices_C = which(comp$groupe == "C")
  
  if (type1[l] == "tte") {
    t_obs1 = extract_tte(comp[indices_T, ], l + 1)[, 1]
    censure1 = extract_tte(comp[indices_T, ], l + 1)[, 2]
    t_obs2 = extract_tte(comp[indices_C, ], l + 1)[, 1]
    censure2 = extract_tte(comp[indices_C, ], l + 1)[, 2]
    
    if (is.na(t_obs1[i]) || is.na(t_obs2[j])) {
      return("non-informative")
    }
    
    diff_tte = t_obs1[i] - t_obs2[j]
    
    if (censure1[i] == 0 && censure2[j] == 0) {
      return(ifelse(diff_tte > threshold, "favorable",
                    ifelse(diff_tte < -threshold, "défavorable", "neutre")))
    } else if (censure1[i] == 1 && censure2[j] == 0) {
      return(ifelse(diff_tte > threshold, "favorable", "non-informative"))
    } else if (censure1[i] == 0 && censure2[j] == 1) {
      return(ifelse(diff_tte > threshold, "non-informative", "défavorable"))
    } else {
      return("non-informative")
    }
  } 
  
  if (type1[l] == "continue") {
    diff = comp[indices_T, l + 1][i] - comp[indices_C, l + 1][j]
    
    if (is.na(diff)) {
      return("non-informative")
    }
    
    return(ifelse(diff > threshold, "favorable",
                  ifelse(diff < -threshold, "défavorable", "neutre")))
  }
  
  if (type1[l] == "binaire") {
    val_T = comp[indices_T, l + 1][i]
    val_C = comp[indices_C, l + 1][j]
    
    if (is.na(val_T) || is.na(val_C)) {
      return("non-informative")
    }
    
    return(ifelse(val_T == 1 & val_C == 0, "favorable",
                  ifelse(val_T == 0 & val_C == 1, "défavorable", "neutre")))
  }
}

affect_crit_strata = function(treatmentdata, controldata, threshold = 0, strata = NULL) {
  n1 = nrow(treatmentdata)
  n2 = nrow(controldata)
  
  ### On regarde le nb d'outcome 
  if (is.null(strata)) {
    L = ncol(treatmentdata)
  } else {
    L = ncol(treatmentdata) - 1
  }
  
  col = colnames(treatmentdata) 
  
  if (ncol(treatmentdata) != ncol(controldata)) {
    stop("Il n'y a pas le même nombre d'outcomes")
  }
  
  
  
  ### On crée les tableaux qui nous seront utile suivant s'il y a de la stratification ou non 
  
  if (is.null(strata)) {
    type1 <- type_variable(treatmentdata, L) 
    groupe <- as.factor(rep(c("T", "C"), c(n1, n2)))
    comp <- data.frame(groupe = groupe, outcome = rbind(treatmentdata, controldata))
  } else {
    type1 <- type_variable(treatmentdata[, -which(colnames(treatmentdata) == "strata")], L)
    
    groupe <- as.factor(rep(c("T", "C"), c(n1, n2)))
    comp <- data.frame(groupe = groupe, 
                       outcome = rbind(treatmentdata[, -which(colnames(treatmentdata) == "strata")], 
                                       controldata[, -which(colnames(controldata) == "strata")]),
                       strata = c(treatmentdata[["strata"]], controldata[["strata"]]))
  }
  
  
  ### On écrase le nom du groupe à chaque fois (problème lors de la parallélisation sinon)
  
  colnames(comp)[which(colnames(comp) == "groupe")] <- "groupe"
  
  outcome_cols <- colnames(comp)[grep("outcome", colnames(comp))]
  for (l in 1:L) {
    if (l <= length(outcome_cols)) {
      colnames(comp)[which(colnames(comp) == outcome_cols[l])] <- paste("Y_", l, sep = "")
    }
  }
  
  if ("strata" %in% colnames(comp)) {
    colnames(comp)[which(colnames(comp) == "strata")] <- "strata"
  }
  
  ### On trie suivant les strates
  
  if (!is.null(strata)) {
    if (!any("strata" %in% colnames(comp))) {
      stop("La colonne strata spécifiée n'existe pas dans les données")
    }
    comp = comp[order(comp[["strata"]]), ]
  }
  
  ### On crée les quantités qui nous intéressent U pour les strates sous forme de liste
  
  matrices_list = list() # les listes de paires
  U_list = list() # les listes de la matrice U_ijk
  TT_list = list() # les listes des quantités T qui sont T=sum_(1:N)U_iD_i où D_i=1 si i est dans le grp tr et 0 sinon
  V_list = list() # les listes des variances de TT
  
  
  
  if (is.null(strata)) {
    U = matrix(NA, ncol = n1, nrow = n2)
    pairs = expand.grid(i = 1:n1, j = 1:n2)
    paire = matrix("", nrow = nrow(pairs), ncol = L)
    
    for (l in 1:L) {
      paire[, l] = mapply(eval_diff, pairs$i, pairs$j, 
                          MoreArgs = list(l = l, type1 = type1, comp = comp, threshold = threshold))
    }
    
    for (idx in 1:nrow(pairs)) {
      i <- pairs$i[idx]
      j <- pairs$j[idx]
      
      ligne_paire <- paire[idx, ]
      
      ligne_paire[is.na(ligne_paire)] <- "non-informative"
      
      res <- ligne_paire[ligne_paire != "non-informative"][1]
      if (!is.na(res)) {
        U[i, j] <- ifelse(res == "favorable", 1, 
                          ifelse(res == "défavorable", -1, 0))
      } else {
        U[i, j] <- 0
      }
    }
    
    V_T=rowsum(U, na.rm = TRUE)
    V_C=colSums(U, na.rm= TRUE) 
    V_f= sum(V_T + V_C)
    TT=sum(V_T)
    V=((n1*n2)/((n1+n2)(n1+n2-1)))*sum(V_f**2)
    Z=TT/sqrt(V)
    matrices_list[["all"]] = paire
    U_list[["all"]] = U
    
  } else {
    for (s in unique(comp$strata)) { 
      comp_s = subset(comp, strata == s) # on crée un data set par strate
      n_T = sum(comp_s$groupe == "T")  
      n_C = sum(comp_s$groupe == "C") 
      U = matrix(NA, ncol = n_T, nrow = n_C) # la matrice U_ijk 
      
      if (n_T > 0 & n_C > 0) {  
        pairs = expand.grid(i = 1:n_T, j = 1:n_C)
        
        indices_T = which(comp_s$groupe == "T")
        indices_C = which(comp_s$groupe == "C")
        
        paire = matrix("", nrow = nrow(pairs), ncol = L)
        
        for (l in 1:L) {
          paire[, l] = sapply(1:nrow(pairs), function(idx) {
            eval_diff(pairs$i[idx], pairs$j[idx], l, type1, comp_s, threshold) # on évalue chaque paire de patient 
          })
        }
        
        
        for (idx in 1:nrow(pairs)) {
          i <- pairs$i[idx]
          j <- pairs$j[idx]
          
          ligne_paire <- paire[idx, ]
          
          ligne_paire[is.na(ligne_paire)] <- "non-informative"
          
          res <- ligne_paire[ligne_paire != "non-informative"][1] # on trie les valeurs noon-informative pour prendre le critère suivant
          
          if (!is.na(res)) {
            U[j, i] <- ifelse(res == "favorable", 1,
                              ifelse(res == "défavorable", -1, 0))
          } else {
            U[j, i] <- 0
          }
        }
        
        V_T = rowSums(U, na.rm = TRUE) # On crée le vecteur V des patients avec le nv. tr.
        V_C = colSums(U, na.rm = TRUE) # On crée le vecteur V des patients avec le tr. de contrôle
        V_f = sum(V_T + V_C)
        TT = sum(V_T) 
        V = ((n_T * n_C) / ((n_T + n_C) * (n_T + n_C - 1))) * sum(V_f^2)
        
        V_list[[paste0("strata_", s)]] = V
        TT_list[[paste0("strata_", s)]] = TT
        
        matrices_list[[paste0("strata_", s)]] = paire
        U_list[[paste0("strata_", s)]] = U
      }
    }
    
    matrices_list[["all"]] = do.call(rbind, matrices_list) # On crée la matrice de paire globale
    #U_list[["all"]] = do.call(rbind, U_list)
    V = do.call(sum, V_list)# lorsque l'on stratifie, les variance s'ajoutent
    TT = do.call(sum, TT_list)

     Z=TT/sqrt(V)
    }
  
  return(list(
    paire_glob = matrices_list[["all"]],
    paire_strata = matrices_list[names(matrices_list) != "all"],
    U_strata = U_list[names(matrices_list) != "all"],
    #U_global = U_list[["all"]],
    Z = Z,
    V = V,
    TT = TT
  ))
}

# rentre en argument une matrice de paire comportant des valeurs de type charactère donné par la fonction affect_crit
# et la colonne de stratification
# donne en sortie le nombre de win, de lose et de tie
calcul_stat_strata = function(paire, strata=NULL) {
  if(is.null(strata)){
 
    eval_ligne = function(ligne) {
      valeur = ligne[ligne != "non-informative"][1]
      if (is.na(valeur)) return(NA) 
      return(valeur)
    }
    
    resultats = apply(paire, 1, eval_ligne)
    
    N_w = sum(resultats, na.rm = TRUE)
    N_l = sum(resultats, na.rm = TRUE)
    N_t = sum(resultats, na.rm = TRUE)
    
  } else {
    
    N_w_list=list()
    N_l_list=list()
    N_t_list=list()
    
    for (s in unique(strata)) { 
      paire_s = paire[[paste0("strata_", s)]]
      
      eval_ligne = function(ligne) {
        valeur = ligne[ligne != "non-informative"][1]
        if (is.na(valeur)) return(NA) 
        return(valeur)
      }
      
      resultats = apply(paire_s, 1, eval_ligne)
      
      N_w_list[[as.character(s)]] = sum(resultats == "favorable", na.rm = TRUE)
      N_l_list[[as.character(s)]] = sum(resultats == "défavorable", na.rm = TRUE)
      N_t_list[[as.character(s)]] = sum(resultats == "neutre", na.rm = TRUE)
    }
    N_w = sum(unlist(N_w_list))
    N_l = sum(unlist(N_l_list))
    N_t = sum(unlist(N_t_list))
    
  }
  
  return(c(N_w, N_l, N_t))
}

# GPC_WO_WR necesiite les package doParallele, parallele et foreach
# rentre en argument treatmentdata le nouveau traitrement, controldata le traitement de contrôle, threshold le seuil, 
#   p.val le test unilatéral ou bilatéral et n_perm le nombre de permutation
# donne en sortie une liste de 3 dataframe avec les résultats de la GPC, des WR et des WO et leur p-valeur,
#    l'intervalle de confiance pour ces 3 valeurs et le nombre de win,lose et tie 
GPC_WO_WR_strata = function(treatmentdata, controldata, threshold = 0, p.val = c("one.sided", "two.sided"), n_perm = 1000, strata=NULL, histo=TRUE) {
  
  n_cores = parallel::detectCores()/2 -2 
  cl = makeCluster(n_cores)
  registerDoParallel(cl)
  
  n1 = nrow(treatmentdata)
  n2 = nrow(controldata)
  col = colnames(treatmentdata)
  
  if (is.null(strata)){
    L = ncol(treatmentdata)
    groupe = as.factor(rep(c("T", "C"), c(n1, n2)))
    comp = rbind(treatmentdata, controldata)
    comp = data.frame(groupe = groupe, outcome = comp)
  } else {
    L = ncol(treatmentdata)-1
    groupe = as.factor(rep(c("T", "C"), c(n1, n2)))
    comp = data.frame(groupe = groupe, 
                      outcome = rbind(treatmentdata[,-which(col=="strata")], controldata[,-which(col=="strata")]),
                      strata = c(treatmentdata$strata, controldata$strata))
  }
  af=affect_crit_strata(treatmentdata = treatmentdata, controldata = controldata, threshold = threshold, strata = strata)
  if(is.null(strata)){
    paire=af$paire_glob
  } else {
    paire=af$paire_strata
  }
  
  stat_init = calcul_stat_strata(paire, strata = strata)
  
  z_WR=af$Z
  
  N_w = stat_init[1]
  N_l = stat_init[2]
  N_t = stat_init[3]
  
  Delta_obs = round((N_w - N_l) / (N_w+N_l+N_t), 3)
  WR_obs = round(N_w/N_l,3)
  WO_obs = round((N_w+0.5*N_t)/(N_l+0.5*N_t),3)
  
  Delta_perm=rep(0, n_perm)
  
  
  Delta_perm_res = foreach(p = 1:n_perm, .combine = rbind, .packages = c("dplyr", "survival"), 
                           .export = c("eval_diff","affect_crit_strata", "calcul_stat_strata", "type_variable", "extract_tte")) %dopar% {
                             set.seed(p)
                             if(is.null(strata)){
                               comp_perm=comp
                               comp_perm$groupe = sample(comp_perm$groupe)
                               
                               compT = subset(comp_perm, groupe == "T")[,-1]
                               compC = subset(comp_perm, groupe == "C")[,-1]
                               
                               paire_perm = affect_crit_strata(compT, compC, threshold,strata)$paire_global
                             }else {
                               comp_perm=comp
                               comp_perm$groupe_permute <- ave(comp_perm$groupe, comp_perm$strata, FUN = function(x) sample(x))
                               sum(comp_perm$groupe_permute=="T" & comp_perm$strata==8)
                               compT = subset(comp_perm, groupe_permute == "T")[,-c(1,L+3)]
                               compC = subset(comp_perm, groupe_permute == "C")[,-c(1,L+3)]
                               
                               paire_perm  = affect_crit_strata(compT, compC, threshold,strata)$paire_strata
                             }
                             
                             stat_perm = calcul_stat_strata(paire_perm, strata = strata)  
                             
                             N_w_perm = stat_perm[1] 
                             N_l_perm = stat_perm[2]  
                             N_t_perm = stat_perm[3]
                             
                             
                             Delta_perm = (N_w_perm - N_l_perm) / (N_w_perm + N_l_perm + N_t_perm)
                             
                             
                             return(unname(Delta_perm))
                           }
  
  Delta_perm = Delta_perm_res
  
  
  quantile_GPC = quantile(Delta_perm, c(0.025, 0.975), na.rm = TRUE)
  CI_GPC = c(Delta_obs + quantile_GPC[1], Delta_obs + quantile_GPC[2])
  
  CI_WR = exp(log(WR_obs) + c(-1.96, 1.96) * sqrt(1/N_w + 1/N_l))
  CI_WO = exp(log(WO_obs) + c(-1.96, 1.96) * sqrt(1/N_w + 1/N_l + 1/(n1*n2)))
  if(histo){
    hist(Delta_perm, breaks = 30, main = "Distribution de Δ sous H0 (permutation)", 
         xlab = "Δ permuté", col = "lightblue", border = "black", xlim=c(-1, 1))
    abline(v = Delta_obs, col = "red", lwd = 2, lty = 2)
    abline(v=CI_GPC[1], col = "green", lwd = 2, lty = 2)
    abline(v=CI_GPC[2], col = "green", lwd = 2, lty = 2)
    abline(v=0, col='black', lwd = 2)
    legend('topright', col=c("green","red","black"), legend = c("95% CI", "Δ_obs", "H0"), lwd=c(2,2,1), lty = c(2,2,1))
  }
  
  #sigma_GPC = 1/2*(quantile_GPC[2]-quantile_GPC[2])/1.96
  s1= sum(Delta_perm >=Delta_obs)
  s2 = sum(abs(Delta_perm)>=abs(Delta_obs))
  
  var_log_WO = (N_w + N_l + N_t) / ((N_w + 0.5 * N_t) * (N_l + 0.5 * N_t))
  SD_log_WO = sqrt(var_log_WO)
  z_WO = log(WO_obs) / SD_log_WO
  
  #p_value_GPC = ifelse(p.val == "one.sided", pnorm(-Delta_obs/sigma_GPC), 2*pnorm(-Delta_obs/sigma_GPC))
  p_value_GPC = ifelse(p.val == "one.sided", s1/n_perm, s2/n_perm)
  p_value_WR = ifelse(p.val == "one.sided", 1-pnorm(z_WR), 2*(1-pnorm(abs(z_WR))))
  p_value_WO = ifelse(p.val == "one.sided", 1-pnorm(z_WO), 2*(1-pnorm(abs(z_WO))))
  
  
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
    p_value_WR < 0.001 ~ "***",
    p_value_WR < 0.01 ~ "**",
    p_value_WR < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  data1 = data.frame(
    Method = c("GPC", "Win Ratio (WR)", "Win Odds (WO)"),
    Estimate = c(Delta_obs, WR_obs, WO_obs),
    Z_score = c(".", z_WR, z_WO),
    P_value = c(p_value_GPC, p_value_WR, p_value_WR),
    Signif. = c(signif_GPC, signif_WR, signif_WO)
  )
  
  data2 = data.frame(
    Method = c("GPC", "WR", "WO"),
    CI_lower = c(CI_GPC[1], CI_WR[1], CI_WO[1]),
    CI_upper = c(CI_GPC[2], CI_WR[2], CI_WO[2])
  )
  
  stopCluster(cl)
  
  return(list(results = data1, confidence_intervals = data2))
}

# treatmentdata
# controldata
# set.seed(4) # 4 donne une p-valeur supéreiur à 1
# strata = sample(rep(c(1,3,5,8), each = 10))
# treatmentdata = as.data.frame(cbind(gener_continue(2.5, 1.5), gener_continue(2, 4), strata))
# controldata = as.data.frame(cbind(gener_continue(2.5, 1.5), gener_continue(2, 4), strata))
# GPC_WO_WR_strata(treatmentdata,controldata, threshold = 0.2, p.val="two.sided", n_perm = 1000, strata=strata, histo = F)
