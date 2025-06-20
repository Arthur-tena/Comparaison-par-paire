n_sim = 2000
n=200


nb_core = parallel::detectCores() - 2
cl = makeCluster(nb_core)
registerDoParallel(cl)

start_time = Sys.time()

results = foreach(s = 1:n_sim, .combine = rbind, .packages = c("dplyr", "survival", "parallel", "foreach", "doParallel", "WINS"), .export = c("win.stat")) %dopar% {
  set.seed(s)
  
  id=1:(2*n) 
  arm=rep(c("T","C"), each=200)
  
  U1 = runif(2*n) 
  U2 = runif(2*n) 
  
  lambda1 = 0.1; k1 = 2
  lambda2 = 0.12; k2 = 1.7
  
  fup_censure1=round(rweibull(2*n, shape = 3, scale =300),3)
  summary(fup_censure1)
  fup_censure2=round(rweibull(2*n, shape = 3, scale =300),3)
  
  Z = ifelse(arm == "T", 1, 0)
  beta = -2
  
  Time1 = round((-log(1 - U1)) / (lambda1 * exp(beta * Z)^(1 / k1)),3)
  summary(Time1)
  summary(Time1[Z==1])
  summary(Time1[Z==0])
  sum(Time1[Z==1]<5)
  
  Time_1 = pmin(Time1,  fup_censure1)
  
  
  delta1=as.numeric(fup_censure1==Time_1)
  sum(delta1==1)
  
  
  Time2 = round((-log(1 - U2)) / (lambda2 * exp(beta * Z)^(1 / k2)),3)
  

  
  Time_2 = pmin(Time2,  fup_censure2)
  
  
  delta2=as.numeric(fup_censure2==Time_2)
  
  
  stratum=sample(rep(c(1,3,5,8), each = 1))
  dataT=data.frame(Y_1 = Time_1[Z==1], Delta_1 = delta1[Z==1], Y_2 = Time_2[Z==1], Delta_2  = delta2[Z==1], stratum = stratum)
  dataC=data.frame(Y_1 = Time_1[Z==0], Delta_1 = delta1[Z==0], Y_2 = Time_2[Z==0], Delta_2 = delta2[Z==0], stratum = stratum)
  
  
  data1=rbind(dataT,dataC)
  
  summT = summary.data.frame(dataT)
  
  min_Y_1T = as.numeric(stringr::str_extract(summT[1,1], "\\d+\\.\\d+"))
  min_Y_2T = as.numeric(stringr::str_extract(summT[1,3], "\\d+\\.\\d+"))
  
  median_Y_1T = as.numeric(stringr::str_extract( summT[3,1], "\\d+\\.\\d+"))
  median_Y_2T = as.numeric(stringr::str_extract( summT[3,3], "\\d+\\.\\d+"))
  
  max_Y_1T =as.numeric(stringr::str_extract(summT[6,1], "\\d+\\.\\d+"))
  max_Y_2T =as.numeric(stringr::str_extract(summT[6,3], "\\d+\\.\\d+"))
  
  summC = summary.data.frame(dataC)
  
  min_Y_1C = as.numeric(stringr::str_extract(summC[1,1], "\\d+\\.\\d+"))
  min_Y_2C = as.numeric(stringr::str_extract(summC[1,3], "\\d+\\.\\d+"))
  
  
  median_Y_1C = as.numeric(stringr::str_extract( summC[3,1], "\\d+\\.\\d+"))
  median_Y_2C = as.numeric(stringr::str_extract( summC[3,3], "\\d+\\.\\d+"))
  
  max_Y_1C =as.numeric(stringr::str_extract(summC[6,1], "\\d+\\.\\d+"))
  max_Y_2C =as.numeric(stringr::str_extract(summC[6,3], "\\d+\\.\\d+"))
  
  censure_rateT1=sum(dataT$Delta_1==1)/n 
  censure_rateC1=sum(dataC$Delta_1==1)/n 
  
  censure_rateT2=sum(dataT$Delta_2==1)/n 
  censure_rateC2=sum(dataC$Delta_2==1)/n 
  
  data=data.frame(id=id, arm = arm, Y_1= data1[,1],Delta_1=data1[,2], Y_2 = data1[,3], Delta_2 = data1[,4], stratum = data1[,5] )
  
  result =
    win.stat(
      data = data,
      ep_type = c("tte", "tte"),
      stratum.weight = "equal",
      tau = c(0,0),
      arm.name = c("T", "C"),
      alpha = 0.05,
      digit = 3,
      pvalue = "two-sided",
      priority = c(1,2),
      summary.print = FALSE
    )
  
  win_edp1=sum(result$summary_ep$Trt_Endpoint1[,2])
  loose_edp1=sum(result$summary_ep$Con_Endpoint1[,2])
  tie_edp1=4*2500-sum(win_edp1+loose_edp1)
  
  win_edp2=sum(result$summary_ep$Trt_Endpoint2[,2])
  loose_edp2=sum(result$summary_ep$Con_Endpoint2[,2])
  tie_edp2=tie_edp1-(loose_edp2+win_edp2)
  
  
  val_GPC = result$Win_statistic$Net_Benefit[1]
  val_WR  = result$Win_statistic$Win_Ratio[1]
  val_WO  = result$Win_statistic$Win_Odds[1]
  
  p_val_GPC = result$p_value[2]
  p_val_WR  = result$p_value[1]
  p_val_WO  = result$p_value[3]
  resultats=c(val_GPC,val_WR,val_WO,win_edp1,loose_edp1,tie_edp1,win_edp2,loose_edp2,tie_edp2,min_Y_1C, min_Y_2C,max_Y_1C, max_Y_2C,median_Y_1C, median_Y_2C,min_Y_1T, min_Y_2T,max_Y_1T, max_Y_2T,median_Y_1T, median_Y_2T, censure_rateC1, censure_rateT1,censure_rateC2, censure_rateT2, p_val_GPC, p_val_WR, p_val_WO)
  
  
  return(unname(resultats))
}
stopCluster(cl)
end_time = Sys.time()

execution_time = end_time - start_time

df_1=results
col_names = c("val_GPC", "val_WR", "val_WO",
              "win_edp1", "loose_edp1", "tie_edp1",
              "win_edp2", "loose_edp2", "tie_edp2",
              "min_Y_1C", "min_Y_2C", "max_Y_1C", "max_Y_2C", "median_Y_1C", "median_Y_2C",
              "min_Y_1T", "min_Y_2T", "max_Y_1T", "max_Y_2T", "median_Y_1T", "median_Y_2T",
              "censure_rateC1", "censure_rateT1", "censure_rateC2", "censure_rateT2",
              "p_val_GPC", "p_val_WR", "p_val_WO")


results_df1 = as.data.frame(df_1)


colnames(results_df1) = col_names

GPC=as.numeric(results_df1[,1])
WR= as.numeric(results_df1[,2])
WO= as.numeric(results_df1[,3])

Win   = round(c(mean(as.numeric(results_df1[,4])), mean(as.numeric(results_df1[,7])), mean(as.numeric(results_df1[,4])) + mean(as.numeric(results_df1[,7]))  ))

Loose = round(c(mean(as.numeric(results_df1[,5])), mean(as.numeric(results_df1[,8])), mean(as.numeric(results_df1[,5]))+mean(as.numeric(results_df1[,8]))  ))

Tie   = round(c(mean(as.numeric(results_df1[,6])), mean(as.numeric(results_df1[,9])), mean(as.numeric(results_df1[,9])) ))



min_Y_1C=mean(as.numeric(results_df1[,15]))
min_Y_2C=mean(as.numeric(results_df1[,11]))

max_Y_1C=mean(as.numeric(results_df1[,12]))
max_Y_2C=mean(as.numeric(results_df1[,13]))

median_Y_1C=mean(as.numeric(results_df1[,14]))
median_Y_2C=mean(as.numeric(results_df1[,15]))

min_Y_1T=mean(as.numeric(results_df1[,16]))
min_Y_2T=mean(as.numeric(results_df1[,15]))

max_Y_1T=mean(as.numeric(results_df1[,18]))
max_Y_2T=mean(as.numeric(results_df1[,19]))

median_Y_1T=mean(as.numeric(results_df1[,20]))
median_Y_2T=mean(as.numeric(results_df1[,21]))


censure_rateC1= mean(as.numeric(results_df1[,22]))
censure_rateT1= mean(as.numeric(results_df1[,23]))

censure_rateC2= mean(as.numeric(results_df1[,24]))
censure_rateT2= mean(as.numeric(results_df1[,25]))

dfC=data.frame(row.names = c("min", "median", "max"), c(min_Y_1C, median_Y_1C, max_Y_1C),c(min_Y_2C, median_Y_2C, max_Y_2C))
colnames(dfC) = c("Y_1_C (tte)", "Y_2_C (tte)")
dfT=data.frame(row.names = c("min", "median", "max"), c(min_Y_1T, median_Y_1T, max_Y_1T),c(min_Y_2T, median_Y_2T, max_Y_2T))
colnames(dfT) = c("Y_1_T (tte)", "Y_2_T (tte)")



p_valGPC = results_df1[,26]
p_valWR = results_df1[,25]
p_valWO = results_df1[,28]

nb_pval_GPC = sum(p_valGPC< 0.05)
nb_pval_WR = sum(p_valWR< 0.05)
nb_pval_WO = sum(p_valWO< 0.05)





quantileGPC_lwr=quantile(GPC, 0.025, na.rm = T)
quantileGPC_upr=quantile(GPC, 0.955, na.rm = T)
étendue_GPC= round(quantileGPC_upr - quantileGPC_lwr,4)

quantileWR_lwr=quantile(WR, 0.025, na.rm = T)
quantileWR_upr=quantile(WR, 0.955, na.rm = T)
étendue_WR= round(quantileWR_upr - quantileWR_lwr,4)

quantileWO_lwr=quantile(WO, 0.025)
quantileWO_upr=quantile(WO, 0.955)
étendue_WO= round(quantileWO_upr - quantileWO_lwr,4)



df3 = data.frame(
  row.names = c("endpoint1", "endpoint2","overall"),
  Win   = Win,
  Loose = Loose,
  Tie   = Tie,
  WR = round(Win/Loose,5),
  WO = round((Win+0.5*Tie)/(Loose+0.5*Tie),5),
  GPC = round((Win-Loose)/(Win+Loose+Tie),5)
)

df4 = data.frame(row.names = c("T", "C"), c(censure_rateT1, censure_rateC1), c(censure_rateT2, censure_rateC2))
colnames(df4) = c("endpoint 1", "endpoint2")

list(Count = df3, value_tte_cont_C = dfC, value_tte_cont_T = dfT, censure = df4, p_val_GPC = paste("probabilité d'avoir des p-valeur < 0.05 pour la GPC: ", nb_pval_GPC/n_sim), p_val_WR =  paste("probabilité d'avoir des p-valeur < 0.05 pour le WR: ", nb_pval_WR/n_sim), p_val_WO = paste("probabilité d'avoir des p-valeur < 0.05 pour le WO: ", nb_pval_WO/n_sim))

vlines = data.frame(
  method = c("GPC", "WR", "WO"),
  intercept = c(0, 1, 1),
  linetype_label = "Hypothèse H0")


val1=unlist(results_df1[,1])
val2=unlist(results_df1[,2])
val3=unlist(results_df1[,3])
val=as.data.frame(unname(cbind(val1,val2,val3)))
colnames(val)=c("val_GPC","val_WR","val_WO")
summ_val3_1=summary(val)

values_long = val %>%
  select(starts_with("val_")) %>%
  pivot_longer(cols = everything(), names_to = "method", values_to = "value") %>%
  mutate(method = recode(method,
                         "val_GPC" = "GPC",
                         "val_WR"  = "WR",
                         "val_WO"  = "WO"))
values_long$value = unlist(values_long$value)
values_long <- values_long %>% filter(!is.infinite(value))


ggplot(values_long, aes(x = value, fill = method)) +
  geom_density(alpha = 0.6, color = "black") +
  geom_vline(data = vlines, aes(xintercept = intercept, linetype = linetype_label),
             color = "black", lwd = 1) +
  theme_minimal() +
  labs(title = "Distribution des statistiques de test",
       x = "Valeur", y = "Densité") +
  scale_linetype_manual(name = "", values = c("Hypothèse H0" = "dashed")) +
  scale_fill_manual(values = c("orange",  "purple", "cyan")) +
  annotate("text", x = 4.5, y = 1.6, label = paste("Étendue GPC : ", étendue_GPC), color = "orange", hjust = 1) +
  annotate("text", x = 4.5, y = 2.35, label = paste("Étendue WR : ", étendue_WR), color = "cyan", hjust = 1) +
  annotate("text", x = 4.5, y = 3.1, label = paste("Étendue WO : ", étendue_WO), color = "purple", hjust = 1) +
  coord_cartesian(xlim = c(-1, 5))
