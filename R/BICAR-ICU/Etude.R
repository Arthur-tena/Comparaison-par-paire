setwd("~/Comparaison-par-paire/R/BICAR-ICU")
library(WINS)
data=readxl::read_xlsx("Export-20250617.xlsx")

data_epur=data[,c(1, 3, 6, 7,8, 62, 69, 76, 78)]
colnames(data_epur)=c("id","arm", "SEPSIS_ON", "AGE_RANDO", "AKIN", "Décès", "Dialyse", "Vaso_free_day", "Ventil_free_day")



#===============================================================================================================
#                    Analyse non-stratifié et sans seuil
#===============================================================================================================



New_data = data.frame(
   id = data_epur$id,
  arm = data_epur$arm,
  death = ifelse(data_epur$Décès == "Non", 0, 1),
  dialysis = ifelse(data_epur$Dialyse == "Non", 0, 1),
  ventil_free_day = data_epur$Ventil_free_day,
  vaso_free_day = data_epur$Vaso_free_day
)


New_data = New_data[!is.na(New_data$vaso_free_day), ]
New_data$id=1:nrow(New_data)
colnames(New_data)=c("id", "arm", "Y_1", "Y_2", "Y_3", "Y_4")


result2 =
  win.stat(
    data = New_data,
    ep_type = c("binary", "binary", "continuous",  "continuous" ),
    stratum.weight = "unstratified",
    tau = c(0,0,0,0),
    arm.name = c("BRAS B", "BRAS A"),
    alpha = 0.05,
    digit = 3,
    pvalue = "two-sided",
    priority = c(1,2,3,4),
    np_direction = c("smaller","smaller","larger","larger"),
    summary.print = FALSE
  )
result2$p_value
result2$

Win=c(result2$summary_ep$Trt_Endpoint1$Count,
      result2$summary_ep$Trt_Endpoint2$Count,
      result2$summary_ep$Trt_Endpoint3$Count,
      result2$summary_ep$Trt_Endpoint4$Count,
      result2$summary_ep$Trt_Endpoint1$Count+result2$summary_ep$Trt_Endpoint2$Count+result2$summary_ep$Trt_Endpoint3$Count+result2$summary_ep$Trt_Endpoint4$Count)

Loose=c(result2$summary_ep$Con_Endpoint1$Count,
        result2$summary_ep$Con_Endpoint2$Count,
        result2$summary_ep$Con_Endpoint3$Count,
        result2$summary_ep$Con_Endpoint4$Count,
        result2$summary_ep$Con_Endpoint1$Count+result2$summary_ep$Con_Endpoint2$Count+result2$summary_ep$Con_Endpoint3$Count+result2$summary_ep$Con_Endpoint4$Count)

tie_edp1=36666-(result2$summary_ep$Trt_Endpoint1$Count+result2$summary_ep$Con_Endpoint1$Count)
tie_edp2 = tie_edp1-(result2$summary_ep$Trt_Endpoint2$Count+result2$summary_ep$Con_Endpoint2$Count)
tie_edp3 = tie_edp2-(result2$summary_ep$Trt_Endpoint3$Count+result2$summary_ep$Con_Endpoint3$Count)
tie_edp4 = tie_edp3-(result2$summary_ep$Trt_Endpoint4$Count+result2$summary_ep$Con_Endpoint4$Count)

Tie=c(tie_edp1, tie_edp2,tie_edp3, tie_edp4, tie_edp4)

df_count = data.frame(
  row.names = c("Décès", "Dialyse","ventil_free_day", "vaso_free_day","overall"),
  Win   = Win,
  Loose = Loose,
  Tie   = Tie
)
df_count

df_pval=data.frame(result2$p_value)
colnames(df_pval)=c("Win Ratio","Net Benefit","Win Odds")
df_pval

df_vals=data.frame(row.names = c("Win Ratio", "Net Benefit", "Win Odds"),
                  CI_lower=c(result2$Win_statistic$Win_Ratio[2],result2$Win_statistic$Net_Benefit[2],result2$Win_statistic$Win_Odds[2]),
                  Valeur=c(result2$Win_statistic$Win_Ratio[1],result2$Win_statistic$Net_Benefit[1],result2$Win_statistic$Win_Odds[1]),
                  CI_upper=c(result2$Win_statistic$Win_Ratio[3],result2$Win_statistic$Net_Benefit[3],result2$Win_statistic$Win_Odds[3]))
df_vals



df_global = data.frame(
  Métrique = c("Win Ratio", "Net Benefit", "Win Odds"),
  Valeur = c(result2$Win_statistic$Win_Ratio[1],
             result2$Win_statistic$Net_Benefit[1],
             result2$Win_statistic$Win_Odds[1]),
  CI_lower = c(result2$Win_statistic$Win_Ratio[2],
               result2$Win_statistic$Net_Benefit[2],
               result2$Win_statistic$Win_Odds[2]),
  CI_upper = c(result2$Win_statistic$Win_Ratio[3],
               result2$Win_statistic$Net_Benefit[3],
               result2$Win_statistic$Win_Odds[3]),
  z_score = unname(result2$z_statistic[c("zstat_WR", "zstat_NB", "zstat_WO")]),
  p_val = unname(result2$p_value[c("pvalue_WR", "pvalue_NB", "pvalue_WO")])
)
df_global

#===============================================================================================================
#                    Analyse stratifié sur l'âge et sans seuil
#===============================================================================================================

data_age = data.frame(
  id = data_epur$id,
  arm = data_epur$arm,
  death = ifelse(data_epur$Décès == "Non", 0, 1),
  dialysis = ifelse(data_epur$Dialyse == "Non", 0, 1),
  ventil_free_day = data_epur$Ventil_free_day,
  vaso_free_day = data_epur$Vaso_free_day,
  stratum=ifelse(data_epur$AGE_RANDO=="< 65 ans", 0, 1)
)


data_age = data_age[!is.na(data_age$vaso_free_day), ]
data_age$id=1:nrow(data_age)
colnames(data_age)=c("id", "arm", "Y_1", "Y_2", "Y_3", "Y_4","stratum")

result3 =
  win.stat(
    data = data_age,
    ep_type = c("binary", "binary", "continuous",  "continuous" ),
    stratum.weight = "equal",
    tau = c(0,0,0,0),
    arm.name = c("BRAS B", "BRAS A"),
    alpha = 0.05,
    digit = 3,
    pvalue = "two-sided",
    priority = c(1,2,3,4),
    np_direction = c("smaller","smaller","larger","larger"),
    summary.print = FALSE
  )


df_vals=data.frame(row.names = c("Win Ratio", "Net Benefit", "Win Odds"),
                   CI_lower=c(result3$Win_statistic$Win_Ratio[2],result3$Win_statistic$Net_Benefit[2],result3$Win_statistic$Win_Odds[2]),
                   Valeur=c(result3$Win_statistic$Win_Ratio[1],result3$Win_statistic$Net_Benefit[1],result3$Win_statistic$Win_Odds[1]),
                   CI_upper=c(result3$Win_statistic$Win_Ratio[3],result3$Win_statistic$Net_Benefit[3],result3$Win_statistic$Win_Odds[3]))
df_vals


df_global = data.frame(
  Métrique = c("Win Ratio", "Net Benefit", "Win Odds"),
  Valeur = c(result3$Win_statistic$Win_Ratio[1],
             result3$Win_statistic$Net_Benefit[1],
             result3$Win_statistic$Win_Odds[1]),
  CI_lower = c(result3$Win_statistic$Win_Ratio[2],
               result3$Win_statistic$Net_Benefit[2],
               result3$Win_statistic$Win_Odds[2]),
  CI_upper = c(result3$Win_statistic$Win_Ratio[3],
               result3$Win_statistic$Net_Benefit[3],
               result3$Win_statistic$Win_Odds[3]),
  z_score = unname(result3$z_statistic[c("zstat_WR", "zstat_NB", "zstat_WO")]),
  p_val = unname(result3$p_value[c("pvalue_WR", "pvalue_NB", "pvalue_WO")])
)
df_global


#===============================================================================================================
#                    Analyse stratifié sur le septis et sans seuil
#===============================================================================================================

data_septis = data.frame(
  id = data_epur$id,
  arm = data_epur$arm,
  death = ifelse(data_epur$Décès == "Non", 0, 1),
  dialysis = ifelse(data_epur$Dialyse == "Non", 0, 1),
  ventil_free_day = data_epur$Ventil_free_day,
  vaso_free_day = data_epur$Vaso_free_day,
  stratum=ifelse(data_epur$SEPSIS_ON=="Oui", 0, 1)
)


data_septis = data_septis[!is.na(data_septis$vaso_free_day), ]
data_septis$id=1:nrow(data_septis)
colnames(data_septis)=c("id", "arm", "Y_1", "Y_2", "Y_3", "Y_4","stratum")

result4 =
  win.stat(
    data = data_septis,
    ep_type = c("binary", "binary", "continuous",  "continuous" ),
    stratum.weight = "equal",
    tau = c(0,0,0,0),
    arm.name = c("BRAS B", "BRAS A"),
    alpha = 0.05,
    digit = 3,
    pvalue = "two-sided",
    priority = c(1,2,3,4),
    np_direction = c("smaller","smaller","larger","larger"),
    summary.print = FALSE
  )


df_vals=data.frame(row.names = c("Win Ratio", "Net Benefit", "Win Odds"),
                   CI_lower=c(result4$Win_statistic$Win_Ratio[2],result4$Win_statistic$Net_Benefit[2],result4$Win_statistic$Win_Odds[2]),
                   Valeur=c(result4$Win_statistic$Win_Ratio[1],result4$Win_statistic$Net_Benefit[1],result4$Win_statistic$Win_Odds[1]),
                   CI_upper=c(result4$Win_statistic$Win_Ratio[3],result4$Win_statistic$Net_Benefit[3],result4$Win_statistic$Win_Odds[3]))
df_vals


df_global = data.frame(
  Métrique = c("Win Ratio", "Net Benefit", "Win Odds"),
  Valeur = c(result4$Win_statistic$Win_Ratio[1],
             result4$Win_statistic$Net_Benefit[1],
             result4$Win_statistic$Win_Odds[1]),
  CI_lower = c(result4$Win_statistic$Win_Ratio[2],
               result4$Win_statistic$Net_Benefit[2],
               result4$Win_statistic$Win_Odds[2]),
  CI_upper = c(result4$Win_statistic$Win_Ratio[3],
               result4$Win_statistic$Net_Benefit[3],
               result4$Win_statistic$Win_Odds[3]),
  z_score = unname(result4$z_statistic[c("zstat_WR", "zstat_NB", "zstat_WO")]),
  p_val = unname(result4$p_value[c("pvalue_WR", "pvalue_NB", "pvalue_WO")])
)
df_global


#===============================================================================================================
#                    Analyse stratifié sur le score AKIN et sans seuil
#===============================================================================================================

data_akin = data.frame(
  id = data_epur$id,
  arm = data_epur$arm,
  death = ifelse(data_epur$Décès == "Non", 0, 1),
  dialysis = ifelse(data_epur$Dialyse == "Non", 0, 1),
  ventil_free_day = data_epur$Ventil_free_day,
  vaso_free_day = data_epur$Vaso_free_day,
  stratum=ifelse(data_epur$AKIN=="0-1", 0, 1)
)


data_akin = data_akin[!is.na(data_akin$vaso_free_day), ]
data_akin$id=1:nrow(data_akin)
colnames(data_akin)=c("id", "arm", "Y_1", "Y_2", "Y_3", "Y_4","stratum")

result5 =
  win.stat(
    data = data_akin,
    ep_type = c("binary", "binary", "continuous",  "continuous" ),
    stratum.weight = "equal",
    tau = c(0,0,0,0),
    arm.name = c("BRAS B", "BRAS A"),
    alpha = 0.05,
    digit = 3,
    pvalue = "two-sided",
    priority = c(1,2,3,4),
    np_direction = c("smaller","smaller","larger","larger"),
    summary.print = FALSE
  )


df_vals=data.frame(row.names = c("Win Ratio", "Net Benefit", "Win Odds"),
                   CI_lower=c(result5$Win_statistic$Win_Ratio[2],result5$Win_statistic$Net_Benefit[2],result5$Win_statistic$Win_Odds[2]),
                   Valeur=c(result5$Win_statistic$Win_Ratio[1],result5$Win_statistic$Net_Benefit[1],result5$Win_statistic$Win_Odds[1]),
                   CI_upper=c(result5$Win_statistic$Win_Ratio[3],result5$Win_statistic$Net_Benefit[3],result5$Win_statistic$Win_Odds[3]))
df_vals


df_global = data.frame(
  Métrique = c("Win Ratio", "Net Benefit", "Win Odds"),
  Valeur = c(result5$Win_statistic$Win_Ratio[1],
             result5$Win_statistic$Net_Benefit[1],
             result5$Win_statistic$Win_Odds[1]),
  CI_lower = c(result5$Win_statistic$Win_Ratio[2],
               result5$Win_statistic$Net_Benefit[2],
               result5$Win_statistic$Win_Odds[2]),
  CI_upper = c(result5$Win_statistic$Win_Ratio[3],
               result5$Win_statistic$Net_Benefit[3],
               result5$Win_statistic$Win_Odds[3]),
  z_score = unname(result5$z_statistic[c("zstat_WR", "zstat_NB", "zstat_WO")]),
  p_val = unname(result5$p_value[c("pvalue_WR", "pvalue_NB", "pvalue_WO")])
)
df_global


#===============================================================================================================
#                    Analyse stratifié sur le score AKIN + l'âge et sans seuil
#===============================================================================================================

data_akin_age = data.frame(
  id = data_epur$id,
  arm = data_epur$arm,
  death = ifelse(data_epur$Décès == "Non", 0, 1),
  dialysis = ifelse(data_epur$Dialyse == "Non", 0, 1),
  ventil_free_day = data_epur$Ventil_free_day,
  vaso_free_day = data_epur$Vaso_free_day,
  stratum=ifelse(data_epur$AKIN=="0-1" & data_epur$AGE_RANDO=="< 65 ans", 0,
                 ifelse(data_epur$AKIN=="2-3" & data_epur$AGE_RANDO=="< 65 ans",1,
                        ifelse(data_epur$AKIN=="0-1" & data_epur$AGE_RANDO==">= 65 ans",2,3)))
)


data_akin_age = data_akin_age[!is.na(data_akin_age$vaso_free_day), ]
data_akin_age$id=1:nrow(data_akin_age)
colnames(data_akin_age)=c("id", "arm", "Y_1", "Y_2", "Y_3", "Y_4","stratum")

result6 =
  win.stat(
    data = data_akin_age,
    ep_type = c("binary", "binary", "continuous",  "continuous" ),
    stratum.weight = "equal",
    tau = c(0,0,0,0),
    arm.name = c("BRAS B", "BRAS A"),
    alpha = 0.05,
    digit = 3,
    pvalue = "two-sided",
    priority = c(1,2,3,4),
    np_direction = c("smaller","smaller","larger","larger"),
    summary.print = FALSE
  )


df_vals=data.frame(row.names = c("Win Ratio", "Net Benefit", "Win Odds"),
                   CI_lower=c(result6$Win_statistic$Win_Ratio[2],result6$Win_statistic$Net_Benefit[2],result6$Win_statistic$Win_Odds[2]),
                   Valeur=c(result6$Win_statistic$Win_Ratio[1],result6$Win_statistic$Net_Benefit[1],result6$Win_statistic$Win_Odds[1]),
                   CI_upper=c(result6$Win_statistic$Win_Ratio[3],result6$Win_statistic$Net_Benefit[3],result6$Win_statistic$Win_Odds[3]))
df_vals


df_global = data.frame(
  Métrique = c("Win Ratio", "Net Benefit", "Win Odds"),
  Valeur = c(result6$Win_statistic$Win_Ratio[1],
             result6$Win_statistic$Net_Benefit[1],
             result6$Win_statistic$Win_Odds[1]),
  CI_lower = c(result6$Win_statistic$Win_Ratio[2],
               result6$Win_statistic$Net_Benefit[2],
               result6$Win_statistic$Win_Odds[2]),
  CI_upper = c(result6$Win_statistic$Win_Ratio[3],
               result6$Win_statistic$Net_Benefit[3],
               result6$Win_statistic$Win_Odds[3]),
  z_score = unname(result6$z_statistic[c("zstat_WR", "zstat_NB", "zstat_WO")]),
  p_val = unname(result6$p_value[c("pvalue_WR", "pvalue_NB", "pvalue_WO")])
)
df_global
