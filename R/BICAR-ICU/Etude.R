setwd("C:/Users/MASTERS TEMA Arthur/Documents/Comparaison-par-paire/R/BICAR-ICU")
library(WINS)
data=readxl::read_xlsx("Export-20250617.xlsx")

data_epur=data[,c(1, 3, 6, 7,8, 62, 69, 76, 78)]
colnames(data_epur)=c("id","arm", "SEPSIS_ON", "AGE_RANDO", "AKIN", "Décès", "Dialyse", "Vaso_free_day", "Ventil_free_day")



#===============================================================================================================
#                    Analyse non-stratifié et sans seuil
#===============================================================================================================



New_data <- data.frame(
   id = data_epur$id,
  arm = data_epur$arm,
  death = ifelse(data_epur$Décès == "Non", 0, 1),
  dialysis = ifelse(data_epur$Dialyse == "Non", 0, 1),
  ventil_free_day = data_epur$Ventil_free_day,
  vaso_free_day = data_epur$Vaso_free_day
)


New_data <- New_data[!is.na(New_data$vaso_free_day), ]
New_data$id=1:nrow(New_data)
colnames(New_data)=c("id", "arm", "Y_1", "Y_2", "Y_3", "Y_4")

names(New_data)
str(New_data)
result =
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
    summary.print = TRUE
  )

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
result2$Win_statistic

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
Tie=c()

df3 = data.frame(
  row.names = c("endpoint1", "endpoint2","endpoint3", "endpoint4","overall"),
  Win   = Win,
  Loose = Loose,
  Tie   = Tie,
  WR = Win/Loose,
  WO = (Win+0.5*Tie)/(Loose+0.5*Tie),
  GPC = (Win-Loose)/(Win+Loose+Tie)
)
c(result2$summary_ep$Trt_Endpoint1, result2$summary_ep$Con_Endpoint1)
c(result2$summary_ep$Trt_Endpoint2, result2$summary_ep$Con_Endpoint2)
c(result2$summary_ep$Trt_Endpoint3, result2$summary_ep$Con_Endpoint3)
c(result2$summary_ep$Trt_Endpoint4, result2$summary_ep$Con_Endpoint4)
