setwd("C:/Users/MASTERS TEMA Arthur/Documents/Comparaison-par-paire/R/BICAR-ICU")
data=readxl::read_xlsx("Export-20250617.xlsx")

data_epur=data[,c(1, 3, 6, 7,8, 62, 69, 76, 78)]




#===============================================================================================================
#                    Analyse non-stratifié et sans seuil
#===============================================================================================================

colnames(data_epur)=c("id","arm", "SEPSIS_ON", "AGE_RANDO", "AKIN", "Décès", "Dialyse", "Vaso_free_day", "Ventil_free_day")

New_data=data.frame(id=data_epur[,1], arm = data_epur[,2], Y_1= data_epur[,6],Y_2=data_epur[,7], Y_3 = data_epur[,9], Y_4 = data_epur[,8] )
New_data[,3]= as.numeric(ifelse("Non", 0, 1))
New_data[,4]= as.numeric(ifelse("Non", 0, 1))
str(New_data)


result =
  win.stat(
    data = New_data,
    ep_type = c("binary", "binary", "continuous",  "continuous" ),
    stratum.weight = "unstratified",
    tau = c(0,0,0,0),
    arm.name = c("Bras B", "Bras A"),
    alpha = 0.05,
    digit = 3,
    pvalue = "two-sided",
    priority = c(1,2,3,4),
    np_direction = c("smaller","smaller","larger","larger"),
    summary.print = TRUE
  )
