id=1:(2*n) 
arm=rep(c("T","C"), each=200)

U = runif(2*n) 

lambdaT = 0.1
kT = 0.1


Z = ifelse(arm == "T", 1, 0)
beta = 0.3

prob_T = 0.65
prob_C = 0.35

mean_T = 3
mean_C = 2
sd_T = 2
sd_C = 2


Time_1 = round((-log(1 - U)) / (lambdaT * exp(beta * Z)^(1 / kT)),3)
fup_censureT = round(rweibull(n, shape = 1, scale =1),3)

Time_T = pmin(Time_1, fup_censureT)
deltaT=as.numeric(fup_censureT==Time_T)

Y_2_T=as.numeric(gener_binom(prob = prob_T)[,1] )-1
Y_2_C=as.numeric(gener_binom(prob = prob_C)[,1] )-1

Y_3_T= gener_continue(mean = mean_T, sd = sd_T)
Y_3_C= gener_continue(mean = mean_C, sd = sd_C)

stratum=sample(rep(c(1,3,5,8), each = 50))
dataT=data.frame(Y_1_trt = Time_T[Z==1], Delta_1_trt = deltaT[Z==1], Y_2_trt = Y_2_T, Y_3_trt = Y_3_T, stratum = stratum)
dataC=data.frame(Y_1_con = Time_T[Z==0], Delta_1_con = deltaT[Z==0], Y_2_con = Y_2_C, Y_3_con = Y_3_C)

trt_con=cbind(pid_trt = 1:200, pid_con = 1:200, dataT, dataC)
colnames(trt_con)=c("pid_trt", "pid_con", "Y_1_trt", "Delta_1_trt", "Y_2_trt", "Y_3_trt", "stratum", "Y_1_con", "Delta_1_con", "Y_2_con", "Y_3_con")
View(trt_con)
priority=c(1,2,3)
tau = c(0,0,0)
np_direction = c("smaller","smaller","smaller")

win.strategy.default<-function(trt_con, priority, tau, np_direction){
  n_ep = length(priority)

  #### Obtain the indicator of the first endpoint for treatment and control
  colname.trt_con = colnames(trt_con)

  ind.delta1.trt = which(colname.trt_con=="Delta_1_trt")
  ind.delta1.con = which(colname.trt_con=="Delta_1_con")

  ind.time1.trt = which(colname.trt_con=="Y_1_trt")
  ind.time1.con = which(colname.trt_con=="Y_1_con")


  win.status0 = NULL

  #### For TTE outcome: Denote the observed survival time as Y_trt and Y_con, and event status
  #### as Delta_trt and Delta_con. There is a win for the treatment group if we have:
  ####                      Delta_con = 1 and Y_trt > Y_con + tau_l,
  #### where tau_l denote the magnitude of difference to determine win/loss/tie.

  #### For continuous outcome: Denote the observed value as Y_trt and Y_con. The event status
  #### are Delta_trt = 1 and Delta_con = 1. There is a win for the treatment group if we have:
  ####                      Delta_con = 1 and Y_trt > Y_con + tau_l,
  #### where tau_l denote the magnitude of difference to determine win/loss/tie.

  #### For binary outcome (0/1): Denote the observed value as Y_trt and Y_con. The event status
  #### are Delta_trt = 1 and Delta_con = 1. There is a win for the treatment group if we have:
  ####                      Delta_con = 1 and Y_trt > Y_con + tau_l,
  #### where tau_l is 0.
  for(l in priority){
    delta_l_trt = trt_con[,ind.delta1.trt+l-1]; delta_l_con = trt_con[,ind.delta1.con+l-1]
    Y_l_trt = trt_con[,ind.time1.trt+l-1]; Y_l_con = trt_con[,ind.time1.con+l-1]

    tau_l = tau[l]
    direction_l = np_direction[l]

    if(direction_l == "larger"){
      win.temp1 = ((delta_l_trt == 1 & delta_l_con == 1 & Y_l_trt > (Y_l_con + tau_l)) |
                     (delta_l_trt < 1 & delta_l_con == 1 & Y_l_trt > (Y_l_con + tau_l)))
      win.temp2 = ((delta_l_trt == 1  & delta_l_con == 1 & Y_l_con > (Y_l_trt + tau_l)) |
                     (delta_l_trt == 1 & delta_l_con < 1 & Y_l_con > (Y_l_trt + tau_l)))
    }else if(direction_l == "smaller"){
      win.temp1 = ((delta_l_trt == 1 & delta_l_con == 1 & Y_l_trt < (Y_l_con - tau_l)) |
                     (delta_l_trt < 1 & delta_l_con == 1 & Y_l_trt < (Y_l_con - tau_l)))
      win.temp2 = ((delta_l_trt == 1 & delta_l_con == 1 & Y_l_con < (Y_l_trt - tau_l)) |
                     (delta_l_trt == 1 & delta_l_con < 1 & Y_l_con < (Y_l_trt - tau_l)))
    }

    win.status0 = cbind(win.status0, win.temp1, win.temp2)
  }

  #### prioritize: once a winner is determined, then all the subsequent is set to zero
  win_status = t(apply(win.status0, 1, func<-function(x){
    if(sum(x)>1){
      temp = x; temp[min((2*min(ceiling(which(x==1)/2))+1),(2*n_ep-1)):(2*n_ep)] = 0
      return(temp)
    }else{
      return(as.numeric(x))
    }
  }))

  colnames(win_status) = paste0(rep(c("Trt","Con"),n_ep),"_Endpoint",rep(priority,each=2))
  win_status = as.data.frame(win_status)

  return(win_status)
}
win.strategy.default(trt_con, tau, priority, np_direction)
