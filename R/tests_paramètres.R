#============================================================================
#---------------------------Modèle de Cox--------------------------
#============================================================================
n <- 200
arm <- sample(rep(c("T","C"), each=n))
Z <- ifelse(arm == "T", 1, 0)
U <- runif(2*n)

# Grilles de paramètres à tester
beta_vec   <- c(5, 4, 3, 2)
lambda_vec <- c(0.1, 0.5, 1, 2)
k_vec      <- c(0.1, 0.5, 1, 2)


prob_T <- 0.7
prob_C <- 0.3
mean_T <- 3
mean_C <- 1.3
sd_T <- 2
sd_C <- 1
stratum <- sample(rep(c(1,3,5,8), each = 50))

for (beta in beta_vec) {
  for (lambdaT in lambda_vec) {
    for (kT in k_vec) {

        # 1. Simuler les temps à l'événement (Weibull pour Cox)
        event_times <- (-log(1-U) / (lambdaT * exp(beta * Z)))^(1/kT)
        
        # 2. Simuler les temps de censure avec shape variable
        fup_censureT <- round(rweibull(2*n,1 , scale=5), 3)
        
        # 3. Temps observé et statut
        obs_times <- pmin(event_times, fup_censureT)
        status <- as.numeric(event_times <= fup_censureT)
        
        # 4. Variables supplémentaires
        Y_2_T <- rbinom(n, 1, prob_T)
        Y_2_C <- rbinom(n, 1, prob_C)
        Y_3_T <- rnorm(n, mean = mean_T, sd = sd_T)
        Y_3_C <- rnorm(n, mean = mean_C, sd = sd_C)
        
        # 5. Dataframes par groupe
        dataT <- data.frame(
          Y_1 = obs_times[Z==1],
          Delta_1 = status[Z==1],
          Y_2 = Y_2_T,
          Y_3 = Y_3_T,
          stratum = stratum
        )
        dataC <- data.frame(
          Y_1 = obs_times[Z==0],
          Delta_1 = status[Z==0],
          Y_2 = Y_2_C,
          Y_3 = Y_3_C,
          stratum = stratum
        )
        
        # 6. Résumés
        cat("\n==========================\n")
        cat("beta =", beta, "| lambda =", lambdaT, "| k =", kT, "\n")
        cat("Taux de censure groupe T :", round(sum(dataT$Delta_1==0)/n,3), "\n")
        cat("Taux de censure groupe C :", round(sum(dataC$Delta_1==0)/n,3), "\n")
        cat("Summary groupe T:\n")
        print(summary(dataT))
        cat("Summary groupe C:\n")
        print(summary(dataC))
        cat("==========================\n")
      }
    }
  
}


#============================================================================
#---------------------------Modèle AFT--------------------------
#============================================================================
set.seed(123)
n <- 200
arm <- sample(rep(c("T","C"), each=n))
Z <- ifelse(arm == "T", 1, 0)
U <- runif(2*n)

# Grilles de paramètres à tester
beta_vec   <- c(2,3)
lambda_vec <- c(1,2,0.1, 0.5, 0.05, 0.01 )
k_vec      <- c( 0.5,0.01,1, 2, 3)

# Temps de censure commun à tous
fup_censureT <- round(rweibull(2*n, shape=2, scale=6), 3)

# Boucle sur toutes les combinaisons
for (beta in beta_vec) {
  for (lambdaT in lambda_vec) {
    for (kT in k_vec) {
      # Simulation des temps à l'événement (modèle AFT Weibull)
      Time_1_AFT <- round((((1/(1-U)-1)*(1/lambdaT))^(1/kT)*exp(Z*beta)), 3)
      Time_T <- pmin(Time_1_AFT, fup_censureT)
      deltaT <- as.numeric(fup_censureT == Time_T)
      
      # Résumés par groupe
      summary_T <- summary(Time_1_AFT[Z==1])
      summary_C <- summary(Time_1_AFT[Z==0])
      #censor_T <- sum(deltaT[Z==1]==0)/n
      #censor_C <- sum(deltaT[Z==0]==0)/n
      
      cat("\n==========================\n")
      cat("beta =", beta, "| lambda =", lambdaT, "| k =", kT, "\n")
      #cat("Taux de censure groupe T :", round(censor_T,3), "\n")
      #cat("Taux de censure groupe C :", round(censor_C,3), "\n")
      cat("Summary groupe T:\n")
      print(summary_T)
      cat("Summary groupe C:\n")
      print(summary_C)
      cat("==========================\n")
    }
  }
}



