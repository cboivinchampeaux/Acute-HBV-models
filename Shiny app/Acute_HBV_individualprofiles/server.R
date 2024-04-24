#Server code for modeling via differential equations for oral and iv PK

library(deSolve)
library(shiny)
library(mlxR)
library(shinydashboard)
library(ggplot2)
library(scales)

#Set the directory where you saved the observations files
patient1_obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_2/Dataset/patient1-observationV.csv")
patient2_obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_2/Dataset/patient2-observationV.csv")
patient3_obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_2/Dataset/patient3-observationV.csv")
# patient4_obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_2/Dataset/patient4-observationV.csv")
patient5_obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_2/Dataset/patient5-observationV.csv")
patient6_obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_2/Dataset/patient6-observationV.csv")
patient7_obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_2/Dataset/patient7-observationV.csv")

shinyServer(function(input, output) {
  output$plot1 <- renderPlot({
    
    parms1_2I <- c(
      r = 0.01,
      tmax = 13.6 * 10 ** 6,
      beta = 0.7 * 10 ** (-10),
      rho1 = 0.071,
      delta_tr = 0.81,
      mu = 5.8 * 10 ** (-4),
      rho2 = 0.18,
      p1 = 34,
      p2 = 376,
      c = 0.67,
      z = 10,
      alpha = 7.3 * 10 ** (-7),
      tau = 15.1,
      de = 0.5
    )
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = 13.6 * 10 ** 6, I1 = 0, I2 = 0, V = 0.33, E = 20)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau, I1, lagvalue(t - tau, 2))
        lagI2 <- ifelse(t < tau, I2, lagvalue(t - tau, 3))
        lagE <- ifelse(t < tau, E, lagvalue(t - tau, 5))
        
        dT <- r * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta * T * V + rho1 * I1
        dI1 <- beta * T * V - (rho1 + delta_tr) * I1 - mu * E * I1 + rho2 * I2
        dI2 <- r * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu * E * I2
        dV <- p1 * I1 + p2 * I2 - c * V
        dE <- z + alpha * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = parms1_2I)
    
    parms1_R <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 0.653 * 10 ** (-10),
      rho = 0.729 * 10 ** (-3),
      mu = 72.0 * 10 ** (-4),
      mu_1 = 150.0 * 10 ** (-6),
      rho_r = 2.0 * 10 ** (-5),
      p = 400,
      c = 0.67,
      z = 10,
      alpha = 2.2 * 10 ** (-7),
      tau = 15.2,
      de = 0.5
    )
    
    yinitR <- c(T = 13.6 * 10 ** 6, I = 0, R = 0, V = 0.33, E = 20)
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,5)
        
        dT <- r*T*(1-((T+I+R)/tmax)) - beta*T*V + rho_r*R
        dI <- r*I*(1-((T+I+R)/tmax)) + beta*T*V - rho*I*E - mu*I*E 
        dR <- rho*I*E + r*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_1*R*E
        dV <- p*I - c*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = parms1_R)
    
    parms1_A <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 0.953 * 10 ** (-10),
      delta = 0.1538,
      km = 10,
      kp = 10 ** (-12),
      cxav = 2.7,
      pa = 0.1 * 10 ** (-5),
      theta = 1000,
      ra = 0.4242,
      amax = 4 * 10 ** 15,
      da = 0.033,
      p = 297,
      c = 0.67
    )
    
    yinitA <- c(T = 13.6 * 10 ** 6, I = 0, XAV = 0, A = 0, V = 0.33)
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- beta*T*V-delta*I
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*(1+theta)+ra*A*(1-A/amax)+km*XAV*(1+theta)-kp*A*V*(1+theta)-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = parms1_A)
    
    parms1_AE <- c(
      r = 1.0,
      tmax = 1.36e7,
      beta = 2.6e-10,
      delta = 0.1,
      mu = 2e-3,
      km = 10,
      kp = 1e-10,
      c = 1.67,
      cxav = 4*1.67, #4*c
      pa = 1e-3,
      theta = 1e3,
      ra = 0.35,
      amax = 4e11,
      da = 0.033,
      p = 266,
      z= 30,
      alpha = 2.2e-7,
      tau = 15.2,
      de = 0.5
    )
    
    yinitAE <- c(T = 1.36e7, I = 0, XAV = 0, A = 0, V = 0.33, E = 20)
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,6)
        
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- r*I*(1-(T+I)/tmax)+beta*T*V-delta*I-mu*I*E
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*theta+ra*A*(1-A/amax)+km*XAV*theta-kp*A*V*theta-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = parms1_AE)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               V= output1[,5])
    finaloutput2 <- data.frame(timesR= output2[,1],
                               V= output2[,5])
    finaloutput3 <- data.frame(timesA= output3[,1],
                               V= output3[,6])
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               V= output4[,6])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9)
    colormanual <- c("Two infected cells populations" = "darkolivegreen", 
                     "Refractory cells population"="darkolivegreen4",
                     "Antibody population"="darkolivegreen3",
                     "Antibody and Effector cells populations" = "darkolivegreen2")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = V,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075)

      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = V, color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075)

      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = V, color = "Antibody population",
                                               linetype="Antibody population"),size=1.075)
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = V, color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=2)
      
      vdisp[4] <- TRUE}
    pl <- pl + geom_point(patient1_obsV,mapping = aes(x = time, y = V), shape=8, color = "darkgreen", size=4)
    pl <- pl + labs(x='Time (Days)') + labs(y='HBV DNA (virions/mL)') +
      labs(color  = "Model", linetype = "Model")

    pl <- pl + theme(legend.position=c(0.1, 0.25), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(0.001, 10^11),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 

    print(pl)

  })
    
    
  
  output$plot2 <- renderPlot({

    w <- input$num5 #Duration of the Model (INPUT)
    
    parms2_2I = c(
      r = 0.01, 
      tmax = 13.6e6,
      beta = 6.7e-10, 
      rho1 = 0.064,
      mu = 9.9e-4, 
      rho2 = 0.06,
      delta_tr = 0.59,
      p1 = 75,
      p2 = 73,
      c = 0.67,
      z = 10,
      alpha = 3.6e-7,
      tau = 19.6,
      de = 0.5 
    )
    
    yinit2I <- c(T = 13.6 * 10 ** 6, I1 = 0, I2 = 0, V = 0.33, E = 20)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau, I1, lagvalue(t - tau, 2))
        lagI2 <- ifelse(t < tau, I2, lagvalue(t - tau, 3))
        lagE <- ifelse(t < tau, E, lagvalue(t - tau, 5))
        
        dT <- r * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta * T * V + rho1 * I1
        dI1 <- beta * T * V - (rho1 + delta_tr) * I1 - mu * E * I1 + rho2 * I2
        dI2 <- r * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu * E * I2
        dV <- p1 * I1 + p2 * I2 - c * V
        dE <- z + alpha * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = parms2_2I)
    
    parms2_R <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 7.81 * 10 ** (-10),
      rho = 1.06 * 10 ** (-3),
      mu = 1.5 * 10 ** (-4),
      mu_1 = 0.2 * 10 ** (-6),
      rho_r = 2.0 * 10 ** (-5),
      p = 57,
      c = 0.67,
      z = 10,
      alpha = 3.2 * 10 ** (-7),
      tau = 19.6,
      de = 0.5
    )
    
    yinitR <- c(T = 13.6 * 10 ** 6, I = 0, R = 0, V = 0.33, E = 20)
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,5)
        
        dT <- r*T*(1-((T+I+R)/tmax)) - beta*T*V + rho_r*R
        dI <- r*I*(1-((T+I+R)/tmax)) + beta*T*V - rho*I*E - mu*I*E 
        dR <- rho*I*E + r*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_1*R*E
        dV <- p*I - c*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = parms2_R)
    
    parms2_A <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 11.1 * 10 ** (-10),
      delta = 0.061,
      km = 10,
      kp = 10 ** (-12),
      cxav = 2.7,
      pa = 0.7 * 10 ** (-5),
      theta = 1000,
      ra = 0.4656,
      amax = 4 * 10 ** 15,
      da = 0.033,
      p = 43.5,
      c = 0.67
    )
    
    yinitA <- c(T = 13.6 * 10 ** 6, I = 0, XAV = 0, A = 0, V = 0.33)
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- beta*T*V-delta*I
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*(1+theta)+ra*A*(1-A/amax)+km*XAV*(1+theta)-kp*A*V*(1+theta)-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = parms2_A)
    
    parms2_AE <- c(
      r = 1.0,
      tmax = 1.36e7,
      beta = 9.24e-10,
      delta = 0.035,
      mu = 2.5e-3,
      km = 10,
      kp = 1e-10,
      cxav = 4*1.67,
      pa = 5.9e-3,
      theta = 1e3,
      ra = 0.29,
      amax = 4e11,
      da = 0.033,
      p = 117.1,
      c = 1.67,
      z= 30,
      alpha = 2.2e-7,
      tau = 19.6,
      de = 0.5
    )
    
    yinitAE <- c(T = 1.36e7, I = 0, XAV = 0, A = 0, V = 0.33, E = 20)
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,6)
        
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- r*I*(1-(T+I)/tmax)+beta*T*V-delta*I-mu*I*E
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*theta+ra*A*(1-A/amax)+km*XAV*theta-kp*A*V*theta-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = parms2_AE)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               V= output1[,5])
    finaloutput2 <- data.frame(timesR= output2[,1],
                               V= output2[,5])
    finaloutput3 <- data.frame(timesA= output3[,1],
                               V= output3[,6])
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               V= output4[,6])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9)
    colormanual <- c("Two infected cells populations" = "darkolivegreen", 
                     "Refractory cells population"="darkolivegreen4",
                     "Antibody population"="darkolivegreen3",
                     "Antibody and Effector cells populations" = "darkolivegreen2")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = V,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075)
      
      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = V, color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075)
      
      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = V, color = "Antibody population",
                                               linetype="Antibody population"),size=1.075)
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = V, color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=2)
      
      vdisp[4] <- TRUE}
    pl <- pl + geom_point(patient2_obsV,mapping = aes(x = time, y = V), shape=8, color = "darkgreen", size=4)
    pl <- pl + labs(x='Time (Days)') + labs(y='HBV DNA (virions/mL)') +
      labs(color  = "Model", linetype = "Model")
    
    pl <- pl + theme(legend.position=c(0.1, 0.25), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(0.001, 10^11),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 
    
    print(pl)
    
  })
  
  output$plot3 <- renderPlot({
    
    parms3_2I = c(
      r = 0.01, 
      tmax = 13.6e6,
      beta = 0.8e-10, 
      rho1 = 0.05,
      mu = 6.9e-4, 
      rho2 = 0.11,
      delta_tr = 0.98,
      p1 = 201,
      p2 = 428,
      c = 0.67,
      z = 10,
      alpha = 4.1e-7,
      tau = 24.9,
      de = 0.5 
    )
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = 13.6 * 10 ** 6, I1 = 0, I2 = 0, V = 0.33, E = 20)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau, I1, lagvalue(t - tau, 2))
        lagI2 <- ifelse(t < tau, I2, lagvalue(t - tau, 3))
        lagE <- ifelse(t < tau, E, lagvalue(t - tau, 5))
        
        dT <- r * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta * T * V + rho1 * I1
        dI1 <- beta * T * V - (rho1 + delta_tr) * I1 - mu * E * I1 + rho2 * I2
        dI2 <- r * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu * E * I2
        dV <- p1 * I1 + p2 * I2 - c * V
        dE <- z + alpha * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = parms3_2I)
    
    parms3_R <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 1.74 * 10 ** (-10),
      rho = 0.311 * 10 ** (-3),
      mu = 6.9 * 10 ** (-4),
      mu_1 = 37.0 * 10 ** (-6),
      rho_r = 2.0 * 10 ** (-5),
      p = 189,
      c = 0.67,
      z = 10,
      alpha = 2.3 * 10 ** (-7),
      tau = 29.0,
      de = 0.5
    )
    
    yinitR <- c(T = 13.6 * 10 ** 6, I = 0, R = 0, V = 0.33, E = 20)
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,5)
        
        dT <- r*T*(1-((T+I+R)/tmax)) - beta*T*V + rho_r*R
        dI <- r*I*(1-((T+I+R)/tmax)) + beta*T*V - rho*I*E - mu*I*E 
        dR <- rho*I*E + r*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_1*R*E
        dV <- p*I - c*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = parms3_R)
    
    parms3_A <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 18.2 * 10 ** (-10),
      delta = 0.018,
      km = 10,
      kp = 10 ** (-12),
      cxav = 2.7,
      pa = 9.0 * 10 ** (-5),
      theta = 1000,
      ra = 0.2672,
      amax = 4 * 10 ** 15,
      da = 0.033,
      p = 12.4,
      c = 0.67
    )
    
    yinitA <- c(T = 13.6 * 10 ** 6, I = 0, XAV = 0, A = 0, V = 0.33)
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- beta*T*V-delta*I
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*(1+theta)+ra*A*(1-A/amax)+km*XAV*(1+theta)-kp*A*V*(1+theta)-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = parms3_A)
    
    parms3_AE <- c(
      r = 1.0,
      tmax = 1.36e7,
      beta = 6.23e-10,
      delta = 0.071,
      mu = 0.6e-3,
      km = 10,
      kp = 1e-10,
      cxav = 4*1.67,
      pa = 13.8e-3,
      theta = 1e3,
      ra = 0.24,
      amax = 4e11,
      da = 0.033,
      p = 159.4,
      c = 1.67,
      z= 30,
      alpha = 2.2e-7,
      tau = 29,
      de = 0.5
    )
    
    yinitAE <- c(T = 1.36e7, I = 0, XAV = 0, A = 0, V = 0.33, E = 20)
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,6)
        
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- r*I*(1-(T+I)/tmax)+beta*T*V-delta*I-mu*I*E
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*theta+ra*A*(1-A/amax)+km*XAV*theta-kp*A*V*theta-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = parms3_AE)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               V= output1[,5])
    finaloutput2 <- data.frame(timesR= output2[,1],
                               V= output2[,5])
    finaloutput3 <- data.frame(timesA= output3[,1],
                               V= output3[,6])
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               V= output4[,6])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9)
    colormanual <- c("Two infected cells populations" = "darkolivegreen", 
                     "Refractory cells population"="darkolivegreen4",
                     "Antibody population"="darkolivegreen3",
                     "Antibody and Effector cells populations" = "darkolivegreen2")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = V,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075)
      
      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = V, color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075)
      
      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = V, color = "Antibody population",
                                               linetype="Antibody population"),size=1.075)
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = V, color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=2)
      
      vdisp[4] <- TRUE}
    pl <- pl + geom_point(patient3_obsV,mapping = aes(x = time, y = V), shape=8, color = "darkgreen", size=4)
    pl <- pl + labs(x='Time (Days)') + labs(y='HBV DNA (virions/mL)') +
      labs(color  = "Model", linetype = "Model")
    
    pl <- pl + theme(legend.position=c(0.1, 0.25), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(0.001, 10^11),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 
    
    print(pl)
    
  })

  
  output$plot5 <- renderPlot({
    
    parms5_2I = c(
      r = 0.01, 
      tmax = 13.6e6,
      beta = 6.7e-10, 
      rho1 = 0.05,
      mu = 5.4e-4, 
      rho2 = 0.051,
      delta_tr = 0.51,
      p1 = 29,
      p2 = 42,
      c = 0.67,
      z = 10,
      alpha = 3.9e-7,
      tau = 20.7,
      de = 0.5 
    )
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = 13.6 * 10 ** 6, I1 = 0, I2 = 0, V = 0.33, E = 20)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau, I1, lagvalue(t - tau, 2))
        lagI2 <- ifelse(t < tau, I2, lagvalue(t - tau, 3))
        lagE <- ifelse(t < tau, E, lagvalue(t - tau, 5))
        
        dT <- r * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta * T * V + rho1 * I1
        dI1 <- beta * T * V - (rho1 + delta_tr) * I1 - mu * E * I1 + rho2 * I2
        dI2 <- r * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu * E * I2
        dV <- p1 * I1 + p2 * I2 - c * V
        dE <- z + alpha * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = parms5_2I)
    
    parms5_R <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 6.49 * 10 ** (-10),
      rho = 0.119 * 10 ** (-3),
      mu = 20.0 * 10 ** (-4),
      mu_1 = 0.4 * 10 ** (-6),
      rho_r = 2.0 * 10 ** (-5),
      p = 39,
      c = 0.67,
      z = 10,
      alpha = 1.2 * 10 ** (-7),
      tau = 16.8,
      de = 0.5
    )
    
    yinitR <- c(T = 13.6 * 10 ** 6, I = 0, R = 0, V = 0.33, E = 20)
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,5)
        
        dT <- r*T*(1-((T+I+R)/tmax)) - beta*T*V + rho_r*R
        dI <- r*I*(1-((T+I+R)/tmax)) + beta*T*V - rho*I*E - mu*I*E 
        dR <- rho*I*E + r*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_1*R*E
        dV <- p*I - c*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = parms5_R)
    
    parms5_A <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 6.1 * 10 ** (-10),
      delta = 0.049,
      km = 10,
      kp = 10 ** (-12),
      cxav = 2.7,
      pa = 10 * 10 ** (-5),
      theta = 1000,
      ra = 0.3578,
      amax = 4 * 10 ** 15,
      da = 0.033,
      p = 42.5,
      c = 0.67
    )
    
    yinitA <- c(T = 13.6 * 10 ** 6, I = 0, XAV = 0, A = 0, V = 0.33)
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- beta*T*V-delta*I
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*(1+theta)+ra*A*(1-A/amax)+km*XAV*(1+theta)-kp*A*V*(1+theta)-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = parms5_A)
    
    parms5_AE <- c(
      r = 1.0,
      tmax = 1.36e7,
      beta = 4.95e-10,
      delta = 0.082,
      mu = 0.6e-3,
      km = 10,
      kp = 1e-10,
      cxav = 4*1.67,
      pa = 0.9e-3,
      theta = 1e3,
      ra = 0.41,
      amax = 4e11,
      da = 0.033,
      p = 186.7,
      c = 1.67,
      z= 30,
      alpha = 2.2e-7,
      tau = 16.8,
      de = 0.5
    )
    
    yinitAE <- c(T = 1.36e7, I = 0, XAV = 0, A = 0, V = 0.33, E = 20)
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,6)
        
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- r*I*(1-(T+I)/tmax)+beta*T*V-delta*I-mu*I*E
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*theta+ra*A*(1-A/amax)+km*XAV*theta-kp*A*V*theta-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = parms5_AE)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               V= output1[,5])
    finaloutput2 <- data.frame(timesR= output2[,1],
                               V= output2[,5])
    finaloutput3 <- data.frame(timesA= output3[,1],
                               V= output3[,6])
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               V= output4[,6])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9)
    colormanual <- c("Two infected cells populations" = "darkolivegreen", 
                     "Refractory cells population"="darkolivegreen4",
                     "Antibody population"="darkolivegreen3",
                     "Antibody and Effector cells populations" = "darkolivegreen2")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = V,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075)
      
      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = V, color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075)
      
      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = V, color = "Antibody population",
                                               linetype="Antibody population"),size=1.075)
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = V, color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=2)
      
      vdisp[4] <- TRUE}
    pl <- pl + geom_point(patient5_obsV,mapping = aes(x = time, y = V), shape=8, color = "darkgreen", size=4)
    pl <- pl + labs(x='Time (Days)') + labs(y='HBV DNA (virions/mL)') +
      labs(color  = "Model", linetype = "Model")
    
    pl <- pl + theme(legend.position=c(0.1, 0.25), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(0.001, 10^11),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 
    
    print(pl)
    
  })
  
  output$plot6 <- renderPlot({
    
    parms6_2I = c(
      r = 0.01, 
      tmax = 13.6e6,
      beta = 0.8e-10, 
      rho1 = 0.05,
      mu = 6.4e-4, 
      rho2 = 0.18,
      delta_tr = 0.22,
      p1 = 232,
      p2 = 232,
      c = 0.67,
      z = 10,
      alpha = 5.2e-7,
      tau = 24.9,
      de = 0.5 
    )
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = 13.6 * 10 ** 6, I1 = 0, I2 = 0, V = 0.33, E = 20)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau, I1, lagvalue(t - tau, 2))
        lagI2 <- ifelse(t < tau, I2, lagvalue(t - tau, 3))
        lagE <- ifelse(t < tau, E, lagvalue(t - tau, 5))
        
        dT <- r * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta * T * V + rho1 * I1
        dI1 <- beta * T * V - (rho1 + delta_tr) * I1 - mu * E * I1 + rho2 * I2
        dI2 <- r * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu * E * I2
        dV <- p1 * I1 + p2 * I2 - c * V
        dE <- z + alpha * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = parms6_2I)
    
    parms6_R <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 2.22 * 10 ** (-10),
      rho = 1.03 * 10 ** (-3),
      mu = 0.003 * 10 ** (-4),
      mu_1 = 180.0 * 10 ** (-6),
      rho_r = 2.0 * 10 ** (-5),
      p = 76,
      c = 0.67,
      z = 10,
      alpha = 4.4 * 10 ** (-7),
      tau = 29.7,
      de = 0.5
    )
    
    yinitR <- c(T = 13.6 * 10 ** 6, I = 0, R = 0, V = 0.33, E = 20)
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,5)
        
        dT <- r*T*(1-((T+I+R)/tmax)) - beta*T*V + rho_r*R
        dI <- r*I*(1-((T+I+R)/tmax)) + beta*T*V - rho*I*E - mu*I*E 
        dR <- rho*I*E + r*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_1*R*E
        dV <- p*I - c*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = parms6_R)
    
    parms6_A <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 5.87 * 10 ** (-10),
      delta = 0.043,
      km = 10,
      kp = 10 ** (-12),
      cxav = 2.7,
      pa = 1.16 * 10 ** (-5),
      theta = 1000,
      ra = 0.29,
      amax = 4 * 10 ** 15,
      da = 0.033,
      p = 28.5,
      c = 0.67
    )
    
    yinitA <- c(T = 13.6 * 10 ** 6, I = 0, XAV = 0, A = 0, V = 0.33)
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- beta*T*V-delta*I
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*(1+theta)+ra*A*(1-A/amax)+km*XAV*(1+theta)-kp*A*V*(1+theta)-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = parms6_A)
    
    parms6_AE <- c(
      r = 1.0,
      tmax = 1.36e7,
      beta = 2.21e-10,
      delta = 0.073,
      mu = 0.5e-3, 
      km = 10,
      kp = 1e-10,
      cxav = 4*1.67,
      pa = 1.5e-3,
      theta = 1e3,
      ra = 0.22,
      amax = 4e11,
      da = 0.033,
      p = 327.3,
      c = 1.67,
      z= 10,
      alpha = 2.2e-7,
      tau = 29.7,
      de = 0.5
    )
    
    yinitAE <- c(T = 1.36e7, I = 0, XAV = 0, A = 0, V = 0.33, E = 20)
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,6)
        
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- r*I*(1-(T+I)/tmax)+beta*T*V-delta*I-mu*I*E
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*theta+ra*A*(1-A/amax)+km*XAV*theta-kp*A*V*theta-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = parms6_AE)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               V= output1[,5])
    finaloutput2 <- data.frame(timesR= output2[,1],
                               V= output2[,5])
    finaloutput3 <- data.frame(timesA= output3[,1],
                               V= output3[,6])
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               V= output4[,6])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9)
    colormanual <- c("Two infected cells populations" = "darkolivegreen", 
                     "Refractory cells population"="darkolivegreen4",
                     "Antibody population"="darkolivegreen3",
                     "Antibody and Effector cells populations" = "darkolivegreen2")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = V,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075)
      
      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = V, color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075)
      
      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = V, color = "Antibody population",
                                               linetype="Antibody population"),size=1.075)
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = V, color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=2)
      
      vdisp[4] <- TRUE}
    pl <- pl + geom_point(patient6_obsV,mapping = aes(x = time, y = V), shape=8, color = "darkgreen", size=4)
    pl <- pl + labs(x='Time (Days)') + labs(y='HBV DNA (virions/mL)') +
      labs(color  = "Model", linetype = "Model")
    
    pl <- pl + theme(legend.position=c(0.1, 0.25), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(0.001, 10^11),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 
    
    print(pl)
    
  })
  
  output$plot7 <- renderPlot({
    
    parms7_2I = c(
      r = 0.01, 
      tmax = 13.6e6,
      beta = 0.9e-10, 
      rho1 = 0.11,
      mu = 5.9e-4, 
      rho2 = 0.05,
      delta_tr = 0.84,
      p1 = 83,
      p2 = 275,
      c = 0.67,
      z = 10,
      alpha = 2.3e-7,
      tau = 29.9,
      de = 0.5 
    )
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = 13.6 * 10 ** 6, I1 = 0, I2 = 0, V = 0.33, E = 20)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau, I1, lagvalue(t - tau, 2))
        lagI2 <- ifelse(t < tau, I2, lagvalue(t - tau, 3))
        lagE <- ifelse(t < tau, E, lagvalue(t - tau, 5))
        
        dT <- r * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta * T * V + rho1 * I1
        dI1 <- beta * T * V - (rho1 + delta_tr) * I1 - mu * E * I1 + rho2 * I2
        dI2 <- r * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu * E * I2
        dV <- p1 * I1 + p2 * I2 - c * V
        dE <- z + alpha * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = parms7_2I)
    
    parms7_R <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 1.22 * 10 ** (-10),
      rho = 0.338 * 10 ** (-3),
      mu = 1.2 * 10 ** (-4),
      mu_1 = 12.7 * 10 ** (-6),
      rho_r = 2.0 * 10 ** (-5),
      p = 164,
      c = 0.67,
      z = 10,
      alpha = 2.0 * 10 ** (-7),
      tau = 33.4,
      de = 0.5
    )
    
    yinitR <- c(T = 13.6 * 10 ** 6, I = 0, R = 0, V = 0.33, E = 20)
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,5)
        
        dT <- r*T*(1-((T+I+R)/tmax)) - beta*T*V + rho_r*R
        dI <- r*I*(1-((T+I+R)/tmax)) + beta*T*V - rho*I*E - mu*I*E 
        dR <- rho*I*E + r*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_1*R*E
        dV <- p*I - c*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = parms7_R)
    
    parms7_A <- c(
      r = 1.0,
      tmax = 13.6 * 10 ** 6,
      beta = 2.0 * 10 ** (-10),
      delta = 0.0000998,
      km = 10,
      kp = 10 ** (-12),
      cxav = 2.7,
      pa = 0.41 * 10 ** (-5),
      theta = 1000,
      ra = 0.4623,
      amax = 4 * 10 ** 15,
      da = 0.033,
      p = 113,
      c = 0.67
    )
    
    yinitA <- c(T = 13.6 * 10 ** 6, I = 0, XAV = 0, A = 0, V = 0.33)
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- beta*T*V-delta*I
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*(1+theta)+ra*A*(1-A/amax)+km*XAV*(1+theta)-kp*A*V*(1+theta)-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = parms7_A)
    
    parms7_AE <- c(
      r = 1.0,
      tmax = 1.36e7,
      beta = 2.6e-10,
      delta = 0.018,
      mu = 0.1e-3,
      km = 10,
      kp = 1e-10,
      cxav = 4*1.67,
      pa = 0.1e-3,
      theta = 1e3,
      ra = 0.35,
      amax = 4e11,
      da = 0.033,
      p = 305,
      c = 1.67,
      z= 30,
      alpha = 2.2e-7,
      tau = 33.4,
      de = 0.5
    )
    
    yinitAE <- c(T = 1.36e7, I = 0, XAV = 0, A = 0, V = 0.33, E = 20)
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau)
          lagI <- I
        else
          lagI <- lagvalue(t - tau,2)
        if (t < tau)
          lagE <- E
        else
          lagE <- lagvalue(t - tau,6)
        
        dT <- r*T*(1-(T+I)/tmax)-beta*T*V
        dI <- r*I*(1-(T+I)/tmax)+beta*T*V-delta*I-mu*I*E
        dXAV <- -km*XAV+kp*A*V-cxav*XAV
        dA <- pa*V*theta+ra*A*(1-A/amax)+km*XAV*theta-kp*A*V*theta-da*A
        dV <- p*I-c*V+km*XAV-kp*A*V
        dE <- z + alpha*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = parms7_AE)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               V= output1[,5])
    finaloutput2 <- data.frame(timesR= output2[,1],
                               V= output2[,5])
    finaloutput3 <- data.frame(timesA= output3[,1],
                               V= output3[,6])
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               V= output4[,6])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9)
    colormanual <- c("Two infected cells populations" = "darkolivegreen", 
                     "Refractory cells population"="darkolivegreen4",
                     "Antibody population"="darkolivegreen3",
                     "Antibody and Effector cells populations" = "darkolivegreen2")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = V,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075)
      
      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = V, color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075)
      
      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = V, color = "Antibody population",
                                               linetype="Antibody population"),size=1.075)
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = V, color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=2)
      
      vdisp[4] <- TRUE}
    pl <- pl + geom_point(patient7_obsV,mapping = aes(x = time, y = V), shape=8, color = "darkgreen", size=4)
    pl <- pl + labs(x='Time (Days)') + labs(y='HBV DNA (virions/mL)') +
      labs(color  = "Model", linetype = "Model")
    
    pl <- pl + theme(legend.position=c(0.1, 0.25), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(0.001, 10^11),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 
    
    print(pl)
    
  })
  
})