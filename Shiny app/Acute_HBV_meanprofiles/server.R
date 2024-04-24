#Server code for modeling via differential equations for oral and iv PK

library(deSolve)
library(shiny)
library(mlxR)
library(shinydashboard)
library(ggplot2)
library(scales)

#Set the directory where you saved the observations files
obsV <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_3/Dataset/HBVDNA_allpatients_Ciupe.csv")
obsALT <- read.csv("C:/Users/boivinchampeauxc/Documents/Projet/Software - literature modeling/Shiny App/Acute_HBV_3/Dataset/ALT_allpatients_Ciupe.csv")

shinyServer(function(input, output) {
  output$plot1 <- renderPlot({
    
    r1 <- input$r1 #Proliferation rate (INPUT)
    beta1 <- input$beta1/10000000000 #Infection rate (INPUT)
    mu1 <- input$mu1/10000 #Cytolytic rate infected cells (INPUT)
    tau1 <-input$tau1 #Activation delay (INPUT)
    alpha1 <-input$alpha1/10000000 #Activation and proliferation of effector cells (INPUT)
    delta_tr <- input$delta_tr #Transition rate constant (conditional INPUT) 
    rho1 <- input$rho1 #Cure rate I1 to T (conditional INPUT)
    rho2 <- input$rho2 #Cure rate I2 to I1 (conditional INPUT)
    pI1 <- input$pI1 #Production rate constant from I1 (conditional INPUT)
    pI2 <- input$pI2 #Production rate constant from I2 (conditional INPUT)
    r2 <- input$r2 #Proliferation rate (INPUT)
    beta2 <- input$beta2/10000000000 #Infection rate (INPUT)
    mu2 <- input$mu2/10000 #Cytolytic rate infected cells (INPUT)
    tau2 <-input$tau2 #Activation delay (INPUT)
    alpha2 <-input$alpha2/10000000 #Activation and proliferation of effector cells (INPUT)
    rho <- input$rho/1000 #Cure rate I to T (conditional INPUT)
    rho_r <- input$rho_r/100000 #Cure rate R to T (conditional INPUT)
    mu_r <- input$mu_r/1000000 #Cytolytic rate constant of refractory cells (conditional INPUT)
    p2 <- input$p2 #Production rate constant from I (conditional INPUT)
    beta3 <- input$beta3/10000000000 #Infection rate (INPUT)
    p3 <- input$p3 #Production rate constant from I (conditional INPUT)
    delta3 <- input$delta3 #Death rate infected cells
    ra3 <- input$ra3 #Proliferation rate antibody
    pa3 <- input$pa3/100000 #Production rate constant antibody
    kp3 <- input$kp3/1000000000000
    cxav3 <- input$cxav3
    beta4 <- input$beta4/10000000000 #Infection rate (INPUT)
    p4 <- input$p4 #Production rate constant from I (conditional INPUT)
    delta4 <- input$delta4 #Death rate infected cells
    ra4 <- input$ra4 #Proliferation rate antibody
    pa4 <- input$pa4/1000 #Production rate constant antibody
    mu4 <- input$mu4/1000 #Cytolytic rate infected cells (INPUT)
    kp4 <- input$kp4/10000000000
    c4 <- input$c4
    cxav4 <- input$cxav4
    alpha4 <-input$alpha4/10000000 #Activation and proliferation of effector cells (INPUT)
    c <- input$c #Virus clearance rate (INPUT)
    z <-input$z #Source of effector cells (INPUT)
    de <-input$de #Effector cells clearance rate (INPUT)
    km <- input$km
    amax3 <- input$amax3*1000000000000000 
    amax4 <- input$amax4*100000000000 
    theta <- input$theta
    da <- input$da 
    tmax <- input$tmax*10000000 #Maximum liver capacity (INPUT)
    
    t <- input$num2 #Initial uninfected cell (INPUT)
    v <- input$num3 #Initial free virus (INPUT)
    e <- input$num4 #Initial effector cells (INPUT)
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = t, I1 = 0, I2 = 0, V = v, E = e)
    yinitR <- c(T = t, I = 0, R = 0, V = v, E = e)
    yinitA <- c(T = t, I = 0, XAV = 0, A = 0, V = v)
    yinitAE <- c(T = t, I = 0, XAV = 0, A = 0, V = v, E = e)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau1, I1, lagvalue(t - tau1, 2))
        lagI2 <- ifelse(t < tau1, I2, lagvalue(t - tau1, 3))
        lagE <- ifelse(t < tau1, E, lagvalue(t - tau1, 5))
        
        dT <- r1 * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta1 * T * V + rho1 * I1
        dI1 <- beta1 * T * V - (rho1 + delta_tr) * I1 - mu1 * E * I1 + rho2 * I2
        dI2 <- r1 * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu1 * E * I2
        dV <- pI1 * I1 + pI2 * I2 - c * V
        dE <- z + alpha1 * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = NULL)
    
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,5)
        
        dT <- r2*T*(1-((T+I+R)/tmax)) - beta2*T*V + rho_r*R
        dI <- r2*I*(1-((T+I+R)/tmax)) + beta2*T*V - rho*I*E - mu2*I*E 
        dR <- rho*I*E + r2*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_r*R*E
        dV <- p2*I - c*V
        dE <- z + alpha2*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = NULL)
    
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r2*T*(1-(T+I)/tmax)-beta3*T*V
        dI <- beta3*T*V-delta3*I
        dXAV <- -km*XAV+kp3*A*V-cxav3*XAV
        dA <- pa3*V*(1+theta)+ra3*A*(1-A/amax3)+km*XAV*(1+theta)-kp3*A*V*(1+theta)-da*A
        dV <- p3*I-c*V+km*XAV-kp3*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = NULL)
    
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,6)
        
        dT <- r2*T*(1-(T+I)/tmax)-beta4*T*V
        dI <- r2*I*(1-(T+I)/tmax)+beta4*T*V-delta4*I-mu4*I*E
        dXAV <- -km*XAV+kp4*A*V-cxav4*XAV
        dA <- pa4*V*theta+ra4*A*(1-A/amax4)+km*XAV*theta-kp4*A*V*theta-da*A
        dV <- p4*I-c4*V+km*XAV-kp4*A*V
        dE <- z + alpha4*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = NULL)
    
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
    
    pl <- pl + geom_point(obsV,mapping = aes(x = time, y = V), shape=4, color = "black", size=1)
    pl <- pl + labs(x='Time (Days)') + labs(y='HBV DNA (virions/mL)') +
      labs(color  = "Model", linetype = "Model")
    pl <- pl + theme(legend.position=c(.60, 0.25), legend.justification=c(0,1), legend.title=element_blank())
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
    
    r1 <- input$r1 #Proliferation rate (INPUT)
    beta1 <- input$beta1/10000000000 #Infection rate (INPUT)
    mu1 <- input$mu1/10000 #Cytolytic rate infected cells (INPUT)
    tau1 <-input$tau1 #Activation delay (INPUT)
    alpha1 <-input$alpha1/10000000 #Activation and proliferation of effector cells (INPUT)
    delta_tr <- input$delta_tr #Transition rate constant (conditional INPUT) 
    rho1 <- input$rho1 #Cure rate I1 to T (conditional INPUT)
    rho2 <- input$rho2 #Cure rate I2 to I1 (conditional INPUT)
    pI1 <- input$pI1 #Production rate constant from I1 (conditional INPUT)
    pI2 <- input$pI2 #Production rate constant from I2 (conditional INPUT)
    r2 <- input$r2 #Proliferation rate (INPUT)
    beta2 <- input$beta2/10000000000 #Infection rate (INPUT)
    mu2 <- input$mu2/10000 #Cytolytic rate infected cells (INPUT)
    tau2 <-input$tau2 #Activation delay (INPUT)
    alpha2 <-input$alpha2/10000000 #Activation and proliferation of effector cells (INPUT)
    rho <- input$rho/1000 #Cure rate I to T (conditional INPUT)
    rho_r <- input$rho_r/100000 #Cure rate R to T (conditional INPUT)
    mu_r <- input$mu_r/1000000 #Cytolytic rate constant of refractory cells (conditional INPUT)
    p2 <- input$p2 #Production rate constant from I (conditional INPUT)
    beta3 <- input$beta3/10000000000 #Infection rate (INPUT)
    p3 <- input$p3 #Production rate constant from I (conditional INPUT)
    delta3 <- input$delta3 #Death rate infected cells
    ra3 <- input$ra3 #Proliferation rate antibody
    pa3 <- input$pa3/100000 #Production rate constant antibody
    kp3 <- input$kp3/1000000000000
    cxav3 <- input$cxav3
    beta4 <- input$beta4/10000000000 #Infection rate (INPUT)
    p4 <- input$p4 #Production rate constant from I (conditional INPUT)
    delta4 <- input$delta4 #Death rate infected cells
    ra4 <- input$ra4 #Proliferation rate antibody
    pa4 <- input$pa4/1000 #Production rate constant antibody
    mu4 <- input$mu4/1000 #Cytolytic rate infected cells (INPUT)
    kp4 <- input$kp4/10000000000
    c4 <- input$c4
    cxav4 <- input$cxav4
    alpha4 <-input$alpha4/10000000 #Activation and proliferation of effector cells (INPUT)
    c <- input$c #Virus clearance rate (INPUT)
    z <-input$z #Source of effector cells (INPUT)
    de <-input$de #Effector cells clearance rate (INPUT)
    km <- input$km
    amax3 <- input$amax3*1000000000000000 
    amax4 <- input$amax4*100000000000 
    theta <- input$theta
    da <- input$da 
    tmax <- input$tmax*10000000 #Maximum liver capacity (INPUT)
    
    t <- input$num2 #Initial uninfected cell (INPUT)
    v <- input$num3 #Initial free virus (INPUT)
    e <- input$num4 #Initial effector cells (INPUT)
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = t, I1 = 0, I2 = 0, V = v, E = e)
    yinitR <- c(T = t, I = 0, R = 0, V = v, E = e)
    yinitA <- c(T = t, I = 0, XAV = 0, A = 0, V = v)
    yinitAE <- c(T = t, I = 0, XAV = 0, A = 0, V = v, E = e)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau1, I1, lagvalue(t - tau1, 2))
        lagI2 <- ifelse(t < tau1, I2, lagvalue(t - tau1, 3))
        lagE <- ifelse(t < tau1, E, lagvalue(t - tau1, 5))
        
        dT <- r1 * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta1 * T * V + rho1 * I1
        dI1 <- beta1 * T * V - (rho1 + delta_tr) * I1 - mu1 * E * I1 + rho2 * I2
        dI2 <- r1 * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu1 * E * I2
        dV <- pI1 * I1 + pI2 * I2 - c * V
        dE <- z + alpha1 * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = NULL)
    
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,5)
        
        dT <- r2*T*(1-((T+I+R)/tmax)) - beta2*T*V + rho_r*R
        dI <- r2*I*(1-((T+I+R)/tmax)) + beta2*T*V - rho*I*E - mu2*I*E 
        dR <- rho*I*E + r2*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_r*R*E
        dV <- p2*I - c*V
        dE <- z + alpha2*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = NULL)
    
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r2*T*(1-(T+I)/tmax)-beta3*T*V
        dI <- beta3*T*V-delta3*I
        dXAV <- -km*XAV+kp3*A*V-cxav3*XAV
        dA <- pa3*V*(1+theta)+ra3*A*(1-A/amax3)+km*XAV*(1+theta)-kp3*A*V*(1+theta)-da*A
        dV <- p3*I-c*V+km*XAV-kp3*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = NULL)
    
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,6)
        
        dT <- r2*T*(1-(T+I)/tmax)-beta4*T*V
        dI <- r2*I*(1-(T+I)/tmax)+beta4*T*V-delta4*I-mu4*I*E
        dXAV <- -km*XAV+kp4*A*V-cxav4*XAV
        dA <- pa4*V*theta+ra4*A*(1-A/amax4)+km*XAV*theta-kp4*A*V*theta-da*A
        dV <- p4*I-c4*V+km*XAV-kp4*A*V
        dE <- z + alpha4*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = NULL)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               E= output1[,6])
    finaloutput2 <- data.frame(timesR= output2[,1],
                               E= output2[,6])
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               E= output4[,7])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9) 
    colormanual <- c("Two infected cells populations" = "deeppink4", 
                     "Refractory cells population"="#FF3399",
                     "Antibody and Effector cells populations" ="pink")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = E,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075) 

      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = E,color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075) 

      vdisp[2] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = E,color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=2) 
      
      vdisp[3] <- TRUE}
    
    pl <- pl + geom_point(obsALT,mapping = aes(x = Time, y = Serum_ALT), shape=4, color = "black", size=1)
    pl <- pl + labs(x='Time (Days)') + labs(y='Effector cells (cell/mL)') +
      labs(color  = "Model", linetype = "Model")
    pl <- pl + labs(x='Time (Days)') + labs(y='Effector cells (cell/mL)') +
      labs(color  = "Model", linetype = "Model")
    pl <- pl + theme(legend.position=c(.60, 0.25), legend.justification=c(0,1), legend.title=element_blank())
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
    
    r1 <- input$r1 #Proliferation rate (INPUT)
    beta1 <- input$beta1/10000000000 #Infection rate (INPUT)
    mu1 <- input$mu1/10000 #Cytolytic rate infected cells (INPUT)
    tau1 <-input$tau1 #Activation delay (INPUT)
    alpha1 <-input$alpha1/10000000 #Activation and proliferation of effector cells (INPUT)
    delta_tr <- input$delta_tr #Transition rate constant (conditional INPUT) 
    rho1 <- input$rho1 #Cure rate I1 to T (conditional INPUT)
    rho2 <- input$rho2 #Cure rate I2 to I1 (conditional INPUT)
    pI1 <- input$pI1 #Production rate constant from I1 (conditional INPUT)
    pI2 <- input$pI2 #Production rate constant from I2 (conditional INPUT)
    r2 <- input$r2 #Proliferation rate (INPUT)
    beta2 <- input$beta2/10000000000 #Infection rate (INPUT)
    mu2 <- input$mu2/10000 #Cytolytic rate infected cells (INPUT)
    tau2 <-input$tau2 #Activation delay (INPUT)
    alpha2 <-input$alpha2/10000000 #Activation and proliferation of effector cells (INPUT)
    rho <- input$rho/1000 #Cure rate I to T (conditional INPUT)
    rho_r <- input$rho_r/100000 #Cure rate R to T (conditional INPUT)
    mu_r <- input$mu_r/1000000 #Cytolytic rate constant of refractory cells (conditional INPUT)
    p2 <- input$p2 #Production rate constant from I (conditional INPUT)
    beta3 <- input$beta3/10000000000 #Infection rate (INPUT)
    p3 <- input$p3 #Production rate constant from I (conditional INPUT)
    delta3 <- input$delta3 #Death rate infected cells
    ra3 <- input$ra3 #Proliferation rate antibody
    pa3 <- input$pa3/100000 #Production rate constant antibody
    kp3 <- input$kp3/1000000000000
    cxav3 <- input$cxav3
    beta4 <- input$beta4/10000000000 #Infection rate (INPUT)
    p4 <- input$p4 #Production rate constant from I (conditional INPUT)
    delta4 <- input$delta4 #Death rate infected cells
    ra4 <- input$ra4 #Proliferation rate antibody
    pa4 <- input$pa4/1000 #Production rate constant antibody
    mu4 <- input$mu4/1000 #Cytolytic rate infected cells (INPUT)
    kp4 <- input$kp4/10000000000
    c4 <- input$c4
    cxav4 <- input$cxav4
    alpha4 <-input$alpha4/10000000 #Activation and proliferation of effector cells (INPUT)
    c <- input$c #Virus clearance rate (INPUT)
    z <-input$z #Source of effector cells (INPUT)
    de <-input$de #Effector cells clearance rate (INPUT)
    km <- input$km
    amax3 <- input$amax3*1000000000000000 
    amax4 <- input$amax4*100000000000 
    theta <- input$theta
    da <- input$da 
    tmax <- input$tmax*10000000 #Maximum liver capacity (INPUT)
    
    t <- input$num2 #Initial uninfected cell (INPUT)
    v <- input$num3 #Initial free virus (INPUT)
    e <- input$num4 #Initial effector cells (INPUT)
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = t, I1 = 0, I2 = 0, V = v, E = e)
    yinitR <- c(T = t, I = 0, R = 0, V = v, E = e)
    yinitA <- c(T = t, I = 0, XAV = 0, A = 0, V = v)
    yinitAE <- c(T = t, I = 0, XAV = 0, A = 0, V = v, E = e)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau1, I1, lagvalue(t - tau1, 2))
        lagI2 <- ifelse(t < tau1, I2, lagvalue(t - tau1, 3))
        lagE <- ifelse(t < tau1, E, lagvalue(t - tau1, 5))
        
        dT <- r1 * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta1 * T * V + rho1 * I1
        dI1 <- beta1 * T * V - (rho1 + delta_tr) * I1 - mu1 * E * I1 + rho2 * I2
        dI2 <- r1 * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu1 * E * I2
        dV <- pI1 * I1 + pI2 * I2 - c * V
        dE <- z + alpha1 * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = NULL)
    
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,5)
        
        dT <- r2*T*(1-((T+I+R)/tmax)) - beta2*T*V + rho_r*R
        dI <- r2*I*(1-((T+I+R)/tmax)) + beta2*T*V - rho*I*E - mu2*I*E 
        dR <- rho*I*E + r2*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_r*R*E
        dV <- p2*I - c*V
        dE <- z + alpha2*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = NULL)
    
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r2*T*(1-(T+I)/tmax)-beta3*T*V
        dI <- beta3*T*V-delta3*I
        dXAV <- -km*XAV+kp3*A*V-cxav3*XAV
        dA <- pa3*V*(1+theta)+ra3*A*(1-A/amax3)+km*XAV*(1+theta)-kp3*A*V*(1+theta)-da*A
        dV <- p3*I-c*V+km*XAV-kp3*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = NULL)
    
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,6)
        
        dT <- r2*T*(1-(T+I)/tmax)-beta4*T*V
        dI <- r2*I*(1-(T+I)/tmax)+beta4*T*V-delta4*I-mu4*I*E
        dXAV <- -km*XAV+kp4*A*V-cxav4*XAV
        dA <- pa4*V*theta+ra4*A*(1-A/amax4)+km*XAV*theta-kp4*A*V*theta-da*A
        dV <- p4*I-c4*V+km*XAV-kp4*A*V
        dE <- z + alpha4*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = NULL)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               I1= output1[,3],
                               I2= output1[,4])
    
    finaloutput2 <- data.frame(timesR= output2[,1],
                               I= output2[,3])
    
    finaloutput3 <- data.frame(timesA= output3[,1],
                               I= output3[,3])
    
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               I= output4[,3])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9) 
    colormanual <- c("Two infected cells populations" = "darkred", 
                     "Refractory cells population"="firebrick3",
                     "Antibody population"="firebrick2",
                     "Antibody and Effector cells populations" = "firebrick1")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = I1,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075) +
        geom_line(data=finaloutput1, aes(x=times2I, y = I2,color = "Two infected cells populations",
                                         linetype = "Two infected cells populations"),size=1.075)

      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = I,color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075) 

      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = I,color = "Antibody population",
                                               linetype="Antibody population"),size=1.075) 
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = I,color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=1.075) 
      
      vdisp[4] <- TRUE}
    pl <- pl + labs(x='Time (Days)') + labs(y='Infected cells (cell/mL)') +
      labs(color  = "Model", linetype = "Model")

    pl <- pl + theme(legend.position=c(.60, 0.25), legend.justification=c(0,1), legend.title=element_blank())
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
  
  output$plot4 <- renderPlot({
    
    r1 <- input$r1 #Proliferation rate (INPUT)
    beta1 <- input$beta1/10000000000 #Infection rate (INPUT)
    mu1 <- input$mu1/10000 #Cytolytic rate infected cells (INPUT)
    tau1 <-input$tau1 #Activation delay (INPUT)
    alpha1 <-input$alpha1/10000000 #Activation and proliferation of effector cells (INPUT)
    delta_tr <- input$delta_tr #Transition rate constant (conditional INPUT) 
    rho1 <- input$rho1 #Cure rate I1 to T (conditional INPUT)
    rho2 <- input$rho2 #Cure rate I2 to I1 (conditional INPUT)
    pI1 <- input$pI1 #Production rate constant from I1 (conditional INPUT)
    pI2 <- input$pI2 #Production rate constant from I2 (conditional INPUT)
    r2 <- input$r2 #Proliferation rate (INPUT)
    beta2 <- input$beta2/10000000000 #Infection rate (INPUT)
    mu2 <- input$mu2/10000 #Cytolytic rate infected cells (INPUT)
    tau2 <-input$tau2 #Activation delay (INPUT)
    alpha2 <-input$alpha2/10000000 #Activation and proliferation of effector cells (INPUT)
    rho <- input$rho/1000 #Cure rate I to T (conditional INPUT)
    rho_r <- input$rho_r/100000 #Cure rate R to T (conditional INPUT)
    mu_r <- input$mu_r/1000000 #Cytolytic rate constant of refractory cells (conditional INPUT)
    p2 <- input$p2 #Production rate constant from I (conditional INPUT)
    beta3 <- input$beta3/10000000000 #Infection rate (INPUT)
    p3 <- input$p3 #Production rate constant from I (conditional INPUT)
    delta3 <- input$delta3 #Death rate infected cells
    ra3 <- input$ra3 #Proliferation rate antibody
    pa3 <- input$pa3/100000 #Production rate constant antibody
    kp3 <- input$kp3/1000000000000
    cxav3 <- input$cxav3
    beta4 <- input$beta4/10000000000 #Infection rate (INPUT)
    p4 <- input$p4 #Production rate constant from I (conditional INPUT)
    delta4 <- input$delta4 #Death rate infected cells
    ra4 <- input$ra4 #Proliferation rate antibody
    pa4 <- input$pa4/1000 #Production rate constant antibody
    mu4 <- input$mu4/1000 #Cytolytic rate infected cells (INPUT)
    kp4 <- input$kp4/10000000000
    c4 <- input$c4
    cxav4 <- input$cxav4
    alpha4 <-input$alpha4/10000000 #Activation and proliferation of effector cells (INPUT)
    c <- input$c #Virus clearance rate (INPUT)
    z <-input$z #Source of effector cells (INPUT)
    de <-input$de #Effector cells clearance rate (INPUT)
    km <- input$km
    amax3 <- input$amax3*1000000000000000 
    amax4 <- input$amax4*100000000000 
    theta <- input$theta
    da <- input$da 
    tmax <- input$tmax*10000000 #Maximum liver capacity (INPUT)
    
    t <- input$num2 #Initial uninfected cell (INPUT)
    v <- input$num3 #Initial free virus (INPUT)
    e <- input$num4 #Initial effector cells (INPUT)
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = t, I1 = 0, I2 = 0, V = v, E = e)
    yinitR <- c(T = t, I = 0, R = 0, V = v, E = e)
    yinitA <- c(T = t, I = 0, XAV = 0, A = 0, V = v)
    yinitAE <- c(T = t, I = 0, XAV = 0, A = 0, V = v, E = e)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau1, I1, lagvalue(t - tau1, 2))
        lagI2 <- ifelse(t < tau1, I2, lagvalue(t - tau1, 3))
        lagE <- ifelse(t < tau1, E, lagvalue(t - tau1, 5))
        
        dT <- r1 * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta1 * T * V + rho1 * I1
        dI1 <- beta1 * T * V - (rho1 + delta_tr) * I1 - mu1 * E * I1 + rho2 * I2
        dI2 <- r1 * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu1 * E * I2
        dV <- pI1 * I1 + pI2 * I2 - c * V
        dE <- z + alpha1 * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = NULL)
    
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,5)
        
        dT <- r2*T*(1-((T+I+R)/tmax)) - beta2*T*V + rho_r*R
        dI <- r2*I*(1-((T+I+R)/tmax)) + beta2*T*V - rho*I*E - mu2*I*E 
        dR <- rho*I*E + r2*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_r*R*E
        dV <- p2*I - c*V
        dE <- z + alpha2*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = NULL)
    
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r2*T*(1-(T+I)/tmax)-beta3*T*V
        dI <- beta3*T*V-delta3*I
        dXAV <- -km*XAV+kp3*A*V-cxav3*XAV
        dA <- pa3*V*(1+theta)+ra3*A*(1-A/amax3)+km*XAV*(1+theta)-kp3*A*V*(1+theta)-da*A
        dV <- p3*I-c*V+km*XAV-kp3*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = NULL)
    
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,6)
        
        dT <- r2*T*(1-(T+I)/tmax)-beta4*T*V
        dI <- r2*I*(1-(T+I)/tmax)+beta4*T*V-delta4*I-mu4*I*E
        dXAV <- -km*XAV+kp4*A*V-cxav4*XAV
        dA <- pa4*V*theta+ra4*A*(1-A/amax4)+km*XAV*theta-kp4*A*V*theta-da*A
        dV <- p4*I-c4*V+km*XAV-kp4*A*V
        dE <- z + alpha4*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = NULL)
    
    finaloutput1 <- data.frame(times2I= output1[,1], 
                               T= output1[,2])
    
    finaloutput2 <- data.frame(timesR= output2[,1],
                               T= output2[,2])
    
    finaloutput3 <- data.frame(timesA= output3[,1],
                               T= output3[,2])
    
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               T= output4[,2])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9) 
    colormanual <- c("Two infected cells populations" = "dodgerblue4", 
                     "Refractory cells population"="dodgerblue3",
                     "Antibody population"="dodgerblue1",
                     "Antibody and Effector cells populations" = "lightblue")
    linemanual <- c("Two infected cells populations" = "solid", 
                    "Refractory cells population"="dotdash",
                    "Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$II == 1){
      pl=pl + geom_line(data=finaloutput1, aes(x=times2I, y = T,color = "Two infected cells populations",
                                               linetype = "Two infected cells populations"),size=1.075) 

      vdisp[1] <- TRUE}
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = T,color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075) 

      vdisp[2] <- TRUE}
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = T,color = "Antibody population",
                                               linetype="Antibody population"),size=1.075) 
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = T,color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=1.075) 
      
      vdisp[4] <- TRUE}
    pl <- pl + labs(x='Time (Days)') + labs(y='Uninfected cells (cell/mL)') +
      labs(color  = "Model", linetype = "Model")

    pl <- pl + theme(legend.position=c(.60, 0.25), legend.justification=c(0,1), legend.title=element_blank())
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
    
    r1 <- input$r1 #Proliferation rate (INPUT)
    beta1 <- input$beta1/10000000000 #Infection rate (INPUT)
    mu1 <- input$mu1/10000 #Cytolytic rate infected cells (INPUT)
    tau1 <-input$tau1 #Activation delay (INPUT)
    alpha1 <-input$alpha1/10000000 #Activation and proliferation of effector cells (INPUT)
    delta_tr <- input$delta_tr #Transition rate constant (conditional INPUT) 
    rho1 <- input$rho1 #Cure rate I1 to T (conditional INPUT)
    rho2 <- input$rho2 #Cure rate I2 to I1 (conditional INPUT)
    pI1 <- input$pI1 #Production rate constant from I1 (conditional INPUT)
    pI2 <- input$pI2 #Production rate constant from I2 (conditional INPUT)
    r2 <- input$r2 #Proliferation rate (INPUT)
    beta2 <- input$beta2/10000000000 #Infection rate (INPUT)
    mu2 <- input$mu2/10000 #Cytolytic rate infected cells (INPUT)
    tau2 <-input$tau2 #Activation delay (INPUT)
    alpha2 <-input$alpha2/10000000 #Activation and proliferation of effector cells (INPUT)
    rho <- input$rho/1000 #Cure rate I to T (conditional INPUT)
    rho_r <- input$rho_r/100000 #Cure rate R to T (conditional INPUT)
    mu_r <- input$mu_r/1000000 #Cytolytic rate constant of refractory cells (conditional INPUT)
    p2 <- input$p2 #Production rate constant from I (conditional INPUT)
    beta3 <- input$beta3/10000000000 #Infection rate (INPUT)
    p3 <- input$p3 #Production rate constant from I (conditional INPUT)
    delta3 <- input$delta3 #Death rate infected cells
    ra3 <- input$ra3 #Proliferation rate antibody
    pa3 <- input$pa3/100000 #Production rate constant antibody
    kp3 <- input$kp3/1000000000000
    cxav3 <- input$cxav3
    beta4 <- input$beta4/10000000000 #Infection rate (INPUT)
    p4 <- input$p4 #Production rate constant from I (conditional INPUT)
    delta4 <- input$delta4 #Death rate infected cells
    ra4 <- input$ra4 #Proliferation rate antibody
    pa4 <- input$pa4/1000 #Production rate constant antibody
    mu4 <- input$mu4/1000 #Cytolytic rate infected cells (INPUT)
    kp4 <- input$kp4/10000000000
    c4 <- input$c4
    cxav4 <- input$cxav4
    alpha4 <-input$alpha4/10000000 #Activation and proliferation of effector cells (INPUT)
    c <- input$c #Virus clearance rate (INPUT)
    z <-input$z #Source of effector cells (INPUT)
    de <-input$de #Effector cells clearance rate (INPUT)
    km <- input$km
    amax3 <- input$amax3*1000000000000000 
    amax4 <- input$amax4*100000000000 
    theta <- input$theta
    da <- input$da 
    tmax <- input$tmax*10000000 #Maximum liver capacity (INPUT)
    
    t <- input$num2 #Initial uninfected cell (INPUT)
    v <- input$num3 #Initial free virus (INPUT)
    e <- input$num4 #Initial effector cells (INPUT)
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = t, I1 = 0, I2 = 0, V = v, E = e)
    yinitR <- c(T = t, I = 0, R = 0, V = v, E = e)
    yinitA <- c(T = t, I = 0, XAV = 0, A = 0, V = v)
    yinitAE <- c(T = t, I = 0, XAV = 0, A = 0, V = v, E = e)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau1, I1, lagvalue(t - tau1, 2))
        lagI2 <- ifelse(t < tau1, I2, lagvalue(t - tau1, 3))
        lagE <- ifelse(t < tau1, E, lagvalue(t - tau1, 5))
        
        dT <- r1 * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta1 * T * V + rho1 * I1
        dI1 <- beta1 * T * V - (rho1 + delta_tr) * I1 - mu1 * E * I1 + rho2 * I2
        dI2 <- r1 * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu1 * E * I2
        dV <- pI1 * I1 + pI2 * I2 - c * V
        dE <- z + alpha1 * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = NULL)
    
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,5)
        
        dT <- r2*T*(1-((T+I+R)/tmax)) - beta2*T*V + rho_r*R
        dI <- r2*I*(1-((T+I+R)/tmax)) + beta2*T*V - rho*I*E - mu2*I*E 
        dR <- rho*I*E + r2*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_r*R*E
        dV <- p2*I - c*V
        dE <- z + alpha2*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = NULL)
    
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r2*T*(1-(T+I)/tmax)-beta3*T*V
        dI <- beta3*T*V-delta3*I
        dXAV <- -km*XAV+kp3*A*V-cxav3*XAV
        dA <- pa3*V*(1+theta)+ra3*A*(1-A/amax3)+km*XAV*(1+theta)-kp3*A*V*(1+theta)-da*A
        dV <- p3*I-c*V+km*XAV-kp3*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = NULL)
    
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,6)
        
        dT <- r2*T*(1-(T+I)/tmax)-beta4*T*V
        dI <- r2*I*(1-(T+I)/tmax)+beta4*T*V-delta4*I-mu4*I*E
        dXAV <- -km*XAV+kp4*A*V-cxav4*XAV
        dA <- pa4*V*theta+ra4*A*(1-A/amax4)+km*XAV*theta-kp4*A*V*theta-da*A
        dV <- p4*I-c4*V+km*XAV-kp4*A*V
        dE <- z + alpha4*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = NULL)
    
    finaloutput2 <- data.frame(timesR= output2[,1],
                               R= output2[,4])
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^0,10^3,10^6,10^9) 
    colormanual <- ("Refractory cells population"="brown")
    linemanual <- ("Refractory cells population"="dotdash")
    if (input$R == 1){
      pl=pl + geom_line(data=finaloutput2, aes(x=timesR, y = R,color = "Refractory cells population",
                                               linetype="Refractory cells population"),size=1.075) 

      vdisp[2] <- TRUE}
    pl <- pl + labs(x='Time (Days)') + labs(y='Refractory cells (cell/mL)') +
      labs(color  = "Model", linetype = "Model")

    pl <- pl + theme(legend.position=c(.60, 0.25), legend.justification=c(0,1), legend.title=element_blank())
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
    
    r1 <- input$r1 #Proliferation rate (INPUT)
    beta1 <- input$beta1/10000000000 #Infection rate (INPUT)
    mu1 <- input$mu1/10000 #Cytolytic rate infected cells (INPUT)
    tau1 <-input$tau1 #Activation delay (INPUT)
    alpha1 <-input$alpha1/10000000 #Activation and proliferation of effector cells (INPUT)
    delta_tr <- input$delta_tr #Transition rate constant (conditional INPUT) 
    rho1 <- input$rho1 #Cure rate I1 to T (conditional INPUT)
    rho2 <- input$rho2 #Cure rate I2 to I1 (conditional INPUT)
    pI1 <- input$pI1 #Production rate constant from I1 (conditional INPUT)
    pI2 <- input$pI2 #Production rate constant from I2 (conditional INPUT)
    r2 <- input$r2 #Proliferation rate (INPUT)
    beta2 <- input$beta2/10000000000 #Infection rate (INPUT)
    mu2 <- input$mu2/10000 #Cytolytic rate infected cells (INPUT)
    tau2 <-input$tau2 #Activation delay (INPUT)
    alpha2 <-input$alpha2/10000000 #Activation and proliferation of effector cells (INPUT)
    rho <- input$rho/1000 #Cure rate I to T (conditional INPUT)
    rho_r <- input$rho_r/100000 #Cure rate R to T (conditional INPUT)
    mu_r <- input$mu_r/1000000 #Cytolytic rate constant of refractory cells (conditional INPUT)
    p2 <- input$p2 #Production rate constant from I (conditional INPUT)
    beta3 <- input$beta3/10000000000 #Infection rate (INPUT)
    p3 <- input$p3 #Production rate constant from I (conditional INPUT)
    delta3 <- input$delta3 #Death rate infected cells
    ra3 <- input$ra3 #Proliferation rate antibody
    pa3 <- input$pa3/100000 #Production rate constant antibody
    kp3 <- input$kp3/1000000000000
    cxav3 <- input$cxav3
    beta4 <- input$beta4/10000000000 #Infection rate (INPUT)
    p4 <- input$p4 #Production rate constant from I (conditional INPUT)
    delta4 <- input$delta4 #Death rate infected cells
    ra4 <- input$ra4 #Proliferation rate antibody
    pa4 <- input$pa4/1000 #Production rate constant antibody
    mu4 <- input$mu4/1000 #Cytolytic rate infected cells (INPUT)
    kp4 <- input$kp4/10000000000
    c4 <- input$c4
    cxav4 <- input$cxav4
    alpha4 <-input$alpha4/10000000 #Activation and proliferation of effector cells (INPUT)
    c <- input$c #Virus clearance rate (INPUT)
    z <-input$z #Source of effector cells (INPUT)
    de <-input$de #Effector cells clearance rate (INPUT)
    km <- input$km
    amax3 <- input$amax3*1000000000000000 
    amax4 <- input$amax4*100000000000 
    theta <- input$theta
    da <- input$da 
    tmax <- input$tmax*10000000 #Maximum liver capacity (INPUT)
    
    t <- input$num2 #Initial uninfected cell (INPUT)
    v <- input$num3 #Initial free virus (INPUT)
    e <- input$num4 #Initial effector cells (INPUT)
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = t, I1 = 0, I2 = 0, V = v, E = e)
    yinitR <- c(T = t, I = 0, R = 0, V = v, E = e)
    yinitA <- c(T = t, I = 0, XAV = 0, A = 0, V = v)
    yinitAE <- c(T = t, I = 0, XAV = 0, A = 0, V = v, E = e)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau1, I1, lagvalue(t - tau1, 2))
        lagI2 <- ifelse(t < tau1, I2, lagvalue(t - tau1, 3))
        lagE <- ifelse(t < tau1, E, lagvalue(t - tau1, 5))
        
        dT <- r1 * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta1 * T * V + rho1 * I1
        dI1 <- beta1 * T * V - (rho1 + delta_tr) * I1 - mu1 * E * I1 + rho2 * I2
        dI2 <- r1 * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu1 * E * I2
        dV <- pI1 * I1 + pI2 * I2 - c * V
        dE <- z + alpha1 * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = NULL)
    
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,5)
        
        dT <- r2*T*(1-((T+I+R)/tmax)) - beta2*T*V + rho_r*R
        dI <- r2*I*(1-((T+I+R)/tmax)) + beta2*T*V - rho*I*E - mu2*I*E 
        dR <- rho*I*E + r2*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_r*R*E
        dV <- p2*I - c*V
        dE <- z + alpha2*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = NULL)
    
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r2*T*(1-(T+I)/tmax)-beta3*T*V
        dI <- beta3*T*V-delta3*I
        dXAV <- -km*XAV+kp3*A*V-cxav3*XAV
        dA <- pa3*V*(1+theta)+ra3*A*(1-A/amax3)+km*XAV*(1+theta)-kp3*A*V*(1+theta)-da*A
        dV <- p3*I-c*V+km*XAV-kp3*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = NULL)
    
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,6)
        
        dT <- r2*T*(1-(T+I)/tmax)-beta4*T*V
        dI <- r2*I*(1-(T+I)/tmax)+beta4*T*V-delta4*I-mu4*I*E
        dXAV <- -km*XAV+kp4*A*V-cxav4*XAV
        dA <- pa4*V*theta+ra4*A*(1-A/amax4)+km*XAV*theta-kp4*A*V*theta-da*A
        dV <- p4*I-c4*V+km*XAV-kp4*A*V
        dE <- z + alpha4*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = NULL)
    
    finaloutput3 <- data.frame(timesA= output3[,1],
                               A= output3[,5]*2.7*10**(-16))
    
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               A= output4[,5]*2.7*10**(-16))
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^(-10),10^(-5),10^0,10^5) 
    colormanual <- c("Antibody population"="goldenrod4",
                    "Antibody and Effector cells populations" = "gold")
    linemanual <- c("Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = A,color = "Antibody population",
                                               linetype="Antibody population"),size=1.075) 
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = A,color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=1.075) 
      
      vdisp[4] <- TRUE}
    pl <- pl + labs(x='Time (Days)') + labs(y='Free Antibody (mg/mL)') +
      labs(color  = "Model", linetype = "Model")
    
    pl <- pl + theme(legend.position=c(.60, 0.2), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(10^(-10), 10^5),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 
    
    print(pl)
    
  })
  
  output$plot7 <- renderPlot({
    
    r1 <- input$r1 #Proliferation rate (INPUT)
    beta1 <- input$beta1/10000000000 #Infection rate (INPUT)
    mu1 <- input$mu1/10000 #Cytolytic rate infected cells (INPUT)
    tau1 <-input$tau1 #Activation delay (INPUT)
    alpha1 <-input$alpha1/10000000 #Activation and proliferation of effector cells (INPUT)
    delta_tr <- input$delta_tr #Transition rate constant (conditional INPUT) 
    rho1 <- input$rho1 #Cure rate I1 to T (conditional INPUT)
    rho2 <- input$rho2 #Cure rate I2 to I1 (conditional INPUT)
    pI1 <- input$pI1 #Production rate constant from I1 (conditional INPUT)
    pI2 <- input$pI2 #Production rate constant from I2 (conditional INPUT)
    r2 <- input$r2 #Proliferation rate (INPUT)
    beta2 <- input$beta2/10000000000 #Infection rate (INPUT)
    mu2 <- input$mu2/10000 #Cytolytic rate infected cells (INPUT)
    tau2 <-input$tau2 #Activation delay (INPUT)
    alpha2 <-input$alpha2/10000000 #Activation and proliferation of effector cells (INPUT)
    rho <- input$rho/1000 #Cure rate I to T (conditional INPUT)
    rho_r <- input$rho_r/100000 #Cure rate R to T (conditional INPUT)
    mu_r <- input$mu_r/1000000 #Cytolytic rate constant of refractory cells (conditional INPUT)
    p2 <- input$p2 #Production rate constant from I (conditional INPUT)
    beta3 <- input$beta3/10000000000 #Infection rate (INPUT)
    p3 <- input$p3 #Production rate constant from I (conditional INPUT)
    delta3 <- input$delta3 #Death rate infected cells
    ra3 <- input$ra3 #Proliferation rate antibody
    pa3 <- input$pa3/100000 #Production rate constant antibody
    kp3 <- input$kp3/1000000000000
    cxav3 <- input$cxav3
    beta4 <- input$beta4/10000000000 #Infection rate (INPUT)
    p4 <- input$p4 #Production rate constant from I (conditional INPUT)
    delta4 <- input$delta4 #Death rate infected cells
    ra4 <- input$ra4 #Proliferation rate antibody
    pa4 <- input$pa4/1000 #Production rate constant antibody
    mu4 <- input$mu4/1000 #Cytolytic rate infected cells (INPUT)
    kp4 <- input$kp4/10000000000
    c4 <- input$c4
    cxav4 <- input$cxav4
    alpha4 <-input$alpha4/10000000 #Activation and proliferation of effector cells (INPUT)
    c <- input$c #Virus clearance rate (INPUT)
    z <-input$z #Source of effector cells (INPUT)
    de <-input$de #Effector cells clearance rate (INPUT)
    km <- input$km
    amax3 <- input$amax3*1000000000000000 
    amax4 <- input$amax4*100000000000 
    theta <- input$theta
    da <- input$da 
    tmax <- input$tmax*10000000 #Maximum liver capacity (INPUT)
    
    t <- input$num2 #Initial uninfected cell (INPUT)
    v <- input$num3 #Initial free virus (INPUT)
    e <- input$num4 #Initial effector cells (INPUT)
    
    w <- input$num5 #Duration of the Model (INPUT)
    
    
    yinit2I <- c(T = t, I1 = 0, I2 = 0, V = v, E = e)
    yinitR <- c(T = t, I = 0, R = 0, V = v, E = e)
    yinitA <- c(T = t, I = 0, XAV = 0, A = 0, V = v)
    yinitAE <- c(T = t, I = 0, XAV = 0, A = 0, V = v, E = e)
    
    Ciupe2I <- function(t, y, parms) {
      with(as.list(c(y, parms)), {
        lagI1 <- ifelse(t < tau1, I1, lagvalue(t - tau1, 2))
        lagI2 <- ifelse(t < tau1, I2, lagvalue(t - tau1, 3))
        lagE <- ifelse(t < tau1, E, lagvalue(t - tau1, 5))
        
        dT <- r1 * (T + I1) * (1 - ((T + I1 + I2) / tmax)) - beta1 * T * V + rho1 * I1
        dI1 <- beta1 * T * V - (rho1 + delta_tr) * I1 - mu1 * E * I1 + rho2 * I2
        dI2 <- r1 * I2 * (1 - ((T + I1 + I2) / tmax)) + delta_tr * I1 - rho2 * I2 - mu1 * E * I2
        dV <- pI1 * I1 + pI2 * I2 - c * V
        dE <- z + alpha1 * (lagI1 + lagI2) * lagE - de * E
        return(list(c(dT, dI1, dI2, dV, dE)))
      })
    } 
    #The differential equation that models the dynamics with 2 infected cells populations
    
    times <- seq(0, w, by = 0.1)
    output1 <- dede(y = yinit2I, times = times, func = Ciupe2I, parms = NULL)
    
    
    CiupeR = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,5)
        
        dT <- r2*T*(1-((T+I+R)/tmax)) - beta2*T*V + rho_r*R
        dI <- r2*I*(1-((T+I+R)/tmax)) + beta2*T*V - rho*I*E - mu2*I*E 
        dR <- rho*I*E + r2*R*(1-((T+I+R)/tmax)) - rho_r*R - mu_r*R*E
        dV <- p2*I - c*V
        dE <- z + alpha2*lagI*lagE - de*E
        return(list(c(dT, dI, dR, dV, dE)))
      })
    }
    #The differential equation that models the dynamics with Refractory cells populations
    
    output2 <- dede(y = yinitR, times = times, func = CiupeR, parms = NULL)
    
    
    CiupeA = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        dT <- r2*T*(1-(T+I)/tmax)-beta3*T*V
        dI <- beta3*T*V-delta3*I
        dXAV <- -km*XAV+kp3*A*V-cxav3*XAV
        dA <- pa3*V*(1+theta)+ra3*A*(1-A/amax3)+km*XAV*(1+theta)-kp3*A*V*(1+theta)-da*A
        dV <- p3*I-c*V+km*XAV-kp3*A*V
        return(list(c(dT, dI, dXAV, dA, dV)))
      })
    }
    
    output3 <- dede(y = yinitA, times = times, func = CiupeA, parms = NULL)
    
    
    CiupeAE = function(t, y, parms) {
      with(as.list(c(y, parms)), {
        if (t < tau2)
          lagI <- I
        else
          lagI <- lagvalue(t - tau2,2)
        if (t < tau2)
          lagE <- E
        else
          lagE <- lagvalue(t - tau2,6)
        
        dT <- r2*T*(1-(T+I)/tmax)-beta4*T*V
        dI <- r2*I*(1-(T+I)/tmax)+beta4*T*V-delta4*I-mu4*I*E
        dXAV <- -km*XAV+kp4*A*V-cxav4*XAV
        dA <- pa4*V*theta+ra4*A*(1-A/amax4)+km*XAV*theta-kp4*A*V*theta-da*A
        dV <- p4*I-c4*V+km*XAV-kp4*A*V
        dE <- z + alpha4*lagI*lagE - de*E
        return(list(c(dT, dI, dXAV, dA, dV, dE)))
      })
    }
    
    output4 <- dede(y = yinitAE, times = times, func = CiupeAE, parms = NULL)
    
    finaloutput3 <- data.frame(timesA= output3[,1],
                               XAV= output3[,4]*2.7*10**(-16))
    
    finaloutput4 <- data.frame(timesAE= output4[,1],
                               XAV= output4[,4]*2.7*10**(-16))
    #Converting both lists' output into a dataframe for ggplot
    
    vdisp <- rep(FALSE,2)
    pl=ggplotmlx()
    brk_y <-c(10^(-15),10^(-10),10^(-5),10^0) 
    colormanual <- c("Antibody population"="chocolate3",
                     "Antibody and Effector cells populations" = "chocolate1")
    linemanual <- c("Antibody population"="dashed",
                    "Antibody and Effector cells populations" = "dotted")
    if (input$A == 1){
      pl=pl + geom_line(data=finaloutput3, aes(x=timesA, y = XAV,color = "Antibody population",
                                               linetype="Antibody population"),size=1.075) 
      
      vdisp[3] <- TRUE}
    if (input$AE == 1){
      pl=pl + geom_line(data=finaloutput4, aes(x=timesAE, y = XAV,color = "Antibody and Effector cells populations",
                                               linetype="Antibody and Effector cells populations"),size=1.075) 
      
      vdisp[4] <- TRUE}
    pl <- pl + labs(x='Time (Days)') + labs(y='Virus-antibody complexes (mg/mL)') +
      labs(color  = "Model", linetype = "Model")
    
    pl <- pl + theme(legend.position=c(.60, 0.2), legend.justification=c(0,1), legend.title=element_blank())
    pl <- pl + scale_y_continuous(
      trans = 'log10', limits = c(10^(-15), 10^0),
      breaks = brk_y,
      labels = trans_format("log10", math_format(10^.x)),
      expand = c(0, 0)
    ) + scale_x_continuous(limits = c(0, 301), expand = c(0,0)) +
      scale_color_manual(values = colormanual) +
      scale_linetype_manual(values = linemanual) 
    
    print(pl)
    
  })
  
})