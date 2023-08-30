
require(deSolve)
require(dde)
require(phaseR)
require(ggplot2)

y <- c(T=300, I=0, VI=5, VU=0) #variables initial

#parameters
parms <- c(0.1# nT
          ,0.63 #n
          ,0.0017 #beta1
          ,0.0028 #beta2
          ,0.053 #u1 (mu1)
          ,0.01 #tau 1
          ,0.01 #tau 2
          ,0.1 #tau 3
          ,0.041 # u2 (mu2)
          ,0.074 # u3 (mu3)
          ,1.23 #b
          ,0.73) # eta


#time range 
times <- seq(0 , 1000 , 1)
R0 = 0
# i , r , a , h , u, d , al , e

model <- function(t, y, parms) {
    
        T <- y[1]
        I <- y[2]
        VI <- y[3]
        VU <- y[4]


        
        nT <- 1.5 - (parms[1]*T) # SHOULD BE MADE T DEPENDENT
        n <- parms[2]
        beta1 <- parms[3]  #Transmission rate from  virus
        beta2 <- parms[4]  #Transmission rate from infected cells
        u1 <- parms[5] # death rate of infected cells
        tau1 <- parms[6] #time delay 1
        tau2 <- parms[7] #time delay 2
        tau3 <- parms[8] #time delay 3
        u2 <- parms[9] #
        u3 <- parms[10]
        b <- parms[11]# virus production rate
        e <- parms[12]

        R0 <<- (exp(-u1 * tau2) * beta2 * 15 / u1) + ((b * exp(-u2 * tau3)*(1-e)* exp(-u1 * tau1) * (1-n) * beta1*15)/(u1 * u3 ))
        
        print(R0)
       
        if (t < tau1) {
            T_lag <- y[1]
        } else {
            T_lag <- lagvalue(t - tau1, 1)
        }
    
        if (t < tau1) {
            VI_lag <- y[3]
        } else {
            VI_lag <- lagvalue(t - tau1, 3)
        }

        if (t < tau2) {
            T_lag2 <- y[1]
        } else {
            T_lag2 <- lagvalue(t - tau2, 1)
        }

        if (t < tau2) {
            I_lag <- y[2]
        } else {
            I_lag <- lagvalue(t - tau2, 2)
        }

        if (t < tau3) {
            I_lag2 <- y[2]
        } else {
            I_lag2 <- lagvalue(t - tau3, 2)
        }



        #expression fr R0
        R0 <<- (exp(-u1 * tau2) * beta2 * 15 / u1) + ((b * exp(-u2 * tau3)*(1-e)* exp(-u1 * tau1) * (1-n) * beta1*15)/(u1 * u3 ))

        
        #model
        dTdt <- nT - ((1 - n) * beta1 * T * VI ) - (beta2 * T * I)    

        dIdt <- (exp(-u1 * tau1) * (1 - n) * beta1 * T_lag * VI_lag) + ((exp(-u1 * tau2)) * beta2 * T_lag2 * I_lag) - (u1 * I)

        dVIdt <- ((1 - e) * b * exp(-u2 * tau3) * I_lag2) - u3 * VI

        dVUdt <- (e * b * exp(-u2 * tau3) * I_lag2) - (u3 * VU)


    list(c(dTdt,dIdt,dVIdt,dVUdt))
}

#out <- dede(y, times = times, func = model,parms)

out <- c(0,0)
for (x in 1:100){
        rval <- (exp(-parms[5]*parms[7]) * parms[4] * 15 / parms[5]) + (parms[11] * exp(-parms[9] * parms[8]) *(1 - (x/100)) *exp(-parms[5]* parms[6])* (1- parms[2]) * parms[3] * 15)/(parms[5] * parms[10])
        out <<- rbind(out,c(x/100,rval))
    }

out = as.data.frame(out)
head(round(out,1))
out <- out[-c(1),]
print(out)
plot(
    x = out$"V1", y = out$"V2", col = "blue", ylab = "R0",type="l",
    xlab = expression(eta), xlim=c(0,1), ylim=c(0,4)) #T
  #lines(x = out$time, y = out$"V2", col = "green") # I
  #lines(x = out$time, y = out$"VI", col = "red")  # VI
  #lines(x = out$time, y = out$"VU", col = "blue") # VU
  #legend(300,600,legend = c("T","I","V_I","V_U"),col=c("black","green","red","blue"),lty=1:2, cex=0.8)

