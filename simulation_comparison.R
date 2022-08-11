library(expm)
library(MASS)
library(matrixcalc)
library(tidyr)

rm(list = ls()) 

theme_paper <- theme(
  panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"),
  panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "white"),
  panel.grid.minor = element_line(size = 0.01, linetype = 'solid',colour = "white"),
  plot.title = element_text(size = 20, hjust=0.5),
  axis.text = element_text(size = 14),
  axis.title=element_text(size=18),
  legend.key.size = unit(0.8, "cm"),
  legend.key.width = unit(0.8,"cm"),
  legend.title = element_text(color = "black", size = 18),
  legend.text = element_text(color = "black", size = 14),
  legend.key = element_rect(fill = NA)
)

### load simulation data 
Nr = 100 
load(here::here("data", "N=100", "y.sim.RData"))
MAE1 <- numeric(5)
MAE2 <- numeric(5)
MAE3 <- numeric(5)
for (K in c(10, 14, 20)){
  model.name = paste(c("fit0.K", as.character(K), "RData"), collapse =".")
  load(here::here("data", "simulation", model.name))
  F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
  load(here::here("data", "simulation", F.name))
  F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
  for (i in 11:20){
    tempt = matrix(F.A%*%Fit0.KF$m.flt[[i+1]], Nr, Nr)[1:100, 96:100] - matrix(y.sim[[i]], Nr, Nr)[1:100, 96:100]
    if (K == 10){
      MAE1[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
    } else if(K == 14) {
      MAE2[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
    } else if(K == 20) {
      MAE3[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
    }
  }
}
MAE4 <- numeric(5)
load(here::here("data", "simulation", "Fit.KF.f.RData"))
load(here::here("data", "simulation", "F.f.Nr.200.K.f.20.RData"))
F.f.A <- cbind(F.f, matrix(0, nrow(F.f), ncol(F.f)))
for (i in 11:20) {
  tempt = matrix(F.f.A%*%fit.KF.f$m.flt[[i+1]], 2*Nr, 2*Nr)[1:100, 96:100] - matrix(y.sim[[i]], Nr, Nr)[1:100, 96:100]
  MAE4[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
}

data = data.frame(
  time = as.factor(c(11:20)),
  NF100 = MAE1,
  NF196 = MAE2,
  NF400 = MAE3,
  F400 = MAE4
) %>% 
  gather("Model", "MAE", 2:5)

ggplot(data, aes(x=time, y=MAE, group=Model, shape=Model, color=Model)) +
  geom_line(aes(linetype=Model)) + 
  geom_point(size = 3, aes(shape = Model)) +
  theme_paper 


