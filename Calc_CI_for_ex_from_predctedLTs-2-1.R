# Computation of CIs for life expectancy from .csv files like "LTabs4-retro2015-19_2020-21m.csv" MMMMMMMMMMMMM_BEGIN
# Vladimir M. Shkolnikov, 07-02-2023


source(LTabC.R) #Including a function for computation of a complete life table 

# Getting full path to the current folder. Working Directory <-- current folder
library(rstudioapi)
Cur_dir = dirname(getSourceEditorContext()$path)
setwd(Cur_dir)
getwd()

Data1 = read.csv("LTabs4-retro2000-19_2020-21m.csv", header=T, sep=",", na=".")

Data = Data1[, c(2:17) ]       # Dropping the first column
Data = Data[Data$YEAR>=2020, ] # Dropping years before 2020-21
Data = Data[Data$Popx!="ISR", ] # Dropping data for Israel

Npops = length(unique(Data$Popx))
ppp = unique(Data$Popx)
Nages = length(unique(Data$x))
Nyears = length(unique(Data$YEAR))
yyy = unique(Data$YEAR)
NumLTs = Npops*Nyears


Niter = 2000
EEE = matrix(nrow = Nages, ncol = Niter, byrow = FALSE, dimnames = NULL)
MMx = as.numeric()
EElo = as.numeric()
EEhi = as.numeric()
Mx = as.numeric()
SMx = as.numeric()

# LTfinal = data.frame(Popx=character(), YEAR=numeric(), Sexx=character(), x=numeric(), 
#                     nx=numeric(), ax=numeric(), mx=numeric(), qx=numeric(), lx=numeric(),
#                     dx=numeric(), Lx=numeric(), Tx=numeric(), ex=numeric(), mx_l=numeric(),
#                     mx_u=numeric(), mx_s=numeric(), EE=numeric(), EElo=numeric(),
#                     EEhi=numeric())


ipop = 1  # Calculations for the first life table
  Nline = (ipop - 1)*Nages+1
  Cur_Country = Data$Popx[Nline]
  Cur_Year = Data$YEAR[Nline]
  Cur_Gender = Data$Sexx[Nline]
  cat("\n","******* ", "  country = ", Cur_Country, "  sex = ", Cur_Gender, "  year = ", Cur_Year, "\n", "\n")
  
  LTcurr0 = Data[Data$Popx==Cur_Country&Data$YEAR==Cur_Year, ]
  Mx =  LTcurr0$mx
  SMx = LTcurr0$mx_s
  
  for (k in 1:Niter)   { # Monte-Carlo simulations 
    for (ia in 1:Nages) { 
      if (SMx[ia]< 0) {SMx[ia]=SMx[ia]=0.0000001 }
      if (abs(SMx[ia])<0.0000001) {SMx[ia]=0.0000001}
      MMx[ia] = rnorm(1, Mx[ia], SMx[ia])
      if (MMx[ia]<0) {MMx[ia] = Mx[ia]}
      if (is.na(MMx[ia])==TRUE)  {cat("-------!!!!!!  NA value of MMx at age group   ", ia, "\n")}
                         } 
    LTsimul = LTabC(MMx, Cur_Gender, Cur_Country)
    for (ia in 1:Nages) { EEE[ia, k] = LTsimul$ex[ia]  }
                       } # Simulations
  EE = rowMeans(EEE)
  for (ia in 1:Nages) { EElo[ia] = quantile(EEE[ia, ], probs=0.025) }
  for (ia in 1:Nages) { EEhi[ia] = quantile(EEE[ia, ], probs=0.975) }
  
  LTcurr0 = cbind(LTcurr0, EE, EElo, EEhi) 
  
  for (ipop in 2:NumLTs) { # Cycle across populations ppppppppppppppppppppppppppppp
    Nline = (ipop - 1)*Nages+1
    Cur_Country = Data$Popx[Nline]
    Cur_Year = Data$YEAR[Nline]
    Cur_Gender = Data$Sexx[Nline]
    cat("\n","******* ", "  country = ", Cur_Country, "  sex = ", Cur_Gender, "  year = ", Cur_Year, "\n", "\n")
    
    LTcurr1 = Data[Data$Popx==Cur_Country&Data$YEAR==Cur_Year, ]
    Mx =  LTcurr1$mx
    SMx = LTcurr1$mx_s
    
    for (k in 1:Niter)   { # Monte-Carlo simulations 
      for (ia in 1:Nages) { 
        if (SMx[ia]< 0) {SMx[ia]=SMx[ia]=0.0000001 }
        if (abs(SMx[ia])<0.0000001) {SMx[ia]=0.0000001}
        MMx[ia] = rnorm(1, Mx[ia], SMx[ia])
        if (MMx[ia]<0) {MMx[ia] = Mx[ia]}
        if (is.na(MMx[ia])==TRUE)  {cat("-------!!!!!!  NA value of MMx at age group   ", ia, "\n")}
                          }
      LTsimul = LTabC(MMx, Cur_Gender, Cur_Country)
      for (ia in 1:Nages) { EEE[ia, k] = LTsimul$ex[ia]  }
                        } # Simulations 
    EE = rowMeans(EEE)
    for (ia in 1:Nages) { EElo[ia] = quantile(EEE[ia, ], probs=0.025) }
    for (ia in 1:Nages) { EEhi[ia] = quantile(EEE[ia, ], probs=0.975) }
    
    LTcurr1 = cbind(LTcurr1, EE, EElo, EEhi)
    
    LTcurr0 = rbind(LTcurr0, LTcurr1)
    
                                               } # pppppppppppppppppppppppppppppppp
  
write.csv(LTcurr0,"LTabs5-retro2000-19_2020-21m.csv")
  
  # MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM_END
