# =========== Predicting mortality by age and sex for target years ==========================================

# Function for computation of a complete life table FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF_BEGIN

LTabC = function(mx, sex, popname)  
  
{
  Nx = length(mx)   # number of ages
  
  Sexx = rep(sex, Nx)  # Column specifying gender
  
  Popx = rep(popname, Nx)   # Column specifying 
  
  x = as.numeric()
  for (i in 1:Nx) {x[i] = i-1 }  # vector of ages
  nx = as.numeric()
  nx = rep(1, Nx)  # widths of age intervals
  
  if (sex=="m"&mx[1]<0.02300) {a0 = 0.1429 - 1.99545*mx[1]}
  if (sex=="m"&mx[1] >= 0.02300&mx[1]<0.08307) {a0 = 0.02832 + 3.26021*mx[1]} 
  if (sex=="m"&mx[1]>=0.08307) {a0 = 0.29915}
  
  if (sex=="f"&mx[1]<0.01724) {a0 = 0.14903 - 2.05527*mx[1]}  
  if (sex=="f"&mx[1] >= 0.01724&mx[1]<0.06891) {a0 = 0.04667 + 3.88089*mx[1]}
  if (sex=="f"&mx[1]>=0.06891) {a0 = 0.31411}   
  
  if (sex=="b"&mx[1]<0.02012) {a0 = 0.145965 - 2.02536*mx[1]}
  if (sex=="b"&mx[1] >= 0.02012&mx[1]<0.07599) {a0 = 0.037495 + 3.57055*mx[1]}
  if (mx[1]>=0.07599) {a0 = 0.30663}   
  
  ax0 = c(a0)
  ax1 = rep(0.5, Nx-1)
  ax = append(ax0, ax1)     # ax values
  
  qx = (nx*mx)/(1+(1-ax)*nx*mx)   # probabilities of dying
  qx[Nx] = 1
  ax[Nx] = 1/mx[Nx]
  
  lx = numeric()
  lx[1] = 100000  # Survival function lx
  for (j in 2:Nx) { lx[j] = lx[j-1]*(1-qx[j-1]) }                  
  
  dx = lx*qx     # Numbers of dying
  
  Lx = lx*nx - dx*(1-ax)*nx   # Person-years liven withing age group x
  Lx[Nx] = lx[Nx]*(1/mx[Nx])
  
  Tx = numeric()   # Person-years lived after age x
  Tx[Nx] = Lx[Nx]
  for (j in seq(from=Nx-1, to=1, by=-1)) { Tx[j] = Tx[j+1]+Lx[j] }
  
  ex = Tx / lx   # Life expectancy at age x
  
  # LTcomplete <- data.frame(cbind(Popx, Sexx, x, nx, ax, mx, qx, lx, dx, Lx, Tx, ex) )
  
  LTcomplete = data.frame(Popx, Sexx, x, nx, ax, mx, qx, lx, dx, Lx, Tx, ex)
  # LTcomplete = as.data.frame(LTcomplete)
  row.names(LTcomplete) = NULL
  
  for (ii in 3:12) { 
    LTcomplete[, ii] = as.numeric(LTcomplete[, ii])
                 }

  return(LTcomplete)
}
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF_END


# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF_BEGIN
# Function for prediction of annual mortality for a prediction period for one sex 
# and one country

Predict_Mort_LC = function (PathData, CNTR, s, Mortfile, PopFile, Yretro1, Yretro2, Ypredict1, Ypredict2, maxAge)
  
 {
  # Parameters of the function:
  # PathData - Path to data files ending with "/" 
  # CNTR - country abbreviation (as in the HMD)
  # s - sex 1-f  2-m
  # Mortfile - name of the mortality file 1x1
  # Popfile - name of the population or pop exposure file 1x1
  # Yretro1 - first year of the retrospective period
  # Yretro2 - last year of the retrospective period
  # Ypredict1 - first year of the target period
  # Ypredict2 - last year of the target period
  # maxAge - maximal age for prediction and all mortality and life table columns

  library(demography)

  options(scipen=999) #Forcing R not to use scientific format when printing numbers

  
  # ---- Reading data files ---------------------------------------------------

  Mortfile = paste(PathData, CNTR, ".", Mortfile, sep="") 
  M_xt = read.csv(Mortfile, header=T, sep="", skip=2, na.strings=".")
  head(M_xt)
  M_xt$Age[M_xt$Age == "110+"] = 110
  M_xt$Age = as.numeric(M_xt$Age)

  Popfile = paste(PathData, CNTR, ".", Popfile, sep="") 
  Pop_xt = read.csv(Popfile, header=T, sep="", skip=2, na.strings=".")
  head(Pop_xt)
  Pop_xt$Age[Pop_xt$Age == "110+"] = 110
  Pop_xt$Age = as.numeric(Pop_xt$Age)

  # FILENAME = Fertfile
  # Fert_xt = read.csv(FILENAME, header=T, sep="", skip=2, na.strings=".")
  # head(Fert_xt)
  # Fert_xt$Age[Fert_xt$Age == "12-"] = 12
  # Fert_xt$Age[Fert_xt$Age == "55+"] = 55
  # Fert_xt$Age = as.numeric(Fert_xt$Age)


  # ---- Creation of demogr object for its use in "Demography" ---------------- 

  dir = paste(PathData, CNTR, "/", sep="")
  MORT_cntr = read.demogdata(Mortfile,Popfile, type="mortality", label=CNTR)

  # ----- L-C mortality modeling for the period Yretro1:Yretro2 ---------------

  lca_mort = lca(MORT_cntr,series=names(MORT_cntr$rate)[s],years=Yretro1:Yretro2,ages=MORT_cntr$age,
             adjust="e0",interpolate=TRUE,max.age=maxAge)

  # ----- Controlling mort fit for the last retrospective year ----------------
  
  Mx_obs = numeric()
  Px_obs = numeric()

  if (s==1) Mx_obs=M_xt$Female[M_xt$Year==Yretro2] else Mx_obs=M_xt$Male[M_xt$Year==Yretro2]
  if (s==1) Px_obs=Pop_xt$Female[Pop_xt$Year==Yretro2] else Px_obs=Pop_xt$Male[Pop_xt$Year==Yretro2]

  # Mx_Obs recalculation for the last (open-ended) age group 
  
  Mx_obs100 = tail(Mx_obs,111-maxAge)
  Px_obs100 = tail(Px_obs,111-maxAge)
  Mx_obs100p = sum(Mx_obs100*Px_obs100)/sum(Px_obs100)
  Px_obs100p = sum(Px_obs100)
  ratio100 = Mx_obs100p/Mx_obs[maxAge+1]
  Mx_obs[maxAge+1]=Mx_obs100p
  Px_obs[maxAge+1]=Px_obs100p
  
  Mx_obs = head(Mx_obs, maxAge+1)

  Mx_fit_lastY = numeric()
  mx_u = numeric()
  mx_l = numeric()
  mx_s = numeric()
  Diffx = numeric()
  Diffx2 = numeric()
  RelDiffx = numeric()
  dxx_obs = numeric()
  dxx_fit = numeric()

  # Model death rates for the last retro-year
  
  for (x_age in 0:maxAge)  {
    i=x_age+1
    Mx_fit_lastY[i]=exp(lca_mort$ax[i]+lca_mort$bx[i]*lca_mort$kt[Yretro2-Yretro1+1])
                         }
    
  # Mx_fit_lastY recalculation for the last (open-ended) age group 
  Mx_fit_lastY[maxAge+1] = Mx_fit_lastY[maxAge+1]*ratio100
  
  # Calculation of complete life tables from observed and model death rates for the last retro-year
  
  if (s==1) LTcomplete_Obs = LTabC(Mx_obs,"f",CNTR) else LTcomplete_Obs = LTabC(Mx_obs,"m",CNTR)
  if (s==1) LTcomplete_Fit = LTabC(Mx_fit_lastY,"f",CNTR) else LTcomplete_Fit = LTabC(Mx_fit_lastY,"m",CNTR)
  
  # Calculation of dx*xav for the observed and fitted life tables
  
  dxx_obs = LTcomplete_Obs$dx*(LTcomplete_Obs$x+LTcomplete_Obs$ax)/100000
  dxx_fit = LTcomplete_Fit$dx*(LTcomplete_Fit$x+LTcomplete_Fit$ax)/100000
  
  cat("  *********** The last retro-year = ", Yretro2, "\n","\n")
  cat("Checking dx*x  fit", "\n")
  cat("  Age          dxx_obs         dxx_LC           Diff       ","\n","\n")
  
  for (x_age in 0:maxAge)   { 
    i = x_age + 1
    Diffx[i] = dxx_obs[i] - dxx_fit[i]
    Diffx2[i] = Diffx[i]^2
    
    cat(sprintf("%5.0f",x_age), "    ")
    cat(sprintf("%10.5f",dxx_obs[i]),"    ")
    cat(sprintf("%10.5f",dxx_fit[i]),"    ")
    cat(sprintf("%10.5f",Diffx[i]),"    ", "\n")
                            }
                           
  cat("***********  Mean RMSD =  ", sprintf("%9.6f", sqrt(mean(Diffx2, na.rm=TRUE))), "\n")
  
  cat("***********   e0_Obs =  ", sprintf("%8.3f", LTcomplete_Obs$ex[1]),
      "       e0_fit =  ", sprintf("%8.3f", LTcomplete_Fit$ex[1]),"\n")

  # ---- Storing mean(Mx-obs), mean(Mx-LC), mean deviations and mean abs deviations
  Mean_dxx_obs = mean(head(dxx_obs, maxAge+1), na.rm=TRUE)
  Mean_dxx_fit = mean(head(dxx_fit, maxAge+1), na.rm=TRUE)
  RMSD = sqrt(mean(Diffx2, na.rm=TRUE))
  e0_obs  = LTcomplete_Obs$ex[1]
  e0_fit  = LTcomplete_Fit$ex[1]
  e15_obs = LTcomplete_Obs$ex[16]
  e15_fit = LTcomplete_Fit$ex[16]
  e60_obs = LTcomplete_Obs$ex[61]
  e60_fit = LTcomplete_Fit$ex[61]
  e80_obs = LTcomplete_Obs$ex[81]
  e80_fit = LTcomplete_Fit$ex[81]
  e90_obs = LTcomplete_Obs$ex[91]
  e90_fit = LTcomplete_Fit$ex[91]
  
  
  # --- Building a data frame for deviations between observed and model 
  # dx*x and e0 values-----------
  
  if (s==1) SEX = "f" else SEX = "m"   
  
  Deviations = data.frame( CNTR, SEX, Yretro1, Yretro2, Mean_dxx_obs,  Mean_dxx_fit, 
                          RMSD, e0_obs, e0_fit, e15_obs, e15_fit, e60_obs, e60_fit,
                          e80_obs, e80_fit, e90_obs, e90_fit )
  
  

  # ----- Forecasting mortality for the period Ypredict1:Ypredict2 ------------

  MortLCpredict = forecast(lca_mort, h=Ypredict2-Ypredict1+1, se=c("innovdrift"), 
                             jumpchoice=c("fit"), level=95)

  yyy=Ypredict1

  if (s==1) LTcompl_Predict0 = LTabC(MortLCpredict$rate$female[, yyy-Ypredict1+1],"f",CNTR) else 
    LTcompl_Predict0 = LTabC(MortLCpredict$rate$male[, yyy-Ypredict1+1],"m",CNTR)
    
  YEAR = rep(yyy, maxAge+1)
  LTcompl_Predict0 = cbind(LTcompl_Predict0, YEAR)
  LTcompl_Predict0 = LTcompl_Predict0[, c(1,13,2:12)] # Reorder columns YEAR-->second col
  mx_l = MortLCpredict$rate$lower[, yyy-Ypredict1+1]
  mx_u = MortLCpredict$rate$upper[, yyy-Ypredict1+1]
  mx_s = (mx_u - mx_l)/3.92
  LTcompl_Predict0 = cbind(LTcompl_Predict0, mx_l, mx_u, mx_s)

  for (yyyy in (Ypredict1+1):Ypredict2) {
  
    if (s==1) LTcompl_Predict1 = LTabC(MortLCpredict$rate$female[, yyyy-Ypredict1+1],"f", CNTR) else 
      LTcompl_Predict1 = LTabC(MortLCpredict$rate$male[, yyyy-Ypredict1+1],"m",CNTR)
  
    YEAR = rep(yyyy, maxAge+1)
    LTcompl_Predict1 = cbind(LTcompl_Predict1, YEAR)
    LTcompl_Predict1 = LTcompl_Predict1[, c(1,13,2:12)] # Reorder columns YEAR-->second col
    mx_l = MortLCpredict$rate$lower[, yyy-Ypredict1+1]
    mx_u = MortLCpredict$rate$upper[, yyy-Ypredict1+1]
    mx_s = (mx_u - mx_l)/3.92
    LTcompl_Predict1 = cbind(LTcompl_Predict1, mx_l, mx_u, mx_s)
  
    LTcompl_Predict0 = rbind(LTcompl_Predict0, LTcompl_Predict1)
                                        }
    LT_PREDICT = list(LTcompl_Predict0, Deviations)  
  
    return(LT_PREDICT) 
    
  }    
   
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF_END                                         }

# MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM_BEGIN


# Defining the current folder (with the present R-code in it) as the working directory

library(rstudioapi)
Cur_dir = dirname(getSourceEditorContext()$path)
setwd(Cur_dir)
getwd()

Mx_Px_data = read.csv("HMD_Mx_Px_data.csv", header=T,sep=",", na.strings=".")
# head(Mx_Px_data, 10)

Npops = length(Mx_Px_data$country)

Country = Mx_Px_data$country[1]  
CNTR= Mx_Px_data$code[1]
s   = Mx_Px_data$S[1]
Yretro1   = Mx_Px_data$yretro1[1]
Yretro2   = Mx_Px_data$yretro2[1]
Ypredict1 = Mx_Px_data$ypredict1[1]
Ypredict2 = Mx_Px_data$ypredict2[1]
Popfile   = Mx_Px_data$popfile[1]
Mortfile  = Mx_Px_data$mortfile[1]
PathData  = Mx_Px_data$pathdata[1]

cat("\n", "*********", "  country = ", Country, "   ", CNTR, "Sex = ", s, "\n", "\n" )
cat("*********", "Yretro1 = ", Yretro1, "Yretro2 = ", Yretro2, "Ypredict1 = ", Ypredict1,
    "Ypredict2 = ", Ypredict2, "\n", "\n") 

LTpredict_and_Devs0 = Predict_Mort_LC(PathData,CNTR,s,Mortfile,PopFile,Yretro1,Yretro2,Ypredict1,Ypredict2,
                                      maxAge=100)
LT0 =  LTpredict_and_Devs0[[1]]
DEV0 = LTpredict_and_Devs0[[2]]

# To do:   Cycle across countries  and formatting of LT and DEV

if (Npops>1) {
  
  for (j in (2:Npops)) {
    
    Country = Mx_Px_data$country[j]  
    CNTR= Mx_Px_data$code[j]
    s = Mx_Px_data$S[j]
    Yretro1 = Mx_Px_data$yretro1[j]
    Yretro2 = Mx_Px_data$yretro2[j]
    Ypredict1 = Mx_Px_data$ypredict1[j]
    Ypredict2 = Mx_Px_data$ypredict2[j]
    Popfile =  Mx_Px_data$popfile[j]
    Mortfile = Mx_Px_data$mortfile[j]
    PathData = Mx_Px_data$pathdata[j]
    
    cat("\n", "*********", "  country = ", Country, "   ", CNTR, "Sex = ", s, "\n", "\n" )
    cat("*********", "Yretro1 = ", Yretro1, "Yretro2 = ", Yretro2, "Ypredict1 = ", Ypredict1,
        "Ypredict2 = ", Ypredict2, "\n", "\n") 
    
    LTpredict_and_Devs0 = Predict_Mort_LC(PathData,CNTR,s,Mortfile,PopFile,Yretro1,Yretro2,Ypredict1,Ypredict2,
                                          maxAge=100)
    LT1 = LTpredict_and_Devs0[[1]]
    DEV1 = LTpredict_and_Devs0[[2]]  
    
    LT0 = rbind(LT0, LT1) 
    DEV0 = rbind(DEV0, DEV1)
                        }  
  
                 }

for (ii in 3:17) { 
  DEV0[, ii] = as.numeric(DEV0[, ii])
}


DEV0$Mean_dxx_obs = round(DEV0$Mean_dxx_obs, digits=4) 
DEV0$Mean_dxx_fit = round(DEV0$Mean_dxx_fit, digits=4)
DEV0$RMSD = round(DEV0$RMSD, digits=7)
DEV0$e0_obs = round(DEV0$e0_obs, digits=3)
DEV0$e0_fit = round(DEV0$e0_fit, digits=3)
DEV0$e15_obs = round(DEV0$e15_obs, digits=3)
DEV0$e15_fit = round(DEV0$e15_fit, digits=3)
DEV0$e60_obs = round(DEV0$e60_obs, digits=3)
DEV0$e60_fit = round(DEV0$e60_fit, digits=3)
DEV0$e80_obs = round(DEV0$e80_obs, digits=3)
DEV0$e80_fit = round(DEV0$e80_fit, digits=3)
DEV0$e90_obs = round(DEV0$e90_obs, digits=3)
DEV0$e90_fit = round(DEV0$e90_fit, digits=3)

write.csv(DEV0,"Dev4-ex-retro2000-19_2019m.csv")


for (ii in 6:16) { 
  LT0[, ii] = as.numeric(LT0[, ii])
                  }
LT0$ax = round(LT0$ax, digits=2) 
LT0$mx = round(LT0$mx, digits=7) 
LT0$qx = round(LT0$qx, digits=7) 
LT0$lx = round(LT0$lx, digits=0)
LT0$dx = round(LT0$dx, digits=0) 
LT0$Lx = round(LT0$Lx, digits=0) 
LT0$Tx = round(LT0$Tx, digits=0) 
LT0$ex = round(LT0$ex, digits=3)
LT0$mx_l = round(LT0$mx_l, digits=7)
LT0$mx_u = round(LT0$mx_u, digits=7)
LT0$mx_s = round(LT0$mx_s, digits=7)

write.csv(LT0,"LTabs4-retro2000-19_2020-21m.csv")

# MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM_END