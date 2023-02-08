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
  if (sex=="b"&mx[1]>=0.07599) {a0 = 0.30663}   
  
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