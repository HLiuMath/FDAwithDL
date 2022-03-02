# 2022.03.02
# local constant (approximation) and local linear (approximation) mean function estimation
# for functional data with detection limit

rm(list = ls())
setwd("C:\\Users\\stahl\\OneDrive - University of Leeds\\research\\Energy")

############# generate dataset
# Nsimul<-400 # number of simulation replications
nn <- 200 # number of subjects

set.seed(1)
# number of observations for dense scenario: N_i are iid from discrete Unif distribution on [nn/2, nn]
N.dense <- c(sample((7*nn/8):nn, nn-1, replace=T), nn)
# number of observations for sparse scenario: N_i are iid from discrete Unif distribution on {3, 4, ..., 10}
N.sparse <- c(sample(3:10, nn-1, replace=T), 10)  
N.max <- max(N.dense, N.sparse)
t.all <- seq(1, 2, length.out = N.max)

# error
sigma.true <- 1
eps <- matrix(rnorm(nn*N.max, sd=sigma.true), nrow = nn, ncol = N.max)
# mu
mu <- function(t) {-0.5+1.5*sin(10*pi*(t+0.5))+4*(t-1)^3}
plot(t.all, mu(t.all), type = "l")
phi <- function(t) {sqrt(2)*sin(2*pi*t)}
xi <- rnorm(nn, 0, 1)
plot((xi%o%phi(t.all))[2,])
plot(phi(t.all), type = "l")
plot(mu(t.all)+eps[1:N.max], type = "l")
y.mat <- mu(t.all)+eps #+ t(xi%o%phi(t.all))
# plot(y.mat[, 1], type = "l")


position <- lapply(1:nn, function(i){sort(sample(1:N.max, N.dense[i], replace = F))})
y.list <- lapply(1:nn, function(i){y.mat[position[[i]], i]})

# pdf(file="R_figures/FigureObserverationSparse.pdf", width=7, height=4)
par(mfrow=c(1,2))
plot(t.all, y=seq(min(y.mat), max(y.mat), length.out = N.max),
     xlab="time", ylab="Obs", type='n')
for(i in 1:20){
  points(t.all[position[[i]]], y.list[[i]], type="l")
}
par(new=TRUE)
plot(t.all, mu(t.all), ylim = c(min(y.mat), max(y.mat)),
     xlab="time", ylab="Obs",  type = "l", col="red")

timepoints <- lapply(1:nn, function(i){t.all[position[[i]]]})

# for each subject: around xx% of observations are censored below by 0
indexcensored <- lapply(1:nn,function(i){
  temp = y.list[[i]]
  idxcensor = which(temp<0)
  idxcensor
})

delta <- lapply(1:nn, function(i){
  temp <- rep(0, N.dense[i])
  temp[indexcensored[[i]]] <- 1
  temp
})

sum(do.call(c,delta))/length(do.call(c,delta))
# 0.4217425 DL for sparse
# 0.4242691 DL for dense

# the observations (contaminated by detection limit)
observed <- lapply(1:nn,function(i){
  temp <- y.list[[i]]
  temp[indexcensored[[i]]] <- 0
  temp
})
plot(t.all, y=seq(min(y.mat), max(y.mat), length.out = N.max),
     xlab="time", ylab="Obs", type='n')
for(i in 1:20){
  points(t.all[position[[i]]], observed[[i]], type="l")
}
par(new=TRUE)
plot(t.all, mu(t.all), ylim = c(min(y.mat), max(y.mat)),
     xlab="time", ylab="Obs",  type = "l", col="red")
par(mfrow=c(1,1))
# dev.off()

#####################
##################### Estimation
# OBS scheme: each observation has the same weight: w_i=1/sum(N_i)
scheme <- "OBS"
w <- c()
if (scheme=="OBS"){
  N.total <- 0
  for (i in 1:nn)
  {
    N.total <- sum(N.total, length(position[[i]]))
  }
  w <- lapply(1:nn, function(i){
    temp <- rep(1/N.total, length(position[[i]]))
    temp
  })
}else{
  # SUBJ scheme: w_i=1/nN_i
  # w.SUBJ <- c()
  w <- lapply(1:nn, function(i){
    temp <- rep(1/(nn*length(position[[i]])), length(position[[i]]))
    temp
  })
}



t_eval=seq(1, 2, length=200)

# ######### constant approximation
bw = 0.0075
start_time <- Sys.time()
xx <- yy.censored <- c()
for (i in 1:nn){
  xx <- c(xx, t.all[position[[i]]])
  yy.censored <- c(yy.censored, observed[[i]])
}
est.censored <- lm(yy.censored ~ xx)
sigma.est <- sqrt(mean((residuals(est.censored))^2))

# calculate S and R

xx <- yy.S0 <- yy.R0 <- c()
for (i in 1:nn){
  xx <- c(xx, t.all[position[[i]]])
  yy.S0 <- c(yy.S0, w[[i]][1]*(1-0.498*delta[[i]]))
  yy.R0 <- c(yy.R0, w[[i]][1]*(-0.8194*delta[[i]]*sigma.est + (1-0.498*delta[[i]])*observed[[i]]))
}
S0 <- ksmooth(xx, yy.S0, kernel="normal",bandwidth=2*bw, x.points = t_eval)$y
R0 <- ksmooth(xx, yy.R0, kernel="normal",bandwidth=2*bw, x.points = t_eval)$y

constant.approx <- R0/S0
end_time <- Sys.time()
constant.approx.time <- end_time-start_time


# ######### linear approximation
bw = 0.01 
start_time <- Sys.time()
xx <- yy.censored <- c()
for (i in 1:nn){
  xx <- c(xx, t.all[position[[i]]])
  yy.censored <- c(yy.censored, observed[[i]])
}
est.censored <- lm(yy.censored ~ xx)
sigma.est <- sqrt(mean((residuals(est.censored))^2))

# calculate S and R
xx <- yy.S0 <- yy.S1 <- yy.S2 <- yy.R0 <- yy.R1 <- c()
for (i in 1:nn){
  xx <- c(xx, t.all[position[[i]]])
  yy.S0 <- c(yy.S0, w[[i]][1]*(1-0.498*delta[[i]]))
  yy.S1 <- c(yy.S1, w[[i]][1]*(1-0.498*delta[[i]])*t.all[position[[i]]]/bw)
  yy.S2 <- c(yy.S2, w[[i]][1]*(1-0.498*delta[[i]])*(t.all[position[[i]]]/bw)^2)
  yy.R0 <- c(yy.R0, w[[i]][1]*(-0.8194*delta[[i]]*sigma.est + (1-0.498*delta[[i]])*observed[[i]]))
  yy.R1 <- c(yy.R1, w[[i]][1]*(-0.8194*delta[[i]]*sigma.est + (1-0.498*delta[[i]])*observed[[i]])*t.all[position[[i]]]/bw)
}
S0 <- ksmooth(xx, yy.S0, kernel="normal",bandwidth=2*bw, x.points = t.all)$y
S1 <- ksmooth(xx, yy.S1, kernel="normal",bandwidth=2*bw, x.points = t.all)$y-S0*t.all/bw
S1.half <- ksmooth(xx, yy.S1, kernel="normal",bandwidth=2*bw, x.points = t.all)$y
S2 <- ksmooth(xx, yy.S2, kernel="normal",bandwidth=2*bw, x.points = t.all)$y-2*S1.half*t.all/bw+S0*(t.all/bw)^2
R0 <- ksmooth(xx, yy.R0, kernel="normal",bandwidth=2*bw, x.points = t.all)$y
R1 <- ksmooth(xx, yy.R1, kernel="normal",bandwidth=2*bw, x.points = t.all)$y-R0*t.all/bw

linear.approx <- (R0*S2-R1*S1)/(S0*S2-S1^2)
end_time <- Sys.time()
linear.approx.time <- end_time-start_time


############################
# show the results
par(mfrow=c(1,1))
plot(t_eval, mu(t_eval), type = "l", ylim=c(-2,5), xlab="time", ylab="y", col=1, lty=1)
par(new=TRUE)
plot(t_eval, constant.approx, type = "l", ylim=c(-2,5), xlab="time", ylab="y",col=2, lty=2)
par(new=TRUE)
plot(t_eval, linear.approx, type = "l", ylim=c(-2,5), xlab="time", ylab="y",col=3, lty=3)
legend(1,5, c("true", "local constant approx", "local linear approx"), lty=1:3, col=1:3)
