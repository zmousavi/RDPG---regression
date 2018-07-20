

rm(list=ls())
n=500
library("rgl")
library("MASS")

graph_size=n
t=seq(0, 1, by=(1/(n-1)))
#theta_true = t

set.seed(88)
#set.seed(453) #nope
#set.seed(8359) #nope only works at n=100
#set.seed(666) #works
#set.seed(88)

theta_true = c(0, 0.5, 1, runif(n=(graph_size-3), min=0, max=1 ))
theta_true = seq(1/(graph_size), 1, by=1/(graph_size))

theta_true_hw = seq(0, 1, by=1*(1/(graph_size-1)))
theta_true_cube = seq(0, 0.5, by=(1/2)*(1/(graph_size-1)))

theta_hw = function(th) c(th^2,2*th*(1-th),(1-th)^2)
#theta_cube = function(th) c((th)^3,(th)^2,(th))

theta_true = theta_true_hw
#theta_true = theta_true_cube


theta_f= theta_hw
#theta_f = theta_cube



X = matrix(0,nrow=n,ncol=3)
for(i in 1:n) X[i,] = theta_f((theta_true[i]))



P = X %*% t(X)



beta_naive = list()
beta_true = list()
beta_adj_naive = list()
beta_adj_true = list()
beta_true_se = list()
beta_naive_se = list()
beta_adj_true_se=list()
beta_adj_naive_se=list()
beta_adj_est = list()
beta_adj_est_se = list()

mc_runs = 20
r=1
for (r in c(1:mc_runs)){
  print(r)
  set.seed(r)
  A = matrix(0,nrow=n,ncol=n)
  for(i in 1:(n-1)) for(j in (i+1):n) A[i,j] = rbinom(1,1,P[i,j])
  A = A + t(A)
  # compare: diagaug
  svdA = svd(A)
  Xhat = svdA$u %*% diag(sqrt(svdA$d))[,1:3]
  

  Xr_noise=Xhat
  
  df=data.frame(rbind(X, Xhat))
  df$type = c(rep(1, (graph_size)), rep(2, (graph_size)))
  with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
  triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")
  
  
  
  
  
  #y =  1 + 1*(X[,1]+0.5*X[,2]) + rnorm(n=graph_size, 0, 1e-5)
  y =  1 + 1*theta_true + rnorm(n=graph_size, 0, 1e-5)
  X_true = theta_true #(X for regression ... not to be cofused with latent possition)
  





#sort and find one endpoint as t=0 

W.Xhat = Xr_noise

Xr_noise_sorted = Xr_noise
W.Xhat = Xr_noise_sorted

A1_hw = c(2*sum(t^4) , 4*sum((t^3)*(1-t)) , 2*sum((t^2)*((1-t)^2)))
A2_hw =c(4*sum(t^3*(1-t)), 8*sum(t^2*(1-t)^2), 4*sum(t*(1-t)^3))
A3_hw = c(2*sum(t^2*(1-t)^2), 4*sum(t*(1-t)^3), 2*sum((1-t)^4))   
A_hw = matrix(c(A1_hw, A2_hw, A3_hw), byrow=T, nrow=3, ncol=3)       



d=W.Xhat[,1]
b1_hw = 2*sum(t^2*d)
b2_hw  = 4*sum(t*(1-t)*d)
b3_hw  =  2*sum(d*(1-t)^2)
b_hw =rbind(b1_hw, b2_hw, b3_hw)
x_hw = solve(A_hw)%*%b_hw


d=W.Xhat[,2]
b1_hw = 2*sum(t^2*d)
b2_hw  = 4*sum(t*(1-t)*d)
b3_hw  =  2*sum(d*(1-t)^2)
b_hw=rbind(b1_hw, b2_hw, b3_hw)
y_hw = solve(A_hw)%*%b_hw


d=W.Xhat[,3]
b1_hw = 2*sum(t^2*d)
b2_hw  = 4*sum(t*(1-t)*d)
b3_hw  =  2*sum(d*(1-t)^2)
b_hw=rbind(b1_hw, b2_hw, b3_hw)
z_hw = solve(A_hw)%*%b_hw

x_fit_hw = x_hw[1]*t^2+2*x_hw[2]*t*(1-t)+x_hw[3]*(1-t)^2
y_fit_hw = y_hw[1]*t^2+2*y_hw[2]*t*(1-t)+y_hw[3]*(1-t)^2
z_fit_hw = z_hw[1]*t^2+2*z_hw[2]*t*(1-t)+z_hw[3]*(1-t)^2

X_fit_hw=cbind(x_fit_hw, y_fit_hw, z_fit_hw)
Q_fit_hw =cbind(x_hw, y_hw, z_hw)


df=data.frame(rbind(X, Xr_noise,  X_fit_hw))
df$type = c(rep(1, (graph_size)), rep(2, (graph_size)), rep(3, (graph_size)))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")


# ###CUBE!###
# 
# Xr_noise_sorted = Xr_noise
# W.Xhat = Xr_noise_sorted
# 
# t=seq(0, 0.5, by=((1/2)*(1/(n-1))))
# 
# A1_cube = c(2*sum((t/1)^6) , 2*sum(((t/1)^5)) , 2*sum(((t/1)^4)))
# A2_cube=c(2*sum((t/1)^5) , 2*sum(((t/1)^4)) , 2*sum(((t/1)^3)))
# A3_cube = c(2*sum((t/1)^4) , 2*sum(((t/1)^3)) , 2*sum(((t/1)^2)))
# A_cube = matrix(c(A1_cube, A2_cube, A3_cube), byrow=T, nrow=3, ncol=3)
# 
# 
# 
# 
# 
# d=W.Xhat[,1]
# b1_cube = 2*sum((t/1)^3*d)
# b2_cube = 2*sum((t/1)^2*d)
# b3_cube =  2*sum((t/1)^1*d)
# b_cube=rbind(b1_cube, b2_cube, b3_cube)
# x_cube = solve(A_cube)%*%b_cube
# 
# d=W.Xhat[,2]
# b1_cube = 2*sum((t/1)^3*d)
# b2_cube = 2*sum((t/1)^2*d)
# b3_cube =  2*sum((t/1)^1*d)
# b_cube=rbind(b1_cube, b2_cube, b3_cube)
# y_cube = solve(A_cube)%*%b_cube
# 
# d=W.Xhat[,3]
# b1_cube = 2*sum((t/1)^3*d)
# b2_cube = 2*sum((t/1)^2*d)
# b3_cube =  2*sum((t/1)^1*d)
# b_cube=rbind(b1_cube, b2_cube, b3_cube)
# z_cube = solve(A_cube)%*%b_cube
# 
# x_fit_cube = x_cube[1]*t^3+1*x_cube[2]*t^2+x_cube[3]*t
# y_fit_cube = y_cube[1]*t^3+1*y_cube[2]*t^2+y_cube[3]*t
# z_fit_cube = z_cube[1]*t^3+1*z_cube[2]*t^2+z_cube[3]*t
# 
# X_fit_cube=cbind(x_fit_cube, y_fit_cube, z_fit_cube)
# Q_fit_cube =cbind(x_cube, y_cube, z_cube)
# 
# 
# df=data.frame(rbind(X, Xr_noise,  X_fit_cube))
# df$type = c(rep(1, (graph_size)), rep(2, (graph_size)), rep(3, (graph_size)))
# with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
# triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")

###PROJECT ONTO S_HAT

min.RSS_hw <- function(data, par) {
  with(data, sum(((x_hw[1]*par^2+2*x_hw[2]*par*(1-par)+x_hw[3]*(1-par)^2)-h[1])^2+(( y_hw[1]*par^2+2*y_hw[2]*par*(1-par)+y_hw[3]*(1-par)^2)-h[2])^2+((z_hw[1]*par^2+2*z_hw[2]*par*(1-par)+z_hw[3]*(1-par)^2)-h[3])^2))}

 min.RSS_cube <- function(data, par) {
   with(data, sum(((x_cube[1]*par^3+1*x_cube[2]*par^2+x_cube[3]*(par)^1)-h[1])^2+((y_cube[1]*par^3+1*y_cube[2]*par^2+y_cube[3]*(par)^1)-h[2])^2+((z_cube[1]*par^3+1*z_cube[2]*par^2+z_cube[3]*(par)^1)-h[3])^2))}

min.RSS = min.RSS_hw
#min.RSS = min.RSS_cube

theta_noise_fit2 <- c()
for (i in 1:nrow(X)){
  dat=data.frame(h=Xr_noise [i,])
  result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
  theta_noise_fit2[i] <- result$par
}

X_naive =  theta_noise_fit2
if (sum((theta_noise_fit2-theta_true)^2) > sum(((1-theta_noise_fit2)-theta_true)^2)){
  X_naive = 1-theta_noise_fit2
  #plot(theta_true, X_naive)
}
plot(theta_true, X_naive)

# theta_sorted=c()
# for(i in (1:graph_size)){
#   theta_sorted[i] = theta_true[index_sorted[i]]
# }
# 
# sigma_u_theta_true_sorted=c()
# for(i in (1:graph_size)){
#   sigma_u_theta_true_sorted[i] =  sigma_u_theta_true[index_sorted[i]]
# }
# 

#sigma_u_theta_true_sorted=sigma_u_theta_true
X_true = theta_true


# plot(theta_sorted, theta_noise_fit_sorted)
# X_naive_sorted =  theta_noise_fit_sorted
#  if (sum((theta_noise_fit_sorted-theta_sorted)^2) > sum(((1-theta_noise_fit_sorted)-theta_sorted)^2)){
#    X_naive_sorted = 1-theta_noise_fit_sorted
#  }
#plot(theta_sorted, X_naive_sorted) #if using X_naive_sorted remeber gotta use y_sorted as well...

X_naive =  theta_noise_fit2
if (sum((theta_noise_fit2-theta_true)^2) > sum(((1-theta_noise_fit2)-theta_true)^2)){
  X_naive = 1-theta_noise_fit2
  #plot(theta_true, X_naive)
}
#plot(theta_true, X_naive)




X_naive_sorted=X_naive
theta_sorted=theta_true

y2_noise = X_naive_sorted
y2=theta_sorted

window_half = 10
y2_padded=c(rep(0, window_half), y2_noise, rep(0, window_half))
sd_i = c()
ptr=0
for (i in (window_half+1):(window_half+length(y2_noise))){
  sd_g = 1
  gaus = dnorm(seq(-0.5,0.5, length.out=(2*window_half+1)), 0, sd_g) 
  y2_pad_conv = gaus*y2_padded[(i-window_half):(i+window_half)]
  ptr=ptr+1
  sd_i[ptr] = sd(y2_pad_conv)
}
sigma_u_theta_estimated_i = sd_i^2




Xz_big = cbind(rep(1, graph_size), X_true)
X_naive_big = cbind(rep(1, graph_size), X_naive)

beta_true[[r]] = solve(t(Xz_big)%*%Xz_big)%*%t(Xz_big)%*% y
beta_naive[[r]] = solve(t(X_naive_big)%*%(X_naive_big))%*%t(X_naive_big)%*% y






Z = cbind(y, X_naive_big)
Sigma_U_theta_est=sigma_u_theta_estimated_i
Q = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  Q = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U_theta_est[i]))) + Q
}
Q_z = (1/graph_size)*Q
beta_adj_est[[r]] = solve(Q_z[2:3,2:3])%*%Q_z[2:3,1]


beta_true_se[[r]] = (beta_true[[r]]-c(1, 1))^2
beta_naive_se[[r]] = (beta_naive[[r]]-c(1, 1))^2
beta_adj_est_se[[r]] = (beta_adj_est[[r]]-c(1, 1))^2



}



b0_true = c()
b1_true = c()
b0_naive = c()
b1_naive = c()
b0_adj_est = c()
b1_adj_est = c()
b0_true_se = c()
b1_true_se = c()
b0_naive_se = c()
b1_naive_se = c()
b0_adj_est_se = c()
b1_adj_est_se = c()
 
for (z in 1:(r-1)){
  
  if (z%%100 ==0){
    print (z)
  } 
  
  b0_true[z] =  beta_true[[z]][1]
  b1_true[z] =  beta_true[[z]][2]
  
  b0_true_se[z] =  beta_true_se[[z]][1]
  b1_true_se[z] =  beta_true_se[[z]][2]
  
  b0_naive[z] = beta_naive[[z]][1]
  b1_naive[z] = beta_naive[[z]][2]
  
  b0_naive_se[z] = beta_naive_se[[z]][1]
  b1_naive_se[z] = beta_naive_se[[z]][2]
  
  
  
  b0_adj_est[z] = beta_adj_est[[z]][1]
  b1_adj_est[z] = beta_adj_est[[z]][2]
  
  b0_adj_est_se[z] = beta_adj_est_se[[z]][1]
  b1_adj_est_se[z] = beta_adj_est_se[[z]][2]
  
  
}



par(mar=c(10,4,3,2))
boxplot(b0_true, b0_naive, b0_adj_est, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs: 10")), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjsuted_Est"))
abline(h = b0_true, col="red", lty=2, lwd=0.5)


par(mar=c(10,4,3,2))
boxplot(b1_true, b1_naive, b1_adj_est, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs: 10")), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjsuted_Est"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(10,4,3,2))
boxplot(b0_true_se, b0_naive_se, b0_adj_est_se, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b0_true))), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjsuted_Est"))
abline(h = b0_true, col="red", lty=2, lwd=0.5)


par(mar=c(10,4,3,2))
boxplot(b1_true_se, b1_naive_se, b1_adj_est_se, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b1_true))), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjsuted_Est"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)
















###########
###########
###########
setwd("~/Desktop/Feb16")
dev.new()
pdf("b0_estimate_rdpg_sim_5.pdf", 7, 5)
boxplot(b0_true, b0_naive,  b0_adj_true, b0_adj_naive, b0_adj_est,  notch=TRUE, 
        #main=bquote(paste("b1 Estimate","\n graph_size:", n , "\n mc_runs:",length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        ylim=c(0.95, 1.05),
        las=1, names=c(expression(hat(beta)[true], hat(beta)[naive], hat(beta)[adj.~sigma], hat(beta)[adj.~hat(sigma)], hat(beta)[adj.~nonpar])))
        #las=2, names=c("True", "Naive", "Adjusted_true", "Adjusted_naive"))
        abline(h = b0_true, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b1_estimate_rdpg_sim_5.pdf", 7, 5)
boxplot(b1_true, b1_naive,  b1_adj_true, b1_adj_naive,  b0_adj_est, notch=TRUE, 
        #main=bquote(paste("b1 Estimate","\n graph_size:", n , "\n mc_runs:",length(b1_true), "\n m:", m)), 
        ylim=c(0.9, 1.05),
        cex.main=0.5,
        las=1, names=c(expression(hat(beta)[true], hat(beta)[naive], hat(beta)[adj.~sigma], hat(beta)[adj.~hat(sigma)],  hat(beta)[adj.~nonpar])))
       abline(h = b1_true, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b0_se_rdpg_sim_5.pdf", 7, 5)
 boxplot(b0_true_se, b0_naive_se,   b0_adj_true_se, b0_adj_naive_se, b0_adj_est_se, notch=TRUE, 
          cex.main=0.5,
         ylim=c(0, 6e-4),
         las=1, names=c(expression(hat(beta)[true], hat(beta)[naive], hat(beta)[adj.~sigma], hat(beta)[adj.~hat(sigma)], hat(beta)[adj.~nonpar])))
abline(h = 0, col="red", lty=2, lwd=0.5)
dev.off()

dev.new()
pdf("b1_se_rdpg_sim_5.pdf", 7, 5)
boxplot(b1_true_se, b1_naive_se,   b1_adj_true_se, b1_adj_naive_se, b1_adj_est_se,notch=TRUE, 
        cex.main=0.5,
        ylim=c(0, 0.004),
        las=1, names=c(expression(hat(beta)[true], hat(beta)[naive], hat(beta)[adj.~sigma], hat(beta)[adj.~hat(sigma)], hat(beta)[adj.~nonpar])))
abline(h = 0, col="red", lty=2, lwd=0.5)
dev.off()

#SE (sign test )
N = sum((b0_adj_true_se - b0_naive_se) !=0)
k = sum(b0_adj_true_se < b0_naive_se)
1 - pbinom(k, size=N, prob=0.5)


#SE (sign test )
N = sum((b1_adj_true_se - b1_naive_se) !=0)
k = sum(b1_adj_true_se < b1_naive_se)
1 - pbinom(k, size=N, prob=0.5)


#SE (sign test )
N = sum((b0_adj_naive_se - b0_naive_se) !=0)
k = sum(b0_adj_naive_se < b0_naive_se)
1 - pbinom(k, size=N, prob=0.5)


#SE (sign test )
N = sum((b1_adj_naive_se - b1_naive_se) !=0)
k = sum(b1_adj_naive_se < b1_naive_se)
1 - pbinom(k, size=N, prob=0.5)



#SE (sign test )
N = sum((b1_adj_naive_se - b1_adj_true_se) !=0)
k = sum(b1_adj_naive_se < b1_adj_true_se)
1 - pbinom(k, size=N, prob=0.5)


#SE (sign test )
N = sum((b0_adj_naive_se - b0_adj_true_se) !=0)
k = sum(b0_adj_naive_se < b0_adj_true_se)
1 - pbinom(k, size=N, prob=0.5)

boxplot(b1_adj_true_se, b1_adj_naive_se, notch=TRUE)



  
df=data.frame(rbind(X[1:2,],  X_noise[1:2,], X_noise_sorted[1:2,] ))
df$type = c(rep(1, graph_size/2),   rep(2, graph_size/2), rep(3, graph_size/2))#, rep(4, graph_size))


df=data.frame(rbind(X[1:2,],  X_noise[1:2,], Xr_noise_sorted[1:2,] ))
df$type = c(rep(1, 2),   rep(2, 2), rep(3, 2))#, rep(4, graph_size))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")



df=data.frame(rbind(X[1:8,],  Xr_noise[1:8,] ))
df=data.frame(Xr_noise)
df$type = c((1:7), rep(8, (graph_size-7)))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")




#X_noise = X_sorted
#t = t_sorted 
#X_noise = Xr_noise



df=data.frame(rbind(X, W.Xhat,  X_fit))
df$type = c(rep(1, graph_size), rep(2, graph_size),  rep(3, graph_size))#, rep(4, graph_size))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")

df=data.frame(rbind(X[1:(graph_size/2),], Xr_noise[1:(graph_size/2),], W.Xhat[1:(graph_size/2),],  X_fit[1:(graph_size/2),]))
df$type = c(rep(1, (graph_size/2)), rep(2, (graph_size/2)),  rep(3, (graph_size/2)), rep(4, (graph_size/2)))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")



df=data.frame(rbind(X, W.Xhat))
#df$type =  rep(c(1:7, rep(8,(graph_size-7))), 2)
df$type =  c(1:3, rep(7,(graph_size-3)), 4:6, rep(8,(graph_size-3)))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")

