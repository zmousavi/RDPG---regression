

rm(list=ls())
n=500
library("rgl")
library("MASS")

graph_size=n
t=seq(0, 1, by=(1/(n-1)))

seed_start=88
noise_increase = 1
#set.seed(44) #works for b1 not b0 
#set.seed(453) #works for 100!!!! not 500!
#set.seed(8359) #works for 100!!!! not 500!
#set.seed(666) #works
#set.seed(88) #works
set.seed(seed_start)
theta_true = seq(0, 1, by=1/(graph_size-1))


theta_hw = function(th) c(th^2,2*th*(1-th),(1-th)^2)
#theta_z = function(th) c(th^3-th^2,2*th*(1-th^2),(1-th)^2)
#theta_cube = function(th) c(th^3,th^2,th)

theta_f = theta_hw
#theta_f= theta_cube

X = matrix(0,nrow=n,ncol=3)
for(i in 1:n) X[i,] = theta_f(theta_true[i])

beta_naive = list()
beta_true = list()
beta_adj_naive = list()
beta_adj_true = list()
beta_true_se = list()
beta_naive_se = list()
beta_adj_true_se=list()
beta_adj_naive_se=list()

beta_adj_true_hw = list()
beta_naive_hw = list()
beta_adj_true_hw_se = list()
beta_naive_hw_se = list()

beta_adj_est = list()
beta_adj_est_se = list()

mc_runs = 20
r=1
for (r in c(1:mc_runs)){
  print(r)
  set.seed(r)
  y =  1 + 1*theta_true + rnorm(n=graph_size, 0, 1e-5)
  X_true = theta_true #(X for regression ... not to be cofused with latent possition)
  
Delta_inv = matrix(c(9, -9, 3, -9, 21, -9, 3, -9, 9), nrow=3)

X_noise = matrix(, ncol=3, nrow=graph_size)
sigma_u_theta_true = c()
sigma_u_theta_true2 = c()
Sigma_X_matrices = list()
for (i in 1:graph_size){
  f = theta_true[i]
  S11 = (1)*((1/126) + (f/15)+(f^2)/42-(f^3)/105-2*(f^4)/35)
  S12 = (1)*((13/1260) + (4*f/105)-((f^2)/30)+(2*(f^3)/105)-((f^4)/70))
  S13 = (1)*((1/180)+(f/105)-((f^2)/70)+(f^3)/105-(f^4)/210)
  S22 = (1)*((1/45)+(4*f/105)-(2*f^2/35)+(4*f^3/105)-(2*f^4/105))
  S23 = (1)*((5/252)+(f/35)-13*(f^2)/210+(4*f^3)/105-(f^4/70))
  S33 = (1)*((2/63)+(f/7)-(73*(f^2)/210)+(5*(f^3)/21) - (2*(f^4)/35))
  S_matrix = matrix(c(S11, S12, S13, S12, S22, S23, S13, S23, S33), nrow=3, byrow=T)
  if(any(diag(S_matrix)<0)){print("yes")}
  Sigma_X = (1/graph_size)*Delta_inv%*%S_matrix%*%Delta_inv*noise_increase
  Sigma_X_matrices[[i]] = Sigma_X 
  sigma_u_theta_true[i] = (1/4)*(1*Sigma_X[1,1]-2*Sigma_X[1,3]+Sigma_X[3,3])
  sigma_u_theta_true2[i] = (1/4)*(4*Sigma_X[1,1]+4*Sigma_X[1,2]+Sigma_X[2,2])
  X_noise[i,1:3] = X[i,1:3]  + mvrnorm(n=1, rep(0, 3), Sigma=Sigma_X)
  
}


#Q=(1/3)*matrix(c(2, 1, 2, -2, 2, 1, 1, 2, -2), nrow=3)
Q=(1/3)*matrix(c(-2, 1, 2, 2, 2, 1, -1, 2, -2), nrow=3)
det(Q)
Q%*%t(Q)
Xr_noise=X_noise%*%Q

# df=data.frame(rbind(X, Xr_noise))
# df$type = c(rep(1, (graph_size)), rep(2, (graph_size)))
# with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
# triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")



# g=Xr_noise
# g_dist = dist(g, diag=T, upper=T)
# g_dist_m = as.matrix(g_dist, nrow=graph_size)
# which(g_dist_m == max(g_dist_m), arr.ind = TRUE)
# end1_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[1]
# end2_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[2]
#
# index_sorted=c()
# l=sort(g_dist_m[end1_index,])
# for (i in 1:length(l)){
# sorted_index=which(g_dist_m[end1_index,] == l[i], arr.ind = TRUE)[1]
# index_sorted[i]=sorted_index
# }
#
# Xr_noise_sorted = matrix(, ncol=3, nrow=graph_size)
# for (i in 1:length(l)){
#   Xr_noise_sorted[i,] = Xr_noise[index_sorted[i],]
# }


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

# df=data.frame(rbind(X, Xr_noise,  X_fit_hw))
# df$type = c(rep(1, (graph_size)), rep(2, (graph_size)), rep(3, (graph_size)))
# with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
# triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")


# ###CUBE!###
# 
# Xr_noise_sorted = Xr_noise
# W.Xhat = Xr_noise_sorted
# 
# A1_cube = c(2*sum(t^6) , 2*sum((t^5)) , 2*sum((t^4)))
# A2_cube=c(2*sum(t^5) , 2*sum((t^4)) , 2*sum((t^3)))
# A3_cube = c(2*sum(t^4) , 2*sum((t^3)) , 2*sum((t^2)))
# A_cube = matrix(c(A1_cube, A2_cube, A3_cube), byrow=T, nrow=3, ncol=3)
# 
# d=W.Xhat[,1]
# b1_cube = 2*sum(t^3*d)
# b2_cube = 2*sum(t^2*d)
# b3_cube =  2*sum(t^1*d)
# b_cube=rbind(b1_cube, b2_cube, b3_cube)
# x_cube = solve(A_cube)%*%b_cube
# 
# d=W.Xhat[,2]
# b1_cube = 2*sum(t^3*d)
# b2_cube = 2*sum(t^2*d)
# b3_cube =  2*sum(t^1*d)
# b_cube=rbind(b1_cube, b2_cube, b3_cube)
# y_cube = solve(A_cube)%*%b_cube
# 
# d=W.Xhat[,3]
# b1_cube = 2*sum(t^3*d)
# b2_cube = 2*sum(t^2*d)
# b3_cube =  2*sum(t^1*d)
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
# min.RSS = min.RSS_cube
 
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
 
 sigma_u_theta_true_sorted=sigma_u_theta_true
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
 plot(theta_true, X_naive)
 

 
 
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
 


 
 
 
 
 plot(theta_sorted, sigma_u_theta_estimated_i, ylim=c(0, 0.02), col="red")
 par(new=T)
 plot(theta_sorted, sigma_u_theta_true_sorted, ylim=c(0, 0.02))  



 
Xz_big = cbind(rep(1, graph_size), X_true)
X_naive_big = cbind(rep(1, graph_size), X_naive)

beta_true[[r]] = solve(t(Xz_big)%*%Xz_big)%*%t(Xz_big)%*% y
beta_naive[[r]] = solve(t(X_naive_big)%*%(X_naive_big))%*%t(X_naive_big)%*% y


Z = cbind(y, X_naive_big)
Sigma_U_theta_true_sorted=sigma_u_theta_true_sorted


M = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  M = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U_theta_true_sorted[i]))) + M
}
M_z = (1/graph_size)*M
beta_adj_true[[r]] = solve(M_z[2:3,2:3])%*%M_z[2:3,1]



Sigma_U_theta_est=sigma_u_theta_estimated_i
Q = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  Q = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U_theta_est[i]))) + Q
}
Q_z = (1/graph_size)*Q
beta_adj_est[[r]] = solve(Q_z[2:3,2:3])%*%Q_z[2:3,1]


beta_true_se[[r]] = (beta_true[[r]]-c(1, 1))^2
beta_adj_true_se[[r]] = (beta_adj_true[[r]]-c(1, 1))^2
beta_naive_se[[r]] = (beta_naive[[r]]-c(1, 1))^2
beta_adj_est_se[[r]] = (beta_adj_est[[r]]-c(1, 1))^2

beta_adj_true
beta_adj_est
beta_naive


}


b0_true_sim = c()
b1_true_sim = c()
b0_naive_sim = c()
b1_naive_sim = c()
b0_adj_true_sim = c()
b1_adj_true_sim = c()
b0_adj_est_sim = c()
b1_adj_est_sim = c()
b0_true_se_sim = c()
b1_true_se_sim = c()
b0_adj_true_se_sim = c()
b1_adj_true_se_sim = c()
b0_adj_est_se_sim = c()
b1_adj_est_se_sim = c()
b0_naive_se_sim = c()
b1_naive_se_sim = c()

for (z in 1:(r-1)){
  
  if (z%%100 ==0){
    print (z)
  } 
  
  b0_true_sim[z] =  beta_true[[z]][1]
  b1_true_sim[z] =  beta_true[[z]][2]
  
  b0_true_se_sim[z] =  beta_true_se[[z]][1]
  b1_true_se_sim[z] =  beta_true_se[[z]][2]
  
  b0_naive_sim[z] = beta_naive[[z]][1]
  b1_naive_sim[z] = beta_naive[[z]][2]
  
  b0_naive_se_sim[z] = beta_naive_se[[z]][1]
  b1_naive_se_sim[z] = beta_naive_se[[z]][2]
  

  b0_adj_true_se_sim[z] = beta_adj_true_se[[z]][1]
  b1_adj_true_se_sim[z] = beta_adj_true_se[[z]][2]
  
  b0_adj_true_sim[z] = beta_adj_true[[z]][1]
  b1_adj_true_sim[z] = beta_adj_true[[z]][2]
  
  b0_adj_est_sim[z] = beta_adj_est[[z]][1]
  b1_adj_est_sim[z] = beta_adj_est[[z]][2]
  
  b0_adj_est_se_sim[z] = beta_adj_est_se[[z]][1]
  b1_adj_est_se_sim[z] = beta_adj_est_se[[z]][2]
  
 
  
}

###########
###########
###########





par(mar=c(10,4,3,2))
boxplot(b0_true_sim, b0_naive_sim,  b0_adj_true_sim,  b0_adj_est_sim, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b0_true_sim), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjsuted_Est"))
abline(h = b0_true_sim, col="red", lty=2, lwd=0.5)

par(mar=c(10,4,3,2))
boxplot(b1_true_sim, b1_naive_sim,  b1_adj_true_sim, b1_adj_est_sim, notch=TRUE, 
        main=bquote(paste("b1 Estimate","\n graph_size:", graph_size, "\n mc_runs:",length(b1_true_sim), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjsuted_Est"))
abline(h = b1_true_sim, col="red", lty=2, lwd=0.5)

par(mar=c(8,4,2,2))
boxplot(b0_true_se_sim, b0_naive_se_sim,  b0_adj_true_se_sim, b0_adj_est_se_sim, notch=TRUE, 
        main=bquote(paste("b0 Square Error", "\n graph_size:", graph_size ,"Sigma(error)=4x" ,"\n mc_runs:", length(b0_true_sim), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true",  "Adjsuted_Est"))
abline(h = b0_true_sim, col="red", lty=2, lwd=0.5)

par(mar=c(8,4,2,2))
boxplot(b1_true_se_sim, b1_naive_se_sim,  b1_adj_true_se_sim, b1_adj_est_se_sim, notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n graph_size:", graph_size ,"Sigma(error)=4x" ,"\n mc_runs:", length(b1_true_sim), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true",  "Adjsuted_Est"))
abline(h = b1_true_sim, col="red", lty=2, lwd=0.5)



par(mar=c(10,4,3,2))
boxplot(b0_true, b0_naive,  b0_adj_true,  b0_adj_est, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", graph_size , "\n mc_runs:",length(b0_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjsuted_Est"))
abline(h = b0_true, col="red", lty=2, lwd=0.5)

par(mar=c(10,4,3,2))
boxplot(b1_true, b1_naive,  b1_adj_true, b1_adj_est, notch=TRUE, 
        main=bquote(paste("b1 Estimate","\n graph_size:", graph_size, "\n mc_runs:",length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjsuted_Est"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(8,4,2,2))
boxplot(b0_true_se, b0_naive_se,  b0_adj_true_se, b0_adj_est_se, notch=TRUE, 
        main=bquote(paste("b0 Square Error", "\n graph_size:", graph_size ,"Sigma(error)=4x" ,"\n mc_runs:", length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true",  "Adjsuted_Est"))
abline(h = b0_true, col="red", lty=2, lwd=0.5)

par(mar=c(8,4,2,2))
boxplot(b1_true_se, b1_naive_se,  b1_adj_true_se, b1_adj_est_se, notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n graph_size:", graph_size ,"Sigma(error)=4x" ,"\n mc_runs:", length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true",  "Adjsuted_Est"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)


 #######
######
#######


par(mar=c(10,4,3,2))
boxplot(b0_true, b0_naive,  b0_adj_true, b0_naive_hw, b0_adj_true_hw, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", 250 , "\n mc_runs:",length(b0_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        ylim=c(0.9, 1.1),
        las=2, names=c("True", "Naive", "Adjusted_true", "Naive_hw" , "Adjusted_true_hw"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(10,4,3,2))
boxplot(b1_true, b1_naive,  b1_adj_true, b1_naive_hw, b1_adj_true_hw, notch=TRUE, 
        main=bquote(paste("b1 Estimate","\n graph_size:", 250, "\n mc_runs:",length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        ylim=c(0.8, 1.2),
        las=2, names=c("True", "Naive", "Adjusted_true", "Naive_hw" , "Adjusted_true_hw"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(10,4,3,2))
boxplot(b0_true, b0_naive,  b0_adj_true, notch=TRUE, 
        main=bquote(paste("b0 Estimate","\n graph_size:", 250 , "\n mc_runs:",length(b0_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(10,4,3,2))
boxplot(b1_true, b1_naive,  b1_adj_true,  notch=TRUE, 
        main=bquote(paste("b1 Estimate","\n graph_size:", 250, "\n mc_runs:",length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)

par(mar=c(8,4,2,2))
boxplot(b0_true_se, b0_naive_se,  b0_adj_true_se, notch=TRUE, 
        main=bquote(paste("b0 Square Error", "\n graph_size:", n ,"Sigma(error)=4x" ,"\n mc_runs:", length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true"))
abline(h = 0, col="red", lty=2, lwd=0.5)



par(mar=c(8,4,2,2))
boxplot(b0_true_se, b0_naive_se,  b0_adj_true_se, b0_naive_hw_se, b0_adj_true_hw_se, notch=TRUE, 
        main=bquote(paste("b0 Square Error", "\n graph_size:", n ,"Sigma(error)=4x" ,"\n mc_runs:", length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Naive_hw" , "Adjusted_true_hw"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(8,4,2,2))
boxplot(b1_true_se, b1_naive_se,   b1_adj_true_se, b1_naive_hw_se, b1_adj_true_hw_se, notch=TRUE, 
        main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"Sigma(error)=4x" ,"\n mc_runs:", length(b1_true), "\n seed:", seed_start)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Naive_hw" , "Adjusted_true_hw"))
abline(h = 0, col="red", lty=2, lwd=0.5)





#setwd("~/Desktop")
#dev.new()
#pdf("b0_estimate_rdpg_sim_wr.pdf", 7, 5)
par(mar=c(8,4,2,2))
boxplot(b0_true, b0_naive,  b0_adj_true, b0_adj_naive,  notch=TRUE, 
        #main=bquote(paste("b1 Estimate","\n graph_size:", n , "\n mc_runs:",length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjusted_naive"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)
#dev.off()

#dev.new()
#pdf("b1_estimate_rdpg_sim_wr.pdf", 7, 5)
par(mar=c(8,4,2,2))
boxplot(b1_true, b1_naive,  b1_adj_true, b1_adj_naive,  notch=TRUE, 
        #main=bquote(paste("b1 Estimate","\n graph_size:", n , "\n mc_runs:",length(b1_true), "\n m:", m)), 
       # ylim=c(0.8, 1.2),
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjusted_naive"))
abline(h = b1_true, col="red", lty=2, lwd=0.5)
#dev.off()

#dev.new()
#pdf("b0_se_rdpg_sim_wr.pdf", 7, 5)
par(mar=c(8,4,2,2))
boxplot(b0_true_se, b0_naive_se,  b0_adj_true_se, b0_adj_naive_se, notch=TRUE, 
        #   main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjusted_naive"))
abline(h = 0, col="red", lty=2, lwd=0.5)
#dev.off()

#dev.new()
#pdf("b1_se_rdpg_sim_wr.pdf", 7, 5)
par(mar=c(8,4,2,2))
boxplot(b1_true_se, b1_naive_se,   b1_adj_true_se, b1_adj_naive_se, notch=TRUE, 
        #   main=bquote(paste("b1 Square Error", "\n graph_size:", n ,"\n mc_runs:", length(b1_true), "\n m:", m)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted_true", "Adjusted_naive"))
abline(h = 0, col="red", lty=2, lwd=0.5)
#dev.off()

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



-
df=data.frame(rbind(X[1:(graph_size/2),], Xr_noise[1:(graph_size/2),], W.Xhat[1:(graph_size/2),],  X_fit[1:(graph_size/2),]))
df$type = c(rep(1, (graph_size/2)), rep(2, (graph_size/2)),  rep(4, (graph_size/2)), rep(3, (graph_size/2)))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")



df=data.frame(rbind(X, W.Xhat))
#df$type =  rep(c(1:7, rep(8,(graph_size-7))), 2)
df$type =  c(1:3, rep(7,(graph_size-3)), 4:6, rep(8,(graph_size-3)))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")


y  = runif(100)
ys = smth.gaussian(y)
ys3 = smth.gaussian(y, window=0.1)
 plot(ys, ylim=c(0.2, 1))
 par(new=T)
 plot(ys3, col="blue", ylim=c(0.2, 1))
 
 
 
 
 
 set.seed(88)
 graph_size=100 #10000
 y  = runif(graph_size)

 sd_inc=seq(0, 10, by=(10/(n-1)))
 sd_inc_rand_square=sort(runif(n=(graph_size), min=0, max=10 ))^2
 sd_inc_rand_quad=sort(runif(n=(graph_size), min=0, max=10 ))^4
 sd_inc_dec=c(seq(0, 5, by=(5/((n/2)-1))), sort(seq(0, 5, by=(5/((n/2)-1))), dec=T))
 sd_inc_rand_u=sort(runif(n=(graph_size), min=0, max=10 )-5)^2
 sd_inc_rand_n= -sort(runif(n=(graph_size), min=0, max=10 )-5)^2+10^2
 
 
 sd_noise = 5*10^-4*sd_inc_rand_u
 set.seed(745)
 y_noise = y + rnorm(n=graph_size, mean=0,  sd=sd_noise)
 

 
 # 10^0 and 10^-1
 # window_half=round(0.01*0.5*length(y)/2)
 # gaus = dnorm(seq(-0.5,0.5, length.out=(2*window_half+1)), 0, 0.2)
 # 
 #10^-2
 # window_half = 20
 #   sd_g = 0.5
 #   if ((i >1000) & (i<=2000)){
 #     sd_g=2}
 #   if ((i >2000) & (i<8000)){
 #     sd_g=3}
 #   if ((i >8000) & (i<=9000)){
 #     sd_g=2}
 

 
 y_noise = X_naive_sorted
 y=theta_sorted
 sd_noise=sqrt(sigma_u_theta_true_sorted)
 plot(theta_sorted, X_naive_sorted)
 
 #sigma_u_theta_fit_sorted
 #sigma_u_theta_estimated_i
 
 # window_half = 10
 # y_padded=c(rep(0, window_half), y_noise, rep(0, window_half))
 # sd_i = c()
 # ptr=0
 # for (i in (window_half+1):(window_half+length(y))){
 #   sd_g = 10
 #   if ((i >0.1*length(y)) & (i<=0.2*length(y))){
 #     sd_g=20}
 #   if ((i >0.2*length(y)) & (i<0.8*length(y))){
 #     sd_g=50}
 #   if ((i >0.8*length(y)) & (i<=0.9*length(y))){
 #     sd_g=20}
 # 
 
 
 window_half=round(0.1*0.5*length(y)/2)
 window_half = 10
 y_padded=c(rep(0, window_half), y_noise, rep(0, window_half))
 sd_i = c()
 ptr=0
 for (i in (window_half+1):(window_half+length(y))){
   sd_g = 1
   gaus = dnorm(seq(-0.5,0.5, length.out=(2*window_half+1)), 0, sd_g) 
   y_pad_conv = gaus*y_padded[(i-window_half):(i+window_half)]
   ptr=ptr+1
   sd_i[ptr] = sd(y_pad_conv)
 }
 
# sd_i = sd_i^2
# sd_noise = sigma_u_theta_true_sorted
 
 plot( sd_i, col="red", ylim=c(min(min(sd_i), min(sd_noise)), max(max(sd_i), max(sd_noise))), xlim=c(0, length(sd_i)))
 par(new=T)
 plot(sd_noise, col="blue",  ylim=c(min(min(sd_i), min(sd_noise)), max(max(sd_i), max(sd_noise))), xlim=c(0, length(sd_noise)))

mean(sd_i-sd_noise)^2
 
 
plot( sd_i, col="red", ylim=c(min(min(sd_i), min(sd_noise)), max(max(sd_i), max(sd_noise))), xlim=c(0, length(sd_i)))
par(new=T)
plot(sd_noise, col="blue",  ylim=c(min(min(sd_i), min(sd_noise)), max(max(sd_i), max(sd_noise))), xlim=c(0, length(sd_noise)))



  set.seed(999)
  X_naive_sorted_sim = y+rnorm(n=length(theta_sorted), mean=0, sd=10*sigma_u_theta_true_sorted)
  X_naive_sorted =  X_naive_sorted_sim
  plot(theta_sorted, X_naive_sorted)
 window_half=round(0.05*0.5*length(X_naive_sorted))
 X_naive_sorted_padded=c(rep(0, window_half), X_naive_sorted, rep(0, window_half))
 
 sigma_u_theta_estimated_i = c()
 
 for (i in (window_half+1):(window_half+length(y))){
   gaus = dnorm(seq(-0.5,0.5, length.out=(2*window_half+1)), 0, .5) 
   X_naive_pad_conv = 1*X_naive_sorted_padded[(i-window_half):(i+window_half)]
   sigma_u_theta_estimated_i[i] = sd(X_naive_pad_conv)
 }
 plot(sigma_u_theta_estimated_i)
 
 plot(sigma_u_theta_estimated_i, ylim=c(0, 0.02))
 par(new=T)
 plot(10*sigma_u_theta_true_sorted, col="blue", ylim=c(0, 0.02))
 
 


 
 #ys = smth.gaussian(y)


 #ys_padded=c(rep(0, window_half), ys, rep(0, window_half))
 
 

 
 
 sd(y_padded) 
 

 