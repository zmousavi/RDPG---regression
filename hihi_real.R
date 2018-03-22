

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


hw = function(th) c(th^2,2*th*(1-th),(1-th)^2)
X = matrix(0,nrow=n,ncol=3)
for(i in 1:n) X[i,] = hw(theta_true[i])
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

mc_runs = 200
r=1
for (r in c(1:mc_runs)){
  print(r)
  
  A = matrix(0,nrow=n,ncol=n)
  for(i in 1:(n-1)) for(j in (i+1):n) A[i,j] = rbinom(1,1,P[i,j])
  A = A + t(A)
  # compare: diagaug
  svdA = svd(A)
  Xhat = svdA$u %*% diag(sqrt(svdA$d))[,1:3]
  
 
  
  y =  1 + 1*(X[,1]+0.5*X[,2]) + rnorm(n=graph_size, 0, 1e-5)
  #y =  1 + 1*theta_true + rnorm(n=graph_size, 0, 1e-5)
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
  Sigma_X = (1/graph_size)*Delta_inv%*%S_matrix%*%Delta_inv
  Sigma_X_matrices[[i]] = Sigma_X
  sigma_u_theta_true[i] = (1/4)*(1*Sigma_X[1,1]-2*Sigma_X[1,3]+Sigma_X[3,3])
  sigma_u_theta_true2[i] = (1/4)*(4*Sigma_X[1,1]+4*Sigma_X[1,2]+Sigma_X[2,2])
}



Xr_noise=Xhat

#sort and find one endpoint as t=0 

g=Xr_noise
g_dist = dist(g, diag=T, upper=T)
g_dist_m = as.matrix(g_dist, nrow=graph_size)
which(g_dist_m == max(g_dist_m), arr.ind = TRUE)
end1_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[1]
end2_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[2]

index_sorted=c()
l=sort(g_dist_m[end1_index,])
for (i in 1:length(l)){
sorted_index=which(g_dist_m[end1_index,] == l[i], arr.ind = TRUE)[1]
# if(sorted_index==end1_index){
#   sorted_index=which(g_dist_m[end1_index,] == l[i], arr.ind = TRUE)[2]
# }
index_sorted[i]=sorted_index
}

Xr_noise_sorted = matrix(, ncol=3, nrow=graph_size)
for (i in 1:length(l)){
  Xr_noise_sorted[i,] = Xr_noise[index_sorted[i],]
}


# t_sorted=c()
# for(i in (1:graph_size)){
#   t_sorted[i] = t[index_sorted[i]]
# }

y_sorted=c()
for(i in (1:graph_size)){
  y_sorted[i] = y[index_sorted[i]]
}


W.Xhat = Xr_noise_sorted

A1 = c(2*sum(t^4) , 4*sum((t^3)*(1-t)) , 2*sum((t^2)*((1-t)^2)))
A2=c(4*sum(t^3*(1-t)), 8*sum(t^2*(1-t)^2), 4*sum(t*(1-t)^3))
A3 = c(2*sum(t^2*(1-t)^2), 4*sum(t*(1-t)^3), 2*sum((1-t)^4))       
A = matrix(c(A1, A2, A3), byrow=T, nrow=3, ncol=3)       

d=W.Xhat[,1]
b1 = 2*sum(t^2*d)
b2 = 4*sum(t*d) - 4*sum(t^2*d)  
b3 =  2*sum(d)- 4*sum(t*d)+ 2*sum(t^2*d) 
b=rbind(b1, b2, b3)
x = solve(A)%*%b

d=W.Xhat[,2]
b1 = 2*sum(t^2*d)
b2 = 4*sum(t*d) - 4*sum(t^2*d)  
b3 =  2*sum(d)- 4*sum(t*d)+ 2*sum(t^2*d)  
b=rbind(b1, b2, b3)
y.f = solve(A)%*%b

d=W.Xhat[,3]
b1 = 2*sum(t^2*d)
b2 = 4*sum(t*d) - 4*sum(t^2*d)  
b3 =  2*sum(d)- 4*sum(t*d)+ 2*sum(t^2*d) 
b=rbind(b1, b2, b3)
z = solve(A)%*%b


x_fit = x[1]*t^2+2*x[2]*t*(1-t)+x[3]*(1-t)^2
y_fit = y.f[1]*t^2+2*y.f[2]*t*(1-t)+y.f[3]*(1-t)^2
z_fit = z[1]*t^2+2*z[2]*t*(1-t)+z[3]*(1-t)^2

X_fit=cbind(x_fit, y_fit, z_fit)
Q_fit =cbind(x, y.f, z)

 min.RSS <- function(data, par) {
   with(data, sum(((x[1]*par^2+2*x[2]*par*(1-par)+x[3]*(1-par)^2)-h[1])^2+(( y.f[1]*par^2+2*y.f[2]*par*(1-par)+y.f[3]*(1-par)^2)-h[2])^2+((z[1]*par^2+2*z[2]*par*(1-par)+z[3]*(1-par)^2)-h[3])^2))}


 theta_noise_fit <- c()
for (i in 1:nrow(X)){
  dat=data.frame(h=Xr_noise_sorted [i,])
  result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
  theta_noise_fit[i] <- result$par
}

 theta_noise_fit2 <- c()
 for (i in 1:nrow(X)){
   dat=data.frame(h=Xr_noise [i,])
   result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
   theta_noise_fit2[i] <- result$par
 }
 
 
 
 
 
 # X_true = theta_true
 # X_naive = t #1-t #**y=y_sorted
 # if (sum((t-theta_true)^2) > sum(((1-t)-theta_true)^2)){
 #   X_naive = 1-t
 # }
 # 
 
##X_naive = Xr_noise[,1]+0.5*Xr_noise[,2]

 
 theta_sorted=c()
 for(i in (1:graph_size)){
   theta_sorted[i] = theta_true[index_sorted[i]]
 }
 
 sigma_u_theta_true_sorted=c()
 for(i in (1:graph_size)){
   sigma_u_theta_true_sorted[i] =  sigma_u_theta_true[index_sorted[i]]
 }
 
 X_true = theta_true
 X_naive =  theta_noise_fit
 
 #plot(theta_sorted, theta_noise_fit)
 
  if (sum((theta_noise_fit-theta_sorted)^2) > sum(((1-theta_noise_fit)-theta_sorted)^2)){
    X_naive = 1-theta_noise_fit
   # plot(theta_sorted, X_naive)
    
  }
  
 
 #plot(theta_true, theta_noise_fit2)
 
 sigma_u_theta_fit = c()
 Sigma_X_matrices_fit = list()
 for (i in 1:graph_size){
   f = X_naive[i]
   S11 = (1)*((1/126) + (f/15)+(f^2)/42-(f^3)/105-2*(f^4)/35)
   S12 = (1)*((13/1260) + (4*f/105)-((f^2)/30)+(2*(f^3)/105)-((f^4)/70))
   S13 = (1)*((1/180)+(f/105)-((f^2)/70)+(f^3)/105-(f^4)/210)
   S22 = (1)*((1/45)+(4*f/105)-(2*f^2/35)+(4*f^3/105)-(2*f^4/105))
   S23 = (1)*((5/252)+(f/35)-13*(f^2)/210+(4*f^3)/105-(f^4/70))
   S33 = (1)*((2/63)+(f/7)-(73*(f^2)/210)+(5*(f^3)/21) - (2*(f^4)/35))
   S_matrix = matrix(c(S11, S12, S13, S12, S22, S23, S13, S23, S33), nrow=3, byrow=T)
   if(any(diag(S_matrix)<0)){print("yes")}
   Sigma_X = (1/graph_size)*Delta_inv%*%S_matrix%*%Delta_inv
   Sigma_X_matrices_fit[[i]] = Sigma_X
   sigma_u_theta_fit[i] = (1/4)*(1*Sigma_X[1,1]-2*Sigma_X[1,3]+Sigma_X[3,3])
 }
 
 
 sigma_u_theta_fit_sorted=c()
 for(i in (1:graph_size)){
   sigma_u_theta_fit_sorted[i] =  sigma_u_theta_fit[index_sorted[i]]
 }
 
 theta_noise_fit_sorted <- c()
 for (i in 1:nrow(X)){
   dat=data.frame(h=Xr_noise_sorted [i,])
   result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
   theta_noise_fit_sorted[i] <- result$par
 }
 
 X_naive_sorted =  theta_noise_fit_sorted
 if (sum((theta_noise_fit_sorted-theta_sorted)^2) > sum(((1-theta_noise_fit_sorted)-theta_sorted)^2)){
   X_naive_sorted = 1-theta_noise_fit_sorted
 }
 
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
 #plot(theta_true, sigma_u_theta_true)
 par(new=T)
 plot(theta_sorted, sigma_u_theta_fit_sorted, ylim=c(0, 0.02), col="blue")
 #plot(theta_true, sigma_u_theta_fit)
 
 
Xz_big = cbind(rep(1, graph_size), X_true)
X_naive_big = cbind(rep(1, graph_size), X_naive)

beta_true[[r]] = solve(t(Xz_big)%*%Xz_big)%*%t(Xz_big)%*% y
#beta_naive[[r]] = solve(t(X_naive_big)%*%(X_naive_big))%*%t(X_naive_big)%*% y
beta_naive[[r]] = solve(t(X_naive_big)%*%(X_naive_big))%*%t(X_naive_big)%*% y_sorted


Z = cbind(y_sorted, X_naive_big)
Sigma_U_theta_true_sorted=sigma_u_theta_true_sorted


M = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  M = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U_theta_true_sorted[i]))) + M
}
M_z = (1/graph_size)*M
beta_adj_true[[r]] = solve(M_z[2:3,2:3])%*%M_z[2:3,1]


Sigma_U_theta_fit=sigma_u_theta_fit
N = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  N = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U_theta_fit[i]))) + N
}
N_z = (1/graph_size)*N
beta_adj_naive[[r]] = solve(N_z[2:3,2:3])%*%N_z[2:3,1]


Sigma_U_theta_est=sigma_u_theta_estimated_i
Q = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  Q = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U_theta_est[i]))) + Q
}
Q_z = (1/graph_size)*Q
beta_adj_est[[r]] = solve(Q_z[2:3,2:3])%*%Q_z[2:3,1]


beta_true_se[[r]] = (beta_true[[r]]-c(1, 1))^2
beta_adj_true_se[[r]] = (beta_adj_true[[r]]-c(1, 1))^2
beta_adj_naive_se[[r]] = (beta_adj_naive[[r]]-c(1, 1))^2
beta_adj_est_se[[r]] = (beta_adj_est[[r]]-c(1, 1))^2
beta_naive_se[[r]] = (beta_naive[[r]]-c(1, 1))^2

}


b0_true = c()
b1_true = c()
b0_naive = c()
b1_naive = c()
b0_adj_true = c()
b1_adj_true = c()
b0_adj_naive = c()
b1_adj_naive = c()
b0_adj_est = c()
b1_adj_est = c()
b0_true_se = c()
b1_true_se = c()
b0_naive_se = c()
b1_naive_se = c()
b0_adj_true_se = c()
b1_adj_true_se = c()
b0_adj_naive_se = c()
b1_adj_naive_se = c()
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
  
  b0_adj_naive_se[z] = beta_adj_naive_se[[z]][1]
  b1_adj_naive_se[z] = beta_adj_naive_se[[z]][2]
  
  b0_adj_naive[z] = beta_adj_naive[[z]][1]
  b1_adj_naive[z] = beta_adj_naive[[z]][2]
  
  b0_adj_true_se[z] = beta_adj_true_se[[z]][1]
  b1_adj_true_se[z] = beta_adj_true_se[[z]][2]
  
  b0_adj_true[z] = beta_adj_true[[z]][1]
  b1_adj_true[z] = beta_adj_true[[z]][2]
  
  
  b0_adj_est[z] = beta_adj_est[[z]][1]
  b1_adj_est[z] = beta_adj_est[[z]][2]
  
  b0_adj_est_se[z] = beta_adj_est_se[[z]][1]
  b1_adj_est_se[z] = beta_adj_est_se[[z]][2]
  
  
}

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

