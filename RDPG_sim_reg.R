
#Idea: 1 RDPG graph with latent positions (X) on H-W curve. 
#Simulate covariance(X_hat-X) entries per Avanti's theorem to generate X_hat.
#Note: we do not generate RDPG graph here.
#compute an estimate for covariance(X_hat-X)

#convert X_hat to theta_hat 
#set up graph regression with theta_hat
#Using Avanti's theorem construct an estimate for covar(X_hat-X). Find var(theta_hat)
#use the estimate to adjust measurement error 


rm(list=ls())

library("rgl")


graph_size =  500
set.seed(752)
theta = c(0, 0.5, 1, runif(n=(graph_size-3), min=0, max=1 ))

theta_true = theta
p1_true = theta_true^2
p2_true = 2*theta_true*(1-theta_true)
p3_true = (1-theta_true)^2
P_true_col = cbind(p1_true, p2_true, p3_true)
X=P_true_col
P_true = P_true_col%*%t(P_true_col)


beta_naive = list()
beta_true = list()
beta_adj = list()
beta_true_se = list()
beta_naive_se = list()
beta_adj_se=list()

mc_runs = 50
r=1
for (r in c(1:mc_runs)){
  
y = 1 + 1*theta_true + rnorm(n=graph_size, 0, 1e-5)
X_true = theta_true #(X for regression ... not to be cofused with latent possition)


#Generate X_hat (Xr_noise):
#Simulate covaraince(X_hat-X) per Avanti
Delta_inv = matrix(c(9, -9, 3, -9, 21, -9, 3, -9, 9), nrow=3)

X_noise = matrix(, ncol=3, nrow=graph_size)
sigma_u_theta_true = c()
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
  sigma_u_theta_true[i] = (1/4)*(1*Sigma_X[1,1]-2*Sigma_X[1,3]+Sigma_X[3,3])
  X_noise[i,1:3] = X[i,1:3]  + mvrnorm(n=1, rep(0, 3), Sigma=Sigma_X)
}

#generate an orthogonal matrix and rotate the X_hat by this matrix (Q)
Q=(1/3)*matrix(c(2, 1, 2, -2, 2, 1, 1, 2, -2), nrow=3)
Xr_noise = X_noise%*%Q

#rotate Xr_noise and project onto HW curve: 
Xr= X%*%Q

#Take Xr_noise onto HW, call this, X_naive 
#Construct Q_est
g=Xr_noise
g_dist = dist(g, diag=T, upper=T)
g_dist_m = as.matrix(g_dist, nrow=graph_size)
which(g_dist_m == max(g_dist_m), arr.ind = TRUE)
end1_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[1]
end2_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[2]
foo = g_dist_m[end1_index,]-g_dist_m[end2_index,]
mid_index_z = which(foo==min(abs(foo)))
if (is.integer(mid_index_z) && length(mid_index_z)==0) {mid_index_z = which(foo==-1*min(abs(foo)))[[1]]}

Xr_noise_pts = rbind(Xr_noise[end1_index,], Xr_noise[mid_index_z,], Xr_noise[end2_index,])

X_pts = rbind(X[1,], X[2,], X[3,])
Xr_pts = rbind(Xr[1,], Xr[2,], Xr[3,])

Q_transpose_est_noise = t(solve(t(X_pts)%*%X_pts)%*%t(X_pts)%*%Xr_noise_pts)

#Rotate X_hat
Xr_noise_simp = Xr_noise%*%Q_transpose_est_noise

#Project onto H-W
min.RSS <- function(data, par) {
  with(data, sum(((par^2-y[1])^2 + (2*par*(1-par)-y[2])^2 + ((1-par)^2-y[3])^2)^2))}

theta_noise_hw <- c()
for (i in 1:nrow(X)){
  dat=data.frame(y=Xr_noise_simp [i,])
  result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
  theta_noise_hw[i] <- result$par
}

#CHEAT:
flip = matrix(c(0, 0, 1, 0, 1, 0, 1, 0, 0), byrow=T, nrow=3)
X_naive = theta_noise_hw
if(sum((theta_noise_hw-theta_true)^2)>sum(((1-theta_noise_hw)-theta_true)^2)){
  X_naive = 1-theta_noise_hw
  theta_noise_hw = 1-theta_noise_hw
  print("1")
  Xr_noise_simp = (Xr_noise_simp) %*% flip
}

#Construct estimate for covaraince (X_hat-X)
sigma_u_theta_est_hw = c()
for (i in 1:graph_size){
  f = theta_noise_hw[i]
  S11 = (1)*((1/126) + (f/15)+(f^2)/42-(f^3)/105-2*(f^4)/35)
  S12 = (1)*((13/1260) + (4*f/105)-((f^2)/30)+(2*(f^3)/105)-((f^4)/70))
  S13 = (1)*((1/180)+(f/105)-((f^2)/70)+(f^3)/105-(f^4)/210)
  S22 = (1)*((1/45)+(4*f/105)-(2*f^2/35)+(4*f^3/105)-(2*f^4/105))
  S23 = (1)*((5/252)+(f/35)-13*(f^2)/210+(4*f^3)/105-(f^4/70)) 
  S33 = (1)*((2/63)+(f/7)-(73*(f^2)/210)+(5*(f^3)/21) - (2*(f^4)/35))
  S_matrix = matrix(c(S11, S12, S13, S12, S22, S23, S13, S23, S33), nrow=3, byrow=T)
  if(any(diag(S_matrix)<0)){print("ALERT")}
  Sigma_X_est_hw = (1/graph_size)*Delta_inv%*%S_matrix%*%Delta_inv
  sigma_u_theta_est_hw[i] = (1/4)*(1*Sigma_X_est_hw[1,1]-2*Sigma_X_est_hw[1,3]+Sigma_X_est_hw[3,3])
}

Xz_big = cbind(rep(1, graph_size), X_true)
X_naive_big = cbind(rep(1, graph_size), X_naive)
beta_true[[r]] = solve(t(Xz_big)%*%Xz_big)%*%t(Xz_big)%*% y
beta_naive[[r]] = solve(t(X_naive_big)%*%(X_naive_big))%*%t(X_naive_big)%*% y


Z = cbind(y, X_naive_big)

Sigma_U_theta=sigma_u_theta_est_hw

M = matrix(0,nrow=3,ncol=3)
for (i in (1:graph_size)){
  M = (Z[i,]%o%Z[i,]-diag(c(0, 0, Sigma_U_theta[i]))) + M
}
M_z = (1/graph_size)*M
beta_adj[[r]] = solve(M_z[2:3,2:3])%*%M_z[2:3,1]

beta_true_se[[r]] = (beta_true[[r]]-c(1, 1))^2
beta_adj_se[[r]] = (beta_adj[[r]]-c(1, 1))^2
beta_naive_se[[r]] = (beta_naive[[r]]-c(1, 1))^2

}

b0_true = c()
b1_true = c()
b0_naive = c()
b1_naive = c()
b0_adj = c()
b1_adj = c()
b0_true_se = c()
b1_true_se = c()
b0_naive_se = c()
b1_naive_se = c()
b0_adj_se = c()
b1_adj_se = c()


for (z in 1:mc_runs){
  
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
 
  b0_adj_se[z] = beta_adj_se[[z]][1]
  b1_adj_se[z] = beta_adj_se[[z]][2]
 
  b0_adj[z] = beta_adj[[z]][1]
  b1_adj[z] = beta_adj[[z]][2]

  
}

###########
###########
###########
par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b0_true, b0_naive,  b0_adj, notch=TRUE, 
        main=bquote(paste("Single Graph Regression", "\n b0 Estimate", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 1, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true, b1_naive,  b1_adj, notch=TRUE, 
        main=bquote(paste("Single Graph Regression", "\n b1 Estimate", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 1, col="red", lty=2, lwd=0.5)

df=data.frame(rbind(X, Xr_noise,  Xr_noise_simp, Xr_noise_hw ))
df$type = c(rep(1, graph_size), rep(2, graph_size),  rep(3, graph_size), rep(4, graph_size))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")


par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b0_true_se, b0_naive_se,  b0_adj_se,notch=TRUE, 
        main=bquote(paste("Single Graph Regression", "\n b0 Square Error", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)

par(mar=c(7.1,4.1,4.1,2.1))
boxplot(b1_true_se, b1_naive_se,  b1_adj_se,notch=TRUE, 
        main=bquote(paste("Single Graph Regression", "\n b1 Square Error", "\n graph_size:", graph_size , "\n mc_runs:", mc_runs)), 
        cex.main=0.5,
        las=2, names=c("True", "Naive", "Adjusted"))
abline(h = 0, col="red", lty=2, lwd=0.5)


#CI 
t.test(b0_true)$conf.int
t.test(b0_naive)$conf.int
t.test(b0_adj)$conf.int

t.test(b1_true)$conf.int
t.test(b1_naive)$conf.int
t.test(b1_adj)$conf.int


#SE (sign test )
N = sum((b0_adj_se - b0_naive_se) !=0)
k = sum(b0_adj_se < b0_naive_se)
1 - pbinom(k, size=N, prob=0.5)

N = sum((b1_adj_se - b1_naive_se) !=0)
k = sum(b1_adj_se < b1_naive_se)
1 - pbinom(k, size=N, prob=0.5)



