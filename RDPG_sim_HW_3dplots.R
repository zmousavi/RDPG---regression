#Check plots for different rotation matrix (W) estimation methods 

#Idea: 1 RDPG graph with latent positions (X) on H-W curve. 
#Simulate covariance(X_hat-X) entries per Avanti's theorem to generate X_hat.
#Note: we do not generate RDPG graph here.
#Rotate the generated X_hat with a rotation matrix, Q

#Find estimate for Q: Q_est (#Trying different methods to find Q_est)
#Rotate back the X_hat
#Project the rotated X_hat to H-W curve
#plot



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

#Simulate covaraince(X_hat-X) per Avanti

Delta_inv = matrix(c(9, -9, 3, -9, 21, -9, 3, -9, 9), nrow=3)

X_noise = matrix(, ncol=3, nrow=graph_size)
sigma_u_theta_true1 = c()
sigma_u_theta_true2 = c()
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
  #sigma_u_theta_true1[i] = (1/4)*(4*Sigma_X[1,1]+4*Sigma_X[1,2]+Sigma_X[2,2])
  sigma_u_theta_true2[i] = (1/4)*(1*Sigma_X[1,1]-2*Sigma_X[1,3]+Sigma_X[3,3])
  X_noise[i,1:3] = X[i,1:3]  + mvrnorm(n=1, rep(0, 3), Sigma=Sigma_X)
}
sigma_u_theta_true = sigma_u_theta_true2

#generate an orthogonal matrix and rotate the X_hat by this matrix (Q)
Q=(1/3)*matrix(c(2, 1, 2, -2, 2, 1, 1, 2, -2), nrow=3)
Xr_noise = X_noise%*%Q
Xr= X%*%Q

#Find the 2 points that are furthest from each other 
#Find the point that is equi-distant from these end points (call this mid_index)
#Try: Find the average of nearest neighbor of point[mid_index] and use this point instead
#Once three X_hat's are found, compute the rotation matrix that maps them to H-W
#Once X_hat is rotated, find closes projection onto H-W curve/
#Cheat: check if X_hat need to be flipped. If so, also change theta to (1-theta)

g=Xr_noise 
g_dist = dist(g, diag=T, upper=T)
g_dist_m = as.matrix(g_dist, nrow=graph_size)  
which(g_dist_m == max(g_dist_m), arr.ind = TRUE) 
end1_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[1]
end2_index = which(g_dist_m == max(g_dist_m), arr.ind = TRUE)[2]
foo = g_dist_m[end1_index,]-g_dist_m[end2_index,]
mid_index_z = which(foo==min(abs(foo)))
if (is.integer(mid_index_z) && length(mid_index_z)==0) {mid_index_z = which(foo==-1*min(abs(foo)))[[1]]} 

ker_minz=c()
minz = sort(abs(foo))[1:50]
for (i in c(1:length(minz))){
  mid_index_minz = which(foo==minz[[i]])
  if (is.integer(mid_index_minz) && length(mid_index_minz)==0) {mid_index_minz = which(foo==-1*minz[[i]])} 
  ker_minz[i]=mid_index_minz[[1]]
}

Xr_ker_min = colMeans(Xr_noise[ker_minz,])

Xr_noise_pts1 = rbind(Xr_noise[end1_index,], Xr_noise[mid_index_z,], Xr_noise[end2_index,])
Xr_noise_pts2 = rbind(Xr_noise[end1_index,], Xr_ker_min, Xr_noise[end2_index,])

X_pts = rbind(X[1,], X[2,], X[3,])
Xr_pts = rbind(Xr[1,], Xr[2,], Xr[3,])

#Check method to find estimate of Q:
Q_transpose_est = t(solve(t(X_pts)%*%X_pts)%*%t(X_pts)%*%Xr_pts)
Xr_simp = Xr%*%Q_transpose_est

#Construct Q_estimate 
#Rotate X_hat 
Q_transpose_est_noise1 = t(solve(t(X_pts)%*%X_pts)%*%t(X_pts)%*%Xr_noise_pts1)
Xr_noise_simp1 = Xr_noise%*%Q_transpose_est_noise1

# #should i do procrustes instead?
# proc=svd(t(Xr_noise_pts1)%*%X_pts)
# Q_transpose_est_proc = proc$u%*%t(proc$v)
# Xr_noise_simp3 = Xr_noise%*%Q_transpose_est_proc


min.RSS <- function(data, par) {
  with(data, sum(((par^2-y[1])^2 + (2*par*(1-par)-y[2])^2 + ((1-par)^2-y[3])^2)^2))}

theta_noise_hw1 <- c()
for (i in 1:nrow(X)){
  dat=data.frame(y=Xr_noise_simp1 [i,])
  result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
  theta_noise_hw1[i] <- result$par
}

#CHEAT:
flip = matrix(c(0, 0, 1, 0, 1, 0, 1, 0, 0), byrow=T, nrow=3)
X_naive1 = theta_noise_hw1
if(sum((theta_noise_hw1-theta_true)^2)>sum(((1-theta_noise_hw1)-theta_true)^2)){
  X_naive1 = 1-theta_noise_hw1
  theta_noise_hw1 = 1-theta_noise_hw1
  print("1")
  Xr_noise_simp1 = (Xr_noise_simp1) %*% flip
}

Xr_noise_hw1 = t(rbind(theta_noise_hw1^2, 2*theta_noise_hw1*(1-theta_noise_hw1), (1-theta_noise_hw1)^2))

# ####
# theta_noise_hw2 <- c()
# for (i in 1:nrow(X)){
#   dat=data.frame(y=Xr_noise_simp2 [i,])
#   result <- optim(par = 0.5, min.RSS, data = dat, method = "L-BFGS-B", lower = 0, upper = 1.0)
#   theta_noise_hw2[i] <- result$par
# }
# 
# 
# #CHEAT:
# X_naive2 = theta_noise_hw2
# if(sum((theta_noise_hw2-theta_true)^2)>sum(((1-theta_noise_hw2)-theta_true)^2)){
#   X_naive2 = 1-theta_noise_hw2
#   theta_noise_hw2 = 1-theta_noise_hw2
#   print("1")
#   Xr_noise_simp2 = (Xr_noise_simp2) %*% flip
# }
# 
# Xr_noise_hw2 = t(rbind(theta_noise_hw2^2, 2*theta_noise_hw2*(1-theta_noise_hw2), (1-theta_noise_hw2)^2))

# #averaging neighbors of mequi-distant point 
# Q_transpose_est_noise2 = t(solve(t(X_pts)%*%X_pts)%*%t(X_pts)%*%Xr_noise_pts2)
# Xr_noise_simp2 = Xr_noise%*%Q_transpose_est_noise2
# 


##3d plot 
df=data.frame(rbind(X, Xr_noise,  Xr_noise_simp1, Xr_noise_hw1 ))
df$type = c(rep(1, graph_size), rep(2, graph_size),  rep(3, graph_size), rep(4, graph_size))
with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")
decorate3d( sub="1")
#text3d(x=1, y=1, z=1, "np_reg" ,col="Black")

# 
# df=data.frame(rbind(X, Xr_noise,  Xr_noise_simp2, Xr_noise_hw2 ))
# df$type = c(rep(1, graph_size), rep(2, graph_size),  rep(3, graph_size), rep(4, graph_size))
# with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
# triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")
# decorate3d( sub="2")


# df=data.frame(rbind(X, Xr_noise_simp1,  Xr_noise_simp2 ))
# df$type = c(rep(1, graph_size), rep(2, graph_size),  rep(3, graph_size))
# with(df, plot3d(df[,1:3],type="s",xlab="",ylab="",zlab="",size=1/2 ,col=as.numeric(type)))
# triangles3d(c(1,0,0),c(0,0,1),c(0,1,0),alpha=.5,col="lightblue")



#Y = rep(1,nrow(Xr_noise))
# Y[mid_index]=2
# rgl(Xr_noise, col=Y)
# plot3d(Xr_noise, col=Y)
#plot(1:10, col=1:10, pch=19, cex=2)

