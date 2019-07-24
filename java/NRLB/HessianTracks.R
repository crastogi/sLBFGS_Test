ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
# Load data
data = read.table("/Users/chaitanya/Documents/GitWorkspaces/slbfgs/java/output/RawPath_2.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
data = data[,1:(ncol(data)-1)]
# Set batchsize
bs = 500 #(no need to normalize for sLBFGS)

# Create index ranges. First four columns are 'info' ranges
info.idx = 1:4
# Fields for path: x_t, g_t, f_it, a_it, mu_k, effGrad, egNorm. We do not know dimensionality d beforehand, so infer it
d = (ncol(data[,5:ncol(data)])-1)/6
start.idx = 5
x_t.idx = start.idx:(start.idx+d-1)
start.idx = start.idx+d
g_t.idx = start.idx:(start.idx+d-1)
start.idx = start.idx+d
f_it.idx = start.idx:(start.idx+d-1)
start.idx = start.idx+d
a_it.idx = start.idx:(start.idx+d-1)
start.idx = start.idx+d
mu_k.idx = start.idx:(start.idx+d-1)
start.idx = start.idx+d
eg.idx = start.idx:(start.idx+d-1)
egNorm.idx = ncol(data)

# Remove svrg update lines
svrg.update.idx = grep(pattern = ">", x=data[,3])
svrg.updates = data[svrg.update.idx,]
pure.traj = data[-svrg.update.idx,]

# Compute function position, g_t norm, f_it norm, a_it norm
xNorm = sqrt(rowSums(pure.traj[,x_t.idx]^2))
gNorm = sqrt(rowSums(pure.traj[,g_t.idx]^2))
f_itNorm = sqrt(rowSums(pure.traj[,f_it.idx]^2))
a_itNorm = sqrt(rowSums((pure.traj[,a_it.idx]/bs)^2))
egNorm = pure.traj[,egNorm.idx]

# Important: diff norm
diffNorm = sqrt(rowSums((pure.traj[,f_it.idx]-pure.traj[,a_it.idx]/bs)^2))
full_precision_diff = sqrt(rowSums((pure.traj[,g_t.idx]-(pure.traj[,f_it.idx]-pure.traj[,a_it.idx]+pure.traj[,mu_k.idx]))^2))

# Hessian update positions and values
temp = pure.traj[, info.idx]
hess.pos = which(!is.na(temp$V4))
prev.idx = 1
hess.update = NULL
while (TRUE) {
  idx = which(!is.na(temp$V4))
  if (length(idx)==0) {
    break
  }
  hess.update = c(hess.update, rep(temp$V4[idx[1]], length(prev.idx:idx[1])))
  prev.idx = idx[1]+1
  temp$V4[idx[1]] = NA
}
# SVRG update positions
temp = data[,info.idx]
svrg.update = NULL
while(TRUE) {
  idx = which(temp$V3==">")
  if (length(idx)==0) {
    break
  }
  svrg.update = c(svrg.update, idx[1]-1)
  temp = temp[-idx[1],]
}
# Likelihood
lik.vals = data[which(data$V3==">"),ncol(data)]
prev.idx = 1
lik = NULL
for (i in 2:length(svrg.update)) {
  lik = c(lik, rep(lik.vals[i-1], length(prev.idx:svrg.update[i])))
  prev.idx = svrg.update[i]+1
}
lik = log((lik - min(lik))*10000+.01)

# Plot trajectory
#matplot(cbind(xNorm, gNorm, f_itNorm, a_itNorm, egNorm), type="l", lty = 2, col = seq_len(5))
#legend("topright", legend = c("x", "g", "f_it", "a_it", "effGrad"),col=seq_len(5),lty=2)

pdf("~/Desktop/Trajectory_Raw.pdf",width=80,height=8) 
matplot(log10(cbind(xNorm, gNorm, f_itNorm, a_itNorm, egNorm)), type="l", lty = 1, col = seq_len(5))
lines(log10(hess.update), lty=1, col=6)
lines(lik, lty=1, col=7)
abline(v = hess.pos, col="red", lty=2)
abline(v = svrg.update, col="black", lty=4)
legend("topright", legend = c("x", "g", "f_it", "a_it", "effGrad", "u_r dist", "lik"),col=seq_len(7),lty=1, ncol=2)
dev.off()

pdf("~/Desktop/Trajectory.pdf",width=80,height=8) 
matplot(log10(cbind(xNorm, gNorm, f_itNorm, a_itNorm, ma(egNorm, 25))), type="l", lty = 1, col = seq_len(5))
lines(log10(hess.update), lty=1, col=6)
lines(lik, lty=1, col=7)
abline(v = hess.pos, col="red", lty=2)
abline(v = svrg.update, col="black", lty=4)
legend("topright", legend = c("x", "g", "f_it", "a_it", "effGrad", "u_r dist", "lik"),col=seq_len(7),lty=1, ncol=2)
dev.off()

pdf("~/Desktop/Trajectory_Raw.pdf",width=80,height=8) 
matplot(log10(cbind(xNorm, gNorm, f_itNorm, a_itNorm, egNorm)), type="l", lty = 1, col = seq_len(5))
lines(log10(hess.update), lty=1, col=6)
lines(lik, lty=1, col=7)
lines(log10(diffNorm), lty=1, col=8)
abline(v = hess.pos, col="red", lty=2)
abline(v = svrg.update, col="black", lty=4)
abline(h = 0, col="black", lty=5)
legend("topright", legend = c("x", "g", "f_it", "a_it", "effGrad", "u_r dist", "lik", "diffNorm"),col=seq_len(8),lty=1, ncol=2)
dev.off()

matplot(cbind(gNorm, f_itNorm, a_itNorm), type="l", lty = 1, col = seq_len(3))
abline(v = hess.update, col="red", lty=3)
abline(v = svrg.update, col="black", lty=4)
legend("topright", legend = c("g", "f_it", "a_it"),col=seq_len(3),lty=1, ncol=2)

# Unused
update.idx = grep(pattern = "NewXBar", x=subset[,1])
pure.traj = subset[-update.idx,]
