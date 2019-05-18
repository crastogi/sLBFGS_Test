# Generate bvn samples
sample.size = 1E4
dat = cbind(rnorm(n = sample.size, sd = .1, mean = 5), rnorm(n = sample.size, sd = .1, mean = -2))
write.table(x = dat, file = "~/Downloads/probLS/bivariate_normal.tsv", row.names = FALSE, col.names = FALSE, sep = "\t")

# Generate branin function samples
g0=5.1/(4*pi^2)
g1=10*(1-1/(8*pi))
g2=5/pi
branin.funcion = function(x) {
  return(((x[2]+2) - g0*(x[1]-3)^2 + g2*(x[1]-3) - 6)^2 + g1*cos(x[1]-3) + 10)
}

burn.in = 1E4
samples = 1E6
output = cbind(rep(0, samples), rep(0,samples))

# init
x0 = runif(n=1, min = -2, max = 13)
y0 = runif(n=1, min = -2, max = 13)
z0 = branin.funcion(c(x0,y0))

for (i in 1:burn.in) {
  # Propose new x,y
  x = min(max(-2, x0+rnorm(1,0,1)), 13)
  y = min(max(-2, y0+rnorm(1,0,1)), 13)
  z = branin.funcion(c(x,y))
  a = z/z0
  if (runif(1)<a) {
    # Accept
    x0 = x;
    y0 = y;
    z0 = z;
  }
}

# record samples 
for (i in 1:samples) {
  # Propose new x,y
  x = min(max(-2, x0+rnorm(1,0,1)), 13)
  y = min(max(-2, y0+rnorm(1,0,1)), 13)
  z = branin.funcion(c(x,y))
  a = z/z0
  if (runif(1)<a) {
    # Accept
    x0 = x;
    y0 = y;
    z0 = z;
  }
  output[i,] = c(x0,y0)
}

# Demonstrate that sampling is accurate
library(hexbin)
df = data.frame(output)
hexbinplot(X2~X1, data=df, trans=log, inv=exp)

# Subset
branin.output = output[sample(1:nrow(output), sample.size),]
df = data.frame(branin.output)
hexbinplot(X2~X1, data=df, trans=log, inv=exp)

# Save
write.table(x = branin.output, file = "~/Downloads/probLS/braninsample.tsv", row.names = FALSE, col.names = FALSE, sep = "\t")
