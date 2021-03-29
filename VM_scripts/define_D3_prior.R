# determining alpha and beta for D3 prior distribution

mean = 0.1019
sd = 0.025 # 0.01 for all species but red maple, white pine, and eastern hemlock right now which has sd of 0.025

alpha = (((1-mean)/sd^2) - (1/mean)) * mean^2
beta = alpha * ((1/mean) - 1)

print(paste('alpha is:',round(alpha,3)))
print(paste('beta is:',round(beta,3)))

x = seq(0,1,0.001)
#y1 = dbeta(x, 19.378,410.3)
y2 = dbeta(x, alpha, beta)

24.351,455.944

plot(NULL, xlim = c(0,0.2), ylim = c(0,50))
lines(x,y1)
lines(x,y2, col='red')

