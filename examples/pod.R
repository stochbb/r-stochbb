library(stochbb)

#
# Processing model
#
d  <- 500
V  <- gamma(5,30)
S1 <- gamma(10,50)
S2 <- gamma(1,120)
M  <- gamma(1,150) 

# response latency control condition
#   is simply R = C + S1 + M
# The resulting distribution of response latecies is an analytic function on (0, infinity).
Rc <- V %+% S1 %+% M

# And for the experimental condition 
#   R = V + min(S1, S2) + M.
# This might be considered as a two parallel proceses S1 (identical to the S1 under control
# condition) and S2 which is delayed by d=300 ms but once started, much faster than S1.
# The resulting distribution is not analytic on (0, infinity)! It is identical to the distribution
# of R under control condition (Rc) on the interval (0, d] and diverges for T>d.
Re <- V %+% minimum(S1, affine(S2, 1, d)) %+% M

Tmin <- 0; Tmax <- 1500; N <- 1000;
t   <- seq(Tmin, Tmax, length.out=N);
Tc  <- array(0, c(N)); Rc$density$eval(Tmin, Tmax, Tc)
Te  <- array(0, c(N)); Re$density$eval(Tmin, Tmax, Te)

# Check with exact sampler
sam <- new(ExactSampler, c(Rc, Re))
res <- array(0, c(10000,2))
sam$sample(res)

# PDFs
plot(t, Te, type="l", lty=2)
lines(t, Tc, lty=1)
lines(density(res[,2]), lty=2)
lines(density(res[,1]), lty=1)
abline(v=d, lty=3)

