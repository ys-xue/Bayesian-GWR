GWRCode <- nimbleCode({
    for (i in 1:S){
        y[1:N, i] ~ dmnorm(mu_y[1:N, i], prec = tau_y0[1:N, 1:N, i])
        mu_y[1:N, i] <- b0[i] + b[i, 1] * x1[1:N] + b[i, 2] * x2[1:N] + b[i, 3] *
            x3[1:N] + b[i, 4] * x4[1:N] + b[i, 5] * x5[1:N]
        tau_y0[1:N, 1:N, i] <- diag(psi_y[i] * exp(-Dist[1:N, i]/lambda))
        for (j in 1:5) {
            b[i, j] ~ dnorm(0, tau = tau_b[ind[j] + 1])
            ind[j] ~ dbern(gamma[j])
            tau_b[1] <- tau_b0
            tau_b[2] <- tau_b0 * c
            gamma[j] ~ dbeta(0.5, 0.5)
        }
        b0[i] ~ dnorm(0, 1)
        psi_y[i] ~ dgamma(1, 1)
    }
    lambda ~ dunif(0, D)
})

GWRData <- list(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5,
                Dist = Dist)

GWRConst <- list(S = S, N = N, D = 50, tau_b0 = 100, c = 0.0001)

GWRInits <- list(b0 = rnorm(GWRConst$S), psi_y = rep(1, GWRConst$S),
                 ind = c(1, 1, 0, 0, 0), gamma = rep(0.5, 5), lambda = 10)

mcmc.out <- nimbleMCMC(code = GWRCode, data = GWRData, constants = GWRConst,
  inits = GWRInits, thin = 10, niter = 50000, nchains = 1,
  monitors = c("b0", "b", "psi_y", "lambda"), setSeed = TRUE)
