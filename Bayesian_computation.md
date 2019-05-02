To obtain the posterior estimates of the proposed GWR model, MCMC algorithms
should be used for sampling from the corresponding posterior distributions of
the parameters. The development of software and computation techniques simplify
the implementation of complex models. A brief introduction of some existing
programs and software can be seen in Ma and Chen (2019). In this paper, we use
the powerful `R` package **nimble** (de Valpine et al., 2017) and briefly introduce
how to implement the proposed Bayesian GWR model with `nimble` 0.7.0 in `R` 3.4.1.
A nimble model contains four parts including the model code, the constants,
the data, and the initial values for MCMC. Syntax of the model code is similar
to the `BUGS` language. For illustration, here we define S as the number of
locations, N as the number of observations, and p = 5 as the dimension of the
covariates vector. For defining the model, the function `nimbleCode()` is used:

```{r}
GWRCode <- nimbleCode({
    for (i in 1:S){
        y[1:N, i] ~ dmnorm(mu_y[1:N, i], prec = tau_y0[1:N, 1:N, i])
        mu_y[1:N, i] <- b0[i] + b[i, 1] * x1[1:N] + b[i, 2] * x2[1:N] + b[i, 3] *
            x3[1:N] + b[i, 4] * x4[1:N] + b[i, 5] * x5[1:N]
        tau_y0[1:N, 1:N, i] <- diag(psi_y[i] * exp(-Dist[1:N, i]/lambda))
```

The above code represents equation (16), where dmnorm denotes the multivariate
normal distribution with mean mu_y and precision tau_y0. In the expression of
tau_y0, psi_y[i] denotes 1/σ^2(s) in (16), and exp(-Dist[1:N, i]/lambda)
denotes the exponential weighting scheme with the distance matrix Dist and the
bandwidth lambda.

For the spike and slab prior of the coefficients beta_j(s) in (17), we can write
the code as

```{r}
        for (j in 1:5) {
            b[i, j] ~ dnorm(0, tau = tau_b[ind[j] + 1])
            ind[j] ~ dbern(gamma[j])
            tau_b[1] <- tau_b0
            tau_b[2] <- tau_b0 * c
            gamma[j] ~ dbeta(0.5, 0.5)
        }
```

The code ind[j] defines the distribution of b[i, j]. When ind[j] = 0, the precision
of b[i, j] is tau b0, otherwise the precision is tau b0 * c instead. The following
code gives the priors of the intercept, the precision of the responses, and
the bandwidth of the weighting function.

```{r}
        b0[i] ~ dnorm(0, 1)
        psi_y[i] ~ dgamma(1, 1)
    }
    lambda ~ dunif(0, D)
})


The second part is to declare the data list for the above model code, which consists
of the response variable, the covariates and the distance matrix:

```{r}
GWRData <- list(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5,
                Dist = Dist)
```

Then we can give the constant list, which contains the constant quantities in the
model code:

```{r}
GWRConst <- list(S = S, N = N, D = 50, tau_b0 = 100, c = 0.0001)
```

And lastly, we can assign some initial values for the parameters:

```{r}
GWRInits <- list(b0 = rnorm(GWRConst$S), psi_y = rep(1, GWRConst$S),
                 ind = c(1, 1, 0, 0, 0), gamma = rep(0.5, 5), lambda = 10)
```

With the above four parts, **nimble** provides an one-line function of MCMC to
invoke the MCMC engine directly, which generally take the code, data, constants,
and initial values as input, and provides a variety of options for executing and
controlling multiple chains, iterations, thinning intervals, etc.

```{r}
mcmc.out <- nimbleMCMC(code = GWRCode, data = GWRData, constants = GWRConst,
  inits = GWRInits, thin = 10, niter = 50000, nchains = 1,
  monitors = c("b0", "b", "psi_y", "lambda"), setSeed = TRUE)
```
