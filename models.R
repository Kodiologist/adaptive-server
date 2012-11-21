library(rstan)

# ------------------------------------------------------------
# Functions for defining models
# ------------------------------------------------------------

empty.ts = data.frame(
    ssr = numeric(0), ssd = numeric(0),
    llr = numeric(0), lld = numeric(0),
    choice = factor(character(0), levels = qw(ss, ll)))

grid.samples = 150

prev.ts = NULL
prev.lhood = NULL

grid.approx.model = function(sample.thetas, prior, choice.p)
# - choice.p should be vectorizable over trials and over model
#   parameters, but need not handle the case of multiple trials
#   and multiple parameter vectors at the same time.
# - The posterior will never be evaluated at the endpoints
#   of each vector in sample.thetas. This means you can begin or
#   end the vectors with asymptotes.
   {orig.sample.thetas = sample.thetas
    sample.thetas = lapply(sample.thetas, function (v) v[-1])
    pn = length(sample.thetas)
    sample.lens = sapply(sample.thetas, length) - 1
    thetas.df = as.matrix(do.call(expand.grid, lapply(sample.thetas,
        function (v) head(v, -1))))
    prior.masses =
        # Prior density
        prior(thetas.df) *
        # Volume
        maprows(f = prod, as.matrix(do.call(expand.grid, lapply(sample.thetas, diff))))
    jitter.in.thetas.df = function(rows)
        do.call(cbind, lapply(1 : pn, function (p)
           {is = expand.grid.idx(sample.lens, rows, p)
            runif(length(is),
                sample.thetas[[p]][is],
                sample.thetas[[p]][is + 1])}))
    unnorm.likelihood = function(ts, thetas)
        maprows(f = prod, maybeparallel.sapply(1 : nrow(ts), function (trial)
           {ps = do.call(choice.p, c(
                list(ts[trial,]),
                lapply(1 : pn, function (p) thetas[,p])))
            if (ts[trial, "choice"] == "ll") ps else 1 - ps}))
    sample.posterior = function(ts, n = grid.samples)
       {lhood =
            if (nrow(ts) == 0)
               1
            else if (identical(ts, prev.ts))
               prev.lhood
            else if (identical(ts[1 : (nrow(ts) - 1),], prev.ts))
               prev.lhood * unnorm.likelihood(ts[nrow(ts),], thetas.df)
            else
               unnorm.likelihood(ts, thetas.df)
        prev.ts <<- ts
        prev.lhood <<- lhood
        posterior = jitter.in.thetas.df(
          # Sample some rows of theta.df, weighted by their
          # posterior masses.
            sample.int(nrow(thetas.df), n, rep = T, prob =
                lhood * prior.masses))
        dimnames(posterior) = list(c(), names(sample.thetas))
        posterior}
    rand.theta = function()
        drop(sample.posterior(empty.ts, n = 1))
    punl(
        sample.thetas = orig.sample.thetas,
        nrow.thetas.df = nrow(thetas.df),
        prior, choice.p, sample.posterior, rand.theta)}

prior.uniform = function (theta) 1

grid.sr.rho = grid.approx.model(
    sample.thetas = list(
        f = seq(-1, 1, len = 200),
        rho = seq(0, 1, len = 200)),
    prior = prior.uniform,
    choice.p = function(ts, f, rho)
       {gamma = 1 / (100 * rho)
        tau = exp(10 * f) * gamma
        ilogit(
            (log(1 + gamma * ts$llr) - log(1 + gamma * ts$ssr))/gamma -
            (log(1 + tau * ts$lld) - log(1 + tau * ts$ssd))/tau)})

gelman.diag.threshold = 1.1

stan.samples = 50
stan.chains = 7

mkdat.stan = function(ts)
    list(
        N_ts = nrow(ts),
        choose_ll = as.integer(ts$choice == "ll"),
        ssr = ts$ssr,
        ssd = ts$ssd,
        llr = ts$llr,
        lld = ts$lld)

stan.choicemodel = function(
        choice_ll_p = NULL, choice_ll_p_logit = NULL,
        choice.p,
        parameters, transformed_parameters = '', prior = '',
        monitor = NULL, init, n.adapt = 150, thin = 1,
        max.mcmc.rounds = 9)
# N.B. choice.p should be vectorizable over either argument, but
# need not handle the case of multiple trials and multiple
# parameter vectors at the same time.
   {std.init = init
    if (is.null(monitor))
        monitor = names(if (mode(init) == "function") init(1) else init[[1]])
    model.str = sprintf(
       'data
           {int <lower = 0> N_ts;
            vector<lower = 0>[N_ts] ssr;
            vector<lower = 0>[N_ts] ssd;
            vector<lower = 0>[N_ts] llr;
            vector<lower = 0>[N_ts] lld;
            int <lower = 0, upper = 1> choose_ll[N_ts];}
        parameters
           {%s}
        transformed parameters
           {%s}
        model
           {%s
            for (t in 1 : N_ts)
              {choose_ll[t] ~ %s;
               /*print("$", ssr[t], " in ", ssd[t], "d vs. ",
                     "$", llr[t], " in ", lld[t], "d: ",
                     choose_ll_p_logit[t]);*/}}',
        parameters,
        transformed_parameters,
        prior,
        (if (is.null(choice_ll_p_logit))
            sprintf("bernoulli(%s)", choice_ll_p) else
            sprintf("bernoulli_logit(%s)", choice_ll_p_logit)))
    sample.posterior = function(ts, init = std.init, debugging = F)
       {current.thin = thin
        current.adapt = n.adapt
        model = cached_stan_model(model.str)
        round = 1
        subround = 1
        repeat
           {capture.output(fit <- sampling(model, refresh = -1,
                data = mkdat.stan(ts),
                init = init,
                chains = stan.chains,
                iter = current.adapt + current.thin * stan.samples,
                warmup = current.adapt,
                thin = current.thin))
            if (fit@mode == 2)
              # An error occurred. If this is a NaN or step-size
              # error, we can probably avoid it by just trying again.
               {if (subround >= 10)
                 # Give up.
                   return(NULL)
                subround = subround + 1
                next}
            if (debugging)
               {ex = rstan::extract(fit, monitor, perm = F)
                for (n in 1 : dim(ex)[[3]])
                    {mat = as.matrix(mapcols(drop(ex[,,n]), function (v)
                        round(quantile(v, c(.025, .25, .5, .75, .975)), 3)))
                     names(dimnames(mat)) = c("", monitor[[n]])
                     print(mat)}}
            rhats = summary(fit, pars = monitor, use_cache = F)$summary[,"Rhat"]
            if (all(rhats < gelman.diag.threshold))
                break
            current.thin = 2 * current.thin
            current.adapt = 2 * current.adapt
            round = round + 1
            subround = 1
            message(sprintf("r%d: Rhats %s; next: thin %d, adapt %d",
                as.integer(round - 1),
                paste(collapse = ", ", round(rhats, 2)),
                as.integer(current.thin),
                as.integer(current.adapt)))
            if (round > max.mcmc.rounds)
               {message("Using this sample anyway")
                break}}
        simplify2array(rstan::extract(fit, monitor, perm = T))}
    precompile = function()
        cached_stan_model(model.str)
    rand.theta = function()
        simplify2array(std.init(1))
    punl(
        monitor, choice.p, model.str,
        init = std.init,
        precompile, sample.posterior, rand.theta)}

rlunif = function(n, min = 1e-10, max = 1)
   exp(runif(n, log(min), log(max)))

# ------------------------------------------------------------
# Some model definitions
# ------------------------------------------------------------

model.rewards = stan.choicemodel(
    choice_ll_p_logit =
           'b_ssr * ssr[t] +
            b_llr * llr[t]',
    choice.p = function(ts, b_ssr, b_llr)
        ilogit(
            b_ssr * ts$ssr +
            b_llr * ts$llr),
    parameters =
       'real <lower = -5, upper = 0> b_ssr;
        real <lower = 0, upper = 5> b_llr;',
        # Implicit uniform priors.
    init = function(n) list(
        b_ssr = -min(5, abs(rnorm(1, 0, .25))),
        b_llr = min(5, abs(rnorm(1, 0, .25)))))

model.rewards.logprior = stan.choicemodel(
    choice_ll_p_logit =
           'b_ssr * ssr[t] +
            b_llr * llr[t]',
    choice.p = function(ts, b_ssr, b_llr)
        ilogit(
            b_ssr * ts$ssr +
            b_llr * ts$llr),
    parameters =
       'real <lower = -5, upper = 0> b_ssr;
        real <lower = 0, upper = 5> b_llr;',
    prior =
      # exp([-10, 1.6]) is about [.000045, 5].
       'log(-b_ssr) ~ uniform(-10, 1.6);
        log(b_llr) ~ uniform(-10, 1.6);',
    init = function(n) list(
        b_ssr = -rlunif(1, .0001, 1.5),
        b_llr = rlunif(1, .0001, 1.5)))

model.diff = stan.choicemodel(
    choice_ll_p_logit =
           'b_rd * (llr[t] - ssr[t]) -
            b_dd * (lld[t] - ssd[t])',
    choice.p = function(ts, b_rd, b_dd)
        ilogit(
            b_rd * (ts$llr - ts$ssr) -
            b_dd * (ts$lld - ts$ssd)),
    parameters =
       'real <lower = 0, upper = 5> b_rd;
        real <lower = 0, upper = 5> b_dd;',
        # Implicit uniform priors.
    init = function(n) list(
        b_rd = runif(1, 0, 5),
        b_dd = runif(1, 0, 5)))

model.diff.delta = stan.choicemodel(
    choice_ll_p =
           'delta * .5 + (1 - delta) * step(
                b_rd * (llr[t] - ssr[t]) -
                b_dd * (lld[t] - ssd[t]))',
    choice.p = function(ts, b_rd, b_dd, delta)
        delta * .5 + (1 - delta) * (
            b_rd * (ts$llr - ts$ssr) >
            b_dd * (ts$lld - ts$ssd)),
    parameters =
       'real <lower = 0, upper = 5> b_rd;
        real <lower = 0, upper = 5> b_dd;
        real <lower = 0, upper = .5> delta;',
        # Implicit uniform priors.
    init = function(n) list(
        b_rd = runif(1, 0, 5),
        b_dd = runif(1, 0, 5),
        delta = runif(1, 0, .5)))

model.diff.rho = stan.choicemodel(
    choice_ll_p_logit =
           'rho_v *
            ((llr[t] - ssr[t]) -
             a * (lld[t] - ssd[t]))',
    parameters =
       'real <lower = -1, upper = 1> f;
        real <lower = 0, upper = 1> rho;',
        # Implicit uniform priors.
    transformed_parameters =
       'real rho_v; real a;
        rho_v <- 10 * rho;
        a <- exp(10 * f);',
    choice.p = function(ts, f, rho)
        ilogit(10 * rho *
            ((ts$llr - ts$ssr) -
            exp(10 * f) * (ts$lld - ts$ssd))),
    init = function(n) list(
        f = runif(1, -1, 1),
        rho = runif(1, 0, 1)))

model.sr = stan.choicemodel(
    choice_ll_p_logit =
           '(log(1 + gamma * llr[t]) - log(1 + gamma * ssr[t]))/gamma -
            (log(1 + tau * lld[t]) - log(1 + tau * ssd[t]))/tau',
    parameters =
       'real <lower = -10, upper = 5> ln_gamma;
        real <lower = -10, upper = 5> ln_tau;',
        # Implicit uniform priors.
    transformed_parameters =
       'real gamma; real tau;
        gamma <- exp(ln_gamma);
        tau <- exp(ln_tau);',
    # The logarithmic scale of the parameters and the upper bound
    # (148 = exp(5)) are pretty much just what I inferred from
    # fitting this model with vaguer priors to the audTemp data.
    #monitor = qw(gamma, tau),
    choice.p = function(ts, ln_gamma, ln_tau)
       {gamma = exp(ln_gamma)
        tau = exp(ln_tau)
        ilogit(
            (log(1 + gamma * ts$llr) - log(1 + gamma * ts$ssr))/gamma - 
            (log(1 + tau * ts$lld) - log(1 + tau * ts$ssd))/tau)},
    init = function(n) list(
        ln_gamma = runif(1, -10, 5),
        ln_tau = runif(1, -10, 5)))

model.sr.rho = stan.choicemodel(
    choice_ll_p_logit =
           '(log(1 + gamma * llr[t]) - log(1 + gamma * ssr[t]))/gamma -
            (log(1 + tau * lld[t]) - log(1 + tau * ssd[t]))/tau',
    parameters =
       'real <lower = -1, upper = 1> f;
        real <lower = 1e-12, upper = 1> rho;',
        # Implicit uniform priors.
    transformed_parameters =
       'real gamma; real tau;
        gamma <- 1/(100 * rho);
        tau <- exp(10 * f) * gamma;',
    choice.p = function(ts, f, rho)
       {gamma = 1 / (100 * rho)
        tau = exp(10 * f) * gamma
        ilogit(
            (log(1 + gamma * ts$llr) - log(1 + gamma * ts$ssr))/gamma -
            (log(1 + tau * ts$lld) - log(1 + tau * ts$ssd))/tau)},
    init = function(n) list(
        f = runif(1, -1, 1),
        rho = runif(1, 1e-12, 1)))

model.fullglm = stan.choicemodel(
    choice_ll_p_logit =
       'b0 +
            b_ssr * ssr[t] +
            b_ssd * ssd[t] +
            b_llr * llr[t] +
            b_lld * lld[t]',
    choice.p = function(ts, b0, b_ssr, b_ssd, b_llr, b_lld)
        ilogit(b0 +
            b_ssr * ts$ssr +
            b_ssd * ts$ssd +
            b_llr * ts$llr +
            b_lld * ts$lld),
    parameters =
       'real <lower = -5, upper = 5> b0;
        real <lower = -5, upper = 0> b_ssr;
        real <lower = 0, upper = 5> b_ssd;
        real <lower = 0, upper = 5> b_llr;
        real <lower = -5, upper = 0> b_lld;',
        # Implicit uniform priors.
    init = function(n) list(
        b0 = min(5, max(-5, rnorm(1, 0, .5))),
        b_ssr = -min(5, abs(rnorm(1, 0, .25))),
        b_ssd = min(5, abs(rnorm(1, 0, .25))),
        b_llr = min(5, abs(rnorm(1, 0, .25))),
        b_lld = -min(5, abs(rnorm(1, 0, .25)))))

model.expk.rho = stan.choicemodel(
# model.exprho reparametrized to be more like model.expk.rho.
    choice_ll_p_logit =
       'rho_v *
          (llr[t] * exp(a * lld[t]) -
           ssr[t] * exp(a * ssd[t]))',
    parameters =
       'real <lower = 0, upper = 1> v30;
        real <lower = 0, upper = 1> rho;',
        # Implicit uniform priors.
    transformed_parameters =
       'real a; real rho_v;
        a <- log(v30)/30;
        rho_v <- 10 * rho;',
    choice.p = function(ts, v30, rho)
       {a = log(v30)/30
        ssv = ts$ssr * exp(a * ts$ssd)
        llv = ts$llr * exp(a * ts$lld)
        rho.v = 10 * rho
        ilogit(rho.v * (llv - ssv))},
    init = function(n) list(
        v30 = runif(1, 0, 1),
        rho = runif(1, 0, 1)))

ghmrho.f = function(disc, curve, rho, ssr, ssd, llr, lld)
   {ssv = ssr * (if (curve == 0) exp(-disc * ssd) else (1 + curve * disc * ssd)^(-1/curve))
    llv = llr * (if (curve == 0) exp(-disc * lld) else (1 + curve * disc * lld)^(-1/curve))
    ilogit(rho * (llv  - ssv))}

model.ghmk.rho = stan.choicemodel(
    choice_ll_p_logit =
       'rho_v *
          (llr[t] * pow(1 + a * lld[t], e) -
           ssr[t] * pow(1 + a * ssd[t], e))',
    parameters =
       'real <lower = .0005, upper = 1 - 1e-12> v30;
        real <lower = 1e-12, upper = .9> scurve;
        real <lower = 0, upper = 1> rho;',
        # Implicit uniform priors.
        #
        # The ranges of v30 and scurve are chosen to be broad
        # without permitting pow(v30, curve), which is used in the
        # computation of 'a' below, to underflow to 0. (v30 =
        # .0005 and scurve = .91, for example, would underflow.)
    transformed_parameters =
       'real curve; real e; real a; real rho_v;
        curve <- 10/(1 - scurve) - 10;
        e <- 1/-curve;
        a <- (1/pow(v30, curve) - 1)/30;
        rho_v <- 10 * rho;',
    choice.p = function(ts, v30, scurve, rho)
       {curve = 10*(1/(1 - scurve) - 1)
        e = 1/-curve
        a = (1/v30^curve - 1)/30
        ssv = ts$ssr * (1 + a * ts$ssd)^e
        llv = ts$llr * (1 + a * ts$lld)^e
        rho.v = 10 * rho
        ilogit(rho.v * (llv - ssv))},
    init = function(n) list(
        v30 = runif(1, .0005, 1 - 1e-12),
        scurve = runif(1, 1e-12, .9),
        rho = runif(1, 0, 1)),
    n.adapt = 250)
