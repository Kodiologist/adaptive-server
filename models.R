library(rstan)

# ------------------------------------------------------------
# Functions for defining models
# ------------------------------------------------------------

gelman.diag.threshold = 1.1

stan.samples = 50
stan.chains = 3

mkdat.stan = function(ts.train, ts.test = NULL)
    list(
        N_ts = nrow(ts.train),
        choose_ll = as.integer(ts.train$choice == "ll"),
        ssr = ts.train$ssr,
        ssd = ts.train$ssd,
        llr = ts.train$llr,
        lld = ts.train$lld,
        N_sim_ts = if (is.null(ts.test)) 1 else nrow(ts.test),
        sim_ssr = if (is.null(ts.test)) 0 else ts.test$ssr,
        sim_ssd = if (is.null(ts.test)) 0 else ts.test$ssd,
        sim_llr = if (is.null(ts.test)) 0 else ts.test$llr,
        sim_lld = if (is.null(ts.test)) 0 else ts.test$lld)

stan.choicemodel = function(
        choice_ll_p = NULL, choice_ll_p_logit = NULL,
        choice.p = NULL,
        parameters, transformed_parameters = '', prior = '',
        params.to.monitor = NULL, init, n.adapt = 150, thin = 1)
   {std.init = init
    if (is.null(params.to.monitor))
        params.to.monitor = names(if (mode(init) == "function") init(1) else init[[1]])
    model.str = sprintf(
       'data
           {int <lower = 0> N_ts;
            vector<lower = 0>[N_ts] ssr;
            vector<lower = 0>[N_ts] ssd;
            vector<lower = 0>[N_ts] llr;
            vector<lower = 0>[N_ts] lld;
            int <lower = 0, upper = 1> choose_ll[N_ts];
            int <lower = 0> N_sim_ts;
            vector<lower = 0>[N_ts] sim_ssr;
            vector<lower = 0>[N_ts] sim_ssd;
            vector<lower = 0>[N_ts] sim_llr;
            vector<lower = 0>[N_ts] sim_lld;}
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
                     choose_ll_p_logit[t]);*/}}
        generated quantities
           {vector[N_sim_ts] sim_choose_ll_p;
            for (t in 1 : N_sim_ts)
               sim_choose_ll_p[t] <- %s;}',
        parameters,
        transformed_parameters,
        prior,
        (if (is.null(choice_ll_p_logit))
            sprintf("bernoulli(%s)", choice_ll_p) else
            sprintf("bernoulli_logit(%s)", choice_ll_p_logit)),
        gsub("\\bssr\\b", "sim_ssr",
            gsub("\\bssd\\b", "sim_ssd",
            gsub("\\bllr\\b", "sim_llr",
            gsub("\\blld\\b", "sim_lld",
            (if (is.null(choice_ll_p_logit))
                choice_ll_p else
                sprintf("inv_logit(%s)", choice_ll_p_logit)))))))
    gendata.model.str = sprintf(
       "data
           {int <lower = 0> N_ts;
            vector[N_ts] ssr;
            vector[N_ts] ssd;
            vector[N_ts] llr;
            vector[N_ts] lld;
            %s}
        parameters {real ignored;} model {0 ~ normal(ignored, 1);}
        generated quantities
           {vector[N_ts] choose_ll_p;
            %s
            for (t in 1 : N_ts)
               choose_ll_p[t] <- inv_logit(%s);}",
        parameters,
        transformed_parameters,
        choice_ll_p_logit)
    gendata = function(ts, ...)
    # Adds true.p and a simulated choice column to the data frame.
       {ts$choice = NA
        if (is.null(choice.p))
           {model = cached_stan_model(gendata.model.str)
            capture.output(fit <- sampling(model, refresh = -1,
                data = c(
                    mkdat.stan(ts)[qw(N_ts, ssr, ssd, llr, lld)],
                    ...),
                chains = 1, iter = 2, warmup = 0, thin = 1))
            p = c(rstan::extract(fit, "choose_ll_p", perm = T)[[1]])}
        else
            p = choice.p(ts, as.numeric(c(...)))
        transform(ts,
            true.p = p,
            choice = logi2factor(rbinom(length(p), 1, p), qw(ss, ll)))};
    sample.posterior = function(ts, raw.posterior = F, init = std.init, debugging = F)
        predict.choices(ts, NULL, raw.posterior, init, debugging)
    predict.choices = function(ts.train, ts.test = NULL, raw.posterior = F, init = std.init, debugging = F)
       {if (!is.null(ts.test)) ts.test$choice = NA
        monitor = params.to.monitor
        if (!is.null(ts.test))
            monitor = c(monitor, "sim_choose_ll_p")
        current.thin = thin
        current.adapt = n.adapt
        model = cached_stan_model(model.str)
        round = 1
        subround = 1
        repeat
           {capture.output(fit <- sampling(model, refresh = -1,
                data = mkdat.stan(ts.train, ts.test),
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
               {ex = rstan::extract(fit, params.to.monitor, perm = F)
                for (n in 1 : dim(ex)[[3]])
                    {mat = as.matrix(mapcols(drop(ex[,,n]), function (v)
                        round(quantile(v, c(.025, .25, .5, .75, .975)), 3)))
                     names(dimnames(mat)) = c("", params.to.monitor[[n]])
                     print(mat)}}
            rhats = summary(fit, pars = params.to.monitor, use_cache = F)$summary[,"Rhat"]
            if (all(rhats < gelman.diag.threshold))
                break
            current.thin = 2 * current.thin
            current.adapt = 2 * current.adapt
            round = round + 1
            subround = 1
            message(sprintf("r%d: Rhats %s; thin %d, adapt %d",
                as.integer(round),
                paste(collapse = ", ", round(rhats, 2)),
                as.integer(current.thin),
                as.integer(current.adapt)))}
        posterior = rstan::extract(fit, monitor, perm = T)
        if (raw.posterior) return(posterior)
        param.post = posterior
        param.post$sim_choose_ll_p = c()
        sim_choose_ll_p = posterior$sim_choose_ll_p
        d = do.call(rbind, lapply(param.post, function (v)
           {q = quantile(v, c(.025, .975))
            data.frame(
                param = "", trial = NA,
                lo = q[1],
                mean = mean(v),
                hi = q[2],
                irng = q[2] - q[1],
                ppos = mean(v > 0))}))
        d$param = factor(params.to.monitor, levels = params.to.monitor)
        row.names(d) = d$param
        if (!is.null(ts.test))
            d = rbind(d, t(sapply(1 : ncol(sim_choose_ll_p), function (t)
               {q = as.vector(quantile(sim_choose_ll_p[,t], c(.025, .975)))
                c(
                    param = NA, trial = t,
                    lo = q[1],
                    mean = mean(sim_choose_ll_p[,t]),
                    hi = q[2],
                    irng = q[2] - q[1],
                    ppos = 1)})))
        d}
    precompile = function()
        cached_stan_model(model.str)
    list(
        params.to.monitor = params.to.monitor,
        choice.p = choice.p,
        model.str = model.str,
        init = std.init,
        gendata = gendata,
        sample.posterior = sample.posterior,
        predict.choices = predict.choices,
        precompile = precompile)}

rlunif = function(n, min = 1e-10, max = 1)
   exp(runif(n, log(min), log(max)))

# ------------------------------------------------------------
# Some model definitions
# ------------------------------------------------------------

model.rewards = stan.choicemodel(
    choice_ll_p_logit =
           'b_ssr * ssr[t] +
            b_llr * llr[t]',
    choice.p = function(ts, theta)
        ilogit(
            theta[1] * ts$ssr +
            theta[2] * ts$llr),
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
    choice.p = function(ts, theta)
        ilogit(
            theta[1] * ts$ssr +
            theta[2] * ts$llr),
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
    choice.p = function(ts, theta)
        ilogit(
            theta[1] * (ts$llr - ts$ssr) -
            theta[2] * (ts$lld - ts$ssd)),
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
    choice.p = function(ts, theta)
        theta[3] * .5 + (1 - theta[3]) * (
            theta[1] * (ts$llr - ts$ssr) >
            theta[2] * (ts$lld - ts$ssd)),
    parameters =
       'real <lower = 0, upper = 5> b_rd;
        real <lower = 0, upper = 5> b_dd;
        real <lower = 0, upper = .5> delta;',
        # Implicit uniform priors.
    init = function(n) list(
        b_rd = runif(1, 0, 5),
        b_dd = runif(1, 0, 5),
        delta = runif(1, 0, .5)))

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
    #params.to.monitor = qw(gamma, tau),
    choice.p = function(ts, theta)
       {gamma = exp(theta[1])
        tau = exp(theta[2])
        ilogit(
            (log(1 + gamma * ts$llr) - log(1 + gamma * ts$ssr))/gamma - 
            (log(1 + tau * ts$lld) - log(1 + tau * ts$ssd))/tau)},
    init = function(n) list(
        ln_gamma = runif(1, -10, 5),
        ln_tau = runif(1, -10, 5)))

model.fullglm = stan.choicemodel(
    choice_ll_p_logit =
       'b0 +
            b_ssr * ssr[t] +
            b_ssd * ssd[t] +
            b_llr * llr[t] +
            b_lld * lld[t]',
    choice.p = function(ts, theta)
        ilogit(theta[1] +
            theta[2] * ts$ssr +
            theta[3] * ts$ssd +
            theta[4] * ts$llr +
            theta[5] * ts$lld),
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
    choice.p = function(ts, theta)
       {vassign(.(v30, rho), theta)
        a = log(v30)/30
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

model.ghmv30.rs.rho1 = stan.choicemodel(
    choice_ll_p_logit =
       'llr[t] * pow(1 + a * lld[t], e) -
        ssr[t] * pow(1 + a * ssd[t], e)',
    parameters =
       'real <lower = 0, upper = 1> v30;
        real <lower = 0, upper = 1> scurve;',
        # Implicit uniform priors.
    transformed_parameters =
       'real curve; real e; real a;
        curve <- 10/(1 - scurve) - 10;
        e <- 1/-curve;
        a <- (1/pow(v30, curve) - 1)/30;',
    init = function(n) list(
        v30 = runif(1, 0, 1),
        scurve = runif(1, 0, 1)),
    n.adapt = 250)

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
    choice.p = function(ts, theta)
       {vassign(.(v30, scurve, rho), theta)
        curve = 10*(1/(1 - scurve) - 1)
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
