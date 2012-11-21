initify = function (model, v)
    `names<-`(as.list(v), model$monitor)

quartet.repeat.interval = 10
  # The adaptive procedures will only be allowed to present the
  # same quartet every repeat.interval trials.

choose.maxpdiff = function(choice.p, theta1, theta2, ts)
# Choose a quartet maximizing the difference in choice
# probability for theta1 and theta2.
   {recent.quartets = with(tail(ts, quartet.repeat.interval - 1),
        paste(ssr, ssd, llr, lld))
    adaptive.quartets = ss(adaptive.quartets, !(
        paste(ssr, ssd, llr, lld) %in% recent.quartets))
    i = which.max(pdiff <- abs(
        if (length(formals(choice.p)) == 2)
            choice.p(adaptive.quartets, theta1) -
            choice.p(adaptive.quartets, theta2)
        else
            do.call(choice.p, c(list(adaptive.quartets), as.list(theta1))) -
            do.call(choice.p, c(list(adaptive.quartets), as.list(theta2)))))
    adaptive.quartets[i,]}

adapt.simultaneous.trials = 50

adapt.simultaneous.grid = function(ts, model)
   {if (nrow(ts))
       {# Sample the posterior.
        post = model$sample.posterior(ts)
        # Pick the two farthest points in the posterior sample
        # to be the new theta1 and theta2.
        rows = post[which.farthest(mapcols(post, scale01)),]
        theta1 = c(rows[1,])
        theta2 = c(rows[2,])}
    else
       {theta1 = c(model$sample.posterior(empty.ts, n = 1))
        theta2 = c(model$sample.posterior(empty.ts, n = 1))}

    # Return the new quartet and a flag saying whether we're done
    # adapting.
    list(
        quartet = choose.maxpdiff(model$choice.p, theta1, theta2, ts),
        final_trial = as.integer(nrow(ts) + 1 == adapt.simultaneous.trials),
        state = list())}

adapt.simultaneous = function(ts, model)
   {if (nrow(ts))
       {# Sample the posterior.
        post = simplify2array(model$sample.posterior(ts))
        stopifnot(!is.null(post))
        # Pick the two farthest points in the posterior sample
        # to be the new theta1 and theta2.
        rows = post[which.farthest(mapcols(post, scale01)),]
        theta1 = c(rows[1,])
        theta2 = c(rows[2,])}
    else
      # Initialize thetas.
       {theta1 = as.numeric(model$init(1))
        theta2 = as.numeric(model$init(2))}

    # Return the new quartet and a flag saying whether we're done
    # adapting.
    list(
        quartet = choose.maxpdiff(model$choice.p, theta1, theta2, ts),
        final_trial = as.integer(nrow(ts) + 1 == adapt.simultaneous.trials),
        state = list())}

adapt.1patatime.trials = 50

adapt.1patatime = function(ts, model, theta1, theta2, postsamp)
   {if (nrow(ts))
       {# Sample the posterior.
        post = simplify2array(model$sample.posterior(ts,
            init = maprows(postsamp, function (v) initify(model, v))))
              # This means the initial values aren't overdispersed,
              # I know. I'm hoping that isn't too much of a problem
              # since we're fitting the same model lots of times
              # per subject.
        stopifnot(!is.null(post))
        postsamp = samprows(post, 3)
        # For a given parameter (which changes each round),
        # set theta1 and theta2 to differ somewhat on that
        # parameter and be the same for all other parameters.
        param = mod1(nrow(ts) + 1, length(theta1))
        theta1 = theta2 = mapcols(post, median)
        theta1[param] = quantile(post[,param], 1/3)
        theta2[param] = quantile(post[,param], 2/3)}
    else
      # Initialize thetas and postsamp.
       {theta1 = as.numeric(model$init(1))
        theta2 = theta1
        theta2[1] = model$init(2)[[1]]
        postsamp = t(sapply(1:3, model$init))}

    # Return the new quartet, a flag saying whether we're done
    # adapting, and the new thetas and postsamp.
    list(
        quartet = choose.maxpdiff(model$choice.p, theta1, theta2, ts),
        final_trial = as.integer(nrow(ts) + 1 == adapt.1patatime.trials),
        state = list(theta1, theta2, postsamp))}

empty.ts = data.frame(
    ssr = numeric(0), ssd = numeric(0),
    llr = numeric(0), lld = numeric(0),
    choice = factor(character(0), levels = qw(ss, ll)))
