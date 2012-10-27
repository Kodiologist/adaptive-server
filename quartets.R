monthlen = 30.4375
  # This approximation is chosen to have a short decimal
  # representation and to be exactly representable in binary
  # floating point. (The fractional part is 7/16.)

adaptive.delays = c(1, 2, 3, 4, 5, 7, 10, 2*7, 3*7,
    monthlen, 6*7, 2*monthlen, 3*monthlen, 4*monthlen)

adaptive.quartets = expand.grid(
    ssr = seq(5, 100, by = 5),
    ssd = c(0, adaptive.delays),
    llr = seq(5, 80, by = 5),
    lld = adaptive.delays)
adaptive.quartets = transform(
    ss(adaptive.quartets, lld > ssd),
    llr = ssr + llr)
#set.seed(200); test.quartets = samprows(adaptive.quartets, 500, F); row.names(test.quartets) = c()
