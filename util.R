ss = subset
char = as.character

qw = function(...)
# qw(a, b, c) → c("a", "b", "c")
   {m = match.call(expand.dots = FALSE)
    char(m[["..."]])}

punl = function(...)
# Like "list", but uses punning to generate default names, so
#   punl(a, b, c = x, d)
# is equivalent to
#   list(a = a, b = b, c = x, d = d)
   {m = match.call(expand.dots = F)
    inp = as.character(m$"...")
    l = list(...)
    names(l) =
        (if (is.null(names(l)))
            inp
         else
            ifelse(names(l) == "", inp, names(l)))
    l}

vassign = function(vars, values, envir = parent.frame())
# vassign(.(a, b), c(1, 2)) is equivalent to {a = 1; b = 2;}.
   {m = match.call()
    vassign.char(char(m$vars)[-1], values, envir)}

vassign.char = function(vars, values, envir = parent.frame())
    for (i in seq_along(vars))
        assign(vars[[i]], values[[i]], envir)

newseed = function ()
# Sets the RNG seed with the current time including some fractional
# seconds.
    set.seed((as.numeric(Sys.time()) * 1e5) %% .Machine$integer.max)

maprows = function(x, f) apply(x, 1, f)
mapcols = function(x, f) apply(x, 2, f)

col.list = function(m)
# Turns a matrix into a list of columns, like as.data.frame.
    lapply(1 : ncol(m), function (i) m[,i])

samprows = function(m, n.rows, replace = F)
    m[sample.int(nrow(m), n.rows, replace),]

pick = function(v, size = 1, ...)
# Like sample, but without surprising behavior when length(v) == 1
# and with size defaulting to 1.
    v[sample.int(length(v), size, ...)]

ilogit = function (v)  1 / (1 + exp(-v))

mod1 = function (v, n)
# mod1(1:20, 6) =>
# 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2
    ((v - 1) %% n) + 1

logi2factor = function(v, labels)
   `levels<-`(factor(as.logical(v), levels = c(F, T)), labels)

scale01 = function (v) {r = range(v); (v - r[1])/(r[2] - r[1])}

dist.idx = function(d, i)
# Given a dist object and indices into its vector, returns
# the matrix coordinates. Adapted from: http://stackoverflow.com/a/12643509
   {n = attr(d, "Size")
    col = ceiling(n - (1 + sqrt(1 + 4 * (n^2 - n - 2*i)))/2)
    cbind(row = n - (2*n - col + 1) * col / 2 + i + col, col)}

which.farthest = function(m, ...)
# Returns the indices of the rows of m with the greatest distance.
   {d = dist(m, ...)
    dist.idx(d, which.max(d))}

standardize.stan.code = function(stan.code)
# Remove comments and reduce whitespace.
    gsub("\\s+", " ",
    sub("\\s+$", "",
    sub("^\\s+", "",
    gsub("/\\*.+?\\*/", "",
    gsub("//.*?(\n|$)", "\\1",
    stan.code)))))

cached_stan_model = function(code, bypass.cache = F)
   {library(R.cache)
    library(rstan)
    dirs = c("Kodi", "stan", "models")
    key = list(code = standardize.stan.code(code))
    if (!bypass.cache)
       {model = loadCache(key, dirs = dirs)
        if (!is.null(model))
            return(model)}
    model = stan_model(model_code = code)
    saveCache(model, key, dirs = dirs)
    model}
