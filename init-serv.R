# -*- R -*-
if (!exists("server.version")) server.version = "unset"

library(RJSONIO)
jsonize = function (x)
  # Required by my hacky Rserve-perl.
    gsub("\n", " ", fixed = T,
        toJSON(x, digits = 20))

options(fork_maybeparallel.sapply = F)
  # Don't do grid approximation in parallel.
source("util.R")
source("quartets.R")
source("models.R")
source("adaptive.R")
library(RSQLite)

# The following would be automatically loaded when necessary,
# but we want to preempt them.
library(R.methodsS3)
library(R.oo)
library(R.utils)
library(digest)
library(R.cache)

if (!exists("db.path")) db.path = "adaptive-server-jobs.sqlite"

db = dbConnect(dbDriver("SQLite"), dbname = db.path)
sql = function (query, ...)
    dbGetPreparedQuery(db, query, bind.data = data.frame(...))

models = list(
    expk.rho = model.expk.rho,
    ghmk.rho = model.ghmk.rho,
    diff.rho = model.diff.rho,
    sr.rho = model.sr.rho)
for (m in models)
    m$precompile()

msg = function (...)
    message("GNQ: ", ...)

# --------------------------------------------------

get_next_quartet = function(modelname, subject, trial, prev_choose_ll = NULL)
   {msg("in get_next_quartet")
    # Create a lockfile.
    my.lockfile = tempfile("adaptive-lock-")
    cat(file = my.lockfile)
    # Try to insert a new row for this job in the database.
    msg("trying to insert (", subject, ", ", trial, ", ", my.lockfile, ")")
    inserted = tryCatch(
        {
            sql('insert into Trials (subject, trial, lockfile) values (?, ?, ?)',
                subject, trial, my.lockfile)
            T},
        error = function (e) F)
    if (inserted) 
       {msg("working")
        # Update the database.
        if (!is.null(prev_choose_ll))
            sql('update Trials set choose_ll = ? where subject = ? and trial = ?',
                prev_choose_ll, subject, trial - 1)
        # Do the job.
        ts = sql('select cast(ssr as real) as ssr, cast(ssd as real) as ssd, cast(llr as real) as llr, cast(lld as real) as lld, choose_ll from Trials where subject = ? and choose_ll not null',
            subject)
          # The casts are necessary because otherwise, if a column
          # begins with an integer, RSQLite will truncate all
          # floats in the column.
        ts$choice = logi2factor(ts$choose_ll, c("ss", "ll"))
        ts$choose_ll = c()
        if (modelname %in% names(models))
            model = models[[modelname]]
        else
            stop(sprintf("unknown model: %s", modelname))
        msg("adapting")
        newseed()
        result = adapt.simultaneous(ts, model)
        msg("done adapting")
        # Save the result to the database and remove our lockfile.
        if (length(result$diagnostics))
            with(result$diagnostics, sql(
                'insert into MCMCDiagnostics values (?, ?, ?, ?, ?)',
                subject, trial, mcmc.round, quit.early,
                paste(sprintf("%.3f", rhats), collapse = ",")))
        quartet = result$quartet
        final_trial = result$final_trial
        with(quartet,
            sql('update Trials set ssr = ?, ssd = ?, llr = ?, lld = ?, final_trial = ?
                    where subject = ? and trial = ?',
                ssr, ssd, llr, lld, final_trial, subject, trial))
        unlink(my.lockfile)}
    else
      # We failed, assumedly because (subject, trial) is already
      # in the database; that is, another process has already
      # started this job.
       {msg("waiting")
        # Delete our own lockfile (since we're not doing this job,
        # it isn't needed).
        unlink(my.lockfile)
        # Wait for the worker's lockfile to be deleted.
        worker.lockfile = sql(
            'select lockfile from Trials where subject = ? and trial = ?',
            subject, trial)[[1]]
        system2('inotifywait', c('-e', 'move_self', worker.lockfile),
            stdout = F, stderr = F)
        msg("done waiting")
        # The worker must be done, so get the result.
        result = sql(
            'select cast(ssr as real) as ssr, cast(ssd as real) as ssd, cast(llr as real) as llr, cast(lld as real) as lld, final_trial
                from Trials where subject = ? and trial = ?',
            subject, trial)
        quartet = result[qw(ssr, ssd, llr, lld)]
        final_trial = result$final_trial}
    punl(quartet, trial, final_trial,
        server_version = server.version)}
