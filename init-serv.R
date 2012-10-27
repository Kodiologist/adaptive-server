# -*- R -*-

library(RJSONIO)
jsonize = function (x)
  # Required by my hacky Rserve-perl.
    gsub("\n", " ", fixed = T,
        toJSON(x, digits = 20))

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
       {sql('insert into Trials (subject, trial, lockfile) values (?, ?, ?)',
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
       ts = sql('select ssr, ssd, llr, lld, choose_ll from Trials where subject = ? and choose_ll not null',
           subject)
       ts$choice = logi2factor(ts$choose_ll, c("ss", "ll"))
       ts$choose_ll = c()
       state = sql('select astate from AdaptiveStates where subject = ?',
           subject)
       state = if (nrow(state)) fromJSON(state[[1]]) else list()
       if (modelname == "expk.rho")
          {f = adapt.simultaneous
           model = model.expk.rho}
       else if (modelname == "ghmk.rho")
          {f = adapt.simultaneous
           model = model.ghmk.rho}
       else if (modelname == "diff")
          {f = adapt.simultaneous
           model = model.diff}
       else if (modelname == "sr")
          {f = adapt.simultaneous
           model = model.sr}
       else
          {stop(sprintf("unknown model: %s", modelname))}
       msg("adapting")
       result = do.call(f, c(list(ts, model), state))
       msg("done adapting")
       # Save the result to the database and remove our lockfile.
       sql('insert or replace into AdaptiveStates (subject, astate) values (?, ?)',
            subject, jsonize(result$state))
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
           'select ssr, ssd, llr, lld, final_trial from Trials where subject = ? and trial = ?',
           subject, trial)
       quartet = result[qw(ssr, ssd, llr, lld)]
       final_trial = result$final_trial}
    list(
        quartet = quartet,
        trial = trial,
        final_trial = final_trial)}
