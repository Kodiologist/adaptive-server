create table Trials
   (subject         integer         not null,
    trial           integer         not null,
    lockfile        text            not null,
    ssr             numeric,
    ssd             numeric,
    llr             numeric,
    lld             numeric,
    final_trial     integer,
    choose_ll       integer,
    primary key (subject, trial));

create table MCMCDiagnostics
   (subject         integer         not null,
    trial           integer         not null,
    mcmc_round      integer         not null,
    quit_early      integer         not null,
    rhats           text            not null,
    primary key (subject, trial));
