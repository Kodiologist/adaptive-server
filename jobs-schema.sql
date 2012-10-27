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

create table AdaptiveStates
   (subject                         integer primary key
        references Trials(subject),
    astate          text            not null);
