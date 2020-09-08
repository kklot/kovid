minmax <- function(x, minx, maxx) {
  if (missing(minx))  minx <- min(x)
  if (missing(maxx))  maxx <- max(x)
  (x  - minx)  / (maxx-minx)
}

today <- function(as_date=TRUE) {
  td = format(Sys.time(), '%Y-%m-%d')
  if(as_date)
    return(as.Date(td))
  td
}

intervention <- function(mixmat=contacts,
  work=0, school=0, other=0, lockdown=0,
  work_essential=0.1, other_essential=0.1) {
  work_m <- school_m <- others_m <- diag(16)
  if (lockdown)
    work <- school <- other <- 1
  if (work)
    work_m <- work_essential * mixmat[['work']]
  if (school)
    school_m <- 0 * mixmat[['school']]
  if (other)
    other <- other_essential * mixmat[['other']]
  mixmat[['home']] + work_m + school_m + others_m
}

prog_id <- function(x, x0, xend)
  x %>% subtract(x0) %>% paste(collapse=":") %>% 
  parse(text=.) %>% eval %>% is_in(1:xend, .)

gen_impulse <- function(epitime, days=c('Sat', 'Sun'), agr=5:11) {
  # 5:11 20-54
  impulse_v <- epitime %>% format('%a') %>% is_in(days)
  impulse_a <- matrix(0, nrow=16, ncol=length(impulse_v))
  for (i in seq_along(impulse_v)) {
    if (impulse_v[i])
      impulse_a[sample(agr, 1), i] <- 1
  }
  list(impulse_v=impulse_v, impulse_a=impulse_a)
}

gen_C <- function(
  date_start    = '2020-01-01', date_end,
  date_work     = c('2020-03-01', '2020-04-01'),
  date_school   = c('2020-03-01', '2020-04-01'),
  date_lockdown = c('2020-04-01', '2020-12-01'),
  mixmat        = contacts,
  essential     = c(.1, .5, .1),
  impulse_days  = NULL,
  agr           = 5:11,
  time_to_full  = 7
) {

  date_start    %<>% as.Date
  if (missing(date_end)) 
    date_end <- format(Sys.time(), '%Y-%m-%d')
  date_end      %<>% as.Date
  date_work     %<>% as.Date
  date_school   %<>% as.Date
  date_lockdown %<>% as.Date
  
  epi_time <- date_end - date_start
  time <- date_start:(date_end-1) %>% as.Date(origin="1970-01-01")
  if (is.null(impulse_days)) 
    impulse_days <- 'K'
  impulse <- gen_impulse(time, impulse_days, agr=agr)

  work_time     = prog_id(date_work, date_start, epi_time)
  school_time   = prog_id(date_school, date_start, epi_time)
  lockdown_time = prog_id(date_lockdown, date_start, epi_time)

  # school has immediate effect
  ess_work <- as.numeric(!work_time)
  ess_work[ess_work==0] <- essential[1]
  ess_lock <- as.numeric(!lockdown_time)
  ess_lock[ess_lock==0] <- essential[2]
  time_need <- time_to_full - 1
  if (time_need > 0) {
    ess_work[which(work_time)[1:time_to_full]] <- 
      exp((log(essential[1])/time_need)*0:time_need)
    ess_lock[which(lockdown_time)[1:time_to_full]] <- 
      exp((log(essential[2])/time_need)*0:time_need)    
  }
  dims <- dim(mixmat[[1]])
  M <- mapply(intervention,
    work           = work_time, 
    school         = school_time, 
    lockdown       = lockdown_time, 
    work_essential = ess_work, 
    other_essential= ess_lock,
    SIMPLIFY=FALSE) %>% unlist %>% 
    array(c(dims[1], dims[2], epi_time))
  list(M=M, epi_time = epi_time, time= time, 
    impulse_v = impulse$impulse_v,
    impulse_a = impulse$impulse_a,
    inputs=as.list(match.call()))
}