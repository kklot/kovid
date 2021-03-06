minmax <- function(x, minx, maxx) {
  if (missing(minx))  minx <- min(x)
  if (missing(maxx))  maxx <- max(x)
  (x  - minx)  / (maxx-minx)
}

intervention <- function(M=mixmat,
  work=FALSE, school=FALSE, other=FALSE, lockdown=FALSE,
  work_essential=0.1, other_essential=0.1) {
  n_agr <- nrow(M[['home']])
  if (lockdown)
    work <- school <- other <- TRUE
  if (work)
    M[['work']] <- work_essential * M[['work']]
  if (school)
    M[['school']] <- 0 * M[['school']]
  if (other)
    M[['others']] <- other_essential * M[['others']] # note the s
  Reduce('+', M)
}

prog_id <- function(x, x0, xend)
  x %>% subtract(x0) %>% paste(collapse=":") %>% 
  eval_text %>% is_in(1:xend, .)

# generate import case by day, randomize age-group
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

#' Generate intervention
#' 
#' Generate intervention
#' 
#' @param date_start when epidemic starts
#' @param date_end when **simulation** ends, default to today
#' @param date_work start and end date of workplace closure
#' @param date_school start and end date of school closure
#' @param date_lockdown start and end date of national lockdown
#' @param mixmat contact mixing matrix
#' @param essential vector of length 2, effect of lockdown on workplace and 
#' other activities matrices, the contact is reduced to contact x essential
#' @param smoothing vector of length 2, shape of the transition effect of 
#' lockdown
#' @param impulse_days simulate imported cases on day of the week, eg. "Sun"
#' @param agr age-group of the imported cases
#' @return x
#' @seealso \code{\link{gen_impulse}} \code{\link{intervention}}
#' @export
gen_C <- function(
  date_start    = '2020-01-01', date_end = today(),
  date_work     = c('2020-03-01', '2020-04-01'),
  date_school   = c('2020-03-01', '2020-04-01'),
  date_lockdown = c('2020-04-01', '2020-12-01'),
  essential     = c(work = 0.01, other=0.01), # maximum control level
  smoothing     = c(set  = 0.5, lift=0.5),
  post_lockdown = c(work = 0.9, other = 0.9, school = 0.9),
  impulse_days  = NULL,
  agr           = 5:11,
  group_size    = group_size
) {

  date_start    %<>% as.Date
  date_end      %<>% as.Date %>% add(7)
  date_work     %<>% as.Date
  date_school   %<>% as.Date
  date_lockdown %<>% as.Date
  
  epi_time <- date_end - date_start
  time <- date_start:(date_end-1) %>% as.Date(origin="1970-01-01")
  # import case by date/random age-group
  if (is.null(impulse_days)) 
    impulse_days <- 'K'
  impulse <- gen_impulse(time, impulse_days, agr=agr)

  # Dates when the interventions occur
  work_time     = prog_id(date_work, date_start, epi_time)
  school_time   = prog_id(date_school, date_start, epi_time)
  lockdown_time = prog_id(date_lockdown, date_start, epi_time)

  # generate effect of intervention:
  # - school has immediate effect whereas other and work 
  # - starts the same time although the effect can be different
  work_acting  = range(which(work_time))+7
  eff_work = double_logistic2(1:epi_time, 
    1, essential[1], post_lockdown[1],
    smoothing[1], smoothing[2],
    work_acting[1], work_acting[2])
  other_acting = range(which(lockdown_time))+7
  eff_other = double_logistic2(1:epi_time, 
    1, essential[2], post_lockdown[2],
    smoothing[1], smoothing[2],
    other_acting[1], other_acting[2])
  school_acting = range(which(school_time))+7
  eff_school = double_logistic2(1:epi_time, 
    1, 0, post_lockdown[3],
    10, 10,
    school_acting[1], school_acting[2])
  eff <- rbind(rep(1, epi_time), eff_work, eff_school,  eff_other)
  # eff <- sweep(eff, 1, unlist(group_size[c('household_size', 'work_size', 'school_size', 'city_size')]), '*')
  list(
    epi_time = epi_time, time= time,
    eff      = eff,
    impulse_v = impulse$impulse_v,
    impulse_a = impulse$impulse_a,
    inputs=as.list(match.call())
  )
}