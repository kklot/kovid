#' @useDynLib kovid ksimc
ksim = function(
  x,
  POP           = vietnampop,
  mixmat        = contacts,
  ip            = gen_C(),
  nInfected     = 1, # per age group
  verbose       = F
) {
  dnames <- c('S', 'E', 'Is', 'Ic', 'R', 'D', 'Q', 'T', 'lambda')
  kinputs <- list(
    beta       = x['beta'],
    durInf     = x['infect'],
    durLat     = x['latent'],
    durQua     = x['quara'],
    durTrt     = x['treat'],
    nInfected  = nInfected,
    epi_time   = ip$epi_time,
    rho        = p_death * x['asymp'],
    p_death    = p_death * x['death'],
    impulse_v  = ip$impulse_v,
    POP        = POP,
    eff        = ip$eff,
    M          = mixmat,
    impulse_a  = ip$impulse_a
  )
  out <- .Call(ksimc, kinputs)
  dimnames(out)[[1]] <- dnames
  return(list(model = out, args  = as.list(match.call())))
}