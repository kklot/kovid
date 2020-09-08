#' @useDynLib kovid ksimc
ksim = function(
  beta          = .1,
  POP           = vietnampop,
  ip            = gen_C(),
  durInf        = 14,
  durLat        = 07,
  durQua        = 21,
  durTrt        = 7,
  p_death       = rep(0, 16),
  ess           =  c(.1, .1),
  rho           = 1.5^(1:16),
  nInfected     = 10,
  verbose       = F,
  version       = "K"
) {
  dnames <- c('S', 'E', 'Is', 'Ic', 'R', 'D', 'Q', 'T', 'lambda')
  if (version=="K") {
    kinputs <-list(
      beta      = beta,
      ess       = ess,
      p_death   = p_death,
      rho       = rho,
      durInf    = durInf,
      durLat    = durLat,
      durQua    = durQua,
      durTrt    = durTrt,
      POP       = as.matrix(POP),
      epi_time  = ip$epi_time,
      M         = ip$M,
      impulse_v = ip$impulse_v,
      impulse_a = ip$impulse_a,
      nInfected = nInfected
    )
    out <- .Call(ksimc, kinputs)
    dimnames(out)[[1]] <- dnames
    return(list(model = out, args  = as.list(match.call())))
  }

  gamma = 1-exp(-1/durInf)
  alpha = 1-exp(-1/durLat)
  Q_tim = 1-exp(-1/21)

  N     = sum(POP$popage)
  N_age = POP$popage
  p_age = POP$propage

  # S = E = Is = Ic = R = lambda
  S = 1; E = 2; Is = 3; Ic = 4; R = 5; D=6; Q=7; lambda = 8
  o <- array(0, c(8, length(p_age), ip$epi_time))
  dimnames(o)[[1]] <- dnames

  o[Ic,,1] = nInfected
  o[S ,,1] = N_age - sum(o[-c(S, lambda),,1])
  for (i in 1:(ip$epi_time-1)) {
    if (impulse_v[i])
      o[E,,i] = impulse_a[,i]
    o[lambda,,i] = beta * ip$M[,,i] %*% ((o[Ic,,i] + o[Is,,i]) / N_age)
    SE   = o[lambda,,i]    * o[S,,i]
    EIc  = alpha *    rho  * o[E,,i]
    EIs  = alpha * (1-rho) * o[E,,i]
    IcQ  = gamma * o[Ic,,i] * (1 - p_death)
    IcD  = gamma * o[Ic,,i] *      p_death
    Q2R  = Q_tim * o[Q,,i]
    IsR  = gamma * o[Is,,i]
    
    o[S ,,i+1] = o[S,  ,i] - SE;
    o[E ,,i+1] = o[E,  ,i] + SE - EIc - EIs;
    o[Ic,,i+1] = o[Ic, ,i]      + EIc - IcQ - IcD;
    o[Is,,i+1] = o[Is, ,i]            + EIs - IsR;
    o[Q ,,i+1] = o[Q,  ,i]            + IcQ - Q2R;
    o[R ,,i+1] = o[R,  ,i]                  + Q2R;
    o[D ,,i+1] = o[D,  ,i]                  + IcD;
    if (verbose) {cat('\rDay', i); flush.console()}
  }
  list(model = o, args  = as.list(match.call()))
}