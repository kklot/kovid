lprior <- function(x) {
  x <- as.vector(x)
  names(x) <- pname
  o <- dbeta(x['beta'], 1, 4, log=T)
  # hist(rbeta(3000, 1, 4))
  
  o <- o + dbeta(x['work'], 1, 10, log=T)
  # hist(rbeta(3000, 1, 10))
  # o <- o + dbeta(x['other'], 1, 10, log=T)
  o <- o + dgamma(x['backdate'], 7, scale=2, log=T)
  # hist(rgamma(3000, 7, scale=2))
  
  o <- o + dlnorm(x['latent'], log(7), 0.1, log=T)
  o <- o + dlnorm(x['infect'], log(14), 0.1, log=T)
  # hist(rlnorm(3000, log(14), .1))

  o <- o + dunif(x['death'], 1, d_scale, log=T)
  o <- o + dunif(x['asymp'], 1, d_scale, log=T)

  o <- o + dlnorm(x['delay'], log(1), 0.5, log=T)
  # hist(rlnorm(3000, log(1), .5))
  o <- o + dlnorm(x['quara'], log(14), 0.1, log=T)
  # hist(rlnorm(3000, log(14), .1))
  o <- o + dlnorm(x['treat'], log(3), 0.2, log=T)
  # hist(rlnorm(3000, log(3), .2))
  # o <- o + dlnorm(x['i0'], log(1), 1, log=T)
  # hist(rlnorm(3000, log(1), 1))
  o
}

sample_prior <- function(n=100) {
  o <- rbeta(n, 1, 4)

  o <- cbind(o, rbeta(n, 1, 10))
  # o <- cbind(o, rbeta(n, 1, 10))
  o <- cbind(o, rgamma(n, 7, scale=2))
  
  o <- cbind(o, rlnorm(n, log(7), 0.1))
  o <- cbind(o, rlnorm(n, log(14), 0.1))

  o <- cbind(o, runif(n, 1, d_scale))
  o <- cbind(o, runif(n, 1, d_scale))

  o <- cbind(o, rlnorm(n, log(1), 0.5))
  o <- cbind(o, rlnorm(n, log(14), 0.1))
  o <- cbind(o, rlnorm(n, log(3), 0.2))
  # o <- cbind(o, rlnorm(n, log(1), 1))
  colnames(o) <- pname
  o
}

find_best <- function(smp) {
  o = apply(smp, 1, function(x) {
    names(x) <- pname
    l_post(x)
  })
  smp[which.max(o), ]
}

ll <- function(x) {
  .x <<- x
  x <- as.vector(x)
  names(x) <- pname
  if (x['delay']>70) return(-Inf)
  ip <- actions %>% 
      inset2("date_start", cases$firstT - round(x['backdate'])) %>% 
      inset2("time_to_full", round(x['delay']) ) %>% 
      inset2("essential", x[c('work', 'work')]) %>% 
      do.call('gen_C', .)
  epi = ksim(
    beta    = x['beta'],
    durLat  = x['latent'],
    durInf  = x['infect'],
    durQua  = x['quara'],
    durTrt  = x['treat'],
    p_death = p_death * x['death'],
    POP     = pop,
    rho     = p_death * x['asymp'],  
    ip      = ip
  )
  o = data.frame(time=ip$time, I=colSums(epi$model['Ic',,])) %>% 
    inner_join(cases$I, by='time') %>% 
    summarise(mean(abs(cumsum(I.x) - cumsum(I.y)))) %>% as.numeric
  o = o + data.frame(time=ip$time, R=colSums(epi$model['R',,])) %>% 
    inner_join(cases$R, by='time') %>% 
    summarise(mean(abs(R.x - R.y)))%>% as.numeric
  o = o + data.frame(time=ip$time, D=colSums(epi$model['D',,])) %>% 
    inner_join(cases$D, by='time') %>% 
    summarise(mean(abs(D.x - D.y)))%>% as.numeric
  if (is.na(o)) return(-Inf)
  -o
}

l_post <- function(x) {
  a <- lprior(x)
  if (!is.finite(a)) return(-Inf)
  a+ll(x)
}
