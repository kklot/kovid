kovid <- R6::R6Class('kovid', portable = FALSE,
  public = list(
    cases=NULL, pop=NULL, actions=NULL, contacts=NULL, group_size=NULL,p_death=NULL,d_scale=NULL,
    initialize = function(inputs) {
      cases      <<- inputs$cases
      pop        <<- as.matrix(inputs$pop)
      actions    <<- inputs$actions
      group_size <<- inputs$group_size
      # group_size <<- lapply(inputs$group_size, '/', sum(pop[, 'size']))
      p_death    <<- inputs$p_death
      d_scale    <<- inputs$d_scale
      # scale_ct   <- lapply(inputs$contacts, function(x) sweep(x, 1, rowSums(x), '/'))
      # contacts   <<- Reduce(function(x, y) abind::abind(x, y, along=3), scale_ct)
      contacts   <<- Reduce(function(x, y) abind::abind(x, y, along=3), inputs$contacts)
    },
    pname = c(    # don't change the name for now
      'beta',     # transmission
      'work',     # reduction of contact on workplace closure
      'post',     # post lockdown contact level
      'backdate', # how long ago the first case infected, starting from the first
                  # day we observed infected cases data
      'latent',   # latent period
      'infect',   # infection rate
      'death',    # death rate
      'asymp',    # % of infected cases being asymptotic
      'smooth1',   # smoothing shape transition lockdown <-> not lockdown
      'smooth2',   # smoothing shape transition lockdown <-> not lockdown
      'quara',    # % being quarrantined
      'treat'     # % treated
    ),
    # Sample priors
    samples_prior = NULL,
    sample_prior = function(n=1000) {
      o <- rgamma(n, 1, 1)
      # o <- rbeta(n, 1, 2)

      o <- cbind(o, rbeta(n, 1, 10))
      o <- cbind(o, rbeta(n, 1, 10))

      o <- cbind(o, rgamma(n, 7, scale=1))
      
      o <- cbind(o, rlnorm(n, log(4), 0.3))
      o <- cbind(o, rlnorm(n, log(7), 0.2))

      o <- cbind(o, runif(n, 1, d_scale))
      o <- cbind(o, runif(n, 1, d_scale))

      o <- cbind(o, rlnorm(n, log(1), 0.3))
      o <- cbind(o, rlnorm(n, log(1), 0.3))

      o <- cbind(o, rlnorm(n, log(14), 0.1))
      o <- cbind(o, rlnorm(n, log(3), 0.2))
      # o <- cbind(o, rlnorm(n, log(1), 1))
      colnames(o) <- pname
      samples_prior <<- o
    },
    lprior = function(x) {
      x <- as.vector(x)
      names(x) <- pname
      # o <- dbeta(x['beta'], 1, 2, log=T)
      # hist(rbeta(3000, 1, 2))
      o <- dgamma(x['beta'], 1, 1, log=T)
      # hist(rgamma(3000, 1, scale=1))
      o <- o + dbeta(x['work'], 1, 10, log=T)
      # hist(rbeta(3000, 1, 10))
      o <- o + dbeta(x['post'], 10, 1, log=T)
      # hist(rbeta(3000, 10, 1))
      o <- o + dgamma(x['backdate'], 7, scale=1, log=T)
      # hist(rgamma(3000, 7, scale=1))
      o <- o + dlnorm(x['latent'], log(4), 0.3, log=T)
      # hist(rlnorm(3000, log(4), .3))
      o <- o + dlnorm(x['infect'], log(7), 0.2, log=T)
      # hist(rlnorm(3000, log(7), .2))
      o <- o + dunif(x['death'], 1, d_scale, log=T)
      o <- o + dunif(x['asymp'], 1, d_scale, log=T)

      o <- o + dlnorm(x['smooth1'], log(1), 0.3, log=T)
      o <- o + dlnorm(x['smooth2'], log(1), 0.3, log=T)
      # hist(rlnorm(3000, log(1), .3))
      o <- o + dlnorm(x['quara'], log(14), 0.1, log=T)
      # hist(rlnorm(3000, log(14), .1))
      o <- o + dlnorm(x['treat'], log(3), 0.2, log=T)
      # hist(rlnorm(3000, log(3), .2))
      # o <- o + dlnorm(x['i0'], log(1), 1, log=T)
      # hist(rlnorm(3000, log(1), 1))
      o
    },
    sim_mod = function(x,...) {
      x <- as.vector(x)
      names(x) <- pname
      ip <- actions %>% 
          inset2("date_start", as.Date(cases$firstT[1]) - round(x['backdate'])) %>% 
          inset2("smoothing", x[c('smooth1', 'smooth2')]) %>% 
          inset2("post_lockdown", x[c('post', 'post', 'post')]) %>% 
          inset2("essential", x[c('work', 'work')]) %>% 
          inset2("group_size", group_size) %>% 
          do.call('gen_C', .)
      epi <- ksim(x, POP = pop, mixmat=contacts, ip = ip, ...)
      list(input=ip, epi=epi)
    },
    last_par = NULL,
    ll = function(x,...) {
      last_par <<- x
      sim <- sim_mod(x,...)      
      # o = tibble(time=sim$input$time, I = colSums(sim$epi$model['Ic',,]) + colSums(sim$epi$model['Is',,])) %>% 
      o = tibble(time=sim$input$time, I = colSums(sim$epi$model['Ic',,]) ) %>% 
        mutate(time=as.character(time)) %>% 
        right_join(cases, by='time')
      # o$I.x %>% plotl
      # o$I.y %>% points(lty=2)
      o = o %>% summarise(sqrt(mean(square(I.x - I.y))))
      # o = o + data.frame(time=sim$input$time, R=colSums(sim$epi$model['R',,])) %>% 
      #   mutate(time=as.character(time)) %>% 
      #   inner_join(cases, by='time') %>% 
      #   summarise(sqrt(mean(square((R.x - R.y)/max(R.y)))))
      # o = o + data.frame(time=sim$input$time, D=colSums(sim$epi$model['D',,])) %>% 
      #   mutate(time=as.character(time)) %>% 
      #   inner_join(cases, by='time') %>% 
      #   summarise(sqrt(mean(square((D.x - D.y)/max(D.y)))))
      if (is.na(o)) return(Inf)
      as.numeric(o)
    },
    best_start = NULL,
    find_best = function(n=10^3) {
      sample_prior()
      o = apply(samples_prior, 1, function(x) {
        names(x) <- pname
        l_post(x)
      })
      best_start <<- samples_prior[which.min(o), ]
    },
    l_post = function(x) {
      a <- -lprior(x)
      if (!is.finite(a)) return(Inf)
      a+ll(x)
    },
    fit = NULL,
    mod = NULL,
    best_par = NULL,
    fitmod = function(start=best_start, method='BFGS', n0=10^3, ctrl=list(trace=1), ...) {
      if (missing(start))
        if (is.null(best_start)) {
          message('Finding starting parameters value...')
          find_best(n0)
          print(best_start)
        }
      else
        best_start <<- start
      message('Fitting started...')
      fit <<- optim(best_start, l_post, method = method, control = ctrl)
      mod <<- sim_mod(fit$par,...)
      best_start <<- fit$par
      message('Fitting done!')
    },
    update = function(...) {
      message('Refit starting from last best parameters...')
      fit <<- optim(fit$par, l_post,...)
      mod <<- sim_mod(fit$par)
      best_start <<- fit$par
      message('Updated!')
    },
    plot_lambda = function() {
      epi <- mod$epi
      input <- mod$input
      lylim=range(epi$model['lambda',,])
      epi$model['lambda',,] %>% reshape2::melt() %>% 
        ggplot() + 
        geom_line(aes(Var2, value, col=age_lab[Var1])) +
        guides(col=FALSE)
    },
    plot_I = function() {
      o <- mod$epi$model %>% {
        data.frame(
          t = mod$input$time,
          I = colSums(.['Ic',,]),
          Q = colSums(.['Q',,]),
          R = colSums(.['R',,]),
          D = colSums(.['D',,])
        )
      }
      p <- o %>% ggplot() +
        geom_point(aes(as.Date(time), I, col='data'), cases) +
        geom_line(aes(t, I, col='Fitted'), alpha=.5) +
        geom_rect(aes(xmin = as.Date(actions$date_work[1]), xmax = as.Date(actions$date_work[2]), ymin = 0, ymax = Inf), fill='lightgray', alpha = 0.05) +
        geom_rect(aes(xmin = as.Date(actions$date_school[1]), xmax = as.Date(actions$date_school[2]), ymin = 0, ymax = Inf), fill='lightblue', alpha = 0.05) +
        geom_rect(aes(xmin = as.Date(actions$date_lockdown[1]), xmax = as.Date(actions$date_lockdown[2]), ymin = 0, ymax = Inf), fill='lightyellow', alpha = 0.05) +
        geom_point(aes(as.Date(time), I, col='data'), cases) +
        geom_line(aes(t, I, col='Fitted'), alpha=.5) +
        labs(col='', title='Incidence') + 
        xlab('Time') + ylab('Counts') +
        geom_vline(xintercept=today(), lty=3)
      print(p)
    }
  ),
  private = list(), 
  lock_objects=FALSE
) 