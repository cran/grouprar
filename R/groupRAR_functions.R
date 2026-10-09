###############################################################################
####################   Complete Randomization   ###############################
###############################################################################

CRDesign = function(k, p, ssn, nsim = 2000, alpha = 0.05,
                    test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k < 2){
    stop("k must be at least 2")
  }
  if(sum((p < 0) | (p > 1)) != 0){ stop("p must be a positive number between 0 and 1!") }
  if(length(p) != k){ stop("Length of parameter vector p must equal to k.") }
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
    obs.outcome = NULL
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # complete randomization
    alloc = sample(1:k, ssn, replace = TRUE, prob = rep(1/k, k))
    for(j in 1:k){
      obs.outcome[which(alloc == j)] = outcome[which(alloc == j), j]
    }
    pwCalc[s] = rar.test(alpha, obs.outcome, alloc, k, test.fun = test.fun)
    failure.rate[s] = mean(obs.outcome == 0)
    group.prop = rbind(group.prop, tabulate(alloc, nbins = k) / length(alloc))
    alloc.seq[s, ] = alloc
  }
  name = "Complete Randomization"
  out = RAR_Output(name, parameter=p, ssn,
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
####################   Randomized Play-the-winner rule   ######################
###############################################################################

RPWRule = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05,
                   test.fun = NULL, typeI = FALSE, seed = NULL){

  # check the accuracy of inputs
  if(k != 2){
    stop("RPWRule is defined for k = 2 only; use WeiUrn for more arms")
  }
  ## check length
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  ## check value
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  if(any(Y0 < 0) | (sum(Y0) <= 0)){
    stop("Y0 must be non-negative with a positive sum")
  }
  if(!is.null(seed)) set.seed(seed)

  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)
  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data(p, ssn)

    # assigned randomly
    Y = Y0
    sample.prob = Y / sum(Y)
    assign.group = NULL
    obs.outcome = NULL
    for(i in 1:ssn){
      # with components 1 or 2
      assign.group[i] = sample(c(1:k), 1, prob = sample.prob) # initial urn
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]] + 1
      }else{
        Y[-assign.group[i]] = Y[-assign.group[i]] + 1
      }
      sample.prob = Y / sum(Y)
    }
    group.prop = rbind(group.prop, tabulate(assign.group, nbins = k) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = rar.test(alpha, obs.outcome, assign.group, k, test.fun = test.fun)
    alloc.seq[s, ] = assign.group
  }
  name = "Randomized Play-the-winner Rule"
  out = RAR_Output(name, parameter=p, ssn,
                   assignment = assign.group, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
############################   Wei's Urn Model   ##############################
###############################################################################

WeiUrn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05,
                  test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k < 2){
    stop("k must be at least 2")
  }
  # check the accuracy of inputs
  ## check length
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  ## check value
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  if(any(Y0 < 0) | (sum(Y0) <= 0)){
    stop("Y0 must be non-negative with a positive sum")
  }
  if(!is.null(seed)) set.seed(seed)

  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)
  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, k = k)
    # assigned randomly
    Y = Y0
    sample.prob = Y / sum(Y)
    assign.group = NULL
    obs.outcome = NULL
    for(i in 1:ssn){
      assign.group[i] = sample(c(1:k), 1, prob = sample.prob) # initial urn
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]] + 1
      }else{
        Y[-assign.group[i]] = Y[-assign.group[i]] + (1/(k-1))
      }
      sample.prob = Y / sum(Y)
    }
    group.prop = rbind(group.prop, tabulate(assign.group, nbins = k) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = rar.test(alpha, obs.outcome, assign.group, k, test.fun = test.fun)
    alloc.seq[s, ] = assign.group
  }
  name = "Wei's Urn"
  out = RAR_Output(name, parameter=p, ssn,
                   assignment = assign.group, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
############################   Polya Urn Model   ##############################
###############################################################################

PolyaUrn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05,
                    test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k < 2){
    stop("k must be at least 2")
  }
  # check the accuracy of inputs
  ## check length
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  ## check value
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  if(any(Y0 < 0) | (sum(Y0) <= 0)){
    stop("Y0 must be non-negative with a positive sum")
  }
  if(!is.null(seed)) set.seed(seed)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL
  alloc.seq = matrix(NA, nsim, ssn)
  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, k = k)
    # assigned randomly
    Y = Y0
    sample.prob = Y / sum(Y)
    assign.group = NULL
    obs.outcome = NULL
    for(i in 1:ssn){
      assign.group[i] = sample(c(1:k), 1, prob = sample.prob) # initial urn
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]] + 1
      }
      sample.prob = Y / sum(Y)
    }

    failure.rate[s] = mean(obs.outcome == 0)
    group.prop = rbind(group.prop, tabulate(assign.group, nbins = k) / ssn)
    # an arm can lose all its patients, then no test is possible
    if(length(unique(assign.group)) < k){
      pwCalc[s] = NA
    }else{
      pwCalc[s] = rar.test(alpha, obs.outcome, assign.group, k, test.fun = test.fun)
    }
    alloc.seq[s, ] = assign.group
  }
  name = "Polya Urn"
  out = RAR_Output(name, parameter=p, ssn,
                   assignment = assign.group, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
#########################  Drop-the-loser Rule   ##############################
###############################################################################

DLRule = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05,
                  test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k < 2){
    stop("k must be at least 2")
  }
  # check the accuracy of inputs
  ## check length
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  ## check value
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  if(any(Y0 < 0)){
    stop("Y0 must be non-negative")
  }
  # add a immigration ball (balls are drawn in proportion to the positive part of the urn)
  Y0 = c(Y0, 1)
  if(!is.null(seed)) set.seed(seed)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # assigned randomly
    Y = Y0
    sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
    assign.group = NULL
    obs.outcome = NULL

    for(i in 1:ssn){
      assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      while(assign.k == (k+1)){
        Y[-assign.k] = Y[-assign.k] + 1
        # sample probability changed here
        sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
        assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      }
      assign.group[i] = assign.k
      obs.outcome[i] = outcome[i, assign.group[i]]
      # rule
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]] # changed
      }else{
        Y[assign.group[i]] = Y[assign.group[i]] - 1
      }
      sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
    }
    group.prop = rbind(group.prop, tabulate(assign.group, nbins = k) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = rar.test(alpha, obs.outcome, assign.group, k, test.fun = test.fun)
    alloc.seq[s, ] = assign.group
  }
  name = "Drop-the-loser Rule"
  out = RAR_Output(name, parameter = p, ssn,
                   assignment = assign.group, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
#####################   Generalized Drop-the-loser Rule   #####################
###############################################################################

GDLRule = function(k, p, ssn, aK, Y0 = NULL, nsim = 2000, alpha = 0.05,
                   test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k < 2){
    stop("k must be at least 2")
  }
  # check the accuracy of inputs
  ## check length
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  ## check value
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  ## check aK
  if(k != length(aK)){
    stop("Length of aK must be equal to k")
  }
  if(any(!is.finite(aK)) | any(aK < 0) | (sum(aK) == 0)){
    stop("aK must be non-negative with a positive sum")
  }
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  if(any(Y0 < 0)){
    stop("Y0 must be non-negative")
  }
  # add a immigration ball (balls are drawn in proportion to the positive part of the urn)
  Y0 = c(Y0, 1)
  if(!is.null(seed)) set.seed(seed)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # assigned randomly
    Y = Y0
    sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
    assign.group = NULL
    obs.outcome = NULL

    for(i in 1:ssn){
      assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      while(assign.k == (k+1)){
        Y[-assign.k] = Y[-assign.k] + aK #add ak here
        # sample probability changed here
        sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
        assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      }
      assign.group[i] = assign.k
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]]
      }else{
        Y[assign.group[i]] = Y[assign.group[i]] - 1
      }
      sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
    }
    group.prop = rbind(group.prop, tabulate(assign.group, nbins = k) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = rar.test(alpha, obs.outcome, assign.group, k, test.fun = test.fun)
    alloc.seq[s, ] = assign.group
  }
  name = "Generalized Drop-the-loser Rule"
  out = RAR_Output(name, parameter = p, ssn,
                   assignment = assign.group, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
####################   Bai, Hu, Shen's Urn Model   ############################
###############################################################################

Bai.Hu.Shen.Urn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05,
                           test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k < 2){
    stop("k must be at least 2")
  }
  # check the accuracy of inputs
  ## check length
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  ## check value
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  if(any(Y0 < 0) | (sum(Y0) <= 0)){
    stop("Y0 must be non-negative with a positive sum")
  }
  if(!is.null(seed)) set.seed(seed)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL
  alloc.seq = matrix(NA, nsim, ssn)
  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # assigned randomly
    Y = Y0
    sample.prob = Y / sum(Y)
    assign.group = NULL
    obs.outcome = NULL
    # successes and patients of each arm so far
    S = rep(0, k)
    N = rep(0, k)
    for(i in 1:ssn){
      assign.group[i] = sample(c(1:k), 1, prob = sample.prob) # initial urn
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]] + 1
      }else{
        # adaptive design 3 of Bai, Hu and Shen (2002): estimated success rates
        # of the previous patients, R = (S + 1) / (N + 1)
        R = (S + 1) / (N + 1)
        M = sum(R)
        Y[-assign.group[i]] = Y[-assign.group[i]] + R[-assign.group[i]] / (M-R[assign.group[i]])
      }
      S[assign.group[i]] = S[assign.group[i]] + obs.outcome[i]
      N[assign.group[i]] = N[assign.group[i]] + 1
      sample.prob = Y / sum(Y)
    }
    group.prop = rbind(group.prop, tabulate(assign.group, nbins = k) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = rar.test(alpha, obs.outcome, assign.group, k, test.fun = test.fun)
    alloc.seq[s, ] = assign.group
  }
  name = "Bai, Hu and Shen's Urn"
  out = RAR_Output(name, parameter = p, ssn,
                   assignment = assign.group, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
####################   Birth and Death Urn Model   ############################
###############################################################################

BirthDeathUrn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05,
                         test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k < 2){
    stop("k must be at least 2")
  }

  # check the accuracy of inputs
  ## check length
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  ## check value
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  if(any(Y0 < 0)){
    stop("Y0 must be non-negative")
  }
  # add a immigration ball (balls are drawn in proportion to the positive part of the urn)
  Y0 = c(Y0, 1)
  if(!is.null(seed)) set.seed(seed)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # assigned randomly
    Y = Y0
    sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
    assign.group = NULL
    obs.outcome = NULL

    for(i in 1:ssn){
      assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      while(assign.k == (k+1)){
        Y[-assign.k] = Y[-assign.k] + 1
        # sample probability changed here
        sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
        assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      }
      assign.group[i] = assign.k
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]] + 1
      }else{
        Y[assign.group[i]] = Y[assign.group[i]] - 1
      }
      sample.prob = pmax(Y, 0) / sum(pmax(Y, 0))
    }
    group.prop = rbind(group.prop, tabulate(assign.group, nbins = k) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = rar.test(alpha, obs.outcome, assign.group, k, test.fun = test.fun)
    alloc.seq[s, ] = assign.group
  }
  name = "Birth and Death Urn"
  out = RAR_Output(name, parameter = p, ssn,
                   assignment = assign.group, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
##########   Hu and Zhang's Doubly biased coin Design (Binary)  ###############
###############################################################################

DBCD_Bin = function(n0 = 20, p, k, ssn, theta0 = NULL, target.alloc = "RPW", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05,
                    allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0, monitor = NULL,
                    test.fun = NULL, typeI = FALSE, seed = NULL){

  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  check.size(n0, k, ssn)
  dsg = rar.design(k, target.alloc, continuous = FALSE, lower.bound, allocation, r, erade.alpha)
  sq = setup.monitor(monitor, k, ssn, n0, alpha, test.fun)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
    obs.outcome = NULL
    alloc.n0 = rep(NA, n0)
    outcome = generate_data_M(p, ssn, mRate, k)

    # initial allocation rule is complete randomization
    x = rep(1:k, n0/k)
    alloc.n0 = sample(x = x)

    # initial estimate for p
    p.hat = NULL
    sample.prob = NULL

    for(j in 1:k){
      obs.outcome[which(alloc.n0 == j)] = outcome[which(alloc.n0 == j), j]
    }
    extra.n = ssn - n0
    alloc = alloc.n0
    look.result = NA

    for(i in 1:extra.n){
      if(is.null(mRate)){
        prop.k = tabulate(alloc, nbins = k) / length(alloc)
        p.hat = calc_theta(p.hat, k, alloc, outcome, theta0)
      }else{
        prop.k = calc_prop(alloc, outcome, mRate, k)
        p.hat = calc_theta_M(p.hat, k, alloc, outcome, theta0)
      }
      sample.prob = rar.prob(prop.k, p.hat, dsg)
      assign.group = sample(x = c(1:k), 1, prob = sample.prob)
      alloc = c(alloc, assign.group)
      obs.outcome = c(obs.outcome, outcome[(n0+i), assign.group])

      # interim and final analyses of a monitored trial
      if(!is.null(sq) && ((n0+i) %in% sq$look.n)){
        j = which(sq$look.n == (n0+i))
        eff = if(is.null(mRate)) 1:(n0+i) else which(outcome[1:(n0+i), k+1] == 0)
        look.result = sq.test(sq, j, obs.outcome[eff], alloc[eff], k, continuous = FALSE)
        if(isTRUE(look.result == 1) | (j == length(sq$look.n))){
          sq$stage[s] = j
          sq$n.stop[s] = n0+i
          break
        }
      }
    }
    alloc.seq[s, seq_along(alloc)] = alloc

    nobs = length(alloc)
    if(!is.null(mRate)){
      obs.outcome = obs.outcome[outcome[1:nobs, k+1] == 0]
      alloc = alloc[outcome[1:nobs, k+1] == 0]
    }
    if(is.null(sq)){
      pwCalc[s] = rar.test(alpha, obs.outcome, alloc, k, test.fun = test.fun)
    }else{
      pwCalc[s] = look.result
    }
    failure.rate[s] = mean(obs.outcome == 0)
    # allocation proportions of all enrolled patients (Zhai et al., 2024)
    group.prop = rbind(group.prop, tabulate(alloc.seq[s, ], nbins = k) / sum(!is.na(alloc.seq[s, ])))
  }
  name = ifelse(allocation == "DBCD", "Hu and Zhang's DBCD", "ERADE")

  out = RAR_Output(name, parameter=p, ssn,
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq, sq = sq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
########   Hu and Zhang's Doubly biased coin Design (delayed+bin)  ############
###############################################################################

dyldDBCD_Bin = function(n0 = 20, p, k, ssn, ent.param, rspT.dist, rspT.param, theta0 = NULL,
                        target.alloc = "RPW", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05,
                        allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0,
                        test.fun = NULL, typeI = FALSE, seed = NULL) {

  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  check.size(n0, k, ssn)
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }
  dsg = rar.design(k, target.alloc, continuous = FALSE, lower.bound, allocation, r, erade.alpha)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)
  duration = NULL
  enrollment = NULL

  for (s in 1:nsim) {
    obs.outcome = NULL
    alloc.n0 = rep(NA, n0)
    outcome = generate_data_M(p, ssn, mRate = NULL, k)

    # entry Time
    entryT = NULL
    diff.entry = rexp(ssn - 1, 1 / ent.param)
    entryT = c(0, cumsum(diff.entry))

    # delayed time
    rspT = responseDist(rspT.dist, rspT.param, k, level = 2, ssn)
    cumRspT = apply(rspT, 2, function(x) x + entryT)

    # Missing Data
    if (!is.null(mRate)) {
      MD = MissingData(mRate, ssn)
      cumRspT[which(MD == 1), ] = Inf
    }

    # initial allocation rule is restricted randomization
    x = rep(1:k, n0 / k)
    alloc.n0 = sample(x = x)

    # initial estimate for p
    sample.prob = NULL
    for (j in 1:k) {
      obs.outcome[which(alloc.n0 == j)] = outcome[which(alloc.n0 == j), j]
    }
    alloc = alloc.n0

    obs.cumRspT = NULL
    # All observed cumulative time so far
    rsp.idx =  obs.outcome[1:n0] + 1
    rspT.idx = 2 * (alloc.n0 - 1) + rsp.idx
    # m[j]: the last patient who entered before the response of patient j is available
    m = NULL
    for (idx in 1:n0) {
      obs.cumRspT[idx] = cumRspT[idx,  rspT.idx[idx]]
      m[idx] = max(which(entryT <= obs.cumRspT[idx]))
    }
    # estimate p with the responses available before patient n0+1 enters
    p.hat = NULL
    avail = m <= n0
    for (t in 1:k) {
      p.hat[t] = (sum(obs.outcome[avail & alloc == t]) + theta0[t]) / (sum(avail & alloc == t) + 1)
    }
    # calculate sample prob
    prop.k = tabulate(alloc, nbins = k) / length(alloc)
    sample.prob = rar.prob(prop.k, p.hat, dsg)
    for (i in (n0 + 1):ssn) {
      alloc[i] = sample(c(1:k), 1, prob = sample.prob)
      obs.outcome[i] = outcome[i, alloc[i]]

      # detect is there anyone's outcome ready
      # time of m-th patient' outcome is available
      asgn.idx = alloc[i]
      rsp.idx =  obs.outcome[i] + 1
      rspT.idx = 2 * (asgn.idx - 1) + rsp.idx
      obs.cumRspT[i] = cumRspT[i, rspT.idx]

      # the outcome could be observed after the m-th entry
      m[i] = max(which(entryT <= obs.cumRspT[i]))

      # allocation probability of patient i+1: estimate from the responses available
      # before patient i+1 enters, with the current allocation proportions
      if (i < ssn) {
        avail = m <= i
        for (t in 1:k) {
          p.hat[t] = (sum(obs.outcome[avail & alloc == t]) + theta0[t]) / (sum(avail & alloc == t) + 1)
        }
        prop.k = tabulate(alloc, nbins = k) / length(alloc)
        sample.prob = rar.prob(prop.k, p.hat, dsg)
      }
    }
    alloc.seq[s, ] = alloc
    # time from the first entry to the last observed response
    duration[s] = if(any(is.finite(obs.cumRspT))) max(obs.cumRspT[is.finite(obs.cumRspT)]) else NA
    enrollment[s] = entryT[ssn]

    # re-adjust obs.outcome and assign.group
    if (!is.null(mRate)) {
      obs.outcome = obs.outcome[MD == 0]
      alloc = alloc[MD == 0]
    }

    pwCalc[s] = rar.test(alpha, obs.outcome, alloc, k, test.fun = test.fun)
    failure.rate[s] = mean(obs.outcome == 0)
    # allocation proportions of all enrolled patients (Zhai et al., 2024)
    group.prop = rbind(group.prop, tabulate(alloc.seq[s, ], nbins = k) / sum(!is.na(alloc.seq[s, ])))
  }
  name = ifelse(is.null(mRate),
                "Delayed DBCD without Missing Data",
                "Delayed DBCD with Missing Data")
  if(allocation == "ERADE") name = sub("DBCD", "ERADE", name, fixed = TRUE)
  out = RAR_Output(name, parameter=p, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq,
                   duration = duration, enrollment = enrollment)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
########   Hu and Zhang's Doubly biased coin Design (Continuous)  #############
###############################################################################

DBCD_Cont = function(n0 = 20, theta, k, ssn, theta0 = NULL, target.alloc = "Neyman", r = 2, nsim = 2000, alpha = 0.05,
                     allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0, monitor = NULL,
                     test.fun = NULL, typeI = FALSE, seed = NULL){

  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(FALSE, TRUE)] < 0) > 0){
    stop("The variance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  if((target.alloc == "ZR") & any(theta[c(TRUE, FALSE)] <= 0)){
    stop("The ZR target requires positive means")
  }
  check.size(n0, k, ssn)
  dsg = rar.design(k, target.alloc, continuous = TRUE, lower.bound, allocation, r, erade.alpha)
  sq = setup.monitor(monitor, k, ssn, n0, alpha, test.fun)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
    obs.outcome = NULL
    alloc.n0 = rep(NA, n0)
    outcome = generate_GaussianRsp_M(theta, ssn, mRate = NULL, k)

    # initial allocation rule is complete randomization
    x = rep(1:k, n0/k)
    alloc.n0 = sample(x = x)

    # initial estimate for p
    theta.hat = NULL
    sample.prob = NULL

    for(j in 1:k){
      obs.outcome[which(alloc.n0 == j)] = outcome[which(alloc.n0 == j), j]
    }

    extra.n = ssn - n0
    alloc = alloc.n0
    look.result = NA

    for(i in 1:extra.n){
      prop.k = tabulate(alloc, nbins = k) / length(alloc)
      theta.hat = calc_thetaGaussian(theta.hat, k, alloc, outcome, theta0)
      sample.prob = rar.prob(prop.k, theta.hat, dsg)
      assign.group = sample(x = c(1:k), 1, prob = sample.prob)
      alloc = c(alloc, assign.group)
      obs.outcome = c(obs.outcome, outcome[(n0+i), assign.group])

      # interim and final analyses of a monitored trial
      if(!is.null(sq) && ((n0+i) %in% sq$look.n)){
        j = which(sq$look.n == (n0+i))
        look.result = sq.test(sq, j, obs.outcome, alloc, k, continuous = TRUE)
        if(isTRUE(look.result == 1) | (j == length(sq$look.n))){
          sq$stage[s] = j
          sq$n.stop[s] = n0+i
          break
        }
      }
    }
    alloc.seq[s, seq_along(alloc)] = alloc

    if(is.null(sq)){
      pwCalc[s] = rar.test(alpha, obs.outcome, alloc, k, continuous = TRUE, test.fun = test.fun)
    }else{
      pwCalc[s] = look.result
    }
    #failure.rate[s] = mean(obs.outcome == 0)
    failure.rate[s] = mean(obs.outcome)
    group.prop = rbind(group.prop, tabulate(alloc, nbins = k) / length(alloc))
  }
  name = ifelse(allocation == "DBCD", "Hu and Zhang's DBCD (Gaussian Response)", "ERADE (Gaussian Response)")

  out = RAR_Output(name, parameter = theta, ssn,
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, continuous = TRUE, alloc.seq = alloc.seq, sq = sq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(theta = null.theta(theta)))
  }
  return(out)
}



###############################################################################
########   Hu and Zhang's Doubly biased coin Design (delayed+cont)  ###########
###############################################################################

dyldDBCD_Cont = function(n0 = 20, theta, k, ssn, ent.param, rspT.dist, rspT.param,
                         target.alloc = "Neyman", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05,
                         allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0,
                         test.fun = NULL, typeI = FALSE, seed = NULL){

  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(FALSE, TRUE)] < 0) > 0){
    stop("The variance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  if((target.alloc == "ZR") & any(theta[c(TRUE, FALSE)] <= 0)){
    stop("The ZR target requires positive means")
  }
  check.size(n0, k, ssn)
  dsg = rar.design(k, target.alloc, continuous = TRUE, lower.bound, allocation, r, erade.alpha)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  theta.hat.set = c()
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)
  duration = NULL
  enrollment = NULL

  for(s in 1:nsim){
    obs.outcome = NULL
    alloc.n0 = rep(NA, n0)
    outcome = generate_GaussianRsp_M(theta, ssn, mRate, k)

    # entry Time
    entryT = NULL
    diff.entry = rexp(ssn-1, 1/ent.param)
    entryT = c(0, cumsum(diff.entry))

    # delayed time
    rspT = responseDist(rspT.dist, rspT.param, k, level = 1, ssn)
    cumRspT = apply(rspT, 2, function(x) x + entryT )

    # Missing Data
    if(!is.null(mRate)){
      cumRspT[which(outcome[, k+1] == 1),] = Inf
    }

    # initial allocation rule is restricted randomization
    x = rep(1:k, ceiling(n0 / k))
    alloc.n0 = sample(x = x)[1:n0]

    # initial estimate for p
    sample.prob = NULL
    for(j in 1:k){
      obs.outcome[which(alloc.n0 == j)] = outcome[which(alloc.n0 == j), j]
    }
    alloc = alloc.n0

    obs.cumRspT = NULL
    # All observed cumulative time so far
    for(idx in 1:n0){
      obs.cumRspT[idx] = cumRspT[idx,  alloc.n0[idx]]
    }

    # estimate p
    theta.hat = NULL
    temp = calc_thetaGaussian_MD(theta.hat, k, alloc, outcome, entryT = entryT[n0+1], obsRspT = obs.cumRspT, mRate = mRate)
    theta.hat = temp[[1]]
    prop.k = temp[[2]]

    # calc sample prob (equal allocation while the target cannot be estimated)
    sample.prob = rar.prob(prop.k, theta.hat, dsg)
    for(i in (n0+1):ssn){

      alloc[i] = sample(c(1:k), 1, prob = sample.prob)
      obs.outcome[i] = outcome[i, alloc[i]]

      # time when the outcome of patient i is available
      obs.cumRspT[i] = cumRspT[i, alloc[i]]

      # allocation probability of patient i+1: estimate from the responses available
      # before patient i+1 enters, with the current allocation proportions
      if(i < ssn){
        temp = calc_thetaGaussian_MD(theta.hat, k, alloc, outcome, entryT[i+1], obsRspT = obs.cumRspT, mRate = mRate)
        theta.hat = temp[[1]]
        prop.k = temp[[2]]
        sample.prob = rar.prob(prop.k, theta.hat, dsg)
      }
    }
    alloc.seq[s, ] = alloc
    # time from the first entry to the last observed response
    duration[s] = if(any(is.finite(obs.cumRspT))) max(obs.cumRspT[is.finite(obs.cumRspT)]) else NA
    enrollment[s] = entryT[ssn]

    # re-adjust obs.outcome and assign.group
    if(!is.null(mRate)){
      effIdx = which(outcome[, k+1] == 0)
      obs.outcome = obs.outcome[effIdx]
      alloc = alloc[effIdx]
    }

    theta.hat.set = rbind(theta.hat.set, theta.hat)
    pwCalc[s] = rar.test(alpha, obs.outcome, alloc, k, continuous = TRUE, test.fun = test.fun)
    #failure.rate[s] = mean(obs.outcome == 0)
    failure.rate[s] = mean(obs.outcome)
    # allocation proportions of all enrolled patients (Zhai et al., 2024)
    group.prop = rbind(group.prop, tabulate(alloc.seq[s, ], nbins = k) / sum(!is.na(alloc.seq[s, ])))
  }
  name = ifelse(is.null(mRate),
                "Hu and Zhang Delayed DBCD without Missing Data",
                "Hu and Zhang Delayed DBCD with Missing Data")
  if(allocation == "ERADE") name = sub("Hu and Zhang Delayed DBCD", "Delayed ERADE", name, fixed = TRUE)

  out = RAR_Output(name, parameter=theta, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, continuous = TRUE, alloc.seq = alloc.seq,
                   duration = duration, enrollment = enrollment)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(theta = null.theta(theta)))
  }
  return(out)
}



###############################################################################
#################   Group doubly biased coin Design (Binary)  ##################
###############################################################################

Group.DBCD_Bin = function(n0 = 20, p, k, gsize.param, ssn, theta0 = NULL, target.alloc = "RPW", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05,
                          allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0, monitor = NULL,
                          test.fun = NULL, typeI = FALSE, seed = NULL){
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  if(!(gsize.param > 0)){
    stop("'gsize.param' must be positive")
  }
  check.size(n0, k, ssn)
  dsg = rar.design(k, target.alloc, continuous = FALSE, lower.bound, allocation, r, erade.alpha)
  sq = setup.monitor(monitor, k, ssn, n0, alpha, test.fun)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
    p.hat = NULL
    sample.prob = NULL
    obs.outcome = NULL
    gsize = NULL
    outcome = generate_data_M(p, ssn, mRate, k)
    alloc = c()

    # initial part
    initial = TRUE
    i = 1
    while (initial) {
      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
      pmtAlloc = sample(rep(1:k, ceiling(gsize[i] / k)))
      alloc = c(alloc, pmtAlloc[1:gsize[i]])
      if (sum(gsize) >= n0) {
        initial = FALSE
      }
      i = i+1
    }

    look.done = 0
    stopped = FALSE
    while(sum(gsize) < ssn){
      # interim analyses of a monitored trial at the end of a group
      if(!is.null(sq)){
        passed = which(sq$look.n[-length(sq$look.n)] <= sum(gsize))
        if(length(passed) > 0 && max(passed) > look.done){
          j = max(passed)
          look.done = j
          if(sq.look.group(sq, j, alloc, outcome, mRate, k, continuous = FALSE) == 1){
            sq$stage[s] = j
            sq$n.stop[s] = length(alloc)
            stopped = TRUE
            break
          }
        }
      }
      # allocation probability
      if(is.null(mRate)){
        p.hat = calc_theta(p.hat, k, alloc, outcome, theta0 = theta0)
      }else{
        p.hat = calc_theta_M(p.hat, k, alloc, outcome, theta0 = theta0)
      }
      prop.k = calc_prop(alloc, outcome, mRate, k)
      sample.prob = rar.prob(prop.k, p.hat, dsg)
      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
      alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = TRUE))
      i = i+1
    }

    alloc = alloc[1:min(length(alloc), ssn)]
    nobs = length(alloc)
    # final analysis of a monitored trial
    if(!is.null(sq) & !stopped){
      j = length(sq$look.n)
      look.result = sq.look.group(sq, j, alloc, outcome, mRate, k, continuous = FALSE)
      sq$stage[s] = j
      sq$n.stop[s] = nobs
    }
    alloc.seq[s, 1:nobs] = alloc
    for(j in 1:k){
      obs.outcome[which(alloc == j)] = outcome[which(alloc == j), j]
    }
    if(!is.null(mRate)){
      obs.outcome = obs.outcome[outcome[1:nobs, k+1] == 0]
      alloc = alloc[outcome[1:nobs, k+1] == 0]
    }

    if(is.null(sq)){
      pwCalc[s] = rar.test(alpha, obs.outcome, alloc, k, test.fun = test.fun)
    }else{
      pwCalc[s] = if(stopped) 1 else look.result
    }

    failure.rate[s] = mean(obs.outcome == 0)
    # allocation proportions of all enrolled patients (Zhai et al., 2024)
    group.prop = rbind(group.prop, tabulate(alloc.seq[s, ], nbins = k) / sum(!is.na(alloc.seq[s, ])))
  }
  name = ifelse(is.null(mRate),
                "Group DBCD with Binary Response (No Missing)",
                "Group DBCD with Binary Response (Random Missing)")
  if(allocation == "ERADE") name = sub("DBCD", "ERADE", name, fixed = TRUE)

  out = RAR_Output(name, parameter=p, ssn = c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate, #shouldn't be reported
                   pwCalc, k, alloc.seq = alloc.seq, sq = sq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
###############   Group doubly biased coin Design (delayed+bin)  ###############
###############################################################################

Group.dyldDBCD_Bin = function(n0 = 20, p, k, ssn, gsize.param, rspT.dist, rspT.param,
                              theta0 = NULL, target.alloc = "RPW",  r = 2, nsim = 2000, eTime = 7,  mRate = NULL, alpha = 0.05,
                              allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0,
                              test.fun = NULL, typeI = FALSE, seed = NULL){

  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  if(!(gsize.param > 0)){
    stop("'gsize.param' must be positive")
  }
  check.size(n0, k, ssn)
  dsg = rar.design(k, target.alloc, continuous = FALSE, lower.bound, allocation, r, erade.alpha)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)
  duration = NULL
  enrollment = NULL

  for(s in 1:nsim){
    p.hat = NULL
    obs.outcome = NULL
    sample.prob = NULL
    entry = NULL
    gsize = NULL
    obsRspT = c()
    alloc = c()
    outcome = generate_data_M(p, ssn, mRate = mRate, k)
    rspT = responseDist(rspT.dist, rspT.param, k, level = 2, ssn)

    # initial part
    initial = TRUE
    i = 1
    entry[i] = 0

    while (initial) {
      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
      #if(gsize[i] == 0) next
      pmtAlloc = sample(rep(1:k, ceiling(gsize[i] / k)))
      alloc = c(alloc, pmtAlloc[1:gsize[i]])

      if (sum(gsize) >= n0) {
        break
      }

      temp = calc_RspT(alloc, outcome, rspT, gsize, entry, continuous = FALSE)
      obs.outcome = c(obs.outcome, temp[,2])
      obsRspT = c(obsRspT, temp[,1])
      i = i+1
      entry[i] = eTime
    }

    while (sum(gsize) < ssn) {

      # calculate the time for delayed response
      temp = calc_RspT(alloc, outcome, rspT, gsize, entry, continuous = FALSE)
      obs.outcome = c(obs.outcome, temp[,2])
      obsRspT = c(obsRspT, temp[,1])

      # update the entry time
      i = i+1
      entry[i] = eTime
      entryT = sum(entry)

      # calculate allocation probability
      temp = calc_theta_MD(p.hat, k, alloc, outcome, theta0, entryT, obsRspT, mRate)
      p.hat = temp[[1]]
      prop.k = temp[[2]]
      sample.prob = rar.prob(prop.k, p.hat, dsg)

      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf) #rpois(1, gsize.param)
      alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = TRUE))
    }

    alloc = alloc[1:ssn]
    temp = calc_RspT(alloc, outcome, rspT, gsize, entry, adjust = TRUE, continuous = FALSE)
    obs.outcome = c(obs.outcome, temp[, 2])
    obsRspT = c(obsRspT, temp[, 1])
    if(is.null(mRate)){effIdx = c(1:ssn)}else{effIdx = which(outcome[, k+1] == 0)}
    alloc.seq[s, ] = alloc
    # time from the first entry to the last observed response
    duration[s] = if(length(effIdx) > 0) max(obsRspT[effIdx]) else NA
    enrollment[s] = sum(entry)

    pwCalc[s] = rar.test(alpha, obs.outcome[effIdx], alloc[effIdx], k, test.fun = test.fun)

    failure.rate[s] = mean(obs.outcome[effIdx] == 0)
    # allocation proportions of all enrolled patients (Zhai et al., 2024)
    group.prop = rbind(group.prop, tabulate(alloc.seq[s, ], nbins = k) / sum(!is.na(alloc.seq[s, ])))
  }

  name = ifelse(is.null(mRate),
                "Group DBCD with Binary Delayed Response (No Missing)",
                "Group DBCD with Binary Delayed Response (Random Missing)")
  if(allocation == "ERADE") name = sub("DBCD", "ERADE", name, fixed = TRUE)

  out = RAR_Output(name, parameter=p, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, alloc.seq = alloc.seq,
                   duration = duration, enrollment = enrollment)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(p = rep(mean(p), k)))
  }
  return(out)
}



###############################################################################
##############   Group doubly biased coin Design (Continuous)  #################
###############################################################################

Group.DBCD_Cont = function(n0 = 20, theta, k, gsize.param, ssn, theta0 = NULL, target.alloc = "Neyman", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05,
                           allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0, monitor = NULL,
                           test.fun = NULL, typeI = FALSE, seed = NULL){
  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(FALSE, TRUE)] < 0) > 0){
    stop("The variance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  if((target.alloc == "ZR") & any(theta[c(TRUE, FALSE)] <= 0)){
    stop("The ZR target requires positive means")
  }
  if(!(gsize.param > 0)){
    stop("'gsize.param' must be positive")
  }
  check.size(n0, k, ssn)
  dsg = rar.design(k, target.alloc, continuous = TRUE, lower.bound, allocation, r, erade.alpha)
  sq = setup.monitor(monitor, k, ssn, n0, alpha, test.fun)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)

  for(s in 1:nsim){
      theta.hat = NULL
      sample.prob = NULL
      obs.outcome = NULL
      gsize = NULL
      outcome = generate_GaussianRsp_M(theta, ssn, mRate, k)
      alloc = c()

      # initial part
      initial = TRUE
      i = 1
      while (initial) {
        gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
        pmtAlloc = sample(rep(1:k, ceiling(gsize[i] / k)))
        alloc = c(alloc, pmtAlloc[1:gsize[i]])
        if (sum(gsize) >= n0) {
          initial = FALSE
        }
        i = i+1
      }

      look.done = 0
      stopped = FALSE
      while(sum(gsize) < ssn){
        # interim analyses of a monitored trial at the end of a group
        if(!is.null(sq)){
          passed = which(sq$look.n[-length(sq$look.n)] <= sum(gsize))
          if(length(passed) > 0 && max(passed) > look.done){
            j = max(passed)
            look.done = j
            if(sq.look.group(sq, j, alloc, outcome, mRate, k, continuous = TRUE) == 1){
              sq$stage[s] = j
              sq$n.stop[s] = length(alloc)
              stopped = TRUE
              break
            }
          }
        }
        # allocation probability
        theta.hat = calc_thetaGaussian(theta.hat, k, alloc, outcome, theta0, mRate)
        prop.k = calc_prop(alloc, outcome, mRate, k)
        sample.prob = rar.prob(prop.k, theta.hat, dsg)
        gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
        alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = TRUE))
        i = i+1
      }

      alloc = alloc[1:min(length(alloc), ssn)]
      nobs = length(alloc)
      # final analysis of a monitored trial
      if(!is.null(sq) & !stopped){
        j = length(sq$look.n)
        look.result = sq.look.group(sq, j, alloc, outcome, mRate, k, continuous = TRUE)
        sq$stage[s] = j
        sq$n.stop[s] = nobs
      }
      alloc.seq[s, 1:nobs] = alloc
      for(j in 1:k){
        obs.outcome[which(alloc == j)] = outcome[which(alloc == j), j]
      }
      if(!is.null(mRate)){
        obs.outcome = obs.outcome[outcome[1:nobs, k+1] == 0]
        alloc = alloc[outcome[1:nobs, k+1] == 0]
      }

      if(is.null(sq)){
        pwCalc[s] = rar.test(alpha, obs.outcome, alloc, k, continuous = TRUE, test.fun = test.fun)
      }else{
        pwCalc[s] = if(stopped) 1 else look.result
      }

      failure.rate[s] = mean(obs.outcome)
      # allocation proportions of all enrolled patients (Zhai et al., 2024)
      group.prop = rbind(group.prop, tabulate(alloc.seq[s, ], nbins = k) / sum(!is.na(alloc.seq[s, ])))
  }

  name = ifelse(is.null(mRate),
                "Group DBCD with Continuous Response (No Missing)",
                "Group DBCD with Continuous Response (Random Missing)")
  if(allocation == "ERADE") name = sub("DBCD", "ERADE", name, fixed = TRUE)

  out = RAR_Output(name, parameter = theta, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, continuous = TRUE, alloc.seq = alloc.seq, sq = sq)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(theta = null.theta(theta)))
  }
  return(out)
}



###############################################################################
#############   Group doubly biased coin Design (delayed+Cont)  ################
###############################################################################

Group.dyldDBCD_Cont = function(n0 = 20, theta, k, ssn, gsize.param, rspT.dist, rspT.param,
                               target.alloc = "Neyman",  r = 2, nsim = 2000, eTime = 7,  mRate = NULL, alpha = 0.05,
                               allocation = "DBCD", erade.alpha = 0.5, lower.bound = 0,
                               test.fun = NULL, typeI = FALSE, seed = NULL){

  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(FALSE, TRUE)] < 0) > 0){
    stop("The variance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participants 'n0' must be a multiple of k")
  }
  if((target.alloc == "ZR") & any(theta[c(TRUE, FALSE)] <= 0)){
    stop("The ZR target requires positive means")
  }
  if(!(gsize.param > 0)){
    stop("'gsize.param' must be positive")
  }
  check.size(n0, k, ssn)
  dsg = rar.design(k, target.alloc, continuous = TRUE, lower.bound, allocation, r, erade.alpha)
  if(!is.null(seed)) set.seed(seed)

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
  alloc.seq = matrix(NA, nsim, ssn)
  duration = NULL
  enrollment = NULL

  for(s in 1:nsim){
    theta.hat = NULL
    obs.outcome = NULL
    sample.prob = NULL
    entry = NULL
    gsize = NULL
    obsRspT = c()
    alloc = c()
    outcome = generate_GaussianRsp_M(theta, ssn, mRate, k)
    rspT = responseDist(rspT.dist, rspT.param, k, level = 1, ssn)

    # initial part
    initial = TRUE
    i = 1
    entry[i] = 0

    while (initial) {
      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
      pmtAlloc = sample(rep(1:k, ceiling(gsize[i] / k)))
      alloc = c(alloc, pmtAlloc[1:gsize[i]])

      if (sum(gsize) >= n0) {
        break
      }

      temp = calc_RspT(alloc, outcome, rspT, gsize, entry)
      obs.outcome = c(obs.outcome, temp[,2])
      obsRspT = c(obsRspT, temp[,1])
      i = i+1
      entry[i] = eTime
    }

    while (sum(gsize) < ssn) {

      # calculate the time for delayed response
      temp = calc_RspT(alloc, outcome, rspT, gsize, entry)
      obs.outcome = c(obs.outcome, temp[,2])
      obsRspT = c(obsRspT, temp[,1])

      # update the entry time
      i = i+1
      entry[i] = eTime
      entryT = sum(entry)

      # calculate allocation probability (equal allocation while the target cannot be estimated)
      temp = calc_thetaGaussian_MD(theta.hat, k, alloc, outcome, entryT, obsRspT, mRate)
      theta.hat = temp[[1]]
      prop.k = temp[[2]]
      sample.prob = rar.prob(prop.k, theta.hat, dsg)

      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf) #rpois(1, gsize.param)
      alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = TRUE))
    }

    alloc = alloc[1:ssn]
    temp = calc_RspT(alloc, outcome, rspT, gsize, entry, adjust = TRUE)
    obs.outcome = c(obs.outcome, temp[, 2])
    obsRspT = c(obsRspT, temp[, 1])
    if(is.null(mRate)){effIdx = c(1:ssn)}else{effIdx = which(outcome[, k+1] == 0)}
    alloc.seq[s, ] = alloc
    # time from the first entry to the last observed response
    duration[s] = if(length(effIdx) > 0) max(obsRspT[effIdx]) else NA
    enrollment[s] = sum(entry)

    pwCalc[s] = rar.test(alpha, obs.outcome[effIdx], alloc[effIdx], k, continuous = TRUE, test.fun = test.fun)
    failure.rate[s] = mean(obs.outcome[effIdx])
    #failure.rate[s] = mean(obs.outcome[effIdx] == 0)
    # allocation proportions of all enrolled patients (Zhai et al., 2024)
    group.prop = rbind(group.prop, tabulate(alloc.seq[s, ], nbins = k) / sum(!is.na(alloc.seq[s, ])))
  }
  name = ifelse(is.null(mRate),
                "Group DBCD with Continuous Delayed Response (No Missing)",
                "Group DBCD with Continuous Delayed Response (Random Missing)")
  if(allocation == "ERADE") name = sub("DBCD", "ERADE", name, fixed = TRUE)
  out = RAR_Output(name, parameter = theta, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                   assignment = alloc, propotion = group.prop,
                   failRate = failure.rate,
                   pwCalc, k, continuous = TRUE, alloc.seq = alloc.seq,
                   duration = duration, enrollment = enrollment)
  if(typeI){
    out[["type I error"]] = typeI.error(match.call(), sys.function(), parent.frame(), list(theta = null.theta(theta)))
  }
  return(out)
}
