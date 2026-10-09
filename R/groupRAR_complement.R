generate_data = function(p, sample.size){
  outcome = matrix(data = NA, nrow = sample.size, ncol = 2)
  for(j in 1:2){
    outcome[,j] = rbinom(sample.size, 1, p[j])
  }
  return(outcome)
}



generate_data_M = function(p, sample.size, mRate=NULL, k){
  outcome = matrix(data = NA, nrow = sample.size, ncol = k)
  for(j in 1:k){
    outcome[,j] = rbinom(sample.size, 1, p[j])
  }
  if(!is.null(mRate)){
    outcome = cbind(outcome, rbinom(sample.size, 1, mRate))
  }
  return(outcome)
}



generate_GaussianRsp = function(p, k, sample.size){
  outcome = matrix(data = NA, nrow = sample.size, ncol = k)
  numPara = length(p)
  for(j in 1:(numPara/2)){
    outcome[,j] = rnorm(sample.size, p[(2*j - 1)], sqrt(p[(2*j)]))
  }
  return(outcome)
}



# hypothesis test (t test, two sided)
ttest.2 = function(alpha, obs.outcome, assign.group){
  x1 = obs.outcome[assign.group == 1]
  x2 = obs.outcome[assign.group == 2]
  n1 = length(x1)
  n2 = length(x2)
  if((n1 == 0) | (n2 == 0)){
    warning("All subjects were assigned to the same group")
    return(NA)
  }
  if(is.na(var(x1)) | is.na(var(x2))){
    return(NA)
  }
  if((var(x1) == 0) & (var(x2) == 0)){
    # no variability in either arm: reject when the means differ
    return(ifelse(mean(x1) != mean(x2), 1, 0))
  }
  t.stats =  (mean(x1) - mean(x2)) / sqrt(var(x1) / n1 + var(x2) / n2)
  df = (var(x1) / n1 + var(x2) / n2)^2 / (var(x1)^2 / (n1^2 * (n1 - 1)) + var(x2)^2 / (n2^2 * (n2 - 1)))
  p.value = 2 * pt(-abs(t.stats), df)
  test.result = ifelse(p.value <= alpha, 1, 0)
  return(test.result)
}



# Wald statistic for H0: all k arms have the same mean (chi-squared with k-1 df)
# continuous = TRUE uses the sample variance of each arm instead of p(1-p)
# returns NA when an arm has no subjects or the covariance matrix is singular
wald.stat = function(obs.outcome, assign.group, k, continuous = FALSE){
  n = tabulate(assign.group, nbins = k)
  if(any(n == 0)) return(NA)
  grp = factor(assign.group, levels = 1:k)
  hat.p = as.vector(tapply(obs.outcome, grp, mean))
  if(continuous){
    hat.var = as.vector(tapply(obs.outcome, grp, var))
  }else{
    hat.var = hat.p * (1-hat.p)
  }
  if(any(!is.finite(hat.var))) return(NA)
  hat.pc = hat.p[-k] - hat.p[k]
  Sigma = diag(hat.var[-k] / n[-k], nrow = k-1) + hat.var[k] / n[k]
  if(!continuous && (rcond(Sigma) < 1e-10)){
    # arms without variability: adjusted variances of Agresti and Caffo (2000)
    tilde.p = (as.vector(tapply(obs.outcome, grp, sum)) + 1) / (n + 2)
    hat.var = tilde.p * (1-tilde.p)
    Sigma = diag(hat.var[-k] / n[-k], nrow = k-1) + hat.var[k] / n[k]
  }
  # singular, also up to rounding error
  if(rcond(Sigma) < 1e-10) return(NA)
  return(as.numeric(t(hat.pc) %*% solve(Sigma) %*% hat.pc))
}



# hypothesis test (chi-squared test)
chisq.test.k = function(alpha, obs.outcome, assign.group, k, continuous = FALSE){
  if(any(tabulate(assign.group, nbins = k) == 0)){
    warning("At least one arm has no subjects")
    return(NA)
  }
  chisq.stats = wald.stat(obs.outcome, assign.group, k, continuous)
  if(is.na(chisq.stats)){
    warning("variance-covariance matrix is exactly singular")
    test.result = NA
  }else{
    df = k-1
    p.value = pchisq(chisq.stats, df, lower.tail=FALSE)
    test.result = ifelse(p.value <= alpha, 1, 0)
  }
  return(test.result)
}



# k-arm optimal allocation (Tymofyeyev, Rosenberger and Hu, 2007, and its continuous analogue)
# minimizes sum(w * rho) / phi(rho) subject to rho >= B, where phi(rho) is the
# noncentrality of the Wald test of equal means, sum(a * (theta - weighted mean)^2)
# with a = rho / v. w = 1 minimizes the total sample size, w = q (binary) the expected
# failures and w = mu (continuous) the total expected response. For k = 2 and B = 0 this
# gives Neyman, RSIHR and the Zhang and Rosenberger (2006) target, respectively.
opt.rho = function(theta, v, w, B = 0){
  k = length(theta)
  if(any(!is.finite(c(theta, v, w))) || any(v <= 0) || any(w < 0)) return(rep(NA, k))
  if(diff(range(theta)) == 0) return(rep(1/k, k))
  if(k == 2){
    rho1 = sqrt(v[1] * w[2]) / (sqrt(v[1] * w[2]) + sqrt(v[2] * w[1]))
    rho1 = min(max(rho1, B), 1 - B)
    return(c(rho1, 1 - rho1))
  }
  obj = function(rho){
    a = rho / v
    phi = sum(a * (theta - sum(a * theta) / sum(a))^2)
    # a finite value keeps optimize() from warning when two arms have equal estimates
    if(phi <= 0) return(.Machine$double.xmax)
    sum(w * rho) / phi
  }
  # all arms except at most two are at the lower bound at the optimum, so search
  # over every pair of arms sharing the remaining mass
  s = 1 - k * B
  best.rho = rep(1/k, k)
  best.val = obj(best.rho)
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      pair.rho = function(t){
        rho = rep(B, k)
        rho[i] = B + s * t
        rho[j] = B + s * (1 - t)
        rho
      }
      cand = c(optimize(function(t) obj(pair.rho(t)), c(0, 1), tol = 1e-10)$minimum, 0, 1)
      for(t in cand){
        val = obj(pair.rho(t))
        if(val < best.val){
          best.val = val
          best.rho = pair.rho(t)
        }
      }
    }
  }
  return(best.rho)
}



# DBCD common
# establish the target distribution
# discrete
target.rho = function(para.set, target.alc, lower.bound = 0){
  if(!(target.alc %in% c("Neyman", "RSIHR", "RPW", "WeisUrn", "OptimalNeyman", "OptimalRSIHR"))){
    stop("Your target allocation is not in the list")
  }else if(target.alc == "Neyman"){
    temp = sqrt(para.set * (1-para.set))
    rho = temp / sum(temp)
  }else if(target.alc == "RSIHR"){
    temp = sqrt(para.set)
    rho = temp / sum(temp)
  }else if(target.alc == "RPW"){
    temp = 1/(1-para.set)
    rho = temp / sum(temp)
  }else if(target.alc == "WeisUrn"){
    temp = 1/(1-para.set)
    rho = temp / sum(temp)
  }else if(target.alc == "OptimalNeyman"){
    rho = opt.rho(para.set, para.set * (1-para.set), rep(1, length(para.set)), lower.bound)
  }else if(target.alc == "OptimalRSIHR"){
    rho = opt.rho(para.set, para.set * (1-para.set), 1-para.set, lower.bound)
  }
  return(rho)
}



target.rho.Ctinuous = function(para.set, target.alc, lower.bound = 0){
  n = length(para.set)
  mu = NULL
  sigma2 = NULL

  for(i in 1:(n/2)){
    mu[i] =  para.set[(2*i-1)]
    sigma2[i] =  para.set[(2*i)]
  }

  if(!(target.alc %in% c("Neyman", "ZR", "DaOptimal", "OptimalNeyman"))){
    stop("Your target allocation is not in the list")
  }else if(target.alc == "Neyman"){
    rho = sqrt(sigma2) / sum(sqrt(sigma2))
  }else if(target.alc == "ZR"){
    # Zhang and Rosenberger (2006), smaller responses are better and means must be positive
    if(any(mu <= 0, na.rm = TRUE)){
      rho = rep(NA, n/2)
    }else if(n/2 == 2){
      # their rule (7): use 1/2 when the optimal allocation favors the worse arm
      r.zr = sqrt(sigma2[1] * mu[2]) / sqrt(sigma2[2] * mu[1])
      if(is.finite(r.zr) && ((mu[1] < mu[2] && r.zr > 1) || (mu[1] > mu[2] && r.zr < 1))){
        rho = opt.rho(mu, sigma2, mu, lower.bound)
      }else{
        rho = c(1/2, 1/2)
      }
    }else{
      # k-arm version: minimize the total expected response for a fixed noncentrality
      rho = opt.rho(mu, sigma2, mu, lower.bound)
    }
  }else if(target.alc == "DaOptimal"){
    temp =( sqrt(sigma2) )^(4/3)
    rho = (temp) / sum(temp)
  }else if(target.alc == "OptimalNeyman"){
    rho = opt.rho(mu, sigma2, rep(1, n/2), lower.bound)
  }
  return(rho)
}



# fall back to equal allocation when the estimated target cannot be computed
# (e.g. no responses yet or a zero variance estimate)
safe.rho = function(rho, k){
  if(length(rho) != k){
    stop("The target allocation has the wrong length")
  }
  if(any(!is.finite(rho)) || any(rho < 0)){
    rho = rep(1/k, k)
  }
  return(rho)
}



# Hu and Zhang (2004) allocation function for k arms
# x: current allocation proportions (length k), y: estimated target allocation (length k)
# each arm is normalized jointly against all k arms
g.func = function(x, y, r){
  x = as.numeric(x)
  y = as.numeric(y)
  if(any(x == 0)){
    # arms without any allocation yet share all the probability
    g.xy = as.numeric(x == 0) / sum(x == 0)
  }else{
    tmp = y * (y / x) ^ r
    g.xy = tmp / sum(tmp)
  }
  return(g.xy)
}



# estimate p (adjusted)
calc_theta = function(p.hat, k, alloc, outcome, theta0){
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }
  for(j in 1:k){
    denom = length(outcome[which(alloc == j), j]) + 1
    numer = sum(outcome[which(alloc == j), j]) + theta0[j]
    p.hat[j] = numer / denom
  }
  return(p.hat)
}



# estimate p (adjusted)
calc_theta_M = function(p.hat, k, alloc, outcome, theta0){
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }
  for(j in 1:k){
    effIndex = intersect(which((alloc == j)), which((outcome[, (k+1)] != 1)))
    denom = length(effIndex) + 1
    numer = sum(outcome[effIndex, j]) + theta0[j]
    p.hat[j] = numer / denom
  }
  return(p.hat)
}



# especially for delayed
calc_theta_MD = function(p.hat, k, alloc, outcome, theta0, entryT, obsRspT, mRate){
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }
  propk = NULL
  for(j in 1:k){
    if(is.null(mRate)) outcome = cbind(outcome, 0)
    effIndex = intersect(which((alloc == j)), intersect(which((outcome[, (k+1)] != 1)), which(entryT >= obsRspT)))
    # calculate theta
    denom = length(effIndex) + 1
    numer = sum(outcome[effIndex, j]) + theta0[j]
    p.hat[j] = numer / denom
    # adjust alloc proportion
    #propk[j] = sum(alloc == j) - length(MDIndex)
    propk[j] = sum(alloc == j)
  }
  return(list(p.hat, propk / sum(propk)))
}



calc_thetaGaussian = function(p.hat, k, alloc, outcome, theta0, mRate = NULL){
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }
  for(j in 1:k){
    effIndex = which(alloc == j)
    if(!is.null(mRate)){
      # exclude missing responses
      effIndex = intersect(effIndex, which(outcome[, (k+1)] != 1))
    }
    mu = mean(outcome[effIndex, j])
    # the variance cannot be estimated from fewer than two responses
    sigma2 = ifelse(length(effIndex) >= 2, mean((outcome[effIndex, j] - mu)^2), NA)
    p.hat[(2*j - 1)] = mu
    p.hat[(2*j)] = sigma2
  }
  return(p.hat)
}



# Delayed
responseDist = function(dist, param, k, level, sample.size){

  if(!(dist %in% c("exponential", "normal", "uniform"))){
    stop("'dist' must be one of the element of 'exponential', 'normal' or 'uniform'")
  }

  npar = ifelse(dist == "exponential", 1, 2) * k * level
  if(length(param) != npar){
    stop("'rspT.param' must have length ", npar, " for this distribution")
  }

  rspT = matrix(data = NA, nrow = sample.size, ncol = (k*level))

  if(dist == "exponential"){
    for(i in 1:(k*level)){
      rspT[,i] = rexp(sample.size, 1/param[i])
    }
  }else if(dist == "normal"){
    for(i in 1:(k*level)){
      rspT[,i] = rnorm(sample.size, param[2*i-1], param[2*i])
    }
    rspT[which(rspT<0)] = 0
  }else if(dist == "uniform"){
    for(i in 1:(k*level)){
      rspT[,i] = runif(sample.size, param[2*i-1], param[2*i])
    }
  }
  return(rspT)
}



# Missing Data
MissingData = function(mRate, sample.size){
  # if  = 1, it means this is a missing data
  MisInfo = rbinom(sample.size, 1, mRate)
  return(MisInfo)
}



# collect the simulation results of a design
RAR_Output = function(name, parameter, ssn,
                      assignment, propotion,
                      failRate,
                      pwCalc, k,
                      continuous = FALSE, alloc.seq = NULL,
                      duration = NULL, enrollment = NULL, sq = NULL){

  trt = paste0("treatment ", LETTERS[1:k])
  propotion = matrix(propotion, ncol = k)
  colnames(propotion) = trt

  # statistics
  Mean.propotion = colMeans(propotion)
  sdprop = apply(propotion, 2, sd)
  Mean.failRate = mean(failRate)
  sdfrt = sd(failRate)
  power = mean(na.omit(pwCalc))

  # parameter
  if(continuous){
    names(parameter) = as.vector(rbind(paste0("mu", LETTERS[1:k]), paste0("sigma2", LETTERS[1:k])))
  }else{
    names(parameter) = paste0('p',  LETTERS[1:k])
  }

  # output
  outputList = list(method = name,
                    "sample size" = ssn,
                    "parameter" = parameter,
                    "propotion" = Mean.propotion,
                    "sd of propotion" = sdprop,
                    "failure rate" = Mean.failRate,
                    "sd of failure rate" = sdfrt,
                    'power' = power,
                    "data: failureRate" = failRate,
                    "data: test" = pwCalc,
                    "data: assignment" = assignment,
                    "data: propotion" = as.data.frame(propotion))
  if(!is.null(alloc.seq)){
    outputList[["data: allocation"]] = alloc.seq
  }
  # delayed designs: time from the first entry to the last observed response
  if(!is.null(duration)){
    outputList[["duration"]] = mean(duration, na.rm = TRUE)
    outputList[["sd of duration"]] = sd(duration, na.rm = TRUE)
    outputList[["enrollment duration"]] = mean(enrollment)
    outputList[["data: duration"]] = duration
    outputList[["data: enrollment"]] = enrollment
  }
  # group sequential monitoring
  if(!is.null(sq)){
    outputList[["boundary"]] = sq$bound
    outputList[["stopping probability"]] = stats::setNames(tabulate(sq$stage, nbins = length(sq$look.n)) / length(sq$stage),
                                                          paste0("look ", seq_along(sq$look.n)))
    outputList[["expected sample size"]] = mean(sq$n.stop)
    outputList[["data: stage"]] = sq$stage
    outputList[["data: sample size"]] = sq$n.stop
  }
  attr(outputList, "response") = ifelse(continuous, "continuous", "binary")
  class(outputList) = "grouprar"
  return(outputList)
}


# allocation proportions used by the allocation function: all enrolled patients, also
# those with missing responses (Zhai et al., 2024, Section 2.3)
calc_prop = function(alloc, outcome, mRate, k){
  return(tabulate(alloc, nbins = k) / length(alloc))
}


calc_RspT = function(alloc, outcome, rspT, gsize, entry, adjust = FALSE, continuous = TRUE){
  numGroup = length(gsize)
  obstime = rep(NA, sum(gsize))
  obscome = rep(NA, sum(gsize))
  if(adjust){
    idx = (sum(gsize[1:numGroup-1])+1) : length(alloc)
    if(length(idx) >= 2){if(idx[1] > idx[2]) idx = length(alloc)}
  }else{
    idx = (sum(gsize[1:numGroup-1])+1) : (sum(gsize[1:numGroup-1])+gsize[numGroup])
  }
  for(i in idx){
    if(continuous){
      obscome[i] = outcome[i, alloc[i]]
      obstime[i] = rspT[i, (alloc[i])] + sum(entry)
    }else{
      obscome[i] = outcome[i, alloc[i]]
      obstime[i] = rspT[i, (2 * alloc[i] + obscome[i] - 1)] + sum(entry)
    }

  }
  result = na.omit(data.frame(obstime, obscome))
  return(result)
}


# for current function it's not apply to the one without missing
calc_thetaGaussian_MD = function(para.hat, k, alloc, outcome, entryT, obsRspT, mRate){
  propk = NULL
  for(j in 1:k){
    if(is.null(mRate)) outcome = cbind(outcome, 0)
    effIndex = intersect(which((alloc == j)), intersect(which((outcome[, (k+1)] != 1)), which(entryT >= obsRspT)))
    # calculate theta
    mu.hat = mean(outcome[effIndex, j])
    sigma2.hat = var(outcome[effIndex, j])
    para.hat[(2*j - 1)] = mu.hat
    para.hat[(2*j)] = sigma2.hat

    # adjust allocation proportion
    propk[j] = sum(alloc == j)
  }
  return(list(para.hat, propk / sum(propk)))
}

# if theta is NULL, allocate to each arms with equal probability


generate_GaussianRsp_M = function(para, ssn, mRate, k){
  outcome = matrix(data = NA, nrow = ssn, ncol = k)
  for(j in 1:k){
    outcome[,j] = rnorm(ssn, para[2*j - 1], sqrt(para[2*j]))
  }
  if(!is.null(mRate)){
    outcome = cbind(outcome, rbinom(ssn, 1, mRate))
  }
  return(outcome)
}

