###############################################################################
########################   Allocation functions   #############################
###############################################################################

# multi-arm ERADE (Hu, Zhang and He, 2009; Alkhnefr, Hu and Zhai, 2025, eq. (1))
# x: current allocation proportions, rho: estimated target allocation, a: degree of randomization
# over-allocated arms get a * rho, under-allocated arms share the remaining probability equally
erade.func = function(x, rho, a){
  x = as.numeric(x)
  rho = as.numeric(rho)
  over = x > rho
  under = x < rho
  prob = rho
  prob[over] = a * rho[over]
  if(any(under)){
    prob[under] = rho[under] + (1 - a) * sum(rho[over]) / sum(under)
  }
  return(prob)
}



# settings shared by the DBCD-family designs
rar.design = function(k, target.alloc, continuous, lower.bound, allocation, r, erade.alpha){
  targets = if(continuous){
    c("Neyman", "ZR", "DaOptimal", "OptimalNeyman")
  }else{
    c("Neyman", "RSIHR", "RPW", "WeisUrn", "OptimalNeyman", "OptimalRSIHR")
  }
  if(!(target.alloc %in% targets)){
    stop("'target.alloc' must be one of ", paste0("'", targets, "'", collapse = ", "))
  }
  if(!(allocation %in% c("DBCD", "ERADE"))){
    stop("'allocation' must be 'DBCD' or 'ERADE'")
  }
  if(r < 0){
    stop("'r' must be non-negative")
  }
  if((erade.alpha <= 0) | (erade.alpha >= 1)){
    stop("'erade.alpha' must be between 0 and 1")
  }
  if((lower.bound < 0) | (lower.bound > 1/k)){
    stop("'lower.bound' must be between 0 and 1/k")
  }
  return(list(k = k, target.alloc = target.alloc, continuous = continuous, lower.bound = lower.bound,
              allocation = allocation, r = r, erade.alpha = erade.alpha))
}



# allocation probabilities for the next patient (or group) given the current allocation
# proportions and parameter estimates
rar.prob = function(prop.k, theta.hat, dsg){
  # no usable allocation proportions yet (e.g. every patient so far has a missing response)
  prop.k = as.numeric(prop.k)
  if(any(!is.finite(prop.k))) prop.k = rep(1/dsg$k, dsg$k)
  if(dsg$continuous){
    rho = target.rho.Ctinuous(theta.hat, dsg$target.alloc, dsg$lower.bound)
  }else{
    rho = target.rho(theta.hat, dsg$target.alloc, dsg$lower.bound)
  }
  rho = safe.rho(rho, dsg$k)
  if(dsg$allocation == "DBCD"){
    prob = g.func(prop.k, rho, dsg$r)
  }else{
    prob = erade.func(prop.k, rho, dsg$erade.alpha)
  }
  return(prob)
}



# the sample size and the initial patients of the DBCD-family designs
check.size = function(n0, k, ssn){
  if(k < 2){
    stop("k must be at least 2")
  }
  if(n0 < k){
    stop("Number of initial participants 'n0' must be at least k")
  }
  if(ssn <= n0){
    stop("Total sample size 'ssn' must be larger than 'n0'")
  }
}



###############################################################################
############################   Hypothesis tests   #############################
###############################################################################

# two-sided test of equal means: t test for 2 arms, chi-squared test otherwise,
# or the p-value returned by a user-supplied test.fun(outcome, assignment)
rar.test = function(alpha, obs.outcome, assign.group, k, continuous = FALSE, test.fun = NULL){
  if(!is.null(test.fun)){
    p.value = test.fun(obs.outcome, assign.group)
    if(!is.numeric(p.value) | (length(p.value) != 1)){
      stop("'test.fun' must return a single p-value")
    }
    return(ifelse(p.value <= alpha, 1, 0))
  }
  if(k == 2){
    test.result = ttest.2(alpha, obs.outcome, assign.group)
  }else{
    test.result = chisq.test.k(alpha, obs.outcome, assign.group, k, continuous)
  }
  return(test.result)
}



# rejection rate under H0 from a second run of the same design in which all arms share
# the same parameter (null.args replaces the true parameters)
typeI.error = function(mc, fun, envir, null.args){
  mc[[1L]] = fun
  for(nm in names(null.args)){
    mc[[nm]] = null.args[[nm]]
  }
  mc$typeI = FALSE
  mc$seed = NULL
  return(eval(mc, envir)[["power"]])
}



# common mean under H0 for continuous designs, keeping the variances of each arm
null.theta = function(theta){
  theta[c(TRUE, FALSE)] = mean(theta[c(TRUE, FALSE)])
  return(theta)
}



###############################################################################
######################   Group sequential monitoring   ########################
###############################################################################

# two-sided alpha spending functions (Lan and DeMets, 1983): each tail spends the one-sided
# function at level alpha/2 (Proschan, Lan and Wittes, 2006), as in Zhu and Hu (2010)
spend.func = function(spend, alpha, t){
  if(spend == "OBF"){
    alpha.t = 2 * (2 * pnorm(qnorm(1 - alpha/4) / sqrt(t), lower.tail = FALSE))
  }else if(spend == "Pocock"){
    alpha.t = alpha * log(1 + (exp(1) - 1) * t)
  }else if(spend == "Linear"){
    alpha.t = alpha * t
  }else{
    stop("'spend' must be one of 'OBF', 'Pocock' or 'Linear'")
  }
  return(alpha.t)
}



sqMonitor = function(t, spend = "OBF"){
  t = sort(unique(c(t, 1)))
  if(any(t <= 0) | any(t > 1)){
    stop("Information times 't' must be in (0, 1]")
  }
  spend.func(spend, 0.05, t)
  monitor = list(t = t, spend = spend)
  class(monitor) = "sqMonitor"
  return(monitor)
}



# two-sided boundaries for the standardized statistic: the probability of first crossing
# the boundary at look j under H0 equals alpha(t_j) - alpha(t_{j-1}). The density of the
# B-value W(t) = Z(t) sqrt(t) on the continuation region is carried from look to look by
# numerical integration (Armitage, McPherson and Rowe, 1969)
sqBoundary = function(t, alpha = 0.05, spend = "OBF"){
  t = sort(unique(c(t, 1)))
  if(any(t <= 0) | any(t > 1)){
    stop("Information times 't' must be in (0, 1]")
  }
  alpha.inc = diff(c(0, spend.func(spend, alpha, t)))
  nlook = length(t)
  bound = rep(Inf, nlook)
  zmax = 40
  ngrid = 2001

  # first look: W(t_1) ~ N(0, t_1)
  if(alpha.inc[1] > 0) bound[1] = qnorm(alpha.inc[1]/2, lower.tail = FALSE)
  cj = min(bound[1], 12) * sqrt(t[1])
  w = seq(-cj, cj, length.out = ngrid)
  f = dnorm(w, sd = sqrt(t[1]))
  for(j in seq_len(nlook)[-1]){
    d = t[j] - t[j-1]
    # Simpson weights on the current grid
    h = w[2] - w[1]
    wt = h / 3 * c(1, rep(c(4, 2), (ngrid - 3) / 2), 4, 1)
    fw = wt * f
    if(alpha.inc[j] > 0){
      # probability of continuing so far and crossing |W(t_j)| >= c now
      cross = function(cc) sum(fw * (pnorm((-cc - w) / sqrt(d)) + pnorm((w - cc) / sqrt(d)))) - alpha.inc[j]
      cj = uniroot(cross, c(0, zmax * sqrt(t[j])), tol = 1e-12)$root
      bound[j] = cj / sqrt(t[j])
    }else{
      cj = 12 * sqrt(t[j])
    }
    if(j < nlook){
      w.new = seq(-min(cj, 12 * sqrt(t[j])), min(cj, 12 * sqrt(t[j])), length.out = ngrid)
      f = as.vector(dnorm(outer(w.new, w, "-"), sd = sqrt(d)) %*% fw)
      w = w.new
    }
  }
  names(bound) = paste0("look ", 1:nlook)
  return(bound)
}



# look sizes and boundaries of a monitored design
setup.monitor = function(monitor, k, ssn, n0, alpha, test.fun){
  if(is.null(monitor)) return(NULL)
  if(!inherits(monitor, "sqMonitor")){
    stop("'monitor' must be created by sqMonitor()")
  }
  if(k != 2){
    stop("Sequential monitoring is only available for k = 2")
  }
  if(!is.null(test.fun)){
    stop("'test.fun' cannot be combined with 'monitor'")
  }
  look.n = ceiling(monitor$t * ssn)
  if(any(duplicated(look.n))){
    stop("Two looks fall on the same sample size")
  }
  if(look.n[1] <= n0){
    stop("The first interim analysis must come after the n0 initial patients")
  }
  return(list(t = monitor$t, look.n = look.n, bound = sqBoundary(monitor$t, alpha, monitor$spend),
              stage = NULL, n.stop = NULL))
}



# Wald test at look j: 1 if the boundary is crossed, 0 otherwise (NA if the statistic
# cannot be computed at the final look)
sq.test = function(sq, j, obs.outcome, assign.group, k, continuous){
  stat = wald.stat(obs.outcome, assign.group, k, continuous)
  if(is.na(stat)){
    return(ifelse(j == length(sq$look.n), NA, 0))
  }
  return(ifelse(sqrt(stat) >= sq$bound[j], 1, 0))
}



# Wald test at look j of a group design, using all patients enrolled so far
# (responses are immediate; missing responses are excluded)
sq.look.group = function(sq, j, alloc, outcome, mRate, k, continuous){
  n = length(alloc)
  y = outcome[cbind(1:n, alloc)]
  eff = if(is.null(mRate)) 1:n else which(outcome[1:n, k+1] == 0)
  return(sq.test(sq, j, y[eff], alloc[eff], k, continuous))
}



###############################################################################
###################   Allocation for an ongoing trial   #######################
###############################################################################

nextAlloc = function(alloc, outcome, k, response = "binary", target.alloc = NULL,
                     allocation = "DBCD", r = 2, erade.alpha = 0.5, lower.bound = 0,
                     theta0 = NULL, size = 1){
  if(!(response %in% c("binary", "continuous"))){
    stop("'response' must be 'binary' or 'continuous'")
  }
  if(length(alloc) != length(outcome)){
    stop("'alloc' and 'outcome' must have the same length")
  }
  if(any(!(alloc %in% 1:k))){
    stop("'alloc' must contain arm labels 1, ..., k")
  }
  continuous = (response == "continuous")
  if(is.null(target.alloc)){
    target.alloc = ifelse(continuous, "Neyman", "RPW")
  }
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }
  dsg = rar.design(k, target.alloc, continuous, lower.bound, allocation, r, erade.alpha)

  # estimates from the observed responses (NA = not observed yet or missing)
  obs = !is.na(outcome)
  theta.hat = NULL
  for(j in 1:k){
    y = outcome[obs & (alloc == j)]
    if(continuous){
      theta.hat[2*j - 1] = mean(y)
      theta.hat[2*j] = ifelse(length(y) >= 2, mean((y - mean(y))^2), NA)
    }else{
      theta.hat[j] = (sum(y) + theta0[j]) / (length(y) + 1)
    }
  }
  if(continuous){
    rho = target.rho.Ctinuous(theta.hat, target.alloc, lower.bound)
  }else{
    rho = target.rho(theta.hat, target.alloc, lower.bound)
  }
  rho = safe.rho(rho, k)
  prop.k = tabulate(alloc, nbins = k) / max(length(alloc), 1)
  if(length(alloc) == 0){
    prob = rep(1/k, k)
  }else{
    prob = rar.prob(prop.k, theta.hat, dsg)
    prob = prob / sum(prob)
  }
  trt = paste0("treatment ", LETTERS[1:k])
  names(prob) = trt
  names(rho) = trt
  if(continuous){
    names(theta.hat) = as.vector(rbind(paste0("mu", LETTERS[1:k]), paste0("sigma2", LETTERS[1:k])))
  }else{
    names(theta.hat) = paste0("p", LETTERS[1:k])
  }
  assignment = if(size > 0) sample.int(k, size = size, replace = TRUE, prob = prob) else integer(0)
  return(list(prob = prob, target = rho, estimate = theta.hat, assignment = assignment))
}



###############################################################################
###########################   print and summary   #############################
###############################################################################

summary.grouprar = function(object, ...){
  continuous = identical(attr(object, "response"), "continuous")
  arms = data.frame(allocation = as.numeric(object[["propotion"]]),
                    sd = as.numeric(object[["sd of propotion"]]),
                    row.names = names(object[["propotion"]]))
  if(continuous){
    arms = cbind(data.frame(mu = object[["parameter"]][c(TRUE, FALSE)],
                            sigma2 = object[["parameter"]][c(FALSE, TRUE)]), arms)
  }else{
    arms = cbind(data.frame(p = object[["parameter"]]), arms)
  }
  rownames(arms) = names(object[["propotion"]])
  stats = c(object[["failure rate"]], object[["sd of failure rate"]], object[["power"]])
  names(stats) = c(ifelse(continuous, "mean response", "failure rate"),
                   ifelse(continuous, "sd of mean response", "sd of failure rate"), "power")
  for(nm in c("type I error", "duration", "sd of duration", "enrollment duration", "expected sample size")){
    if(!is.null(object[[nm]])) stats[nm] = object[[nm]]
  }
  out = list(method = object[["method"]], "sample size" = object[["sample size"]],
             nsim = nrow(object[["data: propotion"]]), arms = arms, stats = stats,
             "stopping probability" = object[["stopping probability"]])
  class(out) = "summary.grouprar"
  return(out)
}



print.summary.grouprar = function(x, digits = 3, ...){
  cat(x$method, "\n")
  ssn = x[["sample size"]]
  cat("Sample size:", paste0(names(ssn), ifelse(is.null(names(ssn)), "", " "), ssn, collapse = ", "),
      "  Simulations:", x$nsim, "\n\n")
  print(round(x$arms, digits))
  cat("\n")
  print(round(x$stats, digits))
  if(!is.null(x[["stopping probability"]])){
    cat("\nStopping probability:\n")
    print(round(x[["stopping probability"]], digits))
  }
  invisible(x)
}



print.grouprar = function(x, ...){
  print(summary(x), ...)
  invisible(x)
}
