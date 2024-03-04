library(stringr)


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
  if((n1 == 0) | (n2 == 0)){warning("All subjects were assigned to the same group")}
  t.stats =  (mean(x1) - mean(x2)) / sqrt(var(x1) / n1 + var(x2) / n2)
  df = (var(x1) / n1 + var(x2) / n2)^2 / (var(x1)^2 / (n1^2 * (n1 - 1)) + var(x2)^2 / (n2^2 * (n2 - 1)))
  # p.value = 2 * (1 - pt(t.stats, df, lower.tail = FALSE))
  p.value = 2 * min(1 - pt(t.stats, df, lower.tail = FALSE), pt(t.stats, df, lower.tail = FALSE))
  test.result = ifelse(p.value <= alpha, 1, 0)
  return(test.result)
}



# hypothesis test (chi-squared test)
chisq.test.k = function(alpha, obs.outcome, assign.group, k){
  hat.p = aggregate(obs.outcome, by = list(type = assign.group), mean)$x
  hat.pc = hat.p[-k] - hat.p[k]
  n = as.vector(table(assign.group))

  Sigma = diag(hat.p[-k] * (1-hat.p[-k]) / n[-k]) + hat.p[k] * (1-hat.p[k]) / n[k]
  if(det(Sigma) == 0){
    warning("variance-covariance matrix is exactly singular")
    test.result = NA
  }else{
    inverse.Sigma = solve(Sigma)
    chisq.stats = t(hat.pc) %*% inverse.Sigma %*% t(t(hat.pc))
    df = k-1
    p.value = pchisq(chisq.stats, df, lower.tail=FALSE)
    test.result = ifelse(p.value <= alpha, 1, 0)
    test.result
  }
  return(test.result)
}



# DBCD common
# establish the target distribution
# discrete
target.rho = function(para.set, target.alc){
  if((target.alc != "Neyman") & (target.alc != "RSIHR") & (target.alc != "RPW") & (target.alc != "WeisUrn")){
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
  }
  return(rho)
}



target.rho.Ctinuous = function(para.set, target.alc){
  n = length(para.set)
  mu = NULL
  sigma2 = NULL

  for(i in 1:(n/2)){
    mu[i] =  para.set[(2*i-1)]
    sigma2[i] =  para.set[(2*i)]
  }

  if((target.alc != "Neyman") & (target.alc != "ZR") & (target.alc != "DaOptimal")){
    stop("Your target allocation is not in the list")
  }else if(target.alc == "Neyman"){
    rho = sqrt(sigma2) / sum(sqrt(sigma2))
  }else if(target.alc == "ZR"){
    temp1 = sqrt(mu)
    temp2 = sqrt(sigma2)
    rho = c(temp1[2] * temp2[1], temp1[1] * temp2[2]) / (temp1[2] * temp2[1] + temp1[1] * temp2[2])
  }else if(target.alc == "DaOptimal"){
    temp =( sqrt(sigma2) )^(4/3)
    rho = (temp) / sum(temp)
  }
  return(rho)
}



g.func = function(x, y, r){
  if((x != 0) & (x != 1)){
    tmp1 = y * (y / x) ^ r
    tmp2 = (1-y) * ((1-y) / (1-x)) ^ r
    g.xy =  tmp1 / (tmp1 + tmp2)
  }else if(x == 0){
    g.xy = 1
  }else if(x == 1){
    g.xy = 0
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



calc_thetaGaussian = function(p.hat, k, alloc, outcome, theta0){
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }
  for(j in 1:k){
    mu = mean(outcome[which(alloc == j), j])
    sigma2 = mean((outcome[which(alloc == j), j] - mu)^2)
    p.hat[(2*j - 1)] = mu
    p.hat[(2*j)] = sigma2
  }
  return(p.hat)
}



# Delayed
responseDist = function(dist, param, k, level, sample.size){

  if(!(dist %in% c("exponential", "normal", "unifrom"))){
    stop("'dist' must be one of the element of 'exponential', 'normal' or 'unifrom'")
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



# monitor
sqMonitor = function(theta, assign.group, t, n, i, alpha.t, n0, N.set){
  ith = i+n0
  if(sum(ith == floor(n*t)) == 1){
    index = which(ith == floor(n*t))
    test.result = brownian.test(theta, assign.group, tk = t[index], alphak = alpha.t[index])
    if(test.result == 1){
      N.set[index] = N.set[index] + 1
      N.set[4] = 1 #stop sign
    }
  }
  return(N.set)
}


# spend function
spend.func = function(spend, alpha, t){
  if(spend == "BFlike"){

  }else if(spend == "Linear"){
    alpha.t = alpha * t
  }else if(spend == "Pococklike"){
    alpha.t = alpha * (log(1 + exp(1) - 1) * t)
  }
  return(alpha.t)
}


# brownian.test two.sided
brownian.test = function(theta, assign.group, tk, alphak){
  diff.mean = theta[1] - theta[2]
  var.1 = theta[1] * (1-theta[1]) / sum(assign.group == 1)
  var.2 = theta[2] * (1-theta[2]) / sum(assign.group == 2)
  zt = diff.mean / sqrt(var.1 + var.2)
  bt = sqrt(tk) * zt
  b = qnorm(1-alphak/2, mean = 0, sd = sqrt(tk))
  test.result = ifelse(abs(bt) > b, 1, 0)
  return(test.result)
}


library(stringr)
RAR_Output = function(name, parameter, ssn,
                      assignment, propotion,
                      failRate,
                      pwCalc, k){

  # statistics
  Mean.propotion = apply(propotion, 2, mean)
  sdprop = sd(propotion[,1])
  #sdprop = sd(sim.prop)
  Mean.failRate = mean(failRate)
  sdfrt = sd(failRate)
  power = mean(na.omit(pwCalc))

  # parameter
  paraNum = stringr::str_c('p',  LETTERS[1:k])
  names(parameter) = paraNum

  # proportion
  trt = stringr::str_c("treatment ", LETTERS[1:k])
  propotion = as.data.frame(propotion)
  colnames(propotion) = trt

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
                    "data: propotion" = propotion)

  # cat(stringr::str_c(trt, round(Mean.propotion, 3), sep = ": ", collapse = "   "), "\n",
  #     sprintf('SD of allocation = %g', round(sdprop, 3)), "\n","\n")
  # cat(sprintf('Failure Rate = %g', round(Mean.failRate, 3)), "\n",
  #     sprintf('SD of Failure Rate = %g', round(sdfrt, 3)), "\n", "\n")
  # cat(sprintf('Power = %g', round(power, 3)), "\n")
  return(invisible(outputList))
}

# RAR_Output_Cont = function(name, parameter, ssn,
#                       assignment, propotion,
#                       failRate,
#                       pwCalc, k){
#
#   # statistics
#   Mean.propotion = apply(propotion, 2, mean)
#   sdprop = sd(propotion[,1])
#   #sdprop = sd(sim.prop)
#   Mean.failRate = mean(failRate)
#   sdfrt = sd(failRate)
#   power = mean(na.omit(pwCalc))
#
#   # parameter
#   paraNum = stringr::str_c('p',  LETTERS[1:k])
#   names(parameter) = paraNum
#
#   # proportion
#   trt = stringr::str_c("treatment ", LETTERS[1:k])
#   propotion = as.data.frame(propotion)
#   colnames(propotion) = trt
#
#   # output
#   outputList = list(method = name,
#                     "sample size" = ssn,
#                     "parameter" = parameter,
#                     "propotion" = Mean.propotion,
#                     "sd of propotion" = sdprop,
#                     "failure rate" = Mean.failRate,
#                     "sd of failure rate" = sdfrt,
#                     'power' = power,
#                     "data: failureRate" = failRate,
#                     "data: test" = pwCalc,
#                     "data: assignment" = assignment,
#                     "data: propotion" = propotion)
#
#   cat(stringr::str_c(trt, round(Mean.propotion, 3), sep = ": ", collapse = "   "), "\n",
#       sprintf('SD of allocation = %g', round(sdprop, 3)), "\n","\n")
#   cat(sprintf('Failure Rate = %g', round(Mean.failRate, 3)), "\n",
#       sprintf('SD of Failure Rate = %g', round(sdfrt, 3)), "\n", "\n")
#   cat(sprintf('Power = %g', round(power, 3)), "\n")
#   return(invisible(outputList))
# }

# checkError = function(){
#   if(!is.integer(n0/k)){
#     stop("n0 has to be divisible by k")
#   }
#   if(is.null(target.alloc)){
#     target.alloc = "RSIHR"
#   }
#   if(is.null(dist)){
#     dist = "normal"
#   }
#   if(is.null(r)){
#     r = 2
#   }
# }


# checkError = function(n0, k){
#   if(!is.integer(n0/k)){
#     stop("n0 has to be divisible by k")
#   }
#   if(is.null(target.alloc)){
#     target.alloc = "RSIHR"
#   }
#   if(is.null(dist)){
#     dist = "normal"
#   }
#   if(is.null(r)){
#     r = 2
#   }
# }


# compPcdure = function(boxLabel, ylab, ...) {
#   data = data.frame(list(...))
#   colnames(data) = 1:ncol(data)
#   n = ncol(data)
#   alldat = tidyr::gather(data, key = 'Varieties', value = "Yield")
#   ggplot(alldat, aes(x = Varieties, y = Yield, fill = Varieties)) +
#     geom_boxplot(
#       outlier.colour = "black",
#       outlier.shape = 16,
#       outlier.size = 2,
#       notch = F
#     ) +
#     stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),
#                  width=0.2)+
#     stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),
#                  width=0.2)+
#     scale_fill_brewer(palette="Pastel2") +
#     labs(x = "Methods", y = "Data") +
#     scale_x_discrete(labels = boxLabel) +
#     theme(legend.position = "none")
# }


calc_prop = function(alloc, outcome, mRate, k){
  prop.k = NULL
  if(is.null(mRate)){
    for(j in 1:k){
      prop.k[j] = sum(alloc == j) / length(alloc)
    }
  }else{
    for(j in 1:k){
      prop.k[j] = sum(alloc[which(outcome[, k+1] == 0)] == j, na.rm = T) / sum(table(alloc[which(outcome[, k+1] == 0)]))
    }
  }
  return(prop.k)
}


calc_RspT = function(alloc, outcome, rspT, gsize, entry, adjust = F, continuous = T){
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
    if(continuous == T){
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

