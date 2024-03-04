###############################################################################
####################   Complete Randomization   ###############################
###############################################################################
#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha An integer between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#' @import stats
#' @description Simulating complete randomization with two-sided hypothesis testing in a clinical trial context.
#' @details Complete randomization: Allocating participants or subjects to different treatment groups in a clinical trial in such a way that each participant has an equal and independent chance of being assigned to any of the treatment groups.
#' @export
#' @return \item{name}{The name of procedure.}
#' @return \item{parameter}{The true parameters used to do the simulations.}
#' @return \item{assignment}{The randomization sequence.}
#' @return \item{propotion}{Average allocation porpotion for each of treatment groups.}
#' @return \item{failRate}{The proportion of individuals who do not achieve the expected outcome in each simulation, on average.}
#' @return \item{pwClac}{The probability of the study to detect a significant difference or effect if it truly exists.}
#' @return \item{k}{Number of arms involved in the trial.}
#' @examples CRDesign(k=3, p = c(0.7, 0.8, 0.6), nsim = 500, ssn = 400)
CRDesign = function(k, p, ssn, nsim = 2000, alpha = 0.05){
  if(sum((p < 0) | (p > 1)) != 0){ stop("p must be a positive number between 0 and 1!") }
  if(length(p) != k){ stop("Length of parameter vector p must equal to k.") }

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

  for(s in 1:nsim){
    obs.outcome = NULL
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # complete randomization
    alloc = sample(1:k, ssn, replace = TRUE, prob = rep(1/k, k))
    for(j in 1:k){
      obs.outcome[which(alloc == j)] = outcome[which(alloc == j), j]
    }
    if(k == 2){
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group = alloc)
    }else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group = alloc, k)
    }
    failure.rate[s] = mean(obs.outcome == 0)
    group.prop = rbind(group.prop, table(alloc) / length(alloc))
  }
  name = "Complete Randomization"
  return(RAR_Output(name, parameter=p, ssn,
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}



###############################################################################
####################   Randomized Play-the-winner rule   ######################
###############################################################################

#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k = 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param Y0 A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1,)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha An integer between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @description Simulating randomized play-the-winner rule with two-sided hypothesis testing in a clinical trial context.
#' @details The Randomized Play-the-Winner Rule allocates future subjects in a clinical trial to treatment groups based on the performance of previously treated subjects. This rule increases the likelihood of future patients being assigned to the better-performing treatment, as determined by the outcomes of previously treated subjects.
#' @export
#'
#' @examples RPWRule(k = 2, p = c(0.7, 0.8), ssn = 400, Y0 = NULL, nsim = 2000, alpha = 0.05)
RPWRule = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05){

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
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
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
    group.prop = rbind(group.prop, table(assign.group) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = ttest.2(alpha, obs.outcome, assign.group)
  }
  name = "Randomized Play-the-winner Rule"
  return(RAR_Output(name, parameter=p, ssn,
                    assignment = assign.group, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}

###############################################################################
############################   Wei's Urn Model   ##############################
###############################################################################

#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k > 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param Y0 A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples WeiUrn(k = 3, p = c(0.7, 0.8, 0.7), ssn = 400, Y0 = NULL, nsim = 2000, alpha = 0.05)
WeiUrn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05){
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

  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()
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
    group.prop = rbind(group.prop, table(assign.group) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    pwCalc[s] = chisq.test.k(alpha, obs.outcome, assign.group, k)
  }
  name = "Wei's Urn"
  return(RAR_Output(name, parameter=p, ssn,
                    assignment = assign.group, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}


###############################################################################
############################   Polya Urn Model   ##############################
###############################################################################

#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param Y0 A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#' @import stringr
#' @export
#'
#' @examples PolyaUrn(k = 3, p = c(0.6, 0.7, 0.5), ssn = 400, Y0 = NULL, nsim = 500, alpha = 0.05)
PolyaUrn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05){
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
  #group.prop = c()
  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL
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
    if(length(unique(assign.group)) < k){
      countA = rep(NA, k)
      for(i in 1:k){countA[i] = sum(assign.group == i) / ssn}
      group.prop = rbind(group.prop, countA)
      pwCalc[s] = NA}
    else{
      prop = table(assign.group) / ssn
      group.prop = rbind(group.prop, table(assign.group) / ssn)
      pwCalc[s] = chisq.test.k(alpha, obs.outcome, assign.group, k)
    }
  }
  name = "Polya Urn"
  return(RAR_Output(name, parameter=p, ssn,
                    assignment = assign.group, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}

###############################################################################
#########################  Drop-the-loser Rule   ##############################
###############################################################################
#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param Y0 A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples DLRule(k = 2, p = c(0.7, 0.8), ssn = 400, Y0 = NULL, nsim = 2000, alpha = 0.05)
DLRule = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05){
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
  # add a immigration ball
  Y0 = c(Y0, 1)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL

  for(s in 1:nsim){
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # assigned randomly
    Y = Y0
    sample.prob = Y / sum(Y)
    assign.group = NULL
    obs.outcome = NULL

    for(i in 1:ssn){
      assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      while(assign.k == (k+1)){
        Y[-assign.k] = Y[-assign.k] + 1
        # sample probability changed here
        sample.prob = Y / sum(Y)
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
      sample.prob = Y / sum(Y)
    }
    group.prop = rbind(group.prop, table(assign.group) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    if(k == 2){
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group)
    }else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group, k)
    }
  }
  name = "Drop-the-loser Rule"
  return(RAR_Output(name, parameter = p, ssn,
                    assignment = assign.group, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}


###############################################################################
#####################   Generalized Drop-the-loser Rule   #####################
###############################################################################

#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param aK a positive vector of length equals to \code{k}. The values specifies when the immigration ball is drawn, the number of each treatment ball added to the urn.
#' @param Y0 A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples GDLRule(k = 3, p = c(0.6, 0.7, 0.6),
#'                   ssn = 400, aK = c(1, 1, 1), Y0 = NULL, nsim = 2000, alpha = 0.05)
GDLRule = function(k, p, ssn, aK, Y0 = NULL, nsim = 2000, alpha = 0.05){
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
  ## check Y0
  if(length(Y0) != k && !is.null(Y0)){
    stop("Length of Y0 must be equal to k")
  }else if(is.null(Y0)){
    # use default Y0
    Y0 = rep(1, k)
  }
  # add a immigration ball
  Y0 = c(Y0, 1)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL

  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # assigned randomly
    Y = Y0
    sample.prob = Y / sum(Y)
    assign.group = NULL
    obs.outcome = NULL

    for(i in 1:ssn){
      assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      while(assign.k == (k+1)){
        Y[-assign.k] = Y[-assign.k] + aK #add ak here
        # sample probability changed here
        sample.prob = Y / sum(Y)
        assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      }
      assign.group[i] = assign.k
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]]
      }else{
        Y[assign.group[i]] = Y[assign.group[i]] - 1
      }
      sample.prob = Y / sum(Y)
    }
    group.prop = rbind(group.prop, table(assign.group) / ssn)
    failure.rate = mean(obs.outcome == 0)
    #pwCalc[s] = chisq.test.k(alpha, obs.outcome, assign.group, k)
    if(k == 2){
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group)
    }else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group, k)
    }
  }
  name = "Generalized Drop-the-loser Rule"
  return(RAR_Output(name, parameter = p, ssn,
                    assignment = assign.group, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}

###############################################################################
####################   Bai, Hu, Shen's Urn Model   ############################
###############################################################################

#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k > 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param Y0 A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples Bai.Hu.Shen.Urn(k = 3, p = c(0.7, 0.8, 0.6), ssn = 500, Y0 = NULL, nsim = 2000, alpha = 0.05)

Bai.Hu.Shen.Urn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05){
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
  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL
  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
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
        M = sum(p)
        Y[-assign.group[i]] = Y[-assign.group[i]] + p[-assign.group[i]] / (M-p[assign.group[i]])
      }
      sample.prob = Y / sum(Y)
    }
    group.prop = rbind(group.prop, table(assign.group) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    if (k == 2) {
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group)
    } else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group, k)
    }
    #pwCalc[s] = chisq.test.k(alpha, obs.outcome, assign.group, k)
  }
  name = "BaiHuShen's Urn"
  return(RAR_Output(name, parameter = p, ssn,
                    assignment = assign.group, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}

###############################################################################
####################   Birth and Death Urn Model   ############################
###############################################################################

#' Title
#'
#' @param k a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param p a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param ssn a positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param Y0 A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples
#' BirthDeathUrn(k = 3, p = c(0.6, 0.7, 0.6), ssn = 400, Y0 = NULL, nsim = 2000, alpha = 0.05)
BirthDeathUrn = function(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05){

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
  # add a immigration ball
  Y0 = c(Y0, 1)

  group.prop = c()
  failure.rate = NULL
  pwCalc = NULL

  for(s in 1:nsim){
    # outcome matrix
    outcome = generate_data_M(p, ssn, mRate=NULL, k)
    # assigned randomly
    Y = Y0
    sample.prob = Y / sum(Y)
    assign.group = NULL
    obs.outcome = NULL

    for(i in 1:ssn){
      assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      while(assign.k == (k+1)){
        Y[-assign.k] = Y[-assign.k] + 1
        # sample probability changed here
        sample.prob = Y / sum(Y)
        assign.k = sample(c(1:(k+1)), 1, prob = sample.prob)
      }
      assign.group[i] = assign.k
      obs.outcome[i] = outcome[i, assign.group[i]]
      if(obs.outcome[i] == 1){
        Y[assign.group[i]] = Y[assign.group[i]] + 1
      }else{
        Y[assign.group[i]] = Y[assign.group[i]] - 1
      }
      sample.prob = Y / sum(Y)
    }
    group.prop = rbind(group.prop, table(assign.group) / ssn)
    failure.rate[s] = mean(obs.outcome == 0)
    if (k == 2) {
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group)
    } else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group, k)
    }
    #pwCalc[s] = chisq.test.k(alpha, obs.outcome, assign.group, k)
  }
  name = "Birth and Death's Urn"
  return(RAR_Output(name, parameter = p, ssn,
                    assignment = assign.group, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}

###############################################################################
##########   Hu and Zhang's Doubly biased coin Design (Binary)  ###############
###############################################################################

#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.
#' @param p A positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param theta0 A vector of length k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length k, with all values set to 0.5.
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}. The default is \code{"RPW"}.
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param mRate a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples DBCD_Bin(n0 = 20, p = c(0.7, 0.8), k = 2, ssn = 300, theta0 = NULL, target.alloc = "RPW", r = 2, nsim = 80, mRate = NULL, alpha = 0.05)
DBCD_Bin = function(n0 = 20, p, k, ssn, theta0 = NULL, target.alloc = "RPW", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05){

  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

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

    p.hat = calc_theta(p.hat, k, alloc.n0, outcome, theta0)
    for(j in 1:k){
      obs.outcome[which(alloc.n0 == j)] = outcome[which(alloc.n0 == j), j]
    }
    extra.n = ssn - n0
    alloc = alloc.n0

    for(i in 1:extra.n){
      prop.k = table(alloc) / length(alloc)
      p.hat = calc_theta(p.hat, k, alloc, outcome, theta0)
      est.rho = target.rho(p.hat, target.alc = target.alloc)
      for(j in 1:k){
        sample.prob[j] = g.func(prop.k[j], est.rho[j], r)
      }
      assign.group = sample(x = c(1:k), 1, prob = sample.prob)
      alloc = c(alloc, assign.group)
      obs.outcome = c(obs.outcome, outcome[(n0+i), assign.group])
    }

    pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group = alloc)
    failure.rate[s] = mean(obs.outcome == 0)
    group.prop = rbind(group.prop, table(alloc) / length(alloc))
  }
  name = "Hu and Zhang's DCBD"

  return(RAR_Output(name, parameter=p, ssn,
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}


###############################################################################
########   Hu and Zhang's Doubly biased coin Design (delayed+bin)  ############
###############################################################################

#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.
#' @param p A positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param ent.param A positive integer. The value specified the parameter for an expoential distribution which determine the time for each participant enter the trial.
#' @param rspT.dist Distribution Type. Specifies the type of distribution that models the time spent for the availability of patient \eqn{i} under treatment \eqn{k}. Acceptable options for this argument include: \code{"exponential"}, \code{"normal"}, and \code{"uniform"}.
#' @param rspT.param A vector. Specifies the parameters required by the distribution that models the time spent for the availability under each treatment. (eg. If there are 3 treatments groups and each of them follows truncated normal distribution with parameter pair (3, 2), (2, 1), (4, 1), repectively. Then the \code{rspT.param = c(3, 2, 2, 1, 4, 1)})
#' @param theta0 A vector of length k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length k, with all values set to 0.5.
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}. The default is \code{"RPW"}.
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param mRate a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples p = c(0.7, 0.8)
#' k = 2
#' ssn = 300
#' ### for enter time and response time simulation
#' ent.param = 0.7
#' rspT.dist = "exponential"
#' rspT.param = c(1, 1, 3, 1)
#' ## Arguments for the deisgn
#' n0 = 20
#' target.alloc = "RSIHR"
#' dyldDBCD_Bin(n0 = n0, p = p, k = k, ssn = ssn, ent.param, rspT.dist, rspT.param, theta0 = NULL, target.alloc, r = 2, nsim = 500, mRate = NULL, alpha = 0.05)
dyldDBCD_Bin = function(n0 = 20, p, k, ssn, ent.param, rspT.dist, rspT.param, theta0 = NULL,
                        target.alloc = "RPW", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05) {

  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }
  if(is.null(theta0)){
    theta0 = rep(0.5, k)
  }

  # setup
  sim.prop = NULL
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

  for (s in 1:nsim) {
    obs.outcome = NULL
    alloc.n0 = rep(NA, n0)
    outcome = generate_data(p, ssn)

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
    extra.n = ssn - n0
    alloc = alloc.n0

    obs.cumRspT = NULL
    # All observed cumulative time so far
    rsp.idx =  obs.outcome[1:n0] + 1
    rspT.idx = 2 * (alloc.n0 - 1) + rsp.idx
    for (idx in 1:n0) {
      obs.cumRspT[idx] = cumRspT[idx,  rspT.idx[idx]]
    }
    # estimate p
    temp = as.data.frame(cbind(alloc[which(entryT[n0 + 1] > obs.cumRspT)],
                               obs.outcome[which(entryT[n0 + 1] > obs.cumRspT)]))
    p.hat = NULL
    for (t in 1:k) {
      p.hat[t] = (sum(temp[which(temp[, 1] == t), 2]) + theta0[t]) / (nrow(temp[which(temp[, 1] == t), ]) + 1)
    }
    # calculate sample prob
    est.rho = target.rho(p.hat, target.alc = target.alloc)
    prop.k = table(alloc) / length(alloc)
    for (j in 1:k) {
      sample.prob[j] = g.func(prop.k[j], est.rho[j], r)
    }
    m = NULL
    for (i in (n0 + 1):ssn) {
      alloc[i] = sample(c(1:k), 1, prob = sample.prob)
      obs.outcome[i] = outcome[i, alloc[i]]

      # detect is there anyone's outcome ready
      # time of m-th patient' outcome is avaiable
      asgn.idx = alloc[i]
      rsp.idx =  obs.outcome[i] + 1
      rspT.idx = 2 * (asgn.idx - 1) + rsp.idx
      obs.cumRspT[i] = cumRspT[i, rspT.idx]

      # the outcome could be observed after the m-th entry
      m[i - n0] = max(which(entryT <= obs.cumRspT[i]))

      # whether adjust the probability
      if (sum(m == i) == 0) {
        next
      } else{
        temp = rbind(temp, cbind(alloc[which(m == i)], obs.outcome[which(m == i)]))
        for (t in 1:k) {
          p.hat[t] = (sum(temp[which(temp[, 1] == t), 2]) + theta0[t]) / (nrow(temp[which(temp[, 1] == t), ]) + 1)
        }
        # calculate sample prob
        est.rho = target.rho(p.hat, target.alloc)
        prop.k = table(alloc) / length(alloc)
        for (j in 1:k) {
          sample.prob[j] = g.func(prop.k[j], est.rho[j], r)
        }
      }
    }
    # re-adjust obs.outcome and assign.group
    if (!is.null(mRate)) {
      obs.outcome = obs.outcome[-which(MD == 1)]
      alloc = alloc[-which(MD == 1)]
    }

    sim.prop[s] = table(alloc)[1] / ssn
    if (k == 2) {
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group = alloc)
    } else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group = alloc, k)
    }
    failure.rate[s] = mean(obs.outcome == 0)
    group.prop = rbind(group.prop, table(alloc) / length(alloc))
  }
  name = ifelse(is.null(mRate),
                "Delayed DBCD without Missing Data",
                "Delayed DBCD with Missing Data")
  return(RAR_Output(name, parameter=p, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}




###############################################################################
########   Hu and Zhang's Doubly biased coin Design (Continuous)  #############
###############################################################################

#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.
#' @param theta A numerical vector of length equal to \code{2k}. These values specify the true parameters for each treatment and are used for generating data in simulations. For example, if \code{k=2}, you should provide two pairs of parameter values, each consisting of the mean and variance, like: \code{theta = c(13, 4.0^2, 15, 2.5^2)}.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param theta0 A vector of length 2k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length 2k, with all parameter pair setted to be (0, 1).
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"ZR"}, \code{"DaOptimal"}. The default is \code{"Neyman"}. The details see Zhang L. and Rosenberger. W (2006).
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples
#' # A simple use!
#' # Define the arguments
#' ## Arguments for generate the simulated data
#' theta = c(13, 4.0^2, 15, 2.5^2)
#' k = 2
#' ssn = 88
#' ### Other arguments
#' target.alloc = "Neyman"
#'
#' DBCD_Cont(n0 = 20, theta, k, ssn, theta0 = NULL, target.alloc = "Neyman", r = 2, nsim = 500, alpha = 0.05)

DBCD_Cont = function(n0 = 20, theta, k, ssn, theta0 = NULL, target.alloc = "Neyman", r = 2, nsim = 2000, alpha = 0.05){

  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(F, T)] < 0) > 0){
    stop("The vairance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

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

    theta.hat = calc_thetaGaussian(theta.hat, k, alloc.n0, outcome, theta0)

    for(j in 1:k){
      obs.outcome[which(alloc.n0 == j)] = outcome[which(alloc.n0 == j), j]
    }

    extra.n = ssn - n0
    alloc = alloc.n0

    for(i in 1:extra.n){
      prop.k = table(alloc) / length(alloc)
      theta.hat = calc_thetaGaussian(theta.hat, k, alloc, outcome, theta0)
      est.rho = target.rho.Ctinuous(theta.hat, target.alc = target.alloc)
      for(j in 1:k){
        sample.prob[j] = g.func(prop.k[j], est.rho[j], r)
      }
      assign.group = sample(x = c(1:k), 1, prob = sample.prob)
      alloc = c(alloc, assign.group)
      obs.outcome = c(obs.outcome, outcome[(n0+i), assign.group])
    }
    pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group = alloc)
    #failure.rate[s] = mean(obs.outcome == 0)
    failure.rate[s] = mean(obs.outcome)
    group.prop = rbind(group.prop, table(alloc) / length(alloc))
  }
  name = "Hu and Zhang's DCBD (Gaussian Response)"

  return(RAR_Output(name, parameter = theta, ssn,
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc,k))
}



###############################################################################
########   Hu and Zhang's Doubly biased coin Design (delayed+cont)  ###########
###############################################################################

#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.
#' @param theta A numerical vector of length equal to \code{2k}. These values specify the true parameters for each treatment and are used for generating data in simulations. For example, if \code{k=2}, you should provide two pairs of parameter values, each consisting of the mean and variance, like: \code{theta = c(13, 4.0^2, 15, 2.5^2)}.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param ent.param A positive integer. The value specified the parameter for an expoential distribution which determine the time for each participant enter the trial.
#' @param rspT.dist Distribution Type. Specifies the type of distribution that models the time spent for the availability of patient \eqn{i} under treatment \eqn{k}. Acceptable options for this argument include: \code{"exponential"}, \code{"normal"}, and \code{"uniform"}.
#' @param rspT.param A vector. Specifies the parameters required by the distribution that models the time spent for the availability under each treatment. (eg. If there are 3 treatments groups and each of them follows truncated normal distribution with parameter pair (3, 2), (2, 1), (4, 1), repectively. Then the \code{rspT.param = c(3, 2, 2, 1, 4, 1)})
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"ZR"}, \code{"DaOptimal"}. The default is \code{"Neyman"}. The details see Zhang L. and Rosenberger. W (2006).
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param mRate a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#' @export
#'
#' @examples
#' # a simple use
#' # Define the arguments
#' ## Arguments for generate the simulated data
#' ### For response simulation
#' theta = c(13, 4.0^2, 15, 2.5^2)
#' k = 2
#' ssn = 88
#'
#' ### for enter time and response time simulation
#' ent.param = 5
#' rspT.param = rep(10, 2)
#' rspT.dist = "exponential"

#' ## Arguments for the deisgn
#' n0 = 10
#' target.alloc = "Neyman"

#' dyldDBCD_Cont(n0 = 10, theta, k, ssn, ent.param, rspT.dist, rspT.param, target.alloc, r = 2, nsim = 500, mRate = 0.2, alpha = 0.05)
dyldDBCD_Cont = function(n0 = 20, theta, k, ssn, ent.param, rspT.dist, rspT.param,
                         target.alloc = "Neyman", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05){

  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(F, T)] < 0) > 0){
    stop("The vairance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }

  # setup
  pwCalc = NULL
  failure.rate = NULL
  theta.hat.set = c()
  group.prop = c()

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
      cumRspT[which(outcome[k+1] == 1),] = Inf
    }

    # initial allocation rule is restricted randomization
    x = rep(1:k, ceiling(n0 / k))
    alloc.n0 = sample(x = x)[1:n0]

    # initial estimate for p
    sample.prob = NULL
    for(j in 1:k){
      obs.outcome[which(alloc.n0 == j)] = outcome[which(alloc.n0 == j), j]
    }
    extra.n = ssn - n0
    alloc = alloc.n0

    obs.cumRspT = NULL
    # All observed cumulative time so far
    for(idx in 1:n0){
      obs.cumRspT[idx] = cumRspT[idx,  alloc.n0[idx]]
    }

    # estimate p
    theta.hat = NULL
    temp = calc_thetaGaussian_MD(theta.hat, k = 2, alloc, outcome, entryT = entryT[n0+1], obsRspT = obs.cumRspT, mRate = mRate)
    theta.hat = temp[[1]]
    prop.k = temp[[2]]

    # calc sample prob
    est.rho = target.rho.Ctinuous(theta.hat, target.alc = target.alloc) #

    # # if no response came out, use 1/k as probability
    if(sum(is.na(est.rho))>0){
      est.rho = rep(1/k, k)
    }

    for(j in 1:k){
      sample.prob[j] = g.func(prop.k[j], est.rho[j], r = 2)
    }
    m = NULL #?????
    for(i in (n0+1):ssn){

      alloc[i] = sample(c(1:k), 1, prob = sample.prob)
      obs.outcome[i] = outcome[i, alloc[i]]

      # detect is there anyone's outcome ready
      # time of m-th patient' outcome is available
      obs.cumRspT[i] = cumRspT[i, alloc[i]]
      #
      # the outcome could be observed after the m-th entry
      m[i-n0] = max(which(entryT <= obs.cumRspT[i]))

      # whether adjust the probability
      if(sum(m == i) == 0){
        next
      }else{
        if(i+1 > ssn) break
        temp = calc_thetaGaussian_MD(theta.hat, k, alloc, outcome, entryT[i+1], obsRspT = obs.cumRspT, mRate = mRate)
        theta.hat = temp[[1]]
        prop.k = temp[[2]]

        # calc sample prob
        est.rho = target.rho.Ctinuous(theta.hat, target.alc = target.alloc)  #
        if(sum(is.na(est.rho))>0){
          est.rho = rep(1/k, k)
        }
        for(j in 1:k){
          sample.prob[j] = g.func(prop.k[j], est.rho[j], r = 2)
        }
      }
    }

    # re-adjust obs.outcome and assign.group
    if(!is.null(mRate)){
      effIdx = which(outcome[, k+1] == 0)
      obs.outcome = obs.outcome[effIdx]
      alloc = alloc[effIdx]
    }

    theta.hat.set = rbind(theta.hat.set, theta.hat)
    if(k == 2){
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group = alloc)
    }else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group = alloc, k)
    }
    #failure.rate[s] = mean(obs.outcome == 0)
    failure.rate[s] = mean(obs.outcome)
    group.prop = rbind(group.prop, table(alloc) / length(alloc))
  }
  name = ifelse(is.null(mRate),
                "Hu and Zhang Delayed DBCD without Missing Data",
                "Hu and Zhang Delayed DBCD with Missing Data")

  return(RAR_Output(name, parameter=theta, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}




###############################################################################
#################   Group doubly biased coin Design (Binary)  ##################
###############################################################################

# lambda is the parameter for poison distribution which control the rate for patients enter the group
#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assigned through restricted randomization for initial parameter estimation.
#' @param p A positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param gsize.param A positive integer. It represents the expected number of people enrolling in a specified interval, assuming that the enrollment rate per unit time follows a Poisson distribution.
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param theta0 A vector of length k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length k, with all values set to 0.5.
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}. The default is \code{"RPW"}.
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param mRate a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples
#' Group.DBCD_Bin(n0 = 20, p = c(0.65, 0.8), k = 2, gsize.param = 5, ssn = 300, theta0 = NULL, target.alloc = "RPW", r = 2, nsim = 500, mRate = NULL, alpha = 0.05)
Group.DBCD_Bin = function(n0 = 20, p, k, gsize.param, ssn, theta0 = NULL, target.alloc = "RPW", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05){
  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }

  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

  for(s in 1:nsim){
    p.hat = NULL
    sample.prob = NULL
    obs.outcome = NULL
    gsize = NULL
    outcome = generate_data_M(p, ssn, mRate, k)
    alloc = c()

    # initial part
    initial = T
    i = 1
    while (initial) {
      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
      pmtAlloc = sample(rep(1:k, ceiling(gsize[i] / k)))
      alloc = c(alloc, pmtAlloc[1:gsize[i]])
      if (sum(gsize) >= n0) {
        initial = F
      }
      i = i+1
    }

    while(sum(gsize) < ssn){
      # allocation probability
      if(is.null(mRate)){
        p.hat = calc_theta(p.hat, k, alloc, outcome, theta0 = theta0)
      }else{
        p.hat = calc_theta_M(p.hat, k, alloc, outcome, theta0 = theta0)
      }
      est.rho = target.rho(p.hat, target.alc = target.alloc)
      prop.k = calc_prop(alloc, outcome, mRate, k)
      for(j in 1:k) {
        sample.prob[j] = g.func(prop.k[j], est.rho[j], r)
      }
      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
      alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = T ))
      i = i+1
    }

    alloc = alloc[1:ssn]
    for(j in 1:k){
      obs.outcome[which(alloc == j)] = outcome[which(alloc == j), j]
    }
    if(!is.null(mRate)){
      obs.outcome = obs.outcome[outcome[, k+1] == 0]
      alloc = alloc[outcome[, k+1] == 0]
    }

    if(k == 2){
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group = alloc)
    }else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group = alloc, k)
    }

    failure.rate[s] = mean(obs.outcome == 0)
    group.prop = rbind(group.prop, table(alloc) / length(alloc))
  }
  name = ifelse(is.null(mRate),
                "Group DBCD with Binary Response (No Missing)",
                "Group DBCD with Binary Response (Random Missing)")

  return(RAR_Output(name, parameter=p, ssn = c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate, #shouldn't be reported
                    pwCalc, k))
}



###############################################################################
###############   Group doubly biased coin Design (delayed+bin)  ###############
###############################################################################

# eTime = enter Time (daily, weekly, biweekly, monthly)
#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.
#' @param p A positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k = 2})
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param gsize.param A positive integer. It represents the expected number of people enrolling in a specified interval, assuming that the enrollment rate per unit time follows a Poisson distribution.
#' @param rspT.dist Distribution Type. Specifies the type of distribution that models the time spent for the availability of patient \eqn{i} under treatment \eqn{k}. Acceptable options for this argument include: \code{"exponential"}, \code{"normal"}, and \code{"uniform"}.
#' @param rspT.param A vector with length \eqn{2k}. Specifies the parameters required by the distribution that models the time spent for the availability under each treatment and each response. (eg. If there are 3 treatments groups with 0 or 1 as response and each of them follows exponential distribution with parameter (3, 2, 3, 3, 1, 2), repectively. Then the \code{rspT.param = c(3, 2, 2, 1, 4, 1)})
#' @param theta0 A vector of length k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length k, with all values set to 0.5.
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}. The default is \code{"RPW"}.
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim A positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param eTime A positive number. The interval time between enrollment of participants in each group. The default is 7.
#' @param mRate a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.
#' @param alpha a numerical value between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#'
#' @export
#'
#' @examples
#' # a simple use
#' # Define the arguments
#' ## Arguments for generate the simulated data
#' ### For response simulation
#' p = c(0.7, 0.5)
#' k = 2
#' ssn = 200
#'
#' ### for enter time and response time simulation
#' eTime = 7
#' rspT.param = rep(10, 4)
#' rspT.dist = "exponential"
#' gsize.param = 5
#'
#' ## Arguments for the deisgn
#' n0 = 10
#' target.alloc = "RPW"
#'
#' Group.dyldDBCD_Bin(n0, p, k, ssn, gsize.param, rspT.dist, rspT.param,
#'                          theta0 = NULL, target.alloc = "RPW",  r = 2, nsim = 150, eTime = 7,  mRate = NULL, alpha = 0.05)
Group.dyldDBCD_Bin = function(n0 = 20, p, k, ssn, gsize.param, rspT.dist, rspT.param,
                              theta0 = NULL, target.alloc = "RPW",  r = 2, nsim = 2000, eTime = 7,  mRate = NULL, alpha = 0.05){

  if(k != length(p)){
    stop("Length of p must be equal to k")
  }
  if(length(which((p<0) | (p>1))) > 0){
    stop("Each components in the vector p is required to be between 0 and 1")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }
  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

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
    initial = T
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

      temp = calc_RspT(alloc, outcome, rspT, gsize, entry)
      obs.outcome = c(obs.outcome, temp[,2])
      obsRspT = c(obsRspT, temp[,1])
      i = i+1
      entry[i] = eTime
    }

    while (sum(gsize) <= ssn) {

      # calculate the time for delayed response
      temp = calc_RspT(alloc, outcome, rspT, gsize, entry)
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
      est.rho = target.rho(p.hat, target.alc = target.alloc)

      for(j in 1:k) {
        sample.prob[j] = g.func(prop.k[j], est.rho[j], r)
      }

      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf) #rpois(1, gsize.param)
      alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = T))
    }

    alloc = alloc[1:ssn]
    obs.outcome = c(obs.outcome, calc_RspT(alloc, outcome, rspT, gsize, entry, adjust = T)[, 2])
    if(is.null(mRate)){effIdx = c(1:ssn)}else{effIdx = which(outcome[, k+1] == 0)}

    if(k == 2){
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome[effIdx], assign.group = alloc[effIdx])
    }else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group = alloc, k)
    }

    failure.rate[s] = mean(obs.outcome[effIdx] == 0)
    group.prop = rbind(group.prop, table(alloc[effIdx]) / length(alloc[effIdx]))
  }

  name = ifelse(is.null(mRate),
                "Group DBCD with Binary Delayed Response (No Missing)",
                "Group DBCD with Binary Delayed Response (Random Missing)")

  return(RAR_Output(name, parameter=p, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))

}


###############################################################################
##############   Group doubly biased coin Design (Continuous)  #################
###############################################################################

#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.
#' @param theta A numerical vector of length equal to \code{2k}. These values specify the true parameters for each treatment and are used for generating data in simulations. For example, if \code{k=2}, you should provide two pairs of parameter values, each consisting of the mean and variance, like: \code{theta = c(13, 4.0^2, 15, 2.5^2)}.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param gsize.param A positive integer. It represents the expected number of people enrolling in a specified interval, assuming that the enrollment rate per unit time follows a Poisson distribution.
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param theta0 A vector of length 2k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length 2k, with all paramter pair setted to be (0, 1).
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"ZR"}, \code{"DaOptimal"}. The default is \code{"Neyman"}. The details see Zhang L. and Rosenberger. W (2006).
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param mRate a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#' @export
#'
#' @examples
#' theta = c(13, 4.0^2, 15, 2.5^2)
#' k = 2
#' gsize.param = 5
#' ssn = 120
#' Group.DBCD_Cont(n0 = 20, theta, k, gsize.param, ssn, target.alloc = "Neyman", r = 2, nsim = 200, mRate = NULL, alpha = 0.05)
Group.DBCD_Cont = function(n0 = 20, theta, k, gsize.param, ssn, theta0 = NULL, target.alloc = "Neyman", r = 2, nsim = 2000, mRate = NULL, alpha = 0.05){
  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(F, T)] < 0) > 0){
    stop("The vairance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }
  # setup
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

  for(s in 1:nsim){
      theta.hat = NULL
      sample.prob = NULL
      obs.outcome = NULL
      gsize = NULL
      outcome = generate_GaussianRsp(theta, k, ssn)
      alloc = c()

      # initial part
      initial = T
      i = 1
      while (initial) {
        gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
        pmtAlloc = sample(rep(1:k, ceiling(gsize[i] / k)))
        alloc = c(alloc, pmtAlloc[1:gsize[i]])
        if (sum(gsize) >= n0) {
          initial = F
        }
        i = i+1
      }

      while(sum(gsize) < ssn){
        # allocation probability
        theta.hat = calc_thetaGaussian(theta.hat, k, alloc, outcome, theta0)
        est.rho = target.rho.Ctinuous(para.set = theta.hat, target.alloc)
        prop.k = calc_prop(alloc, outcome, mRate, k)
        for(j in 1:k) {
          sample.prob[j] = g.func(prop.k[j], est.rho[j], r = r)
        }
        gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf)
        alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = T))
        i = i+1
      }

      alloc = alloc[1:ssn]
      for(j in 1:k){
        obs.outcome[which(alloc == j)] = outcome[which(alloc == j), j]
      }
      if(!is.null(mRate)){
        obs.outcome = obs.outcome[outcome[, k+1] == 0]
        alloc = alloc[outcome[, k+1] == 0]
      }

      if(k == 2){
        pwCalc[s] = ttest.2(alpha = alpha, obs.outcome, assign.group = alloc)
      }else{
        pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome, assign.group = alloc, k)
      }

      failure.rate[s] = mean(obs.outcome)
      group.prop = rbind(group.prop, table(alloc) / length(alloc))
  }

  name = ifelse(is.null(mRate),
                "Group DBCD with Continuous Response (No Missing)",
                "Group DBCD with Continuous Response (Random Missing)")

  return(RAR_Output(name, parameter = theta, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}


###############################################################################
#############   Group doubly biased coin Design (delayed+Cont)  ################
###############################################################################

#' Title
#'
#' @param n0 A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.
#' @param theta A numerical vector of length equal to \code{2k}. These values specify the true parameters for each treatment and are used for generating data in simulations. For example, if \code{k=2}, you should provide two pairs of parameter values, each consisting of the mean and variance, like: \code{theta = c(13, 4.0^2, 15, 2.5^2)}.
#' @param k A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})
#' @param ssn A positive integer. The value specifies the total number of participants involved in each round of the simulation.
#' @param gsize.param A positive integer. It represents the expected number of people enrolling in a specified interval, assuming that the enrollment rate per unit time follows a Poisson distribution.
#' @param rspT.dist Distribution Type. Specifies the type of distribution that models the time spent for the availability of patient \eqn{i} under treatment \eqn{k}. Acceptable options for this argument include: \code{"exponential"}, \code{"normal"}, and \code{"uniform"}.
#' @param rspT.param A vector. Specifies the parameters required by the distribution that models the time spent for the availability under each treatment. (eg. If there are 3 treatments groups and each of them follows truncated normal distribution with parameter pair (3, 2), (2, 1), (4, 1), repectively. Then the \code{rspT.param = c(3, 2, 2, 1, 4, 1)})
#' @param target.alloc Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"ZR"}, \code{"DaOptimal"}. The default is \code{"Neyman"}. The details see Zhang L. and Rosenberger. W (2006).
#' @param r A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.
#' @param nsim a positive integer. The value specifies the total number of simulations, with a default value of 2000.
#' @param eTime A positive number. The interval time between enrollment of participants in each group. The default is 7.
#' @param mRate a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.
#' @param alpha A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.
#' @export
#'
#' @examples
#' k = 2; ssn = 120
#' theta = c(13, 4.0^2, 15, 2.5^2)
#' gsize.param = 5
#' rspT.param = rep(10, 2)
#' rspT.dist = "exponential"
#' Group.dyldDBCD_Cont(n0 = 20, theta, k, ssn, gsize.param, rspT.dist, rspT.param,
#'                     target.alloc = "Neyman",  r = 2, nsim = 500, eTime = 7,  mRate = 0.2, alpha = 0.05)
Group.dyldDBCD_Cont = function(n0 = 20, theta, k, ssn, gsize.param, rspT.dist, rspT.param,
                               target.alloc = "Neyman",  r = 2, nsim = 2000, eTime = 7,  mRate = NULL, alpha = 0.05){

  if((2 * k) != length(theta)){
    stop("Length of theta vector must be equal to 2k")
  }
  if(sum(theta[c(F, T)] < 0) > 0){
    stop("The vairance should be a positive number")
  }
  if((n0 %% k) != 0){
    stop("Number of initial participates 'n0' must be a multiple of k")
  }

  # setup
  sim.prop = NULL
  pwCalc = NULL
  failure.rate = NULL
  group.prop = c()

  for(s in 1:nsim){
    theta.hat = NULL
    obs.outcome = NULL
    sample.prob = NULL
    entry = NULL
    gsize = NULL
    obsRspT = c()
    alloc = c()
    outcome = generate_GaussianRsp_M(theta, ssn, mRate, k)
    #apply(outcome, 2, sd)
    rspT = responseDist(rspT.dist, rspT.param, k, level = 1, ssn)

    # initial part
    initial = T
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

    while (sum(gsize) <= ssn) {

      # calculate the time for delayed response
      temp = calc_RspT(alloc, outcome, rspT, gsize, entry)
      obs.outcome = c(obs.outcome, temp[,2])
      obsRspT = c(obsRspT, temp[,1])

      # update the entry time
      i = i+1
      entry[i] = eTime
      entryT = sum(entry)

      # calculate allocation probability
      temp = calc_thetaGaussian_MD(theta.hat, k, alloc, outcome, entryT, obsRspT, mRate)
      theta.hat = temp[[1]]
      prop.k = temp[[2]]
      # if theta.hat = na ==> est.rho = c(1/2, 1/2)
      est.rho = target.rho.Ctinuous(theta.hat, target.alc = target.alloc)

      # # if no response came out, use 1/k as probability
      if(sum(is.na(est.rho))>0){
        est.rho = rep(1/k, k)
      }

      for(j in 1:k) {
        sample.prob[j] = g.func(prop.k[j], est.rho[j], r)
      }

      gsize[i] = extraDistr::rtpois(1, gsize.param, a = 0, b = Inf) #rpois(1, gsize.param)
      alloc = c(alloc, sample.int(k, size = gsize[i], prob = sample.prob, replace = T))
    }

    alloc = alloc[1:ssn]
    obs.outcome = c(obs.outcome, calc_RspT(alloc, outcome, rspT, gsize, entry, adjust = T)[, 2])
    if(is.null(mRate)){effIdx = c(1:ssn)}else{effIdx = which(outcome[, k+1] == 0)}
    if(k == 2){
      pwCalc[s] = ttest.2(alpha = alpha, obs.outcome[effIdx], assign.group = alloc[effIdx])
    }else{
      pwCalc[s] = chisq.test.k(alpha = alpha, obs.outcome[effIdx], assign.group = alloc[effIdx], k)
    }
    failure.rate[s] = mean(obs.outcome[effIdx])
    #failure.rate[s] = mean(obs.outcome[effIdx] == 0)
    group.prop = rbind(group.prop, table(alloc[effIdx]) / length(alloc[effIdx]))
  }
  name = ifelse(is.null(mRate),
                "Group DBCD with Continuous Delayed Response (No Missing)",
                "Group DBCD with Continuous Delayed Response (Random Missing)")
  return(RAR_Output(name, parameter = theta, ssn =  c("Total" = ssn, "Effective Size" = ifelse(is.null(mRate), ssn, ssn * (1-mRate))),
                    assignment = alloc, propotion = group.prop,
                    failRate = failure.rate,
                    pwCalc, k))
}


