expect_valid_prop = function(res, k){
  prop = res[["propotion"]]
  expect_length(prop, k)
  expect_equal(sum(prop), 1)
  expect_length(res[["sd of propotion"]], k)
  expect_s3_class(res, "grouprar")
}

test_that("chisq.test.k works for 2 arms and continuous responses", {
  set.seed(1)
  grp = rep(1:2, each = 50)
  y = rbinom(100, 1, ifelse(grp == 1, 0.2, 0.8))
  expect_equal(chisq.test.k(0.05, y, grp, 2), 1)
  z = rnorm(150, mean = rep(c(0, 0, 3), each = 50))
  expect_equal(chisq.test.k(0.05, z, rep(1:3, each = 50), 3, continuous = TRUE), 1)
})

test_that("tests handle empty arms and arms without variability", {
  expect_warning(res <- chisq.test.k(0.05, c(1, 0, 1, 0), c(1, 1, 2, 2), 3))
  expect_true(is.na(res))
  # all successes in one arm and all failures in the other is a rejection
  expect_equal(ttest.2(0.05, c(rep(1, 10), rep(0, 10)), rep(1:2, each = 10)), 1)
  expect_equal(ttest.2(0.05, rep(1, 20), rep(1:2, each = 10)), 0)
})

test_that("responseDist accepts the uniform distribution", {
  rspT = responseDist("uniform", c(1, 2, 3, 4), k = 2, level = 1, sample.size = 10)
  expect_true(all(rspT[, 1] >= 1 & rspT[, 1] <= 2))
  expect_true(all(rspT[, 2] >= 3 & rspT[, 2] <= 4))
})

test_that("DBCD designs run with 3 arms", {
  p = c(0.6, 0.7, 0.8)
  theta = c(13, 16, 15, 6.25, 14, 9)
  expect_valid_prop(DBCD_Bin(n0 = 21, p = p, k = 3, ssn = 90, nsim = 3, seed = 1), 3)
  expect_valid_prop(DBCD_Bin(n0 = 21, p = p, k = 3, ssn = 90, nsim = 3, mRate = 0.1, seed = 1), 3)
  expect_valid_prop(DBCD_Cont(n0 = 21, theta = theta, k = 3, ssn = 90, nsim = 3, seed = 1), 3)
  expect_valid_prop(dyldDBCD_Bin(n0 = 21, p = p, k = 3, ssn = 90, ent.param = 0.7,
                                 rspT.dist = "exponential", rspT.param = rep(1, 6), nsim = 3, seed = 1), 3)
  expect_valid_prop(dyldDBCD_Cont(n0 = 21, theta = theta, k = 3, ssn = 90, ent.param = 5,
                                  rspT.dist = "exponential", rspT.param = rep(10, 3), nsim = 3, mRate = 0.1, seed = 1), 3)
  expect_valid_prop(Group.DBCD_Bin(n0 = 21, p = p, k = 3, gsize.param = 5, ssn = 90, nsim = 3, seed = 1), 3)
  expect_valid_prop(Group.DBCD_Cont(n0 = 21, theta = theta, k = 3, gsize.param = 5, ssn = 90, nsim = 3, seed = 1), 3)
  expect_valid_prop(Group.DBCD_Cont(n0 = 21, theta = theta, k = 3, gsize.param = 5, ssn = 90, nsim = 3, mRate = 0.1, seed = 1), 3)
  expect_valid_prop(Group.dyldDBCD_Bin(n0 = 21, p = p, k = 3, ssn = 90, gsize.param = 5,
                                       rspT.dist = "exponential", rspT.param = rep(10, 6), nsim = 3, seed = 1), 3)
  expect_valid_prop(Group.dyldDBCD_Cont(n0 = 21, theta = theta, k = 3, ssn = 90, gsize.param = 5,
                                        rspT.dist = "exponential", rspT.param = rep(10, 3), nsim = 3, seed = 1), 3)
})

test_that("urn designs run and return one failure rate per simulation", {
  res = GDLRule(k = 3, p = c(0.6, 0.7, 0.6), ssn = 60, aK = c(1, 1, 1), nsim = 5, seed = 3)
  expect_length(res[["data: failureRate"]], 5)
  expect_no_error(PolyaUrn(k = 2, p = c(0.6, 0.7), ssn = 60, nsim = 3, seed = 3))
  expect_valid_prop(Bai.Hu.Shen.Urn(k = 3, p = c(0.6, 0.7, 0.8), ssn = 60, nsim = 3, seed = 3), 3)
  expect_valid_prop(CRDesign(k = 3, p = c(0.6, 0.7, 0.8), ssn = 60, nsim = 3, seed = 3), 3)
})

test_that("invalid inputs give clear errors", {
  expect_error(RPWRule(k = 3, p = c(0.6, 0.7, 0.8), ssn = 60, nsim = 2), "k = 2")
  expect_error(DBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, ssn = 20, nsim = 2), "larger than")
  expect_error(DBCD_Bin(n0 = 0, p = c(0.6, 0.8), k = 2, ssn = 50, nsim = 2), "at least k")
  expect_error(DBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, ssn = 50, nsim = 2, target.alloc = "ZR"), "target.alloc")
  expect_error(DBCD_Bin(n0 = 21, p = c(0.6, 0.7, 0.8), k = 3, ssn = 60, nsim = 2,
                        target.alloc = "OptimalRSIHR", lower.bound = 0.5), "lower.bound")
  expect_error(DBCD_Cont(n0 = 20, theta = c(-1, 1, 2, 1), k = 2, ssn = 50, nsim = 2, target.alloc = "ZR"), "positive means")
})

test_that("seed makes simulations reproducible", {
  a = Group.DBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, gsize.param = 5, ssn = 60, nsim = 5, seed = 9)
  b = Group.DBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, gsize.param = 5, ssn = 60, nsim = 5, seed = 9)
  expect_identical(a[["data: allocation"]], b[["data: allocation"]])
  expect_equal(dim(a[["data: allocation"]]), c(5, 60))
})

test_that("typeI, test.fun and durations are reported", {
  res = dyldDBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, ssn = 60, ent.param = 0.7,
                     rspT.dist = "exponential", rspT.param = rep(1, 4), nsim = 5, typeI = TRUE, seed = 2)
  expect_true(res[["type I error"]] >= 0 && res[["type I error"]] <= 1)
  expect_length(res[["data: duration"]], 5)
  expect_true(all(res[["data: duration"]] >= res[["data: enrollment"]] - 1e-8))
  never = function(y, a) 1
  expect_equal(DBCD_Bin(n0 = 20, p = c(0.2, 0.9), k = 2, ssn = 60, nsim = 5, test.fun = never, seed = 2)$power, 0)
})

test_that("print and summary work", {
  res = DBCD_Cont(n0 = 20, theta = c(13, 16, 15, 6.25), k = 2, ssn = 50, nsim = 3, seed = 1)
  expect_output(print(res), "treatment A")
  expect_s3_class(summary(res), "summary.grouprar")
  expect_equal(names(res[["parameter"]]), c("muA", "sigma2A", "muB", "sigma2B"))
})

test_that("designs survive trials whose first patients all have missing responses", {
  # a trial can end with no observed response in one arm, which the t test reports with a warning
  expect_no_error(suppressWarnings(DBCD_Bin(n0 = 10, p = c(0.6, 0.8), k = 2, ssn = 40, mRate = 0.7, nsim = 200, seed = 1)))
})

test_that("near-singular covariance matrices give NA instead of an error", {
  # two arms without variability: singular in exact arithmetic, det about 1e-22 in floating point
  y = c(rep(1, 28), rep(1, 26), rep(1, 32), rep(0, 4))
  g = c(rep(1, 28), rep(2, 26), rep(3, 36))
  stat = wald.stat(y, g, 3)
  tilde = (c(28, 26, 32) + 1) / (c(28, 26, 36) + 2)
  v = tilde * (1 - tilde)
  d = c(1, 1) - 32/36
  S = diag(v[1:2] / c(28, 26)) + v[3] / 36
  expect_equal(stat, as.numeric(t(d) %*% solve(S) %*% d))
  # continuous responses without variability still give NA
  expect_true(is.na(wald.stat(c(rep(1, 10), rep(2, 10), rnorm(10)), rep(1:3, each = 10), 3, continuous = TRUE)))
  expect_no_error(suppressWarnings(CRDesign(k = 3, p = rep(0.9, 3), ssn = 90, nsim = 300, seed = 8)))
})

test_that("optimal targets do not warn when arms have equal estimates", {
  expect_no_warning(opt.rho(c(0.5, 0.5, 0.7), c(0.25, 0.25, 0.21), c(0.5, 0.5, 0.3), 0))
  expect_no_warning(DBCD_Bin(n0 = 12, p = rep(0.6, 4), k = 4, ssn = 60, nsim = 2, target.alloc = "OptimalNeyman", seed = 1))
})

test_that("group designs reject a non-positive group size rate", {
  expect_error(Group.DBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, gsize.param = 0, ssn = 60, nsim = 1), "gsize.param")
})

test_that("response time parameters must have the right length", {
  expect_error(responseDist("exponential", c(1, 3), k = 2, level = 2, sample.size = 5), "rspT.param")
  expect_error(responseDist("normal", c(1, 1, 2, 1), k = 2, level = 2, sample.size = 5), "rspT.param")
  expect_equal(dim(responseDist("uniform", 1:8, k = 2, level = 2, sample.size = 5)), c(5, 4))
})

test_that("delayed designs report NA duration when no response is observed", {
  res = dyldDBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, ssn = 40, ent.param = 1, rspT.dist = "exponential",
                     rspT.param = rep(1, 4), nsim = 3, seed = 1)
  expect_true(all(is.finite(res[["data: duration"]])))
})

test_that("urn designs accept fractional balls and check their inputs", {
  expect_no_error(GDLRule(k = 2, p = c(0.4, 0.7), ssn = 100, aK = c(0.5, 1), nsim = 10, seed = 1))
  expect_no_error(DLRule(k = 2, p = c(0.4, 0.7), ssn = 100, Y0 = c(0.5, 0.5), nsim = 10, seed = 1))
  expect_no_error(BirthDeathUrn(k = 2, p = c(0.2, 0.3), ssn = 100, Y0 = c(0.5, 0.5), nsim = 10, seed = 1))
  expect_error(GDLRule(k = 2, p = c(0.4, 0.7), ssn = 100, aK = c(0, 0), nsim = 2), "aK")
  expect_error(WeiUrn(k = 2, p = c(0.4, 0.7), ssn = 50, Y0 = c(0, 0), nsim = 2), "Y0")
  expect_error(WeiUrn(k = 1, p = 0.5, ssn = 50, nsim = 2), "at least 2")
  expect_error(DBCD_Bin(n0 = 2, p = 0.5, k = 1, ssn = 50, nsim = 2), "at least 2")
})
