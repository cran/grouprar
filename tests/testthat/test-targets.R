test_that("optimal targets reduce to Neyman, RSIHR and ZR for 2 arms", {
  p = c(0.7, 0.4)
  expect_equal(target.rho(p, "OptimalNeyman"), target.rho(p, "Neyman"))
  expect_equal(target.rho(p, "OptimalRSIHR"), target.rho(p, "RSIHR"))
  # Zhang and Rosenberger (2006) formula when it favors the arm with the smaller mean
  theta = c(13, 16, 15, 6.25)
  zr = sqrt(16) * sqrt(15) / (sqrt(16) * sqrt(15) + sqrt(6.25) * sqrt(13))
  expect_equal(target.rho.Ctinuous(theta, "ZR"), c(zr, 1 - zr))
  # rule (7): 1/2 when the formula would favor the worse arm
  expect_equal(target.rho.Ctinuous(c(13, 1, 15, 16), "ZR"), c(0.5, 0.5))
})

test_that("k-arm optimal target matches a brute-force search", {
  obj = function(rho, theta, v, w){
    a = rho / v
    sum(w * rho) / sum(a * (theta - sum(a * theta) / sum(a))^2)
  }
  set.seed(11)
  for(i in 1:10){
    k = sample(3:4, 1)
    p = runif(k, 0.1, 0.9)
    B = runif(1, 0, 0.9 / k)
    v = p * (1-p)
    sol = opt.rho(p, v, 1-p, B)
    expect_equal(sum(sol), 1)
    expect_true(all(sol >= B - 1e-12))
    u = matrix(rexp(3000 * k), ncol = k)
    cand = B + (1 - k * B) * u / rowSums(u)
    best = min(apply(cand, 1, obj, theta = p, v = v, w = 1-p))
    expect_lte(obj(sol, p, v, 1-p), best + 1e-8)
  }
  # the middle arm sits at the lower bound
  expect_equal(opt.rho(c(0.8, 0.6, 0.4), c(0.16, 0.24, 0.24), c(0.2, 0.4, 0.6), 0.1)[2], 0.1)
})

test_that("multi-arm ERADE follows Alkhnefr, Hu and Zhai (2025) and reduces to two arms", {
  a = 0.5
  # two arms: Hu, Zhang and He (2009)
  expect_equal(erade.func(c(0.7, 0.3), c(0.6, 0.4), a), c(a * 0.6, 1 - a * 0.6))
  expect_equal(erade.func(c(0.5, 0.5), c(0.6, 0.4), a), c(1 - a * 0.4, a * 0.4))
  expect_equal(erade.func(c(0.6, 0.4), c(0.6, 0.4), a), c(0.6, 0.4))
  # three arms with equal target: scenario 2 of their Example 1
  rho = rep(1/3, 3)
  expect_equal(erade.func(c(0.5, 0.25, 0.25), rho, a), c(a/3, (1-a)/6 + 1/3, (1-a)/6 + 1/3))
  expect_equal(sum(erade.func(c(0.5, 0.3, 0.2), c(0.3, 0.3, 0.4), 2/3)), 1)
})

test_that("ERADE allocates closer to the target than DBCD", {
  d = DBCD_Bin(n0 = 30, p = c(0.5, 0.7, 0.8), k = 3, ssn = 200, target.alloc = "RSIHR", nsim = 60, seed = 2)
  e = DBCD_Bin(n0 = 30, p = c(0.5, 0.7, 0.8), k = 3, ssn = 200, target.alloc = "RSIHR", nsim = 60, seed = 2,
               allocation = "ERADE")
  expect_lt(mean(e[["sd of propotion"]]), mean(d[["sd of propotion"]]))
})
