test_that("boundaries match Zhu and Hu (2010)", {
  t = c(0.2, 0.5, 1)
  expect_equal(unname(round(sqBoundary(t, 0.05, "OBF"), 3)), c(4.877, 2.963, 1.969))
  expect_equal(unname(round(sqBoundary(t, 0.05, "Linear"), 3)), c(2.576, 2.377, 2.141))
  expect_equal(unname(round(sqBoundary(t, 0.05, "Pocock"), 3)), c(2.438, 2.333, 2.225))
  expect_equal(unname(sqBoundary(1, 0.05)), qnorm(0.975))
})

test_that("monitored designs stop early and report stopping probabilities", {
  m = sqMonitor(c(1/3, 2/3))
  expect_equal(m$t, c(1/3, 2/3, 1))
  res = DBCD_Bin(n0 = 20, p = c(0.3, 0.9), k = 2, ssn = 120, target.alloc = "RSIHR", nsim = 20, monitor = m, seed = 1)
  expect_equal(sum(res[["stopping probability"]]), 1)
  expect_lt(res[["expected sample size"]], 120)
  expect_true(any(is.na(res[["data: allocation"]])))
  g = Group.DBCD_Cont(n0 = 20, theta = c(10, 4, 14, 4), k = 2, gsize.param = 5, ssn = 120, nsim = 10, monitor = m, seed = 1)
  expect_equal(sum(g[["stopping probability"]]), 1)
  expect_error(DBCD_Bin(n0 = 21, p = c(0.6, 0.7, 0.8), k = 3, ssn = 90, nsim = 2, monitor = m), "k = 2")
})

test_that("boundaries are fast and accurate for many looks and tiny alpha", {
  b = sqBoundary(seq(0.1, 1, by = 0.1), 0.05, "OBF")
  expect_length(b, 10)
  expect_true(all(is.finite(b)))
  expect_true(all(diff(b) < 0))
  expect_true(all(is.finite(sqBoundary(c(0.2, 0.5, 1), 1e-12, "Pocock"))))
  expect_true(is.finite(sqBoundary(c(0.1, 0.3, 1), 0.01, "OBF")[1]))
})
