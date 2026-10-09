test_that("nextAlloc reproduces the DBCD allocation probability", {
  set.seed(6)
  alloc = sample(1:3, 30, TRUE)
  y = rbinom(30, 1, c(0.6, 0.7, 0.8)[alloc])
  y[28:30] = NA
  res = nextAlloc(alloc, y, k = 3, target.alloc = "RSIHR", r = 2, size = 4)
  obs = !is.na(y)
  p.hat = sapply(1:3, function(j) (sum(y[obs & alloc == j]) + 0.5) / (sum(obs & alloc == j) + 1))
  expect_equal(unname(res$estimate), p.hat)
  expect_equal(unname(res$prob), g.func(tabulate(alloc, 3) / 30, sqrt(p.hat) / sum(sqrt(p.hat)), 2))
  expect_length(res$assignment, 4)
  expect_true(all(res$assignment %in% 1:3))
})

test_that("nextAlloc handles continuous responses, ERADE and an empty trial", {
  set.seed(1)
  alloc = rep(1:2, 10)
  y = rnorm(20, c(10, 12)[alloc], 2)
  res = nextAlloc(alloc, y, k = 2, response = "continuous", allocation = "ERADE")
  expect_equal(sum(res$prob), 1)
  expect_equal(unname(nextAlloc(integer(0), numeric(0), k = 3)$prob), rep(1/3, 3))
  expect_error(nextAlloc(c(1, 4), c(1, 0), k = 3), "arm labels")
})
