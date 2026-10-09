# Hu and Zhang (2004) k-arm allocation function
hu.zhang = function(x, y, r){
  tmp = y * (y / x) ^ r
  tmp / sum(tmp)
}

test_that("g.func matches the Hu and Zhang k-arm formula for 3 arms", {
  x = rep(1/3, 3)
  y = c(0.5, 0.3, 0.2)
  expect_equal(g.func(x, y, 2), c(0.78125, 0.16875, 0.05))
  expect_equal(g.func(x, y, 2), hu.zhang(x, y, 2))
})

test_that("g.func matches the Hu and Zhang k-arm formula for random inputs", {
  set.seed(2024)
  for(k in 2:5){
    for(r in c(0, 1, 2, 4)){
      x = prop.table(runif(k))
      y = prop.table(runif(k))
      expect_equal(g.func(x, y, r), hu.zhang(x, y, r))
      expect_equal(sum(g.func(x, y, r)), 1)
    }
  }
})

test_that("g.func is unchanged for 2 arms", {
  old.g = function(x, y, r){
    tmp1 = y * (y / x) ^ r
    tmp2 = (1-y) * ((1-y) / (1-x)) ^ r
    tmp1 / (tmp1 + tmp2)
  }
  set.seed(1)
  for(i in 1:20){
    x = runif(1); y = runif(1); r = runif(1, 0, 4)
    expect_equal(g.func(c(x, 1-x), c(y, 1-y), r)[1], old.g(x, y, r))
  }
})

test_that("g.func returns the target when allocation is on target", {
  y = c(0.2, 0.3, 0.5)
  expect_equal(g.func(y, y, 2), y)
})

test_that("g.func gives arms without allocation all the probability", {
  expect_equal(g.func(c(0, 1), c(0.4, 0.6), 2), c(1, 0))
  expect_equal(g.func(c(0, 0.5, 0, 0.5), rep(0.25, 4), 2), c(0.5, 0, 0.5, 0))
})

test_that("g.func accepts a table of allocation proportions", {
  alloc = c(1, 1, 2, 3, 3, 3)
  x = table(alloc) / length(alloc)
  expect_equal(g.func(x, c(0.5, 0.3, 0.2), 2), hu.zhang(as.numeric(x), c(0.5, 0.3, 0.2), 2))
})
