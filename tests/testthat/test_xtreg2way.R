context("xtreg2way output")
library("xtreg2way")
test_tolerance <- 1e-6

hhid <- c("a","b","c","a","b","c" ,"a","b","c" ,"a","b","c" ,"a","b","c")
tid <- c("1","1" ,"1" ,"2","2" ,"3","3","3" ,"4","4","5" ,"5","6","6" ,"6")
w <- rep(1, 15)
#Hard coded some random values that I know the answer to
x1 <- c(45.90712, 54.00554, 47.24599, 54.67337, 50.16003, 52.01525,
        48.25807, 50.59948, 32.07548, 41.16698, 45.31226, 52.86234, 66.36729,
         26.82314, 52.20836)
x2 <- c(48.37376, 47.24714, 63.05863, 56.89668, 52.63793, 39.35706, 49.70856,
        56.22179, 30.49866, 60.17268, 58.81930, 50.87553, 39.65359, 47.51820,
        55.15367) 
y <- c(94.15298, 118.54982, 104.81954, 105.56196, 100.03579, 98.61846, 
       124.25847, 119.16457,  85.04116, 100.39672,  93.22992,  89.24267,
       106.53424,  82.00621, 113.93457)
output <- xtreg2way(y, data.frame(x1,x2), hhid, tid, w, se="2", noise="1")
output2 <- xtreg2way(y, x1, struc=output$struc)

test_that("First run output has correct objects", {
  expect_named(output, c("betaHat", "aVarHat", "y", "X", "struc"))
})
test_that("Second run output has correct objects", {
  expect_named(output2, c("betaHat", "aVarHat"))
})
test_that("Output Run 1 is correct", {
  expect_equal(output$betaHat, matrix(data = c(0.8793373,0.5381357), 
                                      nrow=1, ncol=2),
               tolerance = test_tolerance)
  expect_equal(output$aVarHat, matrix(data=c(0.0420780855,0.0005067944,
                                             0.0005067944,0.0336658460), 
                                      nrow=2, ncol=2),
               tolerance = test_tolerance)
})
test_that("Output Run 2 is correct", {
  expect_equal(output2$betaHat, as.matrix(2.086334),
               tolerance = test_tolerance)
  expect_equal(output2$aVarHat, as.matrix(0.002414187),
               tolerance = test_tolerance)
})