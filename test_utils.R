source("utilities.R")

context("mdy to days")

test_that("incorrect format is not accepted", {
    dates <- c("md10", "10d", "d.1", "who knows", "mdy23")
    sapply(dates, function(x) { 
        expect_error(mdy_to_day(x))
    })
})

ave_month <- 30.5
ave_year <- 366

test_that("months processed", {
    dates <- c("m1", "m10", "m100")
    res <- ave_month * c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("years processed", {
    dates <- c("y1", "y10", "y100")
    res <- ave_year * c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("days processed", {
    dates <- c("d1", "d10", "d100")
    res <- c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("dot is accepted", {
    dates <- c("d1.", "d10.", "d100.")
    res <- c(1, 10, 100)
    expect_equal(mdy_to_day(dates), res)
})

test_that("number dot number is accepted", {
    dates <- c("m1.1", "m10.1", "m100.1")
    res <- ave_month * c(1.1, 10.1, 100.1)
    expect_equal(mdy_to_day(dates), res)
})

test_that("can process mix of mdy", {
    dates <- c("d1", "d10", "d100")
    res <- c(1, 10, 100)
    dates <- c(dates, c("m1.1", "m10.1", "m100.1"))
    res <- c(res, ave_month * c(1.1, 10.1, 100.1))
    dates <- c(dates, c("y1", "y10", "y100"))
    res <- c(res, ave_year * c(1, 10, 100))
    expect_equal(mdy_to_day(dates), res)
})
