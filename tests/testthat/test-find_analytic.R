
test_that("find_analytic calculates correct values", {
  Av <- 0.5
  Mv <- 2.0
  Af <- 2.0
  
  result <- find_analytic(Av, Mv, Af)
  
  expect_length(result, 3)  # Check that the result has 3 elements
  expect_true(all(is.numeric(result)))  # Check that all elements are numeric
  expect_equal(result[1], 1.333333, tolerance = 1e-6)  # Check the first element (G)
  expect_equal(result[2], 2.309401, tolerance = 1e-6)  # Check the second element (Ac)
  expect_equal(result[3], 0.536360, tolerance = 1e-6)  # Check the third element (p)
})


test_that("is_cycling calculates correct cycling params", {
    set.seed(42)
    rad = seq(0,2*pi, by = .1)
    y1 = 4*cos(rad-pi/3)+20
    y2 = 1*cos(rad)+3
    y3 = rnorm(length(rad), 0, 1)

    pb_mat = rbind(y1, y2, y3)

    res = is_cycling(rad, pb_mat, time_units = "rad")
    expect_equal(as.numeric(res$amplitude_bulk[1:2]), c(4, 1))  
    expect_equal(as.numeric(res$acrophase_bulk[1:2]), c(pi/3, 2*pi))  
    expect_gt(as.numeric(res$p_statistic[3]), 0.05)

   
    hr = seq(0,24, by = 1)
    y1 = 4*cos(hr*pi/12-pi/3)+20
    y2 = 1*cos(hr*pi/12)+3
    y3 = rnorm(length(hr), 0, 1)

    pb_mat2 = rbind(y1, y2, y3)

    res2 = is_cycling(hr, pb_mat2)
    expect_equal(as.numeric(res2$amplitude_bulk[1:2]), c(4, 1))  
    expect_equal(as.numeric(res2$acrophase_bulk[1:2]), c(pi/3, 2*pi))  
    expect_gt(as.numeric(res2$p_statistic[3]), 0.05)

})


test_that("var_cycling calculates correct cycling params", {
    x = seq(0,2*pi, by = 0.1)
    pb1 = 4*cos(x-pi)+6
    pb2 = 3*cos(x-pi/3)+7

    var1 = 3*cos(2*(x-pi/2-pi))+pb1
    var2 = cos(2*(x-pi/2-pi/3))+pb2

    pb_mat = rbind(pb1, pb2)
    var_mat = rbind(var1, var2)


    res = var_cycling(x, var_mat, pb_mat, time_units = "rad")

    expect_equal(as.numeric(res$double_freq_acrophase), c(pi/2, 2.61799387), tolerance = 1e-6)  
    expect_equal(as.numeric(res$double_freq_amplitude), c(3, 1))  
    expect_equal(as.numeric(res$fund_amp), c(4, 3))  
    expect_equal(as.numeric(res$fund_acro), c(pi, 1.047197551), tolerance = 1e-6)  


})