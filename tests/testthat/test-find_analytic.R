
test_that("find_analytic calculates correct values", {
  Av <- 0.5
  Af <- 2.0
  Ab <- 2.0
  Mb <- 3.0
  
  result <- find_analytic(Av, Af, Ab, Mb)
  
  expect_length(result, 4)  # Check that the result has 3 elements
  expect_true(all(is.numeric(result)))  # Check that all elements are numeric
  expect_equal(result[1], 3, tolerance = 1e-6)  # Check the first element (G)
  expect_equal(result[2], 2.309401, tolerance = 1e-6)  # Check the second element (Ac)
  expect_equal(result[4], 0.536360, tolerance = 1e-6)  # Check the third element (p)
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


    res = var_cycling(x, var_mat, pb_mat, c(pi, pi/3), time_units = "rad")

    expect_equal(as.numeric(res$double_freq_acrophase), c(pi/2, 2.61799387), tolerance = 1e-6)  
    expect_equal(as.numeric(res$double_freq_amplitude), c(3, 1))  
    expect_equal(as.numeric(res$fund_amp), c(4, 3))  
    expect_equal(as.numeric(res$fund_acro), c(pi, 1.047197551), tolerance = 1e-6)  


})

test_that("var_cycling_fix_acro calculates correct cycling params", {
  x = seq(0,22, by = 2) #hrs 0-22 by 2
  pb1 = 4*cos(pi/12*x-pi)+50 #bulk rhythm1, amp 4, mesor 6, acrophase 12hrs
  pb2 = 3*cos(pi/12*x-pi/3)+50 #bulk rhythm2, amp 3, mesor 7, acrophase 4hrs

  var1 = 3*cos(2*((pi/12)*x-pi/2-pi))+pb1 #double frq rhythm1, acro 6hrs, amp 3
  var2 = cos(2*((pi/12)*x-pi/2-pi/3))+pb2 #double frq rhythm2, acro 10hrs, amp 1

  pb_mat = rbind(pb1, pb2)
  var_mat = rbind(var1, var2)

  plot(x, var1)
  plot(x, var2)

  res = var_cycling_fix_acro(x, var_mat, pb_mat, c(pi, 4*pi/12))

  expect_length(res, 9)
  expect_equal(as.numeric(res$double_freq_amplitude), c(3, 1), tolerance = 1e-6)  
  expect_equal(as.numeric(res$fund_amp), c(4, 3), tolerance = 1e-6)  

})

test_that("acrodiffs calculation works as intended", {
  
  x = seq(0,2*pi, by = 0.1)
  #pseudobulk rhythms
  pb1 = 4*cos(x-0)+6           #bulk rhythm peaks at x = 0
  pb2 = 3*cos(x-pi/3)+7        #bulk rhythm peaks at x = pi/3
  pb3 = 2*cos(x-0)+5           #bulk rhythm peaks at x = 0
  pb4 = 4*cos(x-3*pi/2)+3      #bulk rhythm peaks at x = 3*pi/2


  #double frequency rhythms
  #double freq phase offset:
  dfpo = 0.1
  var1 = 3*cos(2*(x-dfpo-0))+pb1   #peaks double-freq offset after pb wave
  var2 = cos(2*(x-dfpo-pi/3))+pb2
  var3 = cos(2*(x-dfpo-0))+pb3
  var4 = cos(2*(x-dfpo-3*pi/2))+pb4

  pb_mat = rbind(pb1, pb2, pb3, pb4)
  var_mat = rbind(var1, var2, var3, var4)
  rownames(pb_mat) = rownames(var_mat) = c("Gene1", "Gene2", "Gene3", "Gene4")

  pb_cycling = is_cycling(rad, pb_mat, time_units = "rad")
  pb_phases = as.numeric(unname(unlist(pb_cycling$acrophase_bulk)))
  var_cyclers = var_cycling(x, var_mat, pb_mat,pb_phases , time_units = "rad")

  df = var_cyclers %>% 
    left_join(pb_cycling, by = "Gene_Symbols") %>%
    mutate(across(-Gene_Symbols, as.numeric)) %>%
    mutate(acrodiff = (double_freq_acrophase - acrophase_bulk) %% pi)
  df$acrodiff

  expect_equal(as.numeric(df$acrodiff), rep(dfpo, 4), tolerance = 1e-6)  

})


#Commented out because it generates files
# test_that("metacycle code is runnable", {
#   x = seq(0,22, by = 2) #hrs 0-22 by 2
#   pb1 = 4*cos(pi/12*x-pi)+50 #bulk rhythm1, amp 4, mesor 6, acrophase 12hrs
#   pb2 = 3*cos(pi/12*x-pi/3)+50 #bulk rhythm2, amp 3, mesor 7, acrophase 4hrs

#   var1 = 3*cos(2*((pi/12)*x-pi/2-pi))+pb1 #double frq rhythm1, acro 6hrs, amp 3
#   var2 = cos(2*((pi/12)*x-pi/2-pi/3))+pb2 #double frq rhythm2, acro 10hrs, amp 1

#   pb_mat = rbind(pb1, pb2)
#   var_mat = rbind(var1, var2)

#   plot(x, var1)
#   plot(x, var2)
#   circ_temp = tempfile(fileext = ".txt")
#   mean_temp = tempfile(fileext = ".txt")
#   pb_temp = tempfile(fileext = ".txt")
#   dir = tempdir()
#   res = cycling_meta2d(x, var_mat, pb_mat, out_file_circ = circ_temp,
#                      out_file_mean = mean_temp, 
#                      out_file_pb = pb_temp,
#                      out_dir = dir )
#   expect_type(res, "list")
#   expect_length(res, 3)

# })