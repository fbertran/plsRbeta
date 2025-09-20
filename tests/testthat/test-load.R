test_that("package loads and key exports exist", {
  expect_true("plsRbeta" %in% getNamespaceExports("plsRbeta"))
  expect_true("PLS_beta" %in% getNamespaceExports("plsRbeta"))
  expect_true("PLS_beta_kfoldcv_formula" %in% getNamespaceExports("plsRbeta"))
  expect_true("bootplsbeta" %in% getNamespaceExports("plsRbeta"))
  expect_true("simul_data_UniYX_beta" %in% getNamespaceExports("plsRbeta"))
})
