test_that("GasolineYield dataset can be modeled with plsRbeta", {
  skip_if_not_installed("betareg")
  data("GasolineYield", package = "betareg")
  # basic sanity of data
  expect_true(all(is.finite(GasolineYield$yield)))
  expect_true(all(GasolineYield$yield > 0 & GasolineYield$yield < 1))

  fit <- plsRbetamodel(yield ~ ., data = GasolineYield, nt = 2, modele = "pls-beta", verbose = FALSE)
  expect_s3_class(fit, "plsRbetamodel")
  expect_equal(as.integer(fit$computed_nt), 2L)
  expect_true(all(fit$YChapeau > 0 & fit$YChapeau < 1))

  # correlation should be positive on training (weak threshold to be robust)
  expect_gt(cor(GasolineYield$yield, as.numeric(fit$YChapeau)), 0.3)
})
