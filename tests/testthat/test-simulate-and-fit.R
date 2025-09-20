set.seed(123)
  n <- 80
  p <- 4
  dat <- as.data.frame(t(replicate(n, simul_data_UniYX_beta(totdim = p, ncomp = 2, disp = 5, phi0 = 50))))
  dat <- subset(dat, is.finite(Ybeta) & Ybeta > 0 & Ybeta < 1)
  # ensure columns are the expected names
  expect_true(all(c("Ybeta", paste0("X", 1:p)) %in% colnames(dat)))
  
  test_that("plsRbeta default interface fits a simple beta model", {
    fit <- plsRbeta(object = dat$Ybeta, dataX = as.matrix(dat[paste0("X", 1:p)]), nt = 2, modele = "pls-beta", verbose = FALSE)
    expect_s3_class(fit, "plsRbetamodel")
    expect_true(is.numeric(fit$YChapeau))
    expect_true(all(is.finite(fit$YChapeau)))
    expect_true(all(fit$YChapeau > 0 & fit$YChapeau < 1))
    expect_true(is.matrix(fit$Coeffs) || is.numeric(fit$Coeffs))
    expect_true(is.numeric(fit$computed_nt))
    expect_equal(as.integer(fit$computed_nt), 2L)
    expect_s3_class(summary(fit), "summary.plsRbetamodel")
    expect_invisible(print(fit))
    expect_invisible(print(summary(fit)))
  })
