test_that("bootplsbeta works with coefs.plsRbetanp statistic (lightweight)", {
    skip_on_cran()
    set.seed(202)
    n <- 80
    p <- 4
    dat <- as.data.frame(t(replicate(n, simul_data_UniYX_beta(totdim = p, ncomp = 2, disp = 5, phi0 = 50))))
    dat <- subset(dat, is.finite(Ybeta) & Ybeta > 0 & Ybeta < 1)
    fit <- plsRbeta(object = dat$Ybeta, dataX = as.matrix(dat[paste0("X", 1:p)]), nt = 2, modele = "pls-beta", verbose = FALSE)
    b <- bootplsbeta(fit, R = 10, typeboot = "fmodel_np", statistic = coefs.plsRbetanp, verbose = FALSE)
    expect_true(inherits(b, "boot"))
  })
