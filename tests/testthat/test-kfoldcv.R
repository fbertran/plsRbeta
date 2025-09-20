test_that("k-fold CV runs and returns expected structure", {
    set.seed(456)
    n <- 90
    p <- 4
    dat <- as.data.frame(t(replicate(n, simul_data_UniYX_beta(totdim = p, ncomp = 2, disp = 5, phi0 = 50))))
    dat <- subset(dat, is.finite(Ybeta) & Ybeta > 0 & Ybeta < 1)
    cv <- PLS_beta_kfoldcv(dataY = dat$Ybeta, dataX = as.matrix(dat[paste0("X", 1:p)]),
                           nt = 2, K = 3, NK = 1, modele = "pls-beta", verbose = FALSE)
    expect_true(is.list(cv))
    expect_true(is.numeric(cv$call$nt) || is.numeric(cv[["nt"]]))
    expect_true(is.list(cv$results_kfolds) || is.list(cv[["folds"]]))
    if (!is.null(cv$folds)) {
      expect_equal(length(cv$folds), 3L)
    }
  })
