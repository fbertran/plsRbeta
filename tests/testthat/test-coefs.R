test_that("coefs.plsRbeta* return numeric coefficients on simulated data", {
    set.seed(101)
    n <- 120
    p <- 4
    dat <- as.data.frame(t(replicate(n, simul_data_UniYX_beta(totdim = p, ncomp = 2, disp = 5, phi0 = 50))))
    dat <- subset(dat, is.finite(Ybeta) & Ybeta > 0 & Ybeta < 1)
    ind <- seq_len(nrow(dat))
    
    co_std <- coefs.plsRbeta(
      dataset = dat, ind = ind, nt = 2, modele = "pls-beta",
      method = "logistic", maxcoefvalues = 1e6, ifbootfail = rep(NA_real_, p + 1),
      verbose = FALSE
    )
    expect_true(is.numeric(co_std))
    expect_true(all(is.finite(co_std) | is.na(co_std)))
    expect_true(length(co_std) %in% c(p, p + 1))
    
    co_raw <- coefs.plsRbeta.raw(
      dataset = dat, ind = ind, nt = 2, modele = "pls-beta",
      method = "logistic", maxcoefvalues = 1e6, ifbootfail = rep(NA_real_, p + 1),
      verbose = FALSE
    )
    expect_true(is.numeric(co_raw))
    expect_true(all(is.finite(co_raw) | is.na(co_raw)))
    expect_true(length(co_raw) %in% c(p, p + 1))
  })
