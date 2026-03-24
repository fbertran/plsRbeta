test_that("predict.plsRbetamodel handles new beta data from the default interface", {
  set.seed(123)
  n <- 90
  p <- 4
  dat <- as.data.frame(t(replicate(n, simul_data_UniYX_beta(totdim = p, ncomp = 2, disp = 5, phi0 = 50))))
  dat <- subset(dat, is.finite(Ybeta) & Ybeta > 0 & Ybeta < 1)

  train <- dat[1:60, , drop = FALSE]
  test <- dat[61:75, , drop = FALSE]
  xcols <- paste0("X", 1:p)

  fit <- plsRbeta(
    object = train$Ybeta,
    dataX = as.matrix(train[xcols]),
    dataPredictY = as.matrix(test[xcols]),
    nt = 2,
    modele = "pls-beta",
    verbose = FALSE
  )

  pred <- predict(fit, newdata = as.matrix(test[xcols]), verbose = FALSE)
  scores <- predict(fit, newdata = as.matrix(test[xcols]), type = "scores", verbose = FALSE)

  expect_length(pred, nrow(test))
  expect_equal(nrow(scores), nrow(test))
  expect_equal(as.numeric(pred), as.numeric(fit$ValsPredictY))
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict.plsRbetamodel handles new beta data from the formula interface", {
  skip_if_not_installed("betareg")
  data("GasolineYield", package = "betareg")

  fit <- plsRbeta(yield ~ ., data = GasolineYield, nt = 2, modele = "pls-beta", verbose = FALSE)
  newdata <- GasolineYield[1:5, -1, drop = FALSE]

  pred <- predict(fit, newdata = newdata, verbose = FALSE)
  scores <- predict(fit, newdata = newdata, type = "scores", verbose = FALSE)
  manual <- predict(fit$FinalModel, newdata = data.frame(tt = I(scores)), type = "response")

  expect_length(pred, nrow(newdata))
  expect_equal(as.numeric(pred), as.numeric(manual))
  expect_equal(as.numeric(predict(fit, verbose = FALSE)), as.numeric(fit$YChapeau))
  expect_true(all(pred > 0 & pred < 1))
})
