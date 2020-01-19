test_that("logistic regression (lr) model fit has same parameters as glm", {
  # simulated data
  b0 = -.45
  b1 = .1
  x = rnorm(1000, 4, 1)
  y = rbinom(n = 1000, size = 1, prob = 1/(1+exp(-(b0 + b1*x))))
  d = data.frame(y, x)
  sim_fit = lr(y ~ x, data = d)
  glm_sim_fit = glm(y~x, family=binomial)
  expect_true(all(
    round(glm_sim_fit$coefficients,3) == round(sim_fit$coefficients,3)
  ))
})

test_that("logistic regression (lr) predictors, fit to small set of crime data", {

  data(crime, envir = environment())
  crime_sub = crime[sample(1:nrow(crime), nrow(crime)/10),]
  crime_fit = lr(Arrest ~ dist_from_station + Domestic + Year, data = crime_sub)
  preds = predict(crime_fit)
  probs = predict(crime_fit, type="probs")
  vals = predict(crime_fit, type="vals")

  test_preds = predict(crime_fit, newdata=crime[sample(1:nrow(crime), nrow(crime)/100),])
  test_probs = predict(crime_fit, newdata=crime[sample(1:nrow(crime), nrow(crime)/100),], type="probs")
  test_vals = predict(crime_fit, newdata=crime[sample(1:nrow(crime), nrow(crime)/100),], type="vals")

  expect_true(all(c(
    preds == 0 | preds==1,
    probs <=1 & probs >=0,
    test_preds == 0 | test_preds==1,
    test_probs <=1 & test_probs >=0
  )))

})

test_that("logistic regression (lr) cross-validation and log score, ROC", {

  b0 = -.25
  b1 = .2
  x = rnorm(10000, 4, 1)
  y = rbinom(n = 10000, size = 1, prob = 1/(1+exp(-(b0 + b1*x))))
  d = data.frame(y, x)
  sim_fit = lr(y ~ x, data = d)

  c = cv.lr(sim_fit, leave_out=nrow(d)/100+1, metric = "all", verbose = FALSE, seed = 1)
  l = log_score(predict(sim_fit, type="probs"), d$y)
  r = roc.lr(sim_fit, plot = TRUE)
  rtest = roc.lr(sim_fit, newdata = data.frame(
    y = rbinom(n = 100, size = 1, prob = 1/(1+exp(-(b0 + b1*x)))),
    x = rnorm(10000, 4, 1)
  ))
  expect_true(all(c(
    r >= 0 & r <= 1,
    rtest >= 0 & rtest <= 1,
    length(c) == 3,
    l >= 0
  )))
})
