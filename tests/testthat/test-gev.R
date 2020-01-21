
test_that("fit GEV obs",{
  data("rain_max_obs", envir = environment())
  obs.formula.intercept = list( ~1, ~1, ~1)
  obs.links = list("identity", "exponential", "identity")
  obs.fit = gev_fit(formula = obs.formula.intercept, data = rain_max_obs, response = rain_max_obs$rain, rl_n = 100,
                    links = obs.links)
  obs.formula = list( ~ lon + lat + elev + lon:lat + I(lon^2) + I(lat^2),
                      ~ lon + lat + elev,
                      ~ lon + lat)
  obs.fit = gev_fit(formula = obs.formula, data = rain_max_obs, response = rain_max_obs$rain, rl_n = 100,
                    links = obs.links)
  expect_true(length(unlist(obs.fit$par)) == 14)
})

test_that("fit GEV mod",{
  data("rain_max_model", envir = environment())
  mod.formula = list( ~ lon + lat + elev + lon:lat + I(lon^2) + I(lat^2),
                      ~ lon + lat + elev,
                      ~ lon + lat)
  mod.links = list("identity", "exponential", "identity")
  mod.fit = gev_fit(formula = mod.formula, data = rain_max_model, response = rain_max_model$rain, rl_n = 100,
                    links = mod.links)
  expect_true(length(unlist(mod.fit$par)) == 14)
})


test_that("plot/print GEV model", {
  data("rain_max_obs", envir = environment())
  obs.formula = list( ~ lon + lat + elev + lon:lat + I(lon^2) + I(lat^2),
                      ~ lon + lat + elev,
                      ~ lon + lat)
  obs.links = list("identity", "exponential", "identity")
  obs.fit = gev_fit(formula = obs.formula, data = rain_max_obs, response = rain_max_obs$rain, rl_n = 100,
                    links = obs.links)
  plot.gev(obs.fit)
  print(obs.fit)
})




