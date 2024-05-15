## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "##",
  R.options = list(width = 60)
)

## ----packages,warning=FALSE,message=FALSE-----------------
library(tstests)
library(tsdistributions)
library(tsgarch)
library(flextable)
library(xts)
data("spy")
data("arma_forecast")
data("garch_forecast")
spyr <- na.omit(diff(log(spy)))

## ----gmmtest----------------------------------------------
spec <- garch_modelspec(spyr, model = "egarch", order = c(2,1), constant = TRUE,
distribution = "jsu")
mod <- estimate(spec)
skewness <- dskewness("jsu", skew = coef(mod)["skew"], shape = coef(mod)["shape"])
kurtosis <- dkurtosis("jsu", skew = coef(mod)["skew"], shape = coef(mod)["shape"]) + 3
test <- gmm_test(residuals(mod, standardize = TRUE), lags = 2, skewness = skewness, 
                 kurtosis = kurtosis, conf_level = 0.95)
print(test)

## ---------------------------------------------------------
as_flextable(test, include.decision = FALSE, footnote.reference = TRUE)

## ----honglitest-------------------------------------------
spec <- garch_modelspec(spyr[1:1000], model = "egarch", order = c(1,1), constant = TRUE,
distribution = "jsu")
mod <- estimate(spec)
z <- residuals(mod, standardize = TRUE)
p <- pdist("jsu", z, mu = 0, sigma = 1, skew = coef(mod)["skew"], shape = coef(mod)["shape"])
as_flextable(hongli_test(as.numeric(p), lags = 4, conf_level = 0.95), include.decision = T, footnote.reference = TRUE)

## ---------------------------------------------------------
spec <- garch_modelspec(spyr, model = "egarch", order = c(2,1), constant = TRUE,
distribution = "jsu")
mod <- estimate(spec)
test <- nyblom_test(residuals(mod, standardize = TRUE), scores = estfun(mod), 
                    parameter_names = names(coef(mod)), 
                    parameter_symbols = mod$parmatrix[estimate == 1]$symbol)
as_flextable(test, use.symbols = TRUE, footnote.reference = TRUE, include.decision = TRUE)

## ---------------------------------------------------------
test <- signbias_test(residuals(mod), sigma = sigma(mod))
as_flextable(test, use.symbols = TRUE, footnote.reference = TRUE)

## ---------------------------------------------------------
p <- pdist('jsu', q = garch_forecast$actual, mu = garch_forecast$forecast,
sigma = garch_forecast$sigma, skew = garch_forecast$skew,
shape = garch_forecast$shape)
as_flextable(berkowitz_test(p))

## ---------------------------------------------------------
as_flextable(dac_test(arma_forecast$actual, arma_forecast$forecast))

## ----warning=FALSE,message=FALSE--------------------------
spyr <- na.omit(diff(log(spy)))
mod <- arima(as.numeric(spyr[2000:2500]), order = c(15,0,0), transform.pars = TRUE, include.mean = TRUE)
p <- predict(mod, n.ahead = 15)
test <- minzar_test(actual = as.numeric(spyr[2501:2515]), 
                    forecast = as.numeric(p$pred))
as_flextable(test, footnote.reference = T, digits = 2)

## ----warning=FALSE,message=FALSE--------------------------
# the pit data
x <- pdist("jsu", q = garch_forecast$actual, mu = garch_forecast$forecast,
sigma = garch_forecast$sigma, skew = garch_forecast$skew,
shape = garch_forecast$shape)
test <- shortfall_de_test(x, alpha = 0.05, lags = 4)
as_flextable(test, footnote.reference = T) |> width(j = 1:3,width = 1.2)

## ----warning=FALSE,message=FALSE--------------------------
q <- qdist("jsu", p = 0.05, mu = garch_forecast$forecast, sigma = garch_forecast$sigma,
skew = garch_forecast$skew, shape = garch_forecast$shape)
test <- var_cp_test(actual = garch_forecast$actual, forecast = q, alpha = 0.05)
as_flextable(test, footnote.reference = T) |> width(j = 1:3,width = 1.2)

