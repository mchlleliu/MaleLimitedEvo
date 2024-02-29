
#### make some fake data for testing.
{
  set.seed(123)
  n.fake.data.points <- 1000
  x <- rnorm(n.fake.data.points, sd = 3.5)
  y <- rbinom(n, 1, plogis(0.02*(2 + 3*x + 8*x^2 + 0.3*x^3 - 0.2*x^4)))  # Example logistic function with a fourth-order polynomial
  fake.data = data.frame(x = x, y = y)
}

library(lmtest)

# Fit fourth-order polynomial logistic regression
model4 <- glm(y ~ poly(x, degree = 4, raw = TRUE), data = fake.data, family = "binomial")
model3 <- glm(y ~ poly(x, degree = 3, raw = TRUE), data = fake.data, family = "binomial")
model2 <- glm(y ~ poly(x, degree = 2, raw = TRUE), data = fake.data, family = "binomial")
model1 <- glm(y ~ poly(x, degree = 1, raw = TRUE), data = fake.data, family = "binomial")
model0 <- glm(y ~ 1, data = fake.data, family = "binomial")

## the two lines below are just to (crudely) plot what the model predicts the probability is
f2<- function(x) predict(model4, newdata = data.frame(x = x), type = "response")
plot(seq(from = -10, to = 15, by = 0.5), f2(seq(from = -10, to = 15, by = 0.5)))



AIC.vec = c(AIC(model0), AIC(model1),AIC(model2),AIC(model3),AIC(model4))

# requirement "0" for Twin Peaks
#Does model 4 have the best AIC?
TP.requirement0 = AIC.vec[5] == min(AIC.vec)

# requirement 1 for Twin Peaks
# is quartic model significant over null model
lrtest(model0, model4)
TP.requirement1 = lrtest(model0, model4)$`Pr(>Chisq)`[2] < 0.05


# requirement 2 for Twin Peaks
# Is the coefficient negative for the 4th order negative?
TP.requirement2 = model4$coefficients[5] < 0
as.numeric(model$coefficients)

# Find roots of derivative of polynomial to locate maxima and minima
root.locations = polyroot( (1:4)*as.numeric(model4$coefficients)[2:5] )

abs( Im(root.locations) ) < 10^-13  # Are the roots real? check that imaginary part is very small
TP.requirement3a = all(  abs( Im(root.locations) ) < 10^-13  )

## use this function to check if values are all within the range of interest
## the default range is -10, 15
withinBounds<-function(x, lower = -10, upper = 15) (Re(x )> lower) & ( Re(x) < upper)

# requirement 3 for Twin Peaks
# real roots fall within bounds
withinBounds(polyroot( (1:4)*as.numeric(model$coefficients)[2:5] ))

TP.requirement3b = all( withinBounds(Re(root.locations)) )

AllRequirementsMet = all (c(TP.requirement0, TP.requirement1, TP.requirement2, TP.requirement3a, TP.requirement3b))
AllRequirementsMet



