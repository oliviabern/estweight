
<!-- README.md is generated from README.Rmd. Please edit that file -->

# estweight

<!-- badges: start -->

<!-- badges: end -->

The goal of estweight is to estimate propensity weights for convenience
samples using information from a representative sample. The estimated
weights are incorporated into propensity-weighted generalized linear
models and analytic standard error estimates are reported. When using a
logistic propensity weight estimation method, standard errors are
derived using a simultaneous estimating equation approach to account for
uncertianty from the weight estimation process. Otherwise standard
errors are design-based.

Suppose you want to estimate a model η(Y) = Xβ where Y is the response,
X is a matrix containing covariates and η is a link function. However, Y
and X were only collected in a convenience sample. Fitting the model
directly in the convenience sample may lead to bias from the
unrepresentative sampling. Additionally, supppose you have a second
dataset that is representative of the target population you want to
obtain inference about, but it does not contain both Y and X. The
represerntative sample can be used to estimate sampling weights for the
convenience sample adn these derirved weights can be used to weight the
outcome model of interest. This package uses a representative sample to
estimate weights for the convenience sample and weight the outcome model
of interest using the convenience sample.

## Installation

You can install the development version of estweight from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("oliviabern/estweight")
```

## Example

This is a basic example which shows you how to use a representative
sample to estimate propensity weights and fit a weighted GLM for a
biased convenience sample. First we will estimate a sample drawn from
our reference population.

``` r
library(estweight)

set.seed(10403)

# simulate population
n.pop = 50000
X1 = rbinom(n.pop, size = 1, prob = 0.5)
X2 = rbinom(n.pop, size = 1, prob = .3)
X3 = rbinom(n.pop, size = 1, prob = .3)
X4 = rbinom(n.pop, size = 1, prob = .1)
X5 = rbinom(n.pop, size = 1, prob = .05)
X6 = rbinom(n.pop, size = 1, prob = .1)
X7 = rbinom(n.pop, size = 1, prob = .5)
X8 = rnorm(n.pop, mean = 65, sd = 5)

repsample = data.frame(X1, X2, X3, X4,
                       X5, X6, X7, X8)

# create helper function
expit = function(x){exp(x)/(1+exp(x))}
```

Next we generate the known sampling probabilities and corresponding true
propensity weight. We also generate the response Y. Our goal is to
estimate the outcome model of interest, η(Y) = β<sub>0</sub> +
β<sub>1</sub>X<sub>4</sub> + β<sub>2</sub>X<sub>5</sub> +
β<sub>3</sub>X<sub>6</sub>. In this simulation, we know Y for subjects
in both samples, but in practice we assume Y was only collected in the
convenience sample.

``` r
# Calculate probability of being oversampled and true HT weight
sampprob = expit(.5*(repsample$X1*.3 + repsample$X2*.4 + repsample$X3*.8 +
                       repsample$X4*2.7 + repsample$X5*.9 + (repsample$X3)*(repsample$X5)*2 +
                       repsample$X6*.1 + repsample$X6*repsample$X7*1.5 +
                       -.002*repsample$X8^2 + 8))
repsample$sampprob = sampprob
repsample$weight_true = (1/sampprob)/sum(1/sampprob)

# Calculation E(Y) where Y is the response
repsample$Ey = expit(1 + log(2)*(repsample$X4) + -log(3)*(repsample$X5) + log(2.5)*(repsample$X6) +
                        -log(2)*sampprob + log(4)*(repsample$X4)*sampprob +
                        log(4)*(repsample$X5)*sampprob + -log(3)*(repsample$X6)*sampprob)
```

Now, we draw a biased convenience sample where each subject is sampled
according to their sampling probability and a representative sample
where each subject has equal probability of being sampled. Each sample
contains n<sub>samp</sub> = 500 observations.

``` r
n.samp = 500

# draw biased sample
b.samp.id = sample(1:n.pop, n.samp, prob = repsample$sampprob)
b.samp = repsample[b.samp.id, ]

# draw representative sample
r.samp.id = sample(1:n.pop, n.samp)
r.samp = repsample[r.samp.id, ]
```

Next, we simulate response Y for each subject and create an indicator
biased which is 1 for subjects in the biased sample and 0 otherwise.

``` r
# simulate response
b.samp$y = rbinom(n.samp,1,b.samp$Ey)
r.samp$y = rbinom(n.samp,1,r.samp$Ey)

# biased sample indicator
b.samp$biased = 1; r.samp$biased = 0
```

We need to manipulate the data into the form required for the
convGLMfunction. We need to create a stacked dataset that combines the
biased and representative samples. They should be stacked vertically and
contain a unique subject identifier, ID}, in the first column, the
indicator, biased}, should be in the last column, and any covariates
collected in both datasets to be used for matching should be in the
middle columns. All strings and factors should be of class factor}.
Finally, we need to create a response matrix that contains 2 columns,
ID, and the response, y.

``` r
# Combine bisaed and representative samples
Xcomb = data.frame(ID = 1:(2*n.samp),
                   rbind(b.samp, r.samp))

# Dataframe to pass to function
Xfit = Xcomb[,c("ID", paste0("X",1:8), "biased")]

# convert strings and indicators to factors
convert2factor = paste0("X",1:7)
for(i in 1:length(convert2factor)){
  Xfit[,convert2factor[i]] = as.factor(Xfit[,convert2factor[i]])
}

# get response
response = Xcomb[Xcomb$biased==1,c("ID","y")]
```

Use the convGLM function to fit an propensity-weighted generalized
linear model that utilizes logistic regression estimated weights.

``` r
# formulas for final scientific model
form.outcome = as.formula(y ~ X4 + X5 + X6)

# fit weighted model using a logistic propensity weight estimation method
convGLM(data = Xfit, outcome_formula = form.outcome, response = response,
        weight_model = "logistic", outcome_family = "quasibinomial")
#>                  coef analyticSE
#> (Intercept) 0.3695812  0.1110329
#> X41         1.5897785  0.4284810
#> X51         1.1122235  0.4648974
#> X61         0.3366221  0.3162072
```

Use a random forest propensity weight estimation method
instead.

``` r
# fit weighted model using a random forest propensity weight estimation method
convGLM(data = Xfit, outcome_formula = form.outcome, response = response,
        weight_model = "randomForest", outcome_family = "quasibinomial")
#>                   coef analyticSE
#> (Intercept) 0.50048753  0.1477161
#> X41         1.70037738  0.4721582
#> X51         0.63636903  0.5422092
#> X61         0.03364452  0.3847607
```

Now try a covariate balancing propensity score weight estimation
method.

``` r
# fit weighted model using a random forest propensity weight estimation method
convGLM(data = Xfit, outcome_formula = form.outcome, response = response,
        weight_model = "CBPS", outcome_family = "quasibinomial")
#> [1] "Finding ATT with T=0 as the treatment.  Set ATT=1 to find ATT with T=1 as the treatment"
#>                  coef analyticSE
#> (Intercept) 0.3738329  0.1117208
#> X41         1.5960229  0.4198840
#> X51         1.0953650  0.4677181
#> X61         0.3098333  0.3172941
```

Finally, try a entropy balancing weight estimation
method.

``` r
# fit weighted model using a random forest propensity weight estimation method
convGLM(data = Xfit, outcome_formula = form.outcome, response = response,
        weight_model = "entbal", outcome_family = "quasibinomial")
#> final  value -0.049378 
#> converged
#>                  coef analyticSE
#> (Intercept) 0.3619440  0.1122514
#> X41         1.6440069  0.4228829
#> X51         1.1122146  0.4694176
#> X61         0.3191363  0.3193139
```
