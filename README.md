Fixed or random?
================
Anna Quaglieri - Speed Lab meeting March 2019
13/12/2017

-   [Confounded design](#confounded-design)
    -   [Block as a known fixed effect](#block-as-a-known-fixed-effect)
    -   [What if we removed block 2?](#what-if-we-removed-block-2)
    -   [The estimate of the variance of the effect size depends on the residuals](#the-estimate-of-the-variance-of-the-effect-size-depends-on-the-residuals)
    -   [Block as random effect](#block-as-random-effect)
-   [Compare fixed random models](#compare-fixed-random-models)
    -   [Distribution of b\_j at changing sigma2\_b](#distribution-of-b_j-at-changing-sigma2_b)
    -   [sigma2\_b = sigma2\_e](#sigma2_b-sigma2_e)
    -   [sigma2\_b = 10\*sigma2\_e](#sigma2_b-10sigma2_e)
    -   [sigma2\_b = 100\*sigma2\_e](#sigma2_b-100sigma2_e)
    -   [Comments](#comments)
-   [Simulate four batches](#simulate-four-batches)
    -   [High disease effect and low between batches variability](#high-disease-effect-and-low-between-batches-variability)
    -   [Block as random effect](#block-as-random-effect-1)
    -   [High disease effect and high between batches variability](#high-disease-effect-and-high-between-batches-variability)
    -   [Low disease effect and high between batches variability](#low-disease-effect-and-high-between-batches-variability)
    -   [Low disease effect and low between batches variability](#low-disease-effect-and-low-between-batches-variability)
    -   [Comments](#comments-1)

Confounded design
=================

``` r
mu = 0
d = 8
b = 4
```

1.  *E*(*x*<sub>*i*</sub>)=*μ*
2.  *E*(*y*<sub>*i*</sub>)=*μ* + *d*
3.  *E*(*x*<sub>*i*</sub>′) = *μ* + *b*

``` r
set.seed(100)
x <- rnorm(25,mean = mu,sd = 1)
set.seed(100)
y <- rnorm(25,mean = mu + d,sd = 1.7)
set.seed(100)
x_1 <- rnorm(25,mean = mu + b,sd = 1.4)

value <- c(x,y,x_1)
x <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0),times=c(25,25,25))),
                blocks = rep(c(1,2),times=c(50,25)))

ggplot(x,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot() 
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

Block as a known fixed effect
-----------------------------

``` r
summary(lm(value ~ factor(Treatment) + factor(blocks),data=x))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7374 -0.6726 -0.0167  0.5851  3.7436 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.1082     0.1963   0.551    0.583    
    ## factor(Treatment)1   8.0757     0.2776  29.090   <2e-16 ***
    ## factor(blocks)2      4.0433     0.2776  14.565   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9815 on 72 degrees of freedom
    ## Multiple R-squared:  0.9216, Adjusted R-squared:  0.9194 
    ## F-statistic: 423.1 on 2 and 72 DF,  p-value: < 2.2e-16

-   *μ*

``` r
mean(x$value[x$Treatment==0 & x$blocks == 1])
```

    ## [1] 0.108172

-   *disease effect*

$d=\\bar{y} - \\bar{x}$

``` r
mean(x$value[x$Treatment==1 & x$blocks == 1]) - mean(x$value[x$Treatment==0 & x$blocks == 1])
```

    ## [1] 8.07572

-   *batch effect* $b=\\bar{x}' - \\bar{x}$

``` r
mean(x$value[x$Treatment==0 & x$blocks == 2]) - mean(x$value[x$Treatment==0 & x$blocks == 1])
```

    ## [1] 4.043269

What if we removed block 2?
---------------------------

``` r
summary(lm(value ~ factor(Treatment),data=x[x$blocks == 1,]))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment), data = x[x$blocks == 
    ##     1, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7374 -0.5943 -0.0151  0.5554  3.7436 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.1082     0.1960   0.552    0.584    
    ## factor(Treatment)1   8.0757     0.2773  29.128   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9802 on 48 degrees of freedom
    ## Multiple R-squared:  0.9465, Adjusted R-squared:  0.9453 
    ## F-statistic: 848.4 on 1 and 48 DF,  p-value: < 2.2e-16

-   The second batch does not provide any information in the estimate of the disease effect.
-   Is it still worth putting it in?

The estimate of the variance of the effect size depends on the residuals
------------------------------------------------------------------------

-   Estimate of the residual variance after fitting the model.

${\\sigma}^2 = \\frac{sum(residuals^2)}{n - estimated \\; betas}$

``` r
mod1 <- lm(value ~ factor(Treatment) + factor(blocks),data=x)

sigma <- sum(mod1$residuals^2)/(nrow(x)-3)
x2 <- x
x2[,2] <- as.numeric(as.character(x2[,2]))
x2[,3] <- as.numeric(as.character(x2[,3]))
X <- as.matrix(cbind(1,x2[,2:3]))

unscaled_var <- solve(t(X) %*% X)

sqrt(unscaled_var*sigma)
```

    ## Warning in sqrt(unscaled_var * sigma): NaNs produced

    ##                   1 Treatment    blocks
    ## 1         0.4389351       NaN       NaN
    ## Treatment       NaN 0.2776069 0.1962977
    ## blocks          NaN 0.1962977 0.2776069

``` r
summary(mod1)
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = x)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7374 -0.6726 -0.0167  0.5851  3.7436 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.1082     0.1963   0.551    0.583    
    ## factor(Treatment)1   8.0757     0.2776  29.090   <2e-16 ***
    ## factor(blocks)2      4.0433     0.2776  14.565   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9815 on 72 degrees of freedom
    ## Multiple R-squared:  0.9216, Adjusted R-squared:  0.9194 
    ## F-statistic: 423.1 on 2 and 72 DF,  p-value: < 2.2e-16

-   If we only use block 1

``` r
mod1 <- lm(value ~ factor(Treatment),data=x[x$blocks == 1,])
sigma <- sum(mod1$residuals^2)/(nrow(x[x$blocks == 1,])-2)
x2 <- x[x$blocks == 1,]
x2[,2] <- as.numeric(as.character(x2[,2]))
X <- as.matrix(cbind(1,x2[,2]))

unscaled_var <- solve(t(X) %*% X)

sqrt(unscaled_var*sigma)
```

    ## Warning in sqrt(unscaled_var * sigma): NaNs produced

    ##           [,1]      [,2]
    ## [1,] 0.1960459       NaN
    ## [2,]       NaN 0.2772508

``` r
summary(mod1)
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment), data = x[x$blocks == 
    ##     1, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7374 -0.5943 -0.0151  0.5554  3.7436 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.1082     0.1960   0.552    0.584    
    ## factor(Treatment)1   8.0757     0.2773  29.128   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9802 on 48 degrees of freedom
    ## Multiple R-squared:  0.9465, Adjusted R-squared:  0.9453 
    ## F-statistic: 848.4 on 1 and 48 DF,  p-value: < 2.2e-16

-   Which is equivalent to estimating:

*V*(*d**i**s**e**a**s**e* *e**f**f**e**c**t*)=*V*(*m**e**a**n* *i**n* *d**i**s**e**a**s**e*)−*V*(*m**e**a**n* *i**n* *c**o**n**t**r**o**l**s*)

$V(mean \\; in \\; disease) = \\frac{{\\sigma\\; disease}^2}{n}$

``` r
nc <- length(x$value[x$Treatment==0 & x$blocks == 0])
vc <- var(x$value[x$Treatment==0 & x$blocks == 0])/nc

nd <- length(x$value[x$Treatment==1 & x$blocks == 0])
vd <- var(x$value[x$Treatment==1 & x$blocks == 0])/nd

sqrt(vd + vc)
```

    ## [1] NA

-   What if we don't include *blocks*?

``` r
mod1 <- lm(value ~ factor(Treatment),data=x)

summary(mod1)
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment), data = x)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -3.044 -1.678  0.015  1.346  5.105 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          2.1298     0.2738   7.778 3.68e-11 ***
    ## factor(Treatment)1   6.0541     0.4743  12.764  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.936 on 73 degrees of freedom
    ## Multiple R-squared:  0.6906, Adjusted R-squared:  0.6863 
    ## F-statistic: 162.9 on 1 and 73 DF,  p-value: < 2.2e-16

Block as random effect
----------------------

``` r
x$blocks <- factor(x$blocks)
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = x)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: x
    ## 
    ## REML criterion at convergence: 216.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7701 -0.6786 -0.0164  0.5961  3.8142 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 8.1347   2.8521  
    ##  Residual             0.9633   0.9815  
    ## Number of obs: 75, groups:  blocks, 2
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)          2.1298     2.0215   1.054
    ## factor(Treatment)1   8.0662     0.2774  29.073
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.034

Compare fixed random models
===========================

-   Only group one contains both treatments

``` r
library(lme4)
library(ggplot2)
```

Distribution of b\_j at changing sigma2\_b
------------------------------------------

``` r
set.seed(100)
plot(density(rnorm(1000,mean = 0,sd = sqrt(4))),xlab="b_j",ylab="density",main="1k Random sampling from b_j distribution",xlim=c(-20,20))
lines(density(rnorm(1000,mean = 0,sd = sqrt(40))),col=2)
lines(density(rnorm(1000,mean = 0,sd = sqrt(400))),col=4)
legend("topleft",legend=c("sigma2_b = 4","sigma2_b = 40","sigma2_b = 400"),pch = rep(3,16),col=c(1,2,4))
```

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)

sigma2\_b = sigma2\_e
---------------------

``` r
mu <- 0 # mean in controls in batch1
d <- 8 # disease effect
sigma2_b <- 4 
sigma2_e <- 4
```

-   Random *b*<sub>*j**s*</sub>. The random intercepts below are quite close to each other, meaning that the batch effect will be quite small. The random effects could be increased by increasing *s**i**g**m**a*<sub>*b*</sub>

``` r
set.seed(100)
b <- rnorm(2,mean = 0,sd = sqrt(sigma2_b))

b[1]
```

    ## [1] -1.004385

``` r
b[2]
```

    ## [1] 0.2630623

``` r
set.seed(100)
x <- rnorm(25,mean = mu + b[1],sd = sqrt(sigma2_e))
set.seed(100)
y <- rnorm(25,mean = mu + d + b[1],sd = sqrt(sigma2_e))
set.seed(3546)
x_1 <- rnorm(25,mean = mu + b[2],sd = sqrt(sigma2_e))

value <- c(x,y,x_1)
datax <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0),times=c(25,25,25))),
                blocks = rep(c(1,2),times=c(50,25)))

ggplot(datax,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot()
```

![](README_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
datax$Treatment <- factor(datax$Treatment) 
datax$blocks <- factor(datax$blocks) 
table(datax$Treatment,datax$blocks)
```

    ##    
    ##      1  2
    ##   0 25 25
    ##   1 25  0

``` r
summary(lm(value ~ factor(Treatment) + factor(blocks),data=datax))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = datax)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3889 -1.0925 -0.0366  1.1621  4.4042 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -0.7880     0.3120  -2.526  0.01373 *  
    ## factor(Treatment)1   8.0000     0.4412  18.133  < 2e-16 ***
    ## factor(blocks)2      1.3820     0.4412   3.133  0.00251 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.56 on 72 degrees of freedom
    ## Multiple R-squared:  0.8392, Adjusted R-squared:  0.8347 
    ## F-statistic: 187.9 on 2 and 72 DF,  p-value: < 2.2e-16

$\\bar{y} - \\bar{x}$

``` r
mean(datax$value[datax$Treatment == 1 & datax$blocks == 1]) - mean(datax$value[datax$Treatment == 0 & datax$blocks == 1])
```

    ## [1] 8

``` r
datax$blocks <- factor(datax$blocks)
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = datax)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: datax
    ## 
    ## REML criterion at convergence: 281.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.12751 -0.69515 -0.03386  0.74856  2.82362 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 0.8577   0.9261  
    ##  Residual             2.4329   1.5598  
    ## Number of obs: 75, groups:  blocks, 2
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)        -0.09702    0.69102   -0.14
    ## factor(Treatment)1  7.92958    0.43552   18.21
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.162

sigma2\_b = 10\*sigma2\_e
-------------------------

``` r
mu <- 0 # mean in controls in batch1
d <- 8 # disease effect
sigma2_e <- 4
sigma2_b <- 10*sigma2_e
```

``` r
set.seed(100)
b <- rnorm(2,mean = 0,sd = sqrt(sigma2_b))

b[1]
```

    ## [1] -3.176143

``` r
b[2]
```

    ## [1] 0.8318761

``` r
set.seed(100)
x <- rnorm(25,mean = mu + b[1],sd = sqrt(sigma2_e))
set.seed(100)
y <- rnorm(25,mean = mu + d + b[1],sd = sqrt(sigma2_e))
set.seed(3546)
x_1 <- rnorm(25,mean = mu + b[2],sd = sqrt(sigma2_e))

value <- c(x,y,x_1)
datax <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0),times=c(25,25,25))),
                blocks = rep(c(1,2),times=c(50,25)))

ggplot(datax,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot()
```

![](README_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
datax$Treatment <- factor(datax$Treatment) 
datax$blocks <- factor(datax$blocks) 
table(datax$Treatment,datax$blocks)
```

    ##    
    ##      1  2
    ##   0 25 25
    ##   1 25  0

``` r
summary(lm(value ~ factor(Treatment) + factor(blocks),data=datax))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = datax)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3889 -1.0925 -0.0366  1.1621  4.4042 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -2.9598     0.3120  -9.488 2.60e-14 ***
    ## factor(Treatment)1   8.0000     0.4412  18.133  < 2e-16 ***
    ## factor(blocks)2      4.1226     0.4412   9.345 4.78e-14 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.56 on 72 degrees of freedom
    ## Multiple R-squared:  0.8204, Adjusted R-squared:  0.8154 
    ## F-statistic: 164.5 on 2 and 72 DF,  p-value: < 2.2e-16

``` r
datax$blocks <- factor(datax$blocks)
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = datax)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: datax
    ## 
    ## REML criterion at convergence: 283.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.15752 -0.70800 -0.03821  0.74856  2.82362 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 8.401    2.898   
    ##  Residual             2.433    1.560   
    ## Number of obs: 75, groups:  blocks, 2
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)         -0.8985     2.0613  -0.436
    ## factor(Treatment)1   7.9764     0.4405  18.106
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.054

sigma2\_b = 100\*sigma2\_e
--------------------------

``` r
mu <- 0 # mean in controls in batch1
d <- 8 # disease effect
sigma2_e <- 4 
sigma2_b <- 100*sigma2_e
```

-   Distribution of batch effects

``` r
set.seed(100)
b <- rnorm(2,mean = 0,sd = sqrt(sigma2_b))

b[1]
```

    ## [1] -10.04385

``` r
b[2]
```

    ## [1] 2.630623

``` r
set.seed(100)
x <- rnorm(25,mean = mu + b[1],sd = sqrt(sigma2_e))
set.seed(100)
y <- rnorm(25,mean = mu + d + b[1],sd = sqrt(sigma2_e))
set.seed(3546)
x_1 <- rnorm(25,mean = mu + b[2],sd = sqrt(sigma2_e))

value <- c(x,y,x_1)
datax <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0),times=c(25,25,25))),
                blocks = rep(c(1,2),times=c(50,25)))

ggplot(datax,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot()
```

![](README_files/figure-markdown_github/unnamed-chunk-30-1.png)

``` r
datax$Treatment <- factor(datax$Treatment) 
datax$blocks <- factor(datax$blocks) 
table(datax$Treatment,datax$blocks)
```

    ##    
    ##      1  2
    ##   0 25 25
    ##   1 25  0

``` r
summary(lm(value ~ factor(Treatment) + factor(blocks),data=datax))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = datax)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3889 -1.0925 -0.0366  1.1621  4.4042 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -9.8275     0.3120  -31.50   <2e-16 ***
    ## factor(Treatment)1   8.0000     0.4412   18.13   <2e-16 ***
    ## factor(blocks)2     12.7891     0.4412   28.99   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.56 on 72 degrees of freedom
    ## Multiple R-squared:  0.9226, Adjusted R-squared:  0.9204 
    ## F-statistic:   429 on 2 and 72 DF,  p-value: < 2.2e-16

``` r
datax$blocks <- factor(datax$blocks)
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = datax)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: datax
    ## 
    ## REML criterion at convergence: 285.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.16778 -0.70287 -0.02833  0.74856  2.82362 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 81.680   9.038   
    ##  Residual              2.433   1.560   
    ## Number of obs: 75, groups:  blocks, 2
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)         -3.4330     6.3944  -0.537
    ## factor(Treatment)1   7.9924     0.4411  18.119
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.017

Comments
--------

-   In the examples above we increased *s**i**g**m**a*<sub>*b*</sub><sup>2</sup> to allow for higher variability between batches. It is interesting to observe that, as the variability between batches tends to infinity, the effect size estimated by the mixed model approaches the fixed effect estimate.

-   This case is very limited in how accurately we the *σ*<sub>*b*</sub><sup>2</sup> and *σ*<sub>*e*</sub><sup>2</sup> are estimated since we only have 2 observations (2 blocks) from which we have to estimate both of them. Both blocks contain both *σ*<sub>*b*</sub><sup>2</sup> + *σ*<sub>*e*</sub><sup>2</sup>. That's why they are so different in the REML estimate. There are still a lot of issues in the way these two components are estimated (Terry's is trying to figure out the algebra). This could be also a problem in the way the actual std.error estimates of the estimators are estimated and we might doubt a fixed effects model with only two blocks.

-   **But, does this problem becomes better if we have two blocks but across thousands of genes? The assumption done in Limma-Dupl correlations is that the *ρ* is the same across genes. Does that increase the number that we use to make estimates?** - to be better understood maybe by

-   The use of fixed and random really depends on the context in which we are. Say, that in this case we only have two batches, that's all we got, then the fixed effect way is the closest approach that we can use. If we instead think that these two batches are two random observations out of a common distribution then the random effect is the design that you would think of. In the fixed effect we use one block as the baseline and we estimate the batch effect as a delta of difference in the other batches. However, with the random model we assign its own *b*<sub>*j*</sub> to every batch.

Simulate four batches
=====================

![](img/four-batches.png)

-   Identifiability issue raised by Jinjin. The random variable from every batch will have the same sum of *σ*<sub>*b*</sub><sup>2</sup> + *σ*<sub>*e**p**s**i**l**o**n*</sub><sup>2</sup> and that's why their estimation does not happen in the usual way as in the OLS.

High disease effect and low between batches variability
-------------------------------------------------------

``` r
mu <- 0 # mean in controls in batch1
d <- 30 # disease effect
sigma2_e <- 2 
sigma2_b <- 2

set.seed(100)
b <- rnorm(5,mean = 0,sd = sqrt(sigma2_b))

b[1]
```

    ## [1] -0.7102072

``` r
b[2]
```

    ## [1] 0.1860132

``` r
b[3]
```

    ## [1] -0.1116056

``` r
b[4]
```

    ## [1] 1.254103

``` r
b[5]
```

    ## [1] 0.1654224

``` r
set.seed(100)
x <- rnorm(25,mean = mu + b[1],sd = sqrt(sigma2_e))
set.seed(100)
y <- rnorm(25,mean = mu + d + b[1],sd = sqrt(sigma2_e))
set.seed(3546)
x_1 <- rnorm(25,mean = mu + b[2],sd = sqrt(sigma2_e))
set.seed(100)
x_2 <- rnorm(25,mean = mu + b[3],sd = sqrt(sigma2_e))
set.seed(100)
y_2 <- rnorm(25,mean = mu + d + b[3],sd = sqrt(sigma2_e))
set.seed(3546)
x_3 <- rnorm(25,mean = mu + b[4],sd = sqrt(sigma2_e))
set.seed(100)
y_3 <- rnorm(25,mean = mu + d + b[4],sd = sqrt(sigma2_e))
set.seed(100)
x_4 <- rnorm(25,mean = mu + b[5],sd = sqrt(sigma2_e))

value <- c(x,y,x_1,x_2,y_2,x_3,y_3,x_4)
datax <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0,0,1,0,1,0),times=rep(25,8))),
                blocks = rep(c(1,2,3,4,5),times=c(50,25,50,50,25)))

ggplot(datax,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot()
```

![](README_files/figure-markdown_github/unnamed-chunk-35-1.png)

``` r
mod1 <- lm(value ~ factor(Treatment) + factor(blocks),data=datax)

summary(mod1)
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = datax)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.39630 -0.76240 -0.03108  0.82316  3.12778 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -0.5437     0.1749  -3.108  0.00217 ** 
    ## factor(Treatment)1  29.9730     0.1749 171.330  < 2e-16 ***
    ## factor(blocks)2      0.9637     0.2766   3.484  0.00061 ***
    ## factor(blocks)3      0.5986     0.2143   2.794  0.00573 ** 
    ## factor(blocks)4      2.0048     0.2143   9.357  < 2e-16 ***
    ## factor(blocks)5      0.8621     0.2766   3.117  0.00211 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.071 on 194 degrees of freedom
    ## Multiple R-squared:  0.9947, Adjusted R-squared:  0.9946 
    ## F-statistic:  7348 on 5 and 194 DF,  p-value: < 2.2e-16

Block as random effect
----------------------

``` r
datax$blocks <- factor(datax$blocks)
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = datax)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: datax
    ## 
    ## REML criterion at convergence: 609.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.23137 -0.70820 -0.03701  0.75827  2.92730 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 0.5194   0.7207  
    ##  Residual             1.1472   1.0711  
    ## Number of obs: 200, groups:  blocks, 5
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)          0.3421     0.3362   1.018
    ## factor(Treatment)1  29.9715     0.1735 172.766
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.157

High disease effect and high between batches variability
--------------------------------------------------------

``` r
mu <- 0 # mean in controls in batch1
d <- 30 # disease effect
sigma2_e <- 2 
sigma2_b <- 100*sigma2_e

set.seed(100)
b <- rnorm(5,mean = 0,sd = sqrt(sigma2_b))

b[1]
```

    ## [1] -7.102072

``` r
b[2]
```

    ## [1] 1.860132

``` r
b[3]
```

    ## [1] -1.116056

``` r
b[4]
```

    ## [1] 12.54103

``` r
b[5]
```

    ## [1] 1.654224

``` r
set.seed(100)
x <- rnorm(25,mean = mu + b[1],sd = sqrt(sigma2_e))
set.seed(100)
y <- rnorm(25,mean = mu + d + b[1],sd = sqrt(sigma2_e))
set.seed(3546)
x_1 <- rnorm(25,mean = mu + b[2],sd = sqrt(sigma2_e))
set.seed(100)
x_2 <- rnorm(25,mean = mu + b[3],sd = sqrt(sigma2_e))
set.seed(100)
y_2 <- rnorm(25,mean = mu + d + b[3],sd = sqrt(sigma2_e))
set.seed(3546)
x_3 <- rnorm(25,mean = mu + b[4],sd = sqrt(sigma2_e))
set.seed(100)
y_3 <- rnorm(25,mean = mu + d + b[4],sd = sqrt(sigma2_e))
set.seed(100)
x_4 <- rnorm(25,mean = mu + b[5],sd = sqrt(sigma2_e))

value <- c(x,y,x_1,x_2,y_2,x_3,y_3,x_4)
datax <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0,0,1,0,1,0),times=rep(25,8))),
                blocks = rep(c(1,2,3,4,5),times=c(50,25,50,50,25)))

ggplot(datax,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot()
```

![](README_files/figure-markdown_github/unnamed-chunk-39-1.png)

``` r
summary(lm(value ~ factor(Treatment) + factor(blocks),data=datax))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = datax)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.39630 -0.76240 -0.03108  0.82316  3.12778 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -6.9356     0.1749  -39.65   <2e-16 ***
    ## factor(Treatment)1  29.9730     0.1749  171.33   <2e-16 ***
    ## factor(blocks)2      9.0297     0.2766   32.64   <2e-16 ***
    ## factor(blocks)3      5.9860     0.2143   27.94   <2e-16 ***
    ## factor(blocks)4     19.6836     0.2143   91.87   <2e-16 ***
    ## factor(blocks)5      8.7428     0.2766   31.61   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.071 on 194 degrees of freedom
    ## Multiple R-squared:  0.9957, Adjusted R-squared:  0.9956 
    ## F-statistic:  9049 on 5 and 194 DF,  p-value: < 2.2e-16

-   Resisuals after using fixed effects

``` r
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = datax)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: datax
    ## 
    ## REML criterion at convergence: 627.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.23653 -0.71435 -0.03076  0.77297  2.91854 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 50.940   7.137   
    ##  Residual              1.148   1.071   
    ## Number of obs: 200, groups:  blocks, 5
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)          1.7528     3.1933   0.549
    ## factor(Treatment)1  29.9729     0.1749 171.346
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.016

Low disease effect and high between batches variability
-------------------------------------------------------

``` r
mu <- 0 # mean in controls in batch1
d <- 4 # disease effect
sigma2_e <- 2 
sigma2_b <- 100*sigma2_e

set.seed(100)
b <- rnorm(5,mean = 0,sd = sqrt(sigma2_b))

b[1]
```

    ## [1] -7.102072

``` r
b[2]
```

    ## [1] 1.860132

``` r
b[3]
```

    ## [1] -1.116056

``` r
b[4]
```

    ## [1] 12.54103

``` r
b[5]
```

    ## [1] 1.654224

``` r
set.seed(100)
x <- rnorm(25,mean = mu + b[1],sd = sqrt(sigma2_e))
set.seed(100)
y <- rnorm(25,mean = mu + d + b[1],sd = sqrt(sigma2_e))
set.seed(3546)
x_1 <- rnorm(25,mean = mu + b[2],sd = sqrt(sigma2_e))
set.seed(100)
x_2 <- rnorm(25,mean = mu + b[3],sd = sqrt(sigma2_e))
set.seed(100)
y_2 <- rnorm(25,mean = mu + d + b[3],sd = sqrt(sigma2_e))
set.seed(3546)
x_3 <- rnorm(25,mean = mu + b[4],sd = sqrt(sigma2_e))
set.seed(100)
y_3 <- rnorm(25,mean = mu + d + b[4],sd = sqrt(sigma2_e))
set.seed(100)
x_4 <- rnorm(25,mean = mu + b[5],sd = sqrt(sigma2_e))

value <- c(x,y,x_1,x_2,y_2,x_3,y_3,x_4)
datax <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0,0,1,0,1,0),times=rep(25,8))),
                blocks = rep(c(1,2,3,4,5),times=c(50,25,50,50,25)))

ggplot(datax,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot()
```

![](README_files/figure-markdown_github/unnamed-chunk-43-1.png)

``` r
summary(lm(value ~ factor(Treatment) + factor(blocks),data=datax))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = datax)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.39630 -0.76240 -0.03108  0.82316  3.12778 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -6.9356     0.1749  -39.65   <2e-16 ***
    ## factor(Treatment)1   3.9730     0.1749   22.71   <2e-16 ***
    ## factor(blocks)2      9.0297     0.2766   32.64   <2e-16 ***
    ## factor(blocks)3      5.9860     0.2143   27.94   <2e-16 ***
    ## factor(blocks)4     19.6836     0.2143   91.87   <2e-16 ***
    ## factor(blocks)5      8.7428     0.2766   31.61   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.071 on 194 degrees of freedom
    ## Multiple R-squared:  0.9799, Adjusted R-squared:  0.9794 
    ## F-statistic:  1896 on 5 and 194 DF,  p-value: < 2.2e-16

``` r
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = datax)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: datax
    ## 
    ## REML criterion at convergence: 627.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.23653 -0.71435 -0.03076  0.77297  2.91854 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 50.944   7.137   
    ##  Residual              1.148   1.071   
    ## Number of obs: 200, groups:  blocks, 5
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)          1.7528     3.1934   0.549
    ## factor(Treatment)1   3.9729     0.1749  22.712
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.016

Low disease effect and low between batches variability
------------------------------------------------------

``` r
mu <- 0 # mean in controls in batch1
d <- 4 # disease effect
sigma2_e <- 2 
sigma2_b <- sigma2_e

set.seed(100)
b <- rnorm(5,mean = 0,sd = sqrt(sigma2_b))

b[1]
```

    ## [1] -0.7102072

``` r
b[2]
```

    ## [1] 0.1860132

``` r
b[3]
```

    ## [1] -0.1116056

``` r
b[4]
```

    ## [1] 1.254103

``` r
b[5]
```

    ## [1] 0.1654224

``` r
set.seed(100)
x <- rnorm(25,mean = mu + b[1],sd = sqrt(sigma2_e))
set.seed(100)
y <- rnorm(25,mean = mu + d + b[1],sd = sqrt(sigma2_e))
set.seed(3546)
x_1 <- rnorm(25,mean = mu + b[2],sd = sqrt(sigma2_e))
set.seed(100)
x_2 <- rnorm(25,mean = mu + b[3],sd = sqrt(sigma2_e))
set.seed(100)
y_2 <- rnorm(25,mean = mu + d + b[3],sd = sqrt(sigma2_e))
set.seed(3546)
x_3 <- rnorm(25,mean = mu + b[4],sd = sqrt(sigma2_e))
set.seed(100)
y_3 <- rnorm(25,mean = mu + d + b[4],sd = sqrt(sigma2_e))
set.seed(100)
x_4 <- rnorm(25,mean = mu + b[5],sd = sqrt(sigma2_e))

value <- c(x,y,x_1,x_2,y_2,x_3,y_3,x_4)
datax <- data.frame(value = value,
  Treatment = c(rep(c(0,1,0,0,1,0,1,0),times=rep(25,8))),
                blocks = rep(c(1,2,3,4,5),times=c(50,25,50,50,25)))

ggplot(datax,aes(x=factor(blocks),y=value,fill=factor(Treatment))) + geom_boxplot()
```

![](README_files/figure-markdown_github/unnamed-chunk-47-1.png)

``` r
summary(lm(value ~ factor(Treatment) + factor(blocks),data=datax))
```

    ## 
    ## Call:
    ## lm(formula = value ~ factor(Treatment) + factor(blocks), data = datax)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.39630 -0.76240 -0.03108  0.82316  3.12778 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         -0.5437     0.1749  -3.108  0.00217 ** 
    ## factor(Treatment)1   3.9730     0.1749  22.710  < 2e-16 ***
    ## factor(blocks)2      0.9637     0.2766   3.484  0.00061 ***
    ## factor(blocks)3      0.5986     0.2143   2.794  0.00573 ** 
    ## factor(blocks)4      2.0048     0.2143   9.357  < 2e-16 ***
    ## factor(blocks)5      0.8621     0.2766   3.117  0.00211 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.071 on 194 degrees of freedom
    ## Multiple R-squared:  0.7903, Adjusted R-squared:  0.7849 
    ## F-statistic: 146.3 on 5 and 194 DF,  p-value: < 2.2e-16

``` r
fm1 <- lmer(value ~ factor(Treatment) + (1 | blocks), data = datax)
summary(fm1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: value ~ factor(Treatment) + (1 | blocks)
    ##    Data: datax
    ## 
    ## REML criterion at convergence: 609.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.23137 -0.70820 -0.03701  0.75827  2.92730 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  blocks   (Intercept) 0.5194   0.7207  
    ##  Residual             1.1472   1.0711  
    ## Number of obs: 200, groups:  blocks, 5
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)          0.3421     0.3362   1.018
    ## factor(Treatment)1   3.9715     0.1735  22.893
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## fctr(Trtm)1 -0.157

Comments
--------

-   As *σ*<sub>*b*</sub><sup>2</sup> goes to infinity the mixed model estimate of the effect of interest approaches the fixed effect model, and in particular $\\bar{y} - \\bar{x}$.

-   In cases when the variability between batches ( *σ*<sub>*b*</sub><sup>2</sup> ) is not that high with respect to the residual variability within bacthes ( *σ*<sub>*ϵ*</sub> ) then the mixed effect models balances and adjusts better the batch effect and provide more precise estimates of the effect of interest.
