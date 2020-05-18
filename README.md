---
title: "Diss Power"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Appendix Rcode 1: Replication and power analysis for Cronbach's Alpha

Bonnett (2002) formula

$$ n_{0} = [8k/(k-1)][Z_{a/2}/ln(e_{1})]^2+2~~ (1)  $$
$$ e_{1} = (1-LL)/(1-UL)~~(2) $$

Replication of Bonnett (2002) results to ensure I am using formula correctly
```{r}
k = 4
Za = 1.96
e1 = 2
n1_first = (8*k)/(k-1)
n1_second = (Za/log(e1))^2
n1 = (n1_first*n1_second)+2
n1

```
Now test with percision of .7 to .9 as stated in dissertation
```{r}
k = 3
Za = 1.96
e1 = (1-.7)/(1-.9)
e1 = 2
n1_first = (8*k)/(k-1)
n1_second = (Za/log(e1))^2
n1 = (n1_first*n1_second)+2
round(n1,0)
```

Appendix Rcode 2: Create function based on sim.VSS to create polytomous item structure
```{r}
library(DescTools)
library(psych)
library(StatMeasures)

sim_fifth =  function (ncases = 1000, nvariables = 10, nfactors = 1, meanloading = 0.7) 
{
    weight = sqrt(1 - meanloading * meanloading)
    theta = matrix(rnorm(ncases * nfactors), nrow = ncases, ncol = nvariables)
    error = matrix(rnorm(ncases * nvariables), nrow = ncases, 
        ncol = nvariables)
    items = meanloading * theta + weight * error
     items <- apply(items,2, function(x){CutQ(x, breaks = quantile(x, seq(0, 1, by = 0.20)), 
    labels = c(1:5))})
    items = apply(items, 2, function(x){as.numeric(x)})
    return(items)
}


```
Poly does not work with more than eight response options 
```{r}
test_dat = sim_one_fifth(ncases = 500)
test_dat
fa(test_dat, 1, rotate="varimax", cor = "poly", correct = 0)

```
Appendix Rcode 3: Parameter recovery analysis
```{r}
set.seed(123)
n_sample = seq(from = 400, to = 500, by = 20)
efa_recover= function(){
n_sample = n_sample
dat_vss = list()
dat_out = list()
fa_vss = list()
fa_vss_mean = list()
for(i in 1:length(n_sample)){
dat_vss[[i]] = sim_decile(ncases=n_sample[[i]], nvariables=10, nfactors=1, meanloading=.7)
fa_vss[[i]] = fa(dat_vss[[i]], 1, rotate="varimax", cor = "cor")
fa_vss_mean[[i]] = mean(fa_vss[[i]]$loadings)
}
return(fa_vss_mean)
}
### grab each of them sum them then divide by respective n's
reps = 10000
power_recover = replicate(n = reps, efa_recover(), simplify = FALSE)
power_recover
### There are 6 data sets of varying n's in a row and then a replication of that.
power_recover_unlist = round(unlist(power_recover),3)
power_recover_matrix = matrix(power_recover_unlist, ncol = length(n_sample), byrow = TRUE)
colnames(power_recover_matrix) = c(n_sample)
power_recover_df = data.frame(power_recover_matrix) 
recover_mean = round(apply(power_recover_df, 2, mean),3)
recover_mean = data.frame(as.vector(recover_mean))
colnames(recover_mean) = "recover_mean"
recover_sd = round(apply(power_recover_df, 2, sd),3)
recover_sd = data.frame(as.vector(recover_sd))
colnames(recover_sd) = "recover_sd"
recover_sd
recover = data.frame(recover_mean, recover_sd)
rownames(recover) = n_sample
recover
write.csv(recover,"recover.csv")
```
Appendix Rcode 4: Power analysis for EFA
```{r}
set.seed(123)
n_sample = seq(from = 100, to = 200, by = 10)

efa_power= function(){

n_sample = n_sample
tli_out = list()
rmsea = list()
chi_squre_p = list()
dat_vss = list()
dat_out = list()
fa_vss = list()
for(i in 1:length(n_sample)){
dat_vss[[i]] = sim_fifth(ncases=n_sample[[i]], nvariables=10, nfactors=1, meanloading=.7)
fa_vss[[i]] = fa(dat_vss[[i]], 1, rotate="varimax", cor = "poly", correct = FALSE)
tli_out[[i]] = fa_vss[[i]]$TLI
rmsea[[i]] = fa_vss[[i]]$RMSEA[1]
chi_squre_p[[i]] = fa_vss[[i]]$PVAL 
dat_out[[i]] = list(tli_out = tli_out[[i]], rmsea = rmsea[[i]], chi_squre_p = chi_squre_p[[i]])
}
return(dat_out)
}
### grab each of them sum them then divide by respective n's
reps = 50
power_efa = replicate(n = reps, efa_power(), simplify = FALSE)
## First 3 tli's are from the first rep meaning they have different sample size.  There are 3 tli's, because there are 3 samples being tested
## the second set of 3 samples is from the second round.  Stack them.
power_efa_unlist = round(unlist(power_efa),3)
power_efa_matrix = matrix(power_efa_unlist, ncol = 3, byrow = TRUE)
### split data every third row by the number 
colnames(power_efa_matrix) = c("tli", "rmsea", "chi_p")
power_efa_df = data.frame(power_efa_matrix) 
power_efa_df$n = rep(n_sample, reps)
power_efa_df$tli = ifelse(power_efa_df$tli >= .95,1,0)
power_efa_df$rmsea = ifelse(power_efa_df$rmsea <= .05,1,0)
power_efa_df$chi_p = ifelse(power_efa_df$chi_p > .05,1,0)


power_efa_agg = aggregate(power_efa_df[,1:3], by=list(n=power_efa_df$n), FUN=sum)
# divide by number of reps for %
power_efa_agg[,2:4] =  round(power_efa_agg[,2:4]/reps,3)
power_efa_agg
write.csv(power_efa_agg,"power_efa_agg.csv", row.names = FALSE)
```



Appendix Rcode 5: Altered wmwpowd function
```{r}
wmwpowd_paired = function (n, m, distn, distm, sides = "two.sided", alpha = 0.05, 
    nsims = 10000) 
{
    dist1 <- distn
    dist2 <- distm
    n1 = n
    n2 = m
    if (is.numeric(n1) == F | is.numeric(n2) == F) {
        stop("n1 and n2 must be numeric")
    }
    if (is.character(dist1) == F | is.character(dist2) == F | 
        !(sub("[^A-z]+", "", dist1) %in% c("norm", 
            "beta", "cauchy", "f", "gamma", 
            "lnorm", "unif", "weibull", "exp", 
            "chisq", "t", "doublex")) | !(sub("[^A-z]+", 
        "", dist2) %in% c("norm", "beta", "cauchy", 
        "f", "gamma", "lnorm", "unif", 
        "weibull", "exp", "chisq", "t", 
        "doublex"))) {
        stop("distn and distm must be characters in the form of distribution(parmater1) or distribution(paramter1, parameter2).\n         See documentation for details.")
    }
    if (is.numeric(alpha) == F) {
        stop("alpha must be numeric")
    }
    if (is.numeric(nsims) == F) {
        stop("nsims must be numeric")
    }
    if (!(sides %in% c("less", "greater", "two.sided"))) {
        stop("sides must be character value of less, greater, or two.sided")
    }
    if (sides == "two.sided") {
        test_sides <- "Two-sided"
    }
    if (sides %in% c("less", "greater")) {
        test_sides <- "One-sided"
    }
    dist1_func_char <- paste("r", sub("\\(.*", "", 
        dist1), "(", n1, ",", sub(".*\\(", 
        "", dist1), sep = "")
    dist2_func_char <- paste("r", sub("\\(.*", "", 
        dist2), "(", n2, ",", sub(".*\\(", 
        "", dist2), sep = "")
    power_sim_func <- function() {
        wilcox.test(eval(parse(text = dist1_func_char)), eval(parse(text = dist2_func_char)), 
            paired = T, correct = F, alternative = sides, exact = T)$p.value
    }
    pval_vect <- replicate(nsims, power_sim_func())
    empirical_power <- round(sum(pval_vect < alpha)/length(pval_vect), 
        3)
    dist1_func_ppp <- paste("r", sub("\\(.*", "", 
        dist1), "(", 1e+07, ",", sub(".*\\(", 
        "", dist1), sep = "")
    dist2_func_ppp <- paste("r", sub("\\(.*", "", 
        dist2), "(", 1e+07, ",", sub(".*\\(", 
        "", dist2), sep = "")
    p <- round(sum(eval(parse(text = dist1_func_ppp)) < eval(parse(text = dist2_func_ppp)))/1e+07, 
        3)
    wmw_odds <- round(p/(1 - p), 3)
    cat("Supplied distribution 1: ", dist1, "; n = ", 
        n1, "\n", "Supplied distribution 2: ", dist2, 
        "; m = ", n2, "\n\n", "p: ", p, "\n", 
        "WMW odds: ", wmw_odds, "\n", "Number of simulated datasets: ", 
        nsims, "\n", test_sides, " WMW test (alpha = ", 
        alpha, ")\n\n", "Empirical power: ", empirical_power, 
        sep = "")
    output_list <- list(empirical_power = empirical_power, alpha = alpha, 
        test_sides = test_sides, p = p, wmw_odds = wmw_odds, 
        distn = dist1, distm = dist2, n = n1, m = n2)
}
```

Appendix Rcode 6: T-test and Wilcoxon rank sum test power analysis for study three
```{r}
## Knowledge
pwr.t.test(n = 8, d = 1, type = "paired", alternative = "greater")
pwr.t.test(n = 10, d = .9, type = "paired", alternative = "greater")
pwr.t.test(n = 11, d = .8, type = "paired", alternative = "greater")
pwr.t.test(n = 14, d = .7, type = "paired", alternative = "greater")

### Knowledge not normal
#1 effect size
library(wmwpow)
(8-4)/4
wmwpowd_paired(n = 15, m = 15, distn = "norm(8,4)", distm = "norm(4,4)", sides = "greater",
alpha = 0.05, nsims=10000)


(8-4.4)/4
wmwpowd_paired(n = 18, m = 18, distn = "norm(8,4)", distm = "norm(4.4,4)", sides = "greater",
alpha = 0.05, nsims=10000)

(8-4.8)/4
wmwpowd_paired(n = 22, m = 22, distn = "norm(8,4)", distm = "norm(4.8,4)", sides = "greater",
alpha = 0.05, nsims=10000)

(8-5.2)/4
wmwpowd_paired(n = 28, m = 28, distn = "norm(8,4)", distm = "norm(5.2,4)", sides = "greater",
alpha = 0.05, nsims=10000)


```


 

#########################
Extra for dissertation
#########################

Demonstration of item analysis
```{r}
library(sjPlot)
dat_decile =  sim_decile(ncases = 110)
head(dat_decile)
dat_decile = data.frame(dat_decile)
tab_itemscale(dat_decile)
```


