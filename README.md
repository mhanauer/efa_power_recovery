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
fa(test_dat, 1, rotate="varimax", cor = "poly", correct = FALSE)

```

Appendix Rcode 4: Power analysis for EFA
```{r}
set.seed(123)
n_sample = seq(from = 200, to = 300, by = 10)

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
