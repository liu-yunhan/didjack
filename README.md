# didjack
Cluster jackknife (CV3) inference for the Callaway–Sant’Anna DID estimator  <br>
This is a post-estimation program for [did (Brantly Callaway and Pedro H. C. Sant'Anna)
](https://bcallaway11.github.io/did/index.html)
## Installation 
install (or update) using the following command

```r
# install.packages("remotes")   # if needed
remotes::install_github("liu-yunhan/didjack")
library(didjack)
```

## Basic Usage
```r
library(did)
library(didjack)

# Example structure (replace df, variable names with yours)
gt <- att_gt(
  yname       = "log_wage_dm",
  tname       = "year",
  gname       = "first_treat",
  data        = df,
  clustervars = "state",
  panel       = FALSE,
  est_method  = "reg",
  bstrap      = FALSE
)

ag  <- aggte(gt)

res <- didjack(ag)          # default: reuse original cluster variable and confidence level = 0.95
summary(res)                # formatted CV3 results
```
