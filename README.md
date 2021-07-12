What Explains India’s Second Wave of COVID-19 Infections?
================

India experienced a rapid surge in COVID-19 infections during
March-April 2021 that overwhelmed the healthcare system. This paper
shows that the circulation of the Delta variant of the SARS-CoV-2 virus,
amplified by a super-spreader event, likely caused India’s second wave.
We analyze publicly available district-wise data of COVID-19 infections
in India over 2020 and 2021. We use a combination of time-series
analysis, epidemiological modeling, and parameter estimation methods.
The data shows that a super-spreader event (Kumbh Mela festival being
the possible candidate) is probable responsible for circulating the
highly transmissible B.1.617.2 (Delta) variant of the SARS-CoV-2 virus
that caused India’s second wave. Thus, the Indian experience serves as a
cautionary tale that calls for increased genomic sequencing to identify
variants of concern and to regulate super-spreader events, while
vaccination remains the only long-term solution.

The code of the analysis presented here in the R Markdown source file
[README.Rmd](/README.Rmd). Draft of paper is in
[COVID19-India.pdf](/COVID19-India.pdf). The data is in [Data](/Data)
folder. We obtained the data from
[www.covid19india.org](https://www.covid19india.org/).

**Figure 1. Plot (a) shows the number of daily new COVID-19 cases for
all states in India, smoothed via a 21-day moving average.**

![](README_files/figure-gfm/Figure%201%20(a)-1.png)<!-- -->

<!-- #### Tp% Box plot -->

**Figure 2(a) portrays
![T\_{p\\%}](https://latex.codecogs.com/png.latex?T_%7Bp%5C%25%7D "T_{p\%}")
across 394 districts during the first wave.**

![](README_files/figure-gfm/Figure%202(a)-1.png)<!-- -->

**Figure 2(b) portrays the same for the second wave.**
![](README_files/figure-gfm/Figure%202(b)-1.png)<!-- -->

## Origins of the COVID-19 Waves in India

<!-- ### Cross-covariance of time-series -->
<!-- ```{r} -->
<!-- District_list <- model_df %>% -->
<!--   distinct(District_ID) %>% -->
<!--   "$"(District_ID) -->
<!-- nDistrict <- length(District_list) -->
<!-- ``` -->
<!-- ```{r} -->
<!-- get_max_ccf_k <- function(x, y, max_k) { -->
<!-- cbind(ccf(x,y,max_k, plot = FALSE)$acf, ccf(x,y,max_k, plot = FALSE)$lag) %>% -->
<!--   as.data.frame() %>% -->
<!--   as_tibble() %>% -->
<!--   rename(acf = V1, lag = V2) %>% -->
<!--   mutate(acf = abs (acf)) %>% -->
<!--   slice(which.max(acf)) %>% -->
<!--   "$"(lag) -->
<!-- } -->
<!-- ``` -->
<!-- ```{r} -->
<!-- corMat <- matrix(nrow = nDistrict, ncol = nDistrict) -->
<!-- for (i in seq_along(District_list)){ -->
<!--   for (j in seq_along(District_list)){ -->
<!--     #print(paste(i,j)) -->
<!--     ci <- first_wave_df %>% -->
<!--       filter(District_ID == District_list[i]) %>% -->
<!--       "$"(delta) -->
<!--     cj <- first_wave_df %>% -->
<!--       filter(District_ID == District_list[j]) %>% -->
<!--       "$"(delta) -->
<!--     corMat[i,j] <- get_max_ccf_k(ci, cj, 30) -->
<!--   } -->
<!-- } -->
<!-- ``` -->
<!-- ```{r} -->
<!-- corMat_second <- matrix(nrow = nDistrict, ncol = nDistrict) -->
<!-- for (i in seq_along(District_list)){ -->
<!--   for (j in seq_along(District_list)){ -->
<!--     #print(paste(i,j)) -->
<!--     ci <- second_wave_df %>% -->
<!--       filter(District_ID == District_list[i]) %>% -->
<!--       "$"(delta) -->
<!--     cj <- second_wave_df %>% -->
<!--       filter(District_ID == District_list[j]) %>% -->
<!--       "$"(delta) -->
<!--     corMat_second[i,j] <- get_max_ccf_k(ci, cj, 30) -->
<!--   } -->
<!-- } -->
<!-- ``` -->
<!-- ```{r} -->
<!-- heatmap(corMat) -->
<!-- ``` -->
<!-- ```{r} -->
<!-- heatmap(corMat_second) -->
<!-- ``` -->

## Epidemiological Model to Explain India’s Infection Dynamics

We consider an epidemiological diffusion model, estimate its parameters
from the district-wise test results, and demonstrate lack of fit during
the second wave. To this end, consider the susceptible-infected-removed
(SIR) compartmental model, described by

![
\\begin{aligned}
  S^i(t+1) &= S^i(t) - \\beta^i\_{\\text{int}}I^i(t)S^i(t) -  \\beta^i\_{\\text{ext}}I^i\_{\\text{ext}}(t)S^i(t)\\\\
  I^i(t+1) &= I^i(t) + \\beta^i\_{\\text{int}}I^i(t)S^i(t) +  \\beta^i\_{\\text{ext}}I^i\_{\\text{ext}}(t)S^i(t) - \\gamma I^i(t),\\\\
  R^i(t+1) &= R^i(t) + \\gamma I^i(t),    
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%20%20S%5Ei%28t%2B1%29%20%26%3D%20S%5Ei%28t%29%20-%20%5Cbeta%5Ei_%7B%5Ctext%7Bint%7D%7DI%5Ei%28t%29S%5Ei%28t%29%20-%20%20%5Cbeta%5Ei_%7B%5Ctext%7Bext%7D%7DI%5Ei_%7B%5Ctext%7Bext%7D%7D%28t%29S%5Ei%28t%29%5C%5C%0A%20%20I%5Ei%28t%2B1%29%20%26%3D%20I%5Ei%28t%29%20%2B%20%5Cbeta%5Ei_%7B%5Ctext%7Bint%7D%7DI%5Ei%28t%29S%5Ei%28t%29%20%2B%20%20%5Cbeta%5Ei_%7B%5Ctext%7Bext%7D%7DI%5Ei_%7B%5Ctext%7Bext%7D%7D%28t%29S%5Ei%28t%29%20-%20%5Cgamma%20I%5Ei%28t%29%2C%5C%5C%0A%20%20R%5Ei%28t%2B1%29%20%26%3D%20R%5Ei%28t%29%20%2B%20%5Cgamma%20I%5Ei%28t%29%2C%20%20%20%20%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
  S^i(t+1) &= S^i(t) - \beta^i_{\text{int}}I^i(t)S^i(t) -  \beta^i_{\text{ext}}I^i_{\text{ext}}(t)S^i(t)\\
  I^i(t+1) &= I^i(t) + \beta^i_{\text{int}}I^i(t)S^i(t) +  \beta^i_{\text{ext}}I^i_{\text{ext}}(t)S^i(t) - \gamma I^i(t),\\
  R^i(t+1) &= R^i(t) + \gamma I^i(t),    
\end{aligned}
")

where

![I\_{\\text{ext}}^i(t) = \\frac{\\sum\_{j:(j,i) \\in G}I^j(t)P^j}{\\sum\_{j:(j,i)\\in G}P^j}.](https://latex.codecogs.com/png.latex?I_%7B%5Ctext%7Bext%7D%7D%5Ei%28t%29%20%3D%20%5Cfrac%7B%5Csum_%7Bj%3A%28j%2Ci%29%20%5Cin%20G%7DI%5Ej%28t%29P%5Ej%7D%7B%5Csum_%7Bj%3A%28j%2Ci%29%5Cin%20G%7DP%5Ej%7D. "I_{\text{ext}}^i(t) = \frac{\sum_{j:(j,i) \in G}I^j(t)P^j}{\sum_{j:(j,i)\in G}P^j}.")

Please see paper for details.

**Figure 5. Plots (a) and (b) capture the quantiles of ’s and ’s across
districts over 2-week time windows from 05/01/2020 to 05/08/2021. The
regression errors in explaining the emergence of daily new infections
with the estimated β’s is given in (c). Plot (d) shows the quantiles of
mobility variations across districts.** text

<img src="README_figs/README-unnamed-chunk-26-1.png" width="4800" />
