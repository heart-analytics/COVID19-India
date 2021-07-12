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

The source code of the analysis presented here is in the R Markdown file
[README.Rmd](/README.Rmd). Draft of paper is in
[COVID19-India.pdf](/Docs/COVID19-India.pdf). The data is in
[Data](/Data) folder. We obtained the data from
[www.covid19india.org](https://www.covid19india.org/).

**Figure 1(a) shows the number of daily new COVID-19 cases for all
states in India, smoothed via a 21-day moving average.**

![](README_files/figure-gfm/Figure-1(a)-1.png)<!-- -->

**Figure 1(b) illustrates the number of districts with 5 or more
cumulative new COVID-19 cases since the start dates of the two waves.**

![](README_files/figure-gfm/Figure-1(b)-1.png)<!-- -->

<!-- #### Tp% Box plot -->

**Figure 2(a) portrays
![T\_{p\\%}](https://latex.codecogs.com/png.latex?T_%7Bp%5C%25%7D "T_{p\%}")
across 394 districts during the first wave.**

![](README_files/figure-gfm/Figure-2(a)-1.png)<!-- -->

**Figure 2(b) portrays
![T\_{p\\%}](https://latex.codecogs.com/png.latex?T_%7Bp%5C%25%7D "T_{p\%}")
for the second wave.**

![](README_files/figure-gfm/Figure-2(b)-1.png)<!-- -->

## Origins of the COVID-19 Waves in India

### Cross-covariance of time-series

Consider the sequence of new daily COVID-19 infections in district
![i](https://latex.codecogs.com/png.latex?i "i") on day
![t](https://latex.codecogs.com/png.latex?t "t") as
![\\rho^i\_t](https://latex.codecogs.com/png.latex?%5Crho%5Ei_t "\rho^i_t"),
where ![t](https://latex.codecogs.com/png.latex?t "t") is measured as
days since the start dates of each wave. Then, the cross-covariance
function over ![T](https://latex.codecogs.com/png.latex?T "T")
time-periods is described by

![
    \\mathcal{C}(i,j; \\kappa) := \\frac{1}{T}\\sum\_{t=1}^T \\left( \\rho^i\_t - \\langle \\rho^i\\rangle \\right) \\left( \\rho^j\_{t+\\kappa} - \\langle \\rho^j\\rangle \\right),
](https://latex.codecogs.com/png.latex?%0A%20%20%20%20%5Cmathcal%7BC%7D%28i%2Cj%3B%20%5Ckappa%29%20%3A%3D%20%5Cfrac%7B1%7D%7BT%7D%5Csum_%7Bt%3D1%7D%5ET%20%5Cleft%28%20%5Crho%5Ei_t%20-%20%5Clangle%20%5Crho%5Ei%5Crangle%20%5Cright%29%20%5Cleft%28%20%5Crho%5Ej_%7Bt%2B%5Ckappa%7D%20-%20%5Clangle%20%5Crho%5Ej%5Crangle%20%5Cright%29%2C%0A "
    \mathcal{C}(i,j; \kappa) := \frac{1}{T}\sum_{t=1}^T \left( \rho^i_t - \langle \rho^i\rangle \right) \left( \rho^j_{t+\kappa} - \langle \rho^j\rangle \right),
")

where ![\\kappa](https://latex.codecogs.com/png.latex?%5Ckappa "\kappa")
denotes the time-shift of one time-series with respect to the other in
calculating the covariance. Here,
![\\langle \\rho^i\\rangle](https://latex.codecogs.com/png.latex?%5Clangle%20%5Crho%5Ei%5Crangle "\langle \rho^i\rangle")
and
![\\langle \\rho^j\\rangle](https://latex.codecogs.com/png.latex?%5Clangle%20%5Crho%5Ej%5Crangle "\langle \rho^j\rangle")
compute the empirical means of the
![T](https://latex.codecogs.com/png.latex?T "T")-length time series
![\\rho^i](https://latex.codecogs.com/png.latex?%5Crho%5Ei "\rho^i") and
![\\rho^j](https://latex.codecogs.com/png.latex?%5Crho%5Ej "\rho^j"),
respectively. The value of
![\\kappa\_\\star(i, j)](https://latex.codecogs.com/png.latex?%5Ckappa_%5Cstar%28i%2C%20j%29 "\kappa_\star(i, j)")
at which the covariance is maximized between the daily new infections in
district ![i](https://latex.codecogs.com/png.latex?i "i") and that in
district ![j](https://latex.codecogs.com/png.latex?j "j") denotes the
number of days by which the infection pattern in district
![j](https://latex.codecogs.com/png.latex?j "j") roughly lags the
pattern in district ![i](https://latex.codecogs.com/png.latex?i "i"). We
vary ![\\kappa](https://latex.codecogs.com/png.latex?%5Ckappa "\kappa")
in
![\[-30, 30\]](https://latex.codecogs.com/png.latex?%5B-30%2C%2030%5D "[-30, 30]")
with ![T=92](https://latex.codecogs.com/png.latex?T%3D92 "T=92") days
for the emergence of the first wave and
![T=69](https://latex.codecogs.com/png.latex?T%3D69 "T=69") days for the
second wave to compute
![\\kappa\_\\star(i, j)](https://latex.codecogs.com/png.latex?%5Ckappa_%5Cstar%28i%2C%20j%29 "\kappa_\star(i, j)").

A time-lag between the patterns of district
![i](https://latex.codecogs.com/png.latex?i "i") and district
![j](https://latex.codecogs.com/png.latex?j "j") does not imply that
infected people from district
![i](https://latex.codecogs.com/png.latex?i "i") came in direct contact
with people in district ![j](https://latex.codecogs.com/png.latex?j "j")
to drive the spread of COVID-19. However, consistent positive values of
![\\kappa\_\\star(i, j)](https://latex.codecogs.com/png.latex?%5Ckappa_%5Cstar%28i%2C%20j%29 "\kappa_\star(i, j)")’s
for multiple ![j](https://latex.codecogs.com/png.latex?j "j")’s suggests
that district ![i](https://latex.codecogs.com/png.latex?i "i") is an
epicenter of the infection spread.

Figure 3 portrays a bubble plot of
![\\sum\_{j=1}^N\\kappa\_\\star(i, j)](https://latex.codecogs.com/png.latex?%5Csum_%7Bj%3D1%7D%5EN%5Ckappa_%5Cstar%28i%2C%20j%29 "\sum_{j=1}^N\kappa_\star(i, j)")
over those locations ![i](https://latex.codecogs.com/png.latex?i "i")
for which this sum is positive. A larger bubble indicates higher
likelihood of a location being a source of the infection spread. Figure
3 reveals important differences between the likely epicenters of the two
waves.

**Figure 3. Bubble plots of
![\\kappa\_\\star](https://latex.codecogs.com/png.latex?%5Ckappa_%5Cstar "\kappa_\star")
for districts that have
![\\kappa\_\\star &gt;0](https://latex.codecogs.com/png.latex?%5Ckappa_%5Cstar%20%3E0 "\kappa_\star >0")
for the first wave in (a) and the second wave in (b). The size of the
bubble indicates the magnitude of
![\\kappa\_\\star](https://latex.codecogs.com/png.latex?%5Ckappa_%5Cstar "\kappa_\star").
The dotted areas identify Maharashtra and the Ganges belt.**

<img src="Docs/Figure-3.png" width="2134" />

**Figure 4. Variation of average population of districts with cumulative
number of 5 or more infections since the start of the waves as in Figure
1.** Start dates of the two waves are 05/01/2020 and 03/01/2021.

![](README_files/figure-gfm/Figure-4-1.png)<!-- -->

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

![I\_{\\text{ext}}^i(t) = \\frac{\\sum\_{j:(j,i) \\in \\mathcal{G}}I^j(t)P^j}{\\sum\_{j:(j,i)\\in \\mathcal{G}}P^j}.](https://latex.codecogs.com/png.latex?I_%7B%5Ctext%7Bext%7D%7D%5Ei%28t%29%20%3D%20%5Cfrac%7B%5Csum_%7Bj%3A%28j%2Ci%29%20%5Cin%20%5Cmathcal%7BG%7D%7DI%5Ej%28t%29P%5Ej%7D%7B%5Csum_%7Bj%3A%28j%2Ci%29%5Cin%20%5Cmathcal%7BG%7D%7DP%5Ej%7D. "I_{\text{ext}}^i(t) = \frac{\sum_{j:(j,i) \in \mathcal{G}}I^j(t)P^j}{\sum_{j:(j,i)\in \mathcal{G}}P^j}.")

We use two sets of data–the cumulative fraction of COVID-19 cases
![Q^i](https://latex.codecogs.com/png.latex?Q%5Ei "Q^i") on
![t=1](https://latex.codecogs.com/png.latex?t%3D1 "t=1") in district
![i](https://latex.codecogs.com/png.latex?i "i") and the fraction of new
COVID-positive cases
![\\Delta^i(t)](https://latex.codecogs.com/png.latex?%5CDelta%5Ei%28t%29 "\Delta^i(t)")
in district ![i](https://latex.codecogs.com/png.latex?i "i") on days
![t=1,\\ldots,T](https://latex.codecogs.com/png.latex?t%3D1%2C%5Cldots%2CT "t=1,\ldots,T").
Specifically, ![Q^i](https://latex.codecogs.com/png.latex?Q%5Ei "Q^i")’s
yield

![
    S^i(1) = 1 - Q^i, \\ I^i(1) = (1-\\gamma) Q^i, \\ R^i(1) = \\gamma Q^i,
](https://latex.codecogs.com/png.latex?%0A%20%20%20%20S%5Ei%281%29%20%3D%201%20-%20Q%5Ei%2C%20%5C%20I%5Ei%281%29%20%3D%20%281-%5Cgamma%29%20Q%5Ei%2C%20%5C%20R%5Ei%281%29%20%3D%20%5Cgamma%20Q%5Ei%2C%0A "
    S^i(1) = 1 - Q^i, \ I^i(1) = (1-\gamma) Q^i, \ R^i(1) = \gamma Q^i,
")

that are then propagated using
![\\Delta](https://latex.codecogs.com/png.latex?%5CDelta "\Delta")’s via

![
\\begin{aligned}
    S^i(t+1) &= S^i(t) - \\Delta^i(t),
    \\\\
    I^i(t+1) &= I^i(t) + \\Delta^i(t) - \\gamma I^i(t), \\\\
    R^i(t+1) &= R^i(t) + \\gamma I^i(t).
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%20%20%20%20S%5Ei%28t%2B1%29%20%26%3D%20S%5Ei%28t%29%20-%20%5CDelta%5Ei%28t%29%2C%0A%20%20%20%20%5C%5C%0A%20%20%20%20I%5Ei%28t%2B1%29%20%26%3D%20I%5Ei%28t%29%20%2B%20%5CDelta%5Ei%28t%29%20-%20%5Cgamma%20I%5Ei%28t%29%2C%20%5C%5C%0A%20%20%20%20R%5Ei%28t%2B1%29%20%26%3D%20R%5Ei%28t%29%20%2B%20%5Cgamma%20I%5Ei%28t%29.%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
    S^i(t+1) &= S^i(t) - \Delta^i(t),
    \\
    I^i(t+1) &= I^i(t) + \Delta^i(t) - \gamma I^i(t), \\
    R^i(t+1) &= R^i(t) + \gamma I^i(t).
\end{aligned}
")

Regression:

![
\\begin{aligned}
    &\\varphi\\left({\\beta}\[1\], \\ldots, {\\beta}\[26\]\\right)
    \\\\
    &:= \\sum\_{i=1}^N \\sum\_{\\tau=1}^{26} \\sum\_{t=14\\tau - 13}^{14\\tau} \\left(\\Delta^i(t) - \\beta^{i}\_{\\textrm{int}}\[\\tau\] I^i(t) S^i(t) 
    -  \\beta^{i}\_{\\textrm{ext}}\[\\tau\] I^i\_{\\textrm{ext}}(t) S^i(t) \\right)^2
    \\\\
    & \\quad + \\lambda \\sum\_{\\tau=1}^{25}
    \\sum\_{i=1}^{N} 
    \\left\[ \\left({\\beta}^i\_\\textrm{int}\[\\tau+1\]) - {\\beta}^i\_\\textrm{int}\[\\tau\]) \\right)^2 
    + \\left({\\beta}^i\_\\textrm{ext}\[\\tau+1\]) - {\\beta}^i\_\\textrm{ext}\[\\tau\]) \\right)^2
    \\right\]
    \\\\
    & \\quad 
    + \\rho \\sum\_{\\tau=1}^{26} \\sum\_{i=1}^{N} \\beta\_\\textrm{ext}\[\\tau\]^2.
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%20%20%20%20%26%5Cvarphi%5Cleft%28%7B%5Cbeta%7D%5B1%5D%2C%20%5Cldots%2C%20%7B%5Cbeta%7D%5B26%5D%5Cright%29%0A%20%20%20%20%5C%5C%0A%20%20%20%20%26%3A%3D%20%5Csum_%7Bi%3D1%7D%5EN%20%5Csum_%7B%5Ctau%3D1%7D%5E%7B26%7D%20%5Csum_%7Bt%3D14%5Ctau%20-%2013%7D%5E%7B14%5Ctau%7D%20%5Cleft%28%5CDelta%5Ei%28t%29%20-%20%5Cbeta%5E%7Bi%7D_%7B%5Ctextrm%7Bint%7D%7D%5B%5Ctau%5D%20I%5Ei%28t%29%20S%5Ei%28t%29%20%0A%20%20%20%20-%20%20%5Cbeta%5E%7Bi%7D_%7B%5Ctextrm%7Bext%7D%7D%5B%5Ctau%5D%20I%5Ei_%7B%5Ctextrm%7Bext%7D%7D%28t%29%20S%5Ei%28t%29%20%5Cright%29%5E2%0A%20%20%20%20%5C%5C%0A%20%20%20%20%26%20%5Cquad%20%2B%20%5Clambda%20%5Csum_%7B%5Ctau%3D1%7D%5E%7B25%7D%0A%20%20%20%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20%0A%20%20%20%20%5Cleft%5B%20%5Cleft%28%7B%5Cbeta%7D%5Ei_%5Ctextrm%7Bint%7D%5B%5Ctau%2B1%5D%29%20-%20%7B%5Cbeta%7D%5Ei_%5Ctextrm%7Bint%7D%5B%5Ctau%5D%29%20%5Cright%29%5E2%20%0A%20%20%20%20%2B%20%5Cleft%28%7B%5Cbeta%7D%5Ei_%5Ctextrm%7Bext%7D%5B%5Ctau%2B1%5D%29%20-%20%7B%5Cbeta%7D%5Ei_%5Ctextrm%7Bext%7D%5B%5Ctau%5D%29%20%5Cright%29%5E2%0A%20%20%20%20%5Cright%5D%0A%20%20%20%20%5C%5C%0A%20%20%20%20%26%20%5Cquad%20%0A%20%20%20%20%2B%20%5Crho%20%5Csum_%7B%5Ctau%3D1%7D%5E%7B26%7D%20%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20%5Cbeta_%5Ctextrm%7Bext%7D%5B%5Ctau%5D%5E2.%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
    &\varphi\left({\beta}[1], \ldots, {\beta}[26]\right)
    \\
    &:= \sum_{i=1}^N \sum_{\tau=1}^{26} \sum_{t=14\tau - 13}^{14\tau} \left(\Delta^i(t) - \beta^{i}_{\textrm{int}}[\tau] I^i(t) S^i(t) 
    -  \beta^{i}_{\textrm{ext}}[\tau] I^i_{\textrm{ext}}(t) S^i(t) \right)^2
    \\
    & \quad + \lambda \sum_{\tau=1}^{25}
    \sum_{i=1}^{N} 
    \left[ \left({\beta}^i_\textrm{int}[\tau+1]) - {\beta}^i_\textrm{int}[\tau]) \right)^2 
    + \left({\beta}^i_\textrm{ext}[\tau+1]) - {\beta}^i_\textrm{ext}[\tau]) \right)^2
    \right]
    \\
    & \quad 
    + \rho \sum_{\tau=1}^{26} \sum_{i=1}^{N} \beta_\textrm{ext}[\tau]^2.
\end{aligned}
")

Please see the paper for details.

**Figure 5. Plots (a) and (b) capture the quantiles of ’s and ’s across
districts over 2-week time windows from 05/01/2020 to 05/08/2021. The
regression errors in explaining the emergence of daily new infections
with the estimated β’s is given in (c). Plot (d) shows the quantiles of
mobility variations across districts.** text

<img src="Docs/Figure-5.png" width="1885" />
