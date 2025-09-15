# Multifidelity Monte Carlo

This is a MATLAB implementation of the Multifidelity Monte Carlo (MFMC) approach to mean estimation and linear regression described in:

1. Peherstorfer, B., Willcox, K., & Gunzburger, M. (2016).
   [Optimal Model Management for Multifidelity Monte Carlo Estimation.](https://doi.org/10.1137/15M1046472) SIAM Journal on Scientific Computing, 38(5), A3163–A3194.
<details markdown="1">
  <summary>BibTeX</summary>

  ```bibtex
   @article{peherstorfer_optimal_2016,
		title = {Optimal {Model} {Management} for {Multifidelity} {Monte} {Carlo} {Estimation}},
		volume = {38},
		issn = {1064-8275, 1095-7197},
		url = {http://epubs.siam.org/doi/10.1137/15M1046472},
		doi = {10.1137/15M1046472},
		language = {en},
		number = {5},
		urldate = {2025-06-24},
		journal = {SIAM Journal on Scientific Computing},
		author = {Peherstorfer, Benjamin and Willcox, Karen and Gunzburger, Max},
		month = jan,
		year = {2016},
		pages = {A3163--A3194},
	}
```
</details>

2. Qian, E., Kang, D., Sella, V., & Chaudhuri, A. (2025). [Multifidelity linear regression for scientific machine learning from scarce data](https://doi.org/10.3934/fods.2024049). Foundations of Data Science, 7(1), 271–297.

<details markdown="1">
  <summary>BibTeX</summary>

  ```bibtex
   @article{qian_multifidelity_2025,
	title = {Multifidelity linear regression for scientific machine learning from scarce data},
	volume = {7},
	copyright = {http://creativecommons.org/licenses/by/3.0/},
	url = {https://www.aimsciences.org/en/article/doi/10.3934/fods.2024049},
	doi = {10.3934/fods.2024049},
	language = {en},
	number = {1},
	urldate = {2025-06-24},
	journal = {Foundations of Data Science},
	author = {Qian, Elizabeth and Kang, Dayoung and Sella, Vignesh and Chaudhuri, Anirban},
	month = mar,
	year = {2025},
	note = {Publisher: Foundations of Data Science},
	pages = {271--297},
	}
```
</details>

## Summary
The scripts `mfmcExp.m` and `mfmcCDR.m` respectively implements the Multifidelity Monte Carlo method for mean estimation [[1](https://doi.org/10.3934/fods.2024049)] for the numerical experiments performed in [[2](https://doi.org/10.3934/fods.2024049)] while the scripts `mfLinearExp.m` and `mfLinearCDR.m` implements Multifidelity Monte Carlo Method for Linear Regression described in [[2](https://doi.org/10.3934/fods.2024049)] with the same examples. Running these scripts will generate plots found in the `plots` folder.

## Background
The Monte Carlo Methods are a large class of computational algorithms that uses Random Sampling to obtain results. For a single-fidelity (one model) implementation of the Monte Carlo mean estimator, we find that while results are robust, convergence to the solution is costly and requires many samples. 

Monte Carlo Mean Estimator:

<p align="center">
	$$\mathbb{E}[f(X)] \sim \hat{y}_N = \frac{1}{N}\sum_{i=1}^{N}f(x_i) $$
</p>

This project introduces the implementation of the Multifidelity Monte Carlo Method that uses multiple models to evaluate mean values using a control variate method. This method balances accuracy with computational cost, acheiving greater accuracy with the same cost as the single fidelity method.

Multifidelity Monte Carlo Mean Estimator:

<p align="center">
	$$\mathbb{E}[f(X)] \sim \hat{y}_{\text{MFMC}} = \hat{y}^{(1)}_{m_1} + \alpha ( \hat{y}^{(2)}_{m_2} - \hat{y}^{(2)}_{m_1})$$ 
</p>


where $m_1$ and $m_2$ denote sample sizes corresponding to each model and $\alpha$ represents a coefficient optimally chosen to decrease the variance.

## Examples

### Analytical Exponential Example
`mfmcExp.m` compares the performance of single-fidelity Monte Carlo mean estimation with the multifidelity counterpart with evaluating the average value of $\exp(x)$ in the interval $x \in [0, 5]$. 

`mfLinearExp.m` compares the performance of a single-fidelity Monte Carlo approach to linear regression with its multifidelity counterpart with finding a polynomial regression of $\exp(x)$ in the interval $x \in [0, 5]$ using a univariate monomial function basis of degree 4: $X = [1, z, z^2, z^3, z^4]^\top$

### Convection-Diffusion-Reaction Example
The CDR example is from data resulting from a computational combustion simulation. With 5 chemical and temperature inputs, we arrive at our quantity of interest, the maximum temperature of the chamber at steady state.

`mfmcCDR.m` compares the performance of single-fidelity Monte Carlo mean estimation with the multifidelity counterpart with evaluating the average QoI from the CDR Data. 

`mfLinearCDR.m` compares the performance of a single-fidelity Monte Carlo approach to linear regression with its multifidelity counterpart. We use a quadratic polynomial basis with our 5 input variables to estimate our QoI. 

