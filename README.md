# Heston Stochastic Volatility
This is a MATLAB implementation of the Heston stochastic volatiliy model. The model is used to price European options two methods:
  1. Monte Carlo
  2. Integrating the characteristic function using:
    a. Gauss Laguerre quadrature method
    b. Trapzeodial numerical integration

Results.m : 
	- Use this to see the outputs. This will execute 2 Live scripts: MonteCarlo_Live.mlx and Integration_Live.mlx

MonteCarlo.m & MonteCarlo_Live.mlx:
	- Use these to obtain prices using Monte Carlo

Integration.m & Integration_Live.mlx:
	- Use these to obtain prices using integration of the characteristic function

HestonModel_Combined.m:
	- Combines the code for Monte Carlo and Integration

FFTIntegration.m:
	- Implements the numerical methods used for integration

SchoutensCF:
	- Implementation of the characteristic function proposed by Schoutens et al. (2004)

CuiCF:
	- Implementation of the characteristic function proposed by Cui et al. (2017)

GenerateGaussLaguerre.m:
	- Generates 32-points for Gauss Laguerre quadrature integration
  
  
## References:
- Yiran Cui, Sebastian del Ban ̃o Rollin, Guido Germano, “Full and fast calibration of the Heston stochastic volatility model”,    European Journal of Operational Research 263 (2), 625–638, 2017, DOI 10.1016/j.ejor.2017.05.018.

- Albrecher, H., Mayer, P., Schoutens, W., & Tistaert, J. (2007). Thelittle Heston trap. Wilmott Magazine, 2007(January), 83–92. http://www.wilmott.com/pdfs/121218_heston.pdf.

- Schoutens, W., Simons, E., & Tistaert, J. (2004). A perfect calibration! Now what? Wilmott Magazine, 2004(March), 66–78. http://www.wilmott.com/pdfs/070319_schoutens.pdf.

- Rouah, F. D. (2013). The Heston model and its extensions in Matlab and C#. Hoboken, New Jersey: Wiley. doi:10.1002/9781118656471.
