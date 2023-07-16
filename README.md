# Derivative Securities Portfolio Optimization

- Final Project for Financial Economics II at University College Dublin Lochlann Quinn School of Business
  
- Equivalent to Options/Futures/Derivatives Securities at US Schools.
  
- Used Fixed income and interest rate term structure concepts, designing a portfolio of fixed-income instruments immune to Delta and Gamma risks
  
- The portfolio was tested in a simulation provided by the Professor with hypothetical changes to the term structure across time (time series of interest rate term structure).
  
- Portfolio Assets: APAC Bonds
  - Government bonds
  - Corporate bonds
    
- Portfolio Liabilities: Pension Fund with three streams of cashflow:
  - Negative Cashflow: For every monthly increase of 10 bp of the 10-year interest rate, the investors extract $100K
  - Positive Cashflow: For every monthly decrease of 10 bp of the 10-year interest rate, they increment investments by $100K
  - Fixed Extract: Monthly extract of 1% of the total managed APAC bonds
    
- Extracted the yield curve with the same data using the Nelson-Siegel model

- Extracted a second yield curve using bootstrapping

- Elaborated an interactive algorithm that repeated this process iteratively for the simulation

- Calculated the optimal portfolio considering assets, liabilities, and potential future scenarios.

