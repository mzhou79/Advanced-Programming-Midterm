# Block Bootstrap to Estimate the Downfall of an Investment Portfolio

[![Slides and report render](../../actions/workflows/check-qmd-render.yaml/badge.svg)](../../actions/workflows/check-qmd-render.yaml)

**PHS 7045 Advanced Programming — Midterm Project**

## Overview

Stock investors commonly hold portfolios containing multiple exchange-traded funds (ETFs). Understanding potential portfolio losses is important for evaluating investment risk.

This project addresses the question:

**Given a portfolio's composition and historical market behavior, how severe could its future losses become?**

I developed a moving block bootstrap approach to simulate possible future portfolio values using historical daily returns. Unlike ordinary bootstrap sampling, which independently samples observations, moving block bootstrap resamples consecutive trading days. This helps preserve short-term temporal dependence and relationships between assets.

The demonstration uses three ETFs: VOO (S&P 500), QQQ (Nasdaq-100), and SCHD (dividend-focused U.S. equities).

## Project Structure

| File | Description |
|---|---|
| `R/riskpath.R` | Five documented R functions for simulation, risk estimation, and visualization |
| `report.qmd` | Quarto report explaining the methodology and results |
| `report.pdf` | Rendered PDF report |
| `slides.qmd` | Quarto presentation source |
| `slides.html` | Rendered presentation |
| `prices.csv` | Historical ETF price data |
| `README.md` | Project documentation |
| `.github/workflows/check-qmd-render.yaml` | GitHub Actions rendering workflow |



## AI Disclosure

I developed the idea of applying a moving block bootstrap to simulate future portfolio outcomes and described my intended workflow to AI, which helped generate the R implementation.

I reviewed the generated code step by step and checked whether it followed my proposed algorithm. In particular, I studied previously unfamiliar functions such as `stopifnot()`, `unlist()`, and `cumprod()` to understand their roles.

I also used AI to help generate the three visualizations based on my intended graphical designs and to organize the project documentation. I reviewed the generated code and explanations for consistency with the intended analysis.

## Course Submission

This project is associated with [PHS 7045 Advanced Programming — Issue #65](https://github.com/UofUEpiBio/PHS7045-advanced-programming/issues/65).