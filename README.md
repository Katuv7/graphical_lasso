# graphical_lasso

Implementation: Sparse inverse covariance estimation with the graphical lasso from Friedman's paper.

The graphical lasso algorithm is a powerful tool for estimating sparse graphical models, specifically focusing on the precision matrix (\(\Theta\)) of a multivariate Gaussian distribution with mean $\mu$ and covariance matrix $\Sigma$. By incorporating \(L_1\) regularization, it promotes sparsity in \(\Theta\), which is instrumental in revealing the conditional independence structure among variables.

# Sources

1. Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. Statistics and Computing, 18(3), 303-313. https://doi.org/10.1007/s11222-008-9206-1
