# pGMCM: Analysis with penalized Gaussian mixture copula models

The GMCM is a copula mixture that generalizes to any dimension. This package implements a general form of the pGMCM as well as a constrained version. It also implements a general and similarly constrained penalized Gaussian mixture model. The penalization allows for selection of the number of clusters, subject to a user-selected upper bound.

General case, with user-selected upper bound 10:
![triangle](triangle.png)

Constrained case with 9 components, with user-selected upper bound 9:
![nine](nine.jpg)

Constrained with 4 components, with user-selected upper bound 9:
![subnine](subnine.png)

*note: the R version is fully functional but the Julia code is just a test that cannot practically be used.*


