# pGMCM: Analysis with penalized Gaussian mixture copula models

## Dependencies
1.  To use the full extent of this package, you need to download the C++ library LEMON graph library; [download it here.](https://lemon.cs.elte.hu/trac/lemon/wiki/Downloads)

If you do not have them already, you do not need all of LEMON's dependencies to run the pGMCM package. To ease this during configuring LEMON, you can go into the INSTALL file after downloading LEMON and add options

```
-DLEMON_ENABLE_GLPK=NO
-DLEMON_ENABLE_COIN=NO
-DLEMON_ENABLE_ILOG=NO
```

<!---
LEMON citation:
Balázs Dezső, Alpár Jüttner, Péter Kovács. LEMON – an Open Source C++ Graph Template Library. Electronic Notes in Theoretical Computer Science, 264:23-45, 2011. Proceedings of the Second Workshop on Generative Technologies (WGT) 2010.
-->

2.  You need a compiler that has support for C++11, such as
    *   GCC: [see here, for example](https://www.gnu.org/software/gcc/projects/cxx-status.html#cxx11)
    *   clang: [see here](http://clang.llvm.org/cxx_status.html)

## Getting the package
This package is currently only maintained on GitHub. To download just the R package portion of the software, you can open R and do the following commands:

```{r}
library(devtools)
install_github("hillarykoch/pGMCM/R")
library(pGMCM)
```
*note: the R version is fully functional but the Julia code is just a test that cannot practically be used.*

## About the software
<!---
 The GMCM is a copula mixture that generalizes to any dimension. This package implements a general form of the pGMCM as well as a constrained version. It also implements a general and similarly constrained penalized Gaussian mixture model. The penalization allows for selection of the number of clusters, subject to a user-selected upper bound.
-->

General case, with user-selected upper bound 10:
![triangle](triangle.png)

Constrained case with 9 components, with user-selected upper bound 9:
![nine](nine.jpg)

Constrained with 4 components, with user-selected upper bound 9:
![subnine](subnine.png)
