# RationalExtensions

Provides arithmetic support for quadratic field extensions of Q. For example, 

```
using RationalExtensions
(Sqrt(3) + Sqrt(9))//(1+Sqrt(3))
norm(1 - 2Sqrt(3))
conj(1//3 - Sqrt(3))
```

To install, run

```
Pkg.add("RationalExtensions")
```
