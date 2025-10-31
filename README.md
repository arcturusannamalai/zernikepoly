# zernikepoly
Zernike polynomials code from article JTu4A.25 at FiO LS 2025

If you find the work useful provide a citation as follows:
> M. Annamalai, "Techniques to Improve Computation of Zernike Polynomials," Frontiers in Optics - Laser Science, Denver, CO (2025).

Times without Zernike field cache:
```text
Elapsed for 10 = 0.09803891181945801
Computed 66 coefficients.
@Order 10 Relative reconstruction error: 1.841e-15
Elapsed for 12 = 0.14208197593688965
Computed 102 coefficients.
@Order 12 Relative reconstruction error: 1.038e-15
Elapsed for 14 = 0.18423676490783691
Computed 144 coefficients.
@Order 14 Relative reconstruction error: 1.450e-15
Elapsed for 16 = 0.2552919387817383
Computed 192 coefficients.
@Order 16 Relative reconstruction error: 1.086e-15
```

Times with field cache:
```text
Elapsed for 10 = 0.09715414047241211
Computed 66 coefficients.
@Order 10 Relative reconstruction error: 1.841e-15
Elapsed for 12 = 0.18396401405334473
Computed 102 coefficients.
@Order 12 Relative reconstruction error: 1.038e-15
Elapsed for 14 = 0.17792916297912598
Computed 144 coefficients.
@Order 14 Relative reconstruction error: 1.450e-15
Elapsed for 16 = 0.18890714645385742
Computed 192 coefficients.
@Order 16 Relative reconstruction error: 1.086e-15
```

Limits:
```commandline
We can compute n=48 order polynomial in span of 1.066983938217163s on a standard laptop/macbook.
```

The same code is ported to MATLAB (from Python) at zernike.m
