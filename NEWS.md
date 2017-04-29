# eppspectrum 0.1

- Update to interface with `epp v0.3` which implements ANC routine testing (ANC-RT) likelihood in EPP 2017.


# eppspectrum 0.1.3

- EPP model setup occurs in fitmod()
- Allow order of spline model penalty to vary (fp$rtpenord)
- Allow number of knots to vary (fp$numKnots)
- Return spline coefficients in fnCreateParam()

# eppspectrum 0.1.5

- Implemented O-spline model (penalized B-spline) using mgcv to construct