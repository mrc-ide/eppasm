# eppspectrum 0.1

- Update to interface with `epp v0.3` which implements ANC routine testing (ANC-RT) likelihood in EPP 2017.

# eppspectrum 0.1.1

- Add option to use .ep5 file inputs instead of .DP file when preparing fit from PJNZ file. (use_ep5=FALSE by default.)
- update prepare_spec_fit() to not use popadjust if EPP projection has a single population.

# eppspectrum 0.1.2
- Parse Spectrum births in read_specdp_demog_param().
- Add births from Spectrum to specfp object.
- For epp regions, distribute number of births based on 15-44 population size and GFR.