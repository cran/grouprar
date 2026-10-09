# grouprar 0.2.0

## New features

* Multi-arm efficient randomized adaptive design (ERADE): the DBCD-family
  designs gain `allocation = "ERADE"` and `erade.alpha` (Hu, Zhang and He,
  2009, for two arms; Alkhnefr, Hu and Zhai, 2025, for k arms).
* k-arm optimal target allocations with a lower bound `lower.bound` on each
  proportion: `"OptimalNeyman"` and `"OptimalRSIHR"` for binary responses
  (Tymofyeyev, Rosenberger and Hu, 2007) and `"OptimalNeyman"` for continuous
  responses. `"ZR"` now works for any number of arms and, for two arms,
  follows rule (7) of Zhang and Rosenberger (2006).
* Group sequential monitoring for two-arm designs (`DBCD_Bin()`,
  `DBCD_Cont()`, `Group.DBCD_Bin()`, `Group.DBCD_Cont()`) through the new
  `monitor` argument, `sqMonitor()` and `sqBoundary()`, with O'Brien-Fleming,
  Pocock and linear alpha spending (Zhu and Hu, 2010). Results include the
  stopping probability at each look and the expected sample size.
* `nextAlloc()` computes the allocation probabilities (and optionally the
  assignments) of the next patient or group in an ongoing trial.
* All designs gain `test.fun` (user-supplied test returning a p-value),
  `typeI` (also estimates the type I error under the null) and `seed`.
* Results are objects of class `"grouprar"` with `print()` and `summary()`
  methods. They now include the allocation sequence of every simulated trial
  (`data: allocation`) and the SD of the allocation proportion of every arm
  (`sd of propotion`, previously arm 1 only). Delayed designs also report the
  trial duration and the enrollment duration.

## Bug fixes

* The Hu and Zhang (2004) allocation function behind all DBCD designs now
  normalizes each arm jointly against all k arms. Previously each arm was
  compared with the pooled remaining arms, which agrees with Hu and Zhang
  (2004) only for two arms. Results for two-arm designs are unchanged. Thanks
  to Stina Zetterstrom, David S. Robertson and Sofía S. Villar for the report
  (StinaZet/RAR-software-bugs#7).
* `dyldDBCD_Bin()` updated its estimates with the responses of the wrong
  patients and ignored late responses of the initial patients. It also failed
  for more than two arms.
* `Bai.Hu.Shen.Urn()` updated the urn with the true success rates. It now
  implements the proposed design of Bai, Hu and Shen (2002), which uses the
  estimated success rates.
* `dyldDBCD_Cont()` ignored the `r` argument, failed for more than two arms,
  did not mark missing responses as unobserved, and did not update the
  allocation when a response of one of the initial patients arrived.
* `dyldDBCD_Bin()` and `dyldDBCD_Cont()` now recompute the allocation
  probability for every patient with the current allocation proportions
  (previously only when a new response arrived).
* `Group.dyldDBCD_Bin()` and `Group.dyldDBCD_Cont()` drew an extra group when
  the group sizes added up exactly to `ssn`.
* `Group.dyldDBCD_Bin()` used the wrong column of the response time matrix, so
  response times did not depend on the observed response as documented.
* `DBCD_Bin()` ignored `mRate`, and `Group.DBCD_Cont()` failed when `mRate` was
  set.
* With missing responses, the allocation proportions used by the allocation
  function now count all enrolled patients, as in Section 2.3 of Zhai et al.
  (2024). `Group.DBCD_Bin()` and `Group.DBCD_Cont()` previously counted only
  patients with observed responses. The reported allocation proportions
  (`propotion`) of all designs with `mRate` also count all enrolled patients
  now; the failure rate and the test still use the observed responses only.
* `DBCD_Bin()` and `DBCD_Cont()` now use the chi-squared test for more than two
  arms, and the chi-squared test for continuous responses uses the sample
  variance of each arm.
* The tests no longer fail for two arms (`PolyaUrn()`, `WeiUrn()`), for arms
  without patients, or for arms without variability. For binary responses with
  more than two arms, arms without variability no longer make the test NA:
  the Wald test then uses the adjusted variances of Agresti and Caffo (2000).
* `GDLRule()` returns one failure rate per simulation, and `GDLRule()`,
  `DLRule()` and `BirthDeathUrn()` no longer fail when `aK` or `Y0` are not whole
  numbers (balls are drawn in proportion to the positive part of the urn).
* `rspT.dist = "uniform"` was rejected because of a typo.
* Allocation proportions are no longer misaligned when an arm has no patients.
* Equal allocation is used while the target cannot be estimated (e.g. fewer
  than two responses in an arm of a continuous design).
* `RPWRule()` stops with an error for `k != 2`, the DBCD-family designs
  check that `n0 >= k` and `ssn > n0`, and the group designs check that
  `gsize.param > 0` (a rate of 0 made the initial phase loop forever). All
  designs check `k >= 2`, `Y0`, `aK` and the length of `rspT.param`.

## Documentation

* Corrected return values, argument descriptions, references and examples in
  all help pages. Added a citation file.

## Other changes

* Removed the unused `ggplot2`, `gridExtra`, `methods` and `tidyr` imports, and
  the dependency on `stringr`, whose current version needs R >= 4.1.
* Added unit tests and a GitHub Actions workflow for R CMD check.
