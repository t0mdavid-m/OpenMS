---
status: accepted
---

# An MS3 spectrum is evidence for its parent's PrSM, scored by a chance-corrected term outside the learned model

An MS3 spectrum is never searched or identified on its own: its masses count only as **MS3 matches** of the candidates of its **Parent
MS2**, through each candidate's own **Precursor identity**, for targets and decoys alike. They enter the score as one separate term,
−log10 of the chance of that many matches, calibrated on the same MS3 spectrum with mass-shifted sub-ions under the same matching
rules. The term is added to the learned model's choice score and confidence log-odds with a weight (default 1), and every PrSM stays in
the one FDR pool. We decided this on 2026-11-09, when we had one MS3 run (Caco-2 histone, 1,114 MS3, 56 % deconvolved to nothing).

## Considered Options

- **Standalone MS3 PrSMs.** An MS3 precursor is a 2–10 kDa fragment, and at a 2 % bond yield the spectrum is weak evidence in a human
  database on its own. It would also put a second, correlated entry into the FDR.
- **Folding the projected sub-ions into the MS2 fragment list.** The learned model would then see mass counts, intensity shares and
  coverage that mean something different from what it was trained on, and its calibration would drift without any sign of it. This would
  also drop the internal sub-ions.
- **Retraining the model with MS3 features.** One run is far too little data for a model of this size.
- **A separate FDR class for MS3-linked PrSMs.** It would hold a few hundred PrSMs, too few decoys for stable q-values.

## Consequences

- FLASHIda's `scan_commands.tsv`, read through FLASHDeconv's `-scan_commands` and a new FLASHAnnotator option, gives the exact parent
  link and the MS3 kind: `R` is a **Characterization MS3**, `Q` a **Quantification MS3**, which is linked but never scored. With a log
  given, an MS3 the log cannot link is dropped and counted. Without a log, parents come from the m/z rule and every MS3 counts as
  characterization. On the Caco-2 run the two links agree on 1,111 of 1,114 MS3s. The other 3 had two m/z candidates, and the m/z rule
  picked the wrong one.

- An MS3 can only reorder candidates the MS2 has already retrieved. It never adds one.
- A candidate that cannot explain the MS3 precursor gets no evidence from that MS3, and no penalty either. The precursor was chosen
  under FLASHIda's live winner.
- The MS3 term needs the raw MS3 spectrum (`-in_raw`): it is matched with FLASHAnnotator's envelope placement, never from deconvolved
  MS3 masses.
- The term is used as a likelihood ratio, but whether it behaves like one has not been measured. It is accepted on the coverage and
  localisation it adds to the 1 % PrSMs.
