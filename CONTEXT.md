# FLASHTnT / FLASHAnnotator search

Top-down proteoform search on deconvolved spectra. This glossary covers the search side; the FLASHIda real-time engine's terms
live in the parent workspace's `CONTEXT.md`.

## Language

**MS3 spectrum**:
A spectrum of the fragments of one MS2 fragment ion (the **MS3 precursor**), acquired after its **Parent MS2**. It is evidence for the
parent's PrSM and never an identification of its own.
_Avoid_: MS3 PrSM, MS3 identification.

**Quantification MS3**:
An **MS3 spectrum** acquired to read isobaric reporter ions from one labelled fragment, over the reporter m/z range only. It holds no
**sub-ions**, so it is never evidence for a PrSM, though it still has a **Parent MS2**.
_Avoid_: MS3 (unqualified, when the quantification kind is meant); TMT MS3 (the kind is defined by its purpose, not the label).

**Characterization MS3**:
An **MS3 spectrum** acquired over a fragment m/z range to sequence its **MS3 precursor**. Only this kind supplies **MS3 matches**.

**Parent MS2**:
The MS2 spectrum whose fragment ion an **MS3 spectrum** isolated. Its PrSM is the one the MS3 spectrum's masses support.
_Avoid_: precursor scan (for an MS3 that names the MS1 too), previous MS2 (the parent is usually not the MS2 just before it).

**MS3 precursor**:
The MS2 fragment ion an **MS3 spectrum** isolated: the first precursor msconvert writes for it, with its m/z and charge. Its mass is a
fragment mass, not a proteoform mass.
_Avoid_: precursor (unqualified, for an MS3).

**Sub-ion**:
A fragment of the **MS3 precursor**, as seen in an **MS3 spectrum**. Read in the proteoform's frame, a sub-ion of the precursor's own
terminus is an ordinary terminal fragment of the proteoform, and one of the other terminus is an internal fragment.
_Avoid_: MS3 fragment (ambiguous between the precursor and its sub-ions).

**Precursor identity**:
Which fragment of a candidate proteoform (ion type and index, e.g. y70) the **MS3 precursor** is, found by matching its mass against that
candidate's fragments. It depends on the candidate, so one **MS3 spectrum** can have a different identity, or none, under each candidate.
FLASHIda's label in the scan description is the identity under the live winner and is only compared against.
_Avoid_: MS3 label, target fragment (FLASHIda's word for what it chose to isolate).

**MS3 match**:
A **sub-ion** explained by a candidate proteoform, through that candidate's identity for the **MS3 precursor**. It belongs to the
candidate, not to the spectrum: targets and decoys get MS3 matches by the same rule, and they count in the candidate's score.
