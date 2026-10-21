# Step 14: Author Review Memo – Unresolved Claims & Strategic Pre-Submission Decisions

**To**: Co-Authors & Principal Investigators  
**From**: Pre-Submission Research Reproducibility Auditor  
**Date**: 2026-10-21  
**Deliverable File**: [`reports/14_unresolved_claims_author_review.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis-02/reports/14_unresolved_claims_author_review.md)  

---

## Purpose

This memorandum summarizes the four key strategic decisions and scientific qualifications that must be aligned upon by co-authors prior to formal journal submission:

### Decision 1: Reporting of Multiplicity-Adjusted p-values ($p_{\text{adj}} = 0.0602$)
- **Background**: Comparing ESM-2 against Random Forest ($p=0.0201$) and Logistic Regression ($p=0.0294$) yields nominally significant raw p-values. However, applying Holm-Bonferroni correction across the 4 planned comparisons adjusts these to $p_{\text{adj}} = 0.0602$, which misses the standard $\alpha = 0.05$ cutoff.
- **Author Decision**: We strongly recommend reporting the Holm-adjusted p-values ($p_{\text{adj}} = 0.0602$) alongside the strictly positive 95% bootstrap confidence intervals for $\Delta \text{Macro-F1}$ ([+0.0915, +0.4626] and [+0.0540, +0.3976]). This demonstrates methodological rigor and prevents reviewer rejection for uncorrected multiplicity.

### Decision 2: Framing the HSV-2 $N=7$ Ortholog Transfer Experiment
- **Background**: The manuscript draft previously described the 7/7 correct classification of HSV-2 genes as "independent external validation." However, all 7 families are core orthologs shared with HSV-1 genes present during training.
- **Author Decision**: The text has been revised to describe this as an "ortholog sequence transfer diagnostic between sister species." Authors must confirm that no claims of "generalization to novel viral families" remain in the discussion.

### Decision 3: Disclosure of Low-Dimensional PCA Generalization Failure
- **Background**: Step 12 initially reported 89.1% accuracy for an 18-component PCA model. Cross-validation auditing in Step 12.1 and Step 13 proved this was in-sample overfitting, with held-out Macro-F1 collapsing to 43.88%.
- **Author Decision**: Retain this negative result in the Results section. It provides a strong biological insight: viral temporal timing discriminative signals are distributed across non-linear embedding space and cannot be compressed into low-rank linear components.

### Decision 4: Acknowledging the External Validation Gap
- **Background**: Zero external clinical datasets exist in the repository.
- **Author Decision**: Authors must endorse the limitation statement indicating that external validation is prospectively planned (e.g. mass spectrometry profiling of clinical isolates under metabolic arrest), rather than claiming completed external clinical validation.
