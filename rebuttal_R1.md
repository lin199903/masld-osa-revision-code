# Reviewer #1 Rebuttal Draft (Skeleton)

## Response 0 (Clarification on Model Version and Terminology)

We thank the reviewer for the careful reading. We first clarify a version-alignment issue before addressing individual points.

1. The current submitted manuscript is based on a **final locked 10-gene immune-hypoxia signature**, not an earlier draft Venn-derived model.  
Evidence in current manuscript:
- Main text abstract/results: `submission_pack/main_revised.tex:102`
- Workflow and methods wording: `main_revised.tex:149`, `main_revised.tex:175-187`
- Main performance figure caption: `main_revised.tex:284-285` (`\label{fig:ml_comprehensive}`)
- Supplementary validation table: `submission_pack/supp_revised.tex:404-405` (`\label{tab:orthogonal_validation}`)

2. We apologise for any confusion caused by version or figure-mapping differences between earlier draft and current submission versions; we clarify the precise terminology below. The current submission uses a final locked **10-gene** model throughout.

3. To avoid confusion in the revision, we will:
- (i) add an explicit sentence in Methods/Results that the reported score is from the **final locked 10-gene model output**;
- (ii) move feature-importance summary (currently emphasized in Supplementary Fig. S9 SHAP panel) to a main-text panel or a prominent in-text pointer;
- (iii) add disease negative-control stress testing (ALD/AH cohorts) as specificity evidence.

4. Consistency gate for this rebuttal:
- We keep the manuscript score definition aligned to the locked 10-gene model output.
- Any Phase A analysis that is definition-incompatible (e.g., exploratory legacy-score direction audit) is treated as **candidate supplementary analysis** and not used to overwrite current main-text claims.
- The reviewer-boundary analyses in Supplementary S10/S15 use an auxiliary marker-set score for stress-testing only and are explicitly separated from the manuscript primary model claim.

5. Clarification on why we do not reframe the manuscript around a 5-marker set:
- We agree that forcing a 5-marker narrative would weaken biological and methodological coherence of the current submission.
- The manuscript's central model is the locked 10-gene SVM-RFE+LDA signature; this is the only score used for primary conclusions.
- The auxiliary marker-set analyses are retained solely to bound over-interpretation in reviewer-requested negative-control and descriptive clinical-utility checks.

---

## Overall Positioning

In this revision, we have substantially strengthened the manuscript with respect to machine-learning rigor, MR sensitivity analyses, feature attribution, and calibration. At the same time, we have further addressed the issues of disease specificity, clinical utility benchmarking, and physiological hypoxia linkage; however, owing to the constraints of publicly available data, these three aspects should still be considered **partially addressed rather than fully resolved**.

---

## Major Comment 1: Mechanistic overlap vs disease specificity

### (a) We thank the reviewer...
We agree that specificity versus broad inflammatory liver biology must be stress-tested with explicit disease negative controls.

### (b) Direct answer
**Status: Partially addressed with strengthened negative-control evidence.**  
We agree that disease specificity is a critical issue. In the revised manuscript, we strengthened the negative-control analysis without overclaiming disease-specificity. Specifically, we include one analyzable case-control negative-control cohort (`GSE28619`) for discrimination-oriented evaluation, and one disease-only cohort (`GSE103580`) as descriptive evidence only. Because `GSE103580` does not contain suitable control samples, it is not used for case-control discrimination metrics and is presented as descriptive context for non-target inflammatory liver disease background. Accordingly, the revised data package provides stronger boundary-setting evidence against overinterpretation, but we do not claim that disease specificity has been fully established.

### (c) Existing evidence pointer
- Supplementary revision asset for specificity status and distributions:  
  `submission_pack/supp_revised.tex:501-507` (`\label{fig:rev1_negcontrol}`; approx. PDF p14)
- Specificity discrimination table (target cohorts + analyzable ALD/AH negative control):  
  `submission_pack/rev1_tables/Table_NegControl_IHS_Performance_Comparison.tex` (Table S10 input).
- Descriptive-only distribution summary (single-class cohort, no AUC claim):  
  `submission_pack/rev1_tables/Table_NegControl_DescriptiveOnly_Distribution.tex`.
- Machine-readable cohort status table (`analyzable` vs `descriptive_only_no_control`):  
  `submission_pack/rev1_tables/Table_NegControl_Analysis_Status.tex`.
- Pipeline status table documenting `GSE103580` as descriptive-only no-control cohort:  
  `rev1_major_revision/outputs/tables/negative_control_analysis_status.csv`.

### (d) What we will add in revision
- Keep `GSE28619` as the primary case-control ALD/AH negative-control discriminator.
- Retain `GSE103580` as disease-only descriptive evidence (no AUC claim due to absent controls) and state this explicitly in S10 subtitle/caption.
- Keep conservative interpretation language: current evidence does not establish disease specificity versus non-MASLD inflammatory liver disease.

### (e) Manuscript changes
- Supplementary figure/table block around `supp_revised.tex:501-519` (`\label{fig:rev1_negcontrol}` + `Table_NegControl_IHS_Performance_Comparison` input).
- Methods and discussion limitation text around `main_revised.tex:175-187` and `main_revised.tex:361-365`.

---

## Major Comment 2: Machine-learning pipeline rigor / leakage prevention

### (a) We thank the reviewer...
We agree leakage control and split logic must be explicit and auditable.

### (b) Direct answer
This is **addressed** in current text. The manuscript now contains an explicit anti-leakage Methods paragraph (nested learning subsection, revised manuscript Methods §"Strictly nested learning") that states: all feature selection and hyperparameter tuning were confined to inner folds, the external cohort was used only for post-training robustness prioritisation, and the hyperparameter grid is detailed in Supplementary Table S3. An explicit common-gene attrition sentence (26/29 genes retained) and the IHS score-definition lock sentence have been added to the same subsection.

### (c) Existing evidence pointer
- Strictly nested and anti-leakage description: `main_revised.tex:175-183` (approx. PDF p5)
- Hyperparameter grid reference: `main_revised.tex:177`; Supplementary table: `supp_revised.tex:357-361` (`\label{tab:ml_grid}`, approx. PDF p13)
- External cohort isolation statement: `main_revised.tex:183`
- Calibration/PRC presentation: `main_revised.tex:267-285` (`\label{fig:ml_comprehensive}`, approx. PDF p9)

### (d) What we will add in revision
- Add explicit sentence-level split checklist (feature selection/tuning strictly inner-fold; external used only post-training robustness ranking).
- Add one concise line quantifying cross-platform setting and preprocessing harmonization summary.

### (e) Manuscript changes
- Methods subsection around `main_revised.tex:175-187`.
- Optional short note in Results around `main_revised.tex:244-251`.

---

## Major Comment 3: Mendelian randomization validity (Steiger/MR-PRESSO/robust methods)

### (a) We thank the reviewer...
We agree MR assumptions and pleiotropy sensitivity must be presented explicitly.

### (b) Direct answer
This is **already addressed**. Specifically: Steiger directionality (prop.\ correct = 100\%, $P=2.38\times10^{-7}$) is reported in the main text MR paragraph; MR-PRESSO outlier detection (1 outlier SNP removed) is in Supplementary Table S6; MR-Egger intercept ($P=0.349$, no directional pleiotropy) and weighted-median are in Supplementary Table S4/S6. All six sensitivity methods are consolidated in the Supplementary MR sensitivity summary figure (\texttt{fig:rev1\_mr\_summary}). No further reanalysis is needed; we have verified these cross-references are intact.

### (c) Existing evidence pointer
- Main text MR sensitivity narrative: `main_revised.tex:291-300` (Steiger/MR-PRESSO stated at line 298; approx. PDF p10-12)
- Main MR figure: `main_revised.tex:305-325` (`\label{fig:mr_bidirectional}`)
- Supplementary MR sensitivity table: `supp_revised.tex:418-432` (`\label{tab:mr_sensitivity}`, approx. PDF p14)
- Supplementary MR sensitivity summary figure (includes weighted median): `supp_revised.tex:510-514` (`\label{fig:rev1_mr_summary}`)

### (d) What we will add in revision
- Add one sentence in Results/Limitations explicitly tying exclusion-restriction risk to high heterogeneity and pleiotropy sensitivity interpretation.

### (e) Manuscript changes
- Results MR paragraph around `main_revised.tex:291-301`.
- Limitations around `main_revised.tex:352-361`.

---

## Major Comment 4: IHS cutoff and DCA clinical utility

### (a) We thank the reviewer...
We agree that clinical interpretability requires explicit threshold and net-benefit context.

### (b) Direct answer
**Status: Partially addressed.**  
We appreciate the reviewer’s point regarding clinical utility and threshold-oriented interpretation. In the revised manuscript, we added calibration assessment and decision-curve analysis to better characterize score behavior in the external cohort. However, because extractable `FIB-4` and `NFS` variables were not available in the public transcriptomic cohorts used here, we were unable to perform a formal incremental net-benefit comparison against these established clinical fibrosis scores. We therefore present DCA strictly as a descriptive analysis and do not use it to define a clinical cutoff or to claim added net benefit over `FIB-4/NFS`.

### (c) Existing evidence pointer
- Main text DCA interpretation caveat: `main_revised.tex:250` (descriptive-only wording)
- Supplementary DCA figure: `supp_revised.tex:270-289` (`\label{fig:s9_dca_shap}`, approx. PDF p11)
- Supplementary locked-threshold metrics table: `supp_revised.tex:485-488` (`\label{tab:S9_threshold_metrics_exactci}`, approx. PDF p14)
- Supplementary threshold-level incremental net-benefit table (IHS vs SLC2A1-only):  
  `submission_pack/rev1_tables/Table_DCA_Incremental_vs_SLC2A1_KeyThresholds.tex`.

### (d) What we will add in revision
- Keep explicit score-definition lock sentence: manuscript "IHS/score" refers to locked 10-gene model output probability/risk score.
- Keep DCA as descriptive-only and explicitly state that no formal incremental comparison versus `FIB-4/NFS` can be performed with currently available public data.

### (e) Manuscript changes
- Methods model-definition area `main_revised.tex:175-187`.
- Supplementary threshold table note area `supp_revised.tex:485-499`.

---

## Figure/Table Comment 1: Figure 2 Venn and feature importance

### (a) We thank the reviewer...
We agree model contribution transparency is essential.

### (b) Direct answer
This is **already addressed**. SHAP feature attribution is already provided in Supplementary Fig. S9, with explicit in-text/caption linkage from the main manuscript.

### (c) Existing evidence pointer
- SHAP panel in Supplementary Fig. S9: `supp_revised.tex:283-289` (`\label{fig:s9_dca_shap}`)
- Main text references SHAP in Fig. 4 caption: `main_revised.tex:284-285`

### (d) What we will add in revision
- Keep explicit cross-reference language so readers can directly locate the SHAP-based attribution evidence.

### (e) Manuscript changes
- Main figure block `main_revised.tex:260-286` (`\label{fig:ml_comprehensive}`)
- Results paragraph around `main_revised.tex:219-223`.

---

## Figure/Table Comment 2: Calibration with ROC

### (a) We thank the reviewer...
We agree ROC without calibration is insufficient for clinical interpretation.

### (b) Direct answer
This is **already addressed** in current submission.

### (c) Existing evidence pointer
- ROC + calibration panel (Fig. 4b): `main_revised.tex:267-270`, caption at `main_revised.tex:284-285` (`\label{fig:ml_comprehensive}`, approx. PDF p9)

### (d) What we will add in revision
- Add one sentence in Results explicitly naming calibration panel as response to reviewer concern.

### (e) Manuscript changes
- Results around `main_revised.tex:244-251`.

---

## Figure/Table Comment 3: Table 1 nocturnal desaturation parameters

### (a) We thank the reviewer...
We agree hypoxia-physiology variables should be tabulated explicitly.

### (b) Direct answer
**Status: Partially addressed.**  
We agree that direct linkage between the score and physiological nocturnal hypoxia metrics would materially strengthen the biological interpretation. In the revised supplement, we added a structured metadata audit summarizing which cohorts contain extractable hypoxia-related variables, including `AHI`, `RDI`, and oxygen saturation descriptors where available. However, most public transcriptomic datasets lacked sample-level `ODI`, `T90`, or `SpO2` measurements suitable for robust correlation testing. Therefore, although data availability has now been explicitly audited and documented, we could not perform a complete patient-level correlation analysis between the score and physiological nocturnal desaturation parameters in the current revision.

### (c) Existing evidence pointer
- Supplementary hypoxia table is included via: `supp_revised.tex:519` (`\input{rev1_tables/Table1_Supplement_Hypoxia_Parameters.tex}`)
- Table label in input file: `submission_pack/rev1_tables/Table1_Supplement_Hypoxia_Parameters.tex` (`\label{tab:rev1_hypoxia}`)
- Supplementary metadata completeness audit table (variable-level coverage):  
  `submission_pack/rev1_tables/Table1_Hypoxia_Metadata_Completeness.tex`.

### (d) What we will add in revision
- Keep unavailable fields as NA with explicit limitation note (`ODI/T90/SpO2` sparsity).
- Keep structured metadata-audit table in supplement as transparent data-availability evidence.

### (e) Manuscript changes
- Supplementary reviewer-assets block around `supp_revised.tex:517-519`.
- Methods/limitations note around `main_revised.tex:361-365`.

---

## Minor Comment 1: MASLD terminology consistency

### (a) We thank the reviewer...
We agree terminology consistency is important.

### (b) Direct answer
This is **partially addressed** (MASLD-first wording exists), but a final global pass is needed.

### (c) Existing evidence pointer
- MASLD-first terminology used repeatedly (e.g., `main_revised.tex:333`).
- Legacy NAFLD-coded contexts still appear in direction labels and cohort naming (e.g., `main_revised.tex:336`, `supp_revised.tex:464`).

### (d) What we will add in revision
- Global harmonization pass: MASLD preferred term; keep "NAFLD-coded" only where needed for dataset/GWAS naming.

### (e) Manuscript changes
- Main discussion and limitation paragraphs (`main_revised.tex:333-374`)
- Supplementary tables with direction labels (`supp_revised.tex:458-466`)

---

## Minor Comment 2: Reproducibility metadata (versions/repository)

### (a) We thank the reviewer...
We agree reproducibility metadata should be explicit.

### (b) Direct answer
This is **partially addressed** (Zenodo package and supplementary data links already provided), but software-version specificity can be strengthened.

### (c) Existing evidence pointer
- Data/code availability statement: `main_revised.tex:402`
- Supplementary code/data availability: `supp_revised.tex:475-476`

### (d) What we will add in revision
- Add R/session/package version summary reference and explicit pipeline entry path in supplement/reproducibility note.

### (e) Manuscript changes
- Data/code availability block in main text (`main_revised.tex:402`)
- Supplementary reproducibility block (`supp_revised.tex:475-476`)

---

## Minor Comment 3: Common-gene count after integration

### (a) We thank the reviewer...
We agree explicit attrition counts improve interpretability.

### (b) Direct answer
This is **not fully addressed** and needs a concise quantitative statement in Methods/Results.

### (c) Existing evidence pointer
- Current text reports candidate distillation qualitatively (`main_revised.tex:219-223`) but does not explicitly report "common genes across all merged cohorts" count in one clear sentence.

### (d) What we will add in revision
- Add one sentence with exact common-gene attrition numbers and where this filtering occurred.

### (e) Manuscript changes
- Methods near transcriptomic processing and nested modeling (`main_revised.tex:172-177`)
- Results near signature derivation (`main_revised.tex:219-223`)

---

## Final Status Summary

- **Addressed:** Major Comment 2, Major Comment 3, Figure/Table Comment 1, Figure/Table Comment 2.
- **Partially addressed:** Major Comment 1, Major Comment 4, Figure/Table Comment 3.
- Remaining unresolved scope is explicitly bounded to data-availability constraints (negative-control control-sample availability and sample-level hypoxia physiology metadata completeness).

---

## Conflict Note (Phase A vs Current Manuscript)

- Phase A candidate output indicates direction flip when an exploratory legacy IHS raw-score is audited in training OOF (`IHS_direction_decision_rev11_phaseA.csv`).
- Current manuscript's primary score claim is locked 10-gene model output (Fig. 4 / Supp S9 context).
- Therefore, Phase A output is currently treated as **supplementary candidate analysis** and is **not used to overwrite** current main-text effect-size claims before definition harmonization.

