# Sensityping Metrics — Step-by-Step Vignette

**Sensityping Metrics** is a command-line evaluation framework for assessing **genotype-based antimicrobial susceptibility predictions** against phenotypic treatment outcomes.

It is designed as a **companion tool** to SensiTyper, providing rigorous performance evaluation for genomic AMR prediction pipelines for gonococcal (and similar) surveillance datasets.

---

## Overview

The script runs in **two complementary modes** that answer fundamentally different evaluation questions:

| Question | Analysis mode |
|--------|--------------|
| How accurate is the prediction for each antibiotic? | `predicted_vs_treatment` |
| How well does the final treatment recommendation perform in practice? | `first_line_vs_treatment` |

These cannot be answered in a single analysis without mixing incompatible assumptions. Both modes are typically required for a complete evaluation.

---

## Requirements

The script is standalone and requires only standard scientific Python:

- python >= 3.6
- numpy
- pandas
- scipy
- scikit-learn
- plotly

---

## Input format

The input is a **tab-delimited table** containing, per isolate:

- Genomic susceptibility predictions from SensiTyper or equivalent (`*_predicted`)
- Observed outcome from a gold-standard method (e.g. MIC testing) (`*_treatment`)
- Optional recommendation columns from SensiTyper (`*_recommend`)

### Example input table (`example_data.tab`)

```text
isolate_id	CRO_predicted	CRO_treatment	AZM_predicted	AZM_treatment	CRO+AZM_recommend	CRO+AZM_treatment
WHO_A	S	S	S	S	YES	YES
WHO_B	S	S	R	R	NO	NO
WHO_C	S	S	S	S	YES	YES
WHO_D	R	R	S	S	NO	NO
WHO_E	S	S	R	R	NO	NO
WHO_F	S	S	S	S	YES	YES
WHO_Q	R	R	R	R	NO	NO
WHO_Z	R	R	S	S	NO	NO
```

Dual-therapy tokens (`CRO+AZM`, `AZM+SPC`, etc.) are supported **only if present as explicit columns**.

---

## Step 1 — Model-centric evaluation (`predicted_vs_treatment`)

This mode evaluates **each antibiotic independently**, asking:

> *Does the genomic predictor correctly classify susceptibility for this drug?*

### Command

```bash
python sensityping_metrics.py \
  -i example_data.tab \
  -o example_predict.tab \
  -d example_plots \
  --analysis_type predicted_vs_treatment \
  --ci_flag \
  --ci_method hybrid \
  --ci_level 0.95 \
  --n_boot 2000 \
  --seed 1 \
  --ssd_flag \
  --ssd_width 0.05 \
  --ssd_mode observed \
  --radar_flag \
  --radar_metrics PPV,one_minus_FDR,coverage_fraction
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-i` | Input tab-delimited file |
| `-o` | Output metrics summary table |
| `-d` | Output directory for plots and extracted ID files |
| `--analysis_type` | `predicted_vs_treatment` for model-centric evaluation |
| `--ci_flag` | Enable confidence interval calculation |
| `--ci_method` | CI method: `wilson`, `bootstrap`, or `hybrid` (recommended) |
| `--ci_level` | CI level (default: 0.95) |
| `--n_boot` | Number of bootstrap replicates for composite metrics |
| `--seed` | Random seed for reproducibility |
| `--ssd_flag` | Enable sample size diagnostics |
| `--ssd_width` | Target CI half-width for SSD calculation |
| `--ssd_mode` | SSD mode: `observed` (use observed prevalence) or `worst` |
| `--radar_flag` | Generate radar plot |
| `--radar_metrics` | Comma-separated metric names for radar axes |

### Expected output (console)

```text
Metrics for CRO (Prediction vs. Treatment):
Total Isolates Analyzed: 28
note: OK
tp: 23
tn: 5
fp: 0
fn: 0
assigned_count: 23
accuracy: 1.0
concordance: 1.0
sensitivity: 1.0
specificity: 1.0
PPV: 1.0
NPV: 1.0
f1_score: 1.0
mcc: 1.0
auc: 1.0
balanced_accuracy: 1.0
cost_sensitive_error_rate: 0.0
kappa: 1.0
ME_rate: 0.0
VME_rate: 0.0
ME_denom: 23
VME_denom: 5
one_minus_ME: 1.0
one_minus_VME: 1.0
coverage_fraction: 0.82
one_minus_FDR: 1.0
PPV_ci_lower: 0.86
PPV_ci_upper: 1.0
sensitivity_ci_lower: 0.85
sensitivity_ci_upper: 1.0
ssd_PPV: 42
ssd_sensitivity: 38

Metrics for AZM (Prediction vs. Treatment):
Total Isolates Analyzed: 28
...
```

### Interpreting model-centric metrics

- **sensitivity** — ability to detect resistant isolates (avoids very major errors)
- **specificity** — ability to correctly classify susceptible isolates
- **PPV** — when predicting resistance, probability that resistance is real (positive predictive value)
- **NPV** — when predicting susceptibility, probability that the isolate is truly susceptible
- **ME_rate** — major error rate (resistant predicted as susceptible); highest-consequence error
- **VME_rate** — very major error rate (susceptible predicted as resistant)
- **coverage_fraction** — fraction of isolates with an unambiguous prediction
- **ssd_*** — minimum sample size needed to achieve the target CI half-width for that metric

---

## Step 2 — Clinical workflow evaluation (`first_line_vs_treatment`)

This mode evaluates **what matters clinically**:

> *Did the recommended first-line treatment work?*

Key characteristics:

- Treatments are evaluated **sequentially** following the priority order
- Once a regimen is assigned, downstream alternatives are not evaluated
- **Ordering matters** — reflects guideline-based decision logic

You must explicitly define:

- The **order** of treatment consideration (`--order`)
- Which recommendation tokens to extract from the input table (`--id_extraction`)

### Command

```bash
python sensityping_metrics.py \
  -i example_data.tab \
  -o example_firstline.tab \
  -d example_out \
  --analysis_type first_line_vs_treatment \
  --order 'CRO+AZM,CRO,AZM,AZM+SPC,SPC,ZOL' \
  --id_extraction CRO+AZM,CRO,AZM,AZM+SPC,SPC,ZOL \
  --ci_flag \
  --ci_method hybrid \
  --ci_level 0.95 \
  --n_boot 2000 \
  --seed 1 \
  --ssd_flag \
  --ssd_width 0.05 \
  --ssd_mode observed \
  --radar_flag \
  --radar_metrics PPV,one_minus_FDR,coverage_fraction
```

### Arguments specific to this mode

| Argument | Description |
|----------|-------------|
| `--order` | Comma-separated treatment priority order (quoted string) |
| `--id_extraction` | Comma-separated list of regimens to extract isolate IDs for |

### Expected output (console)

```text
Metrics for CRO+AZM (First Line vs. Treatment):
Total Isolates Analyzed: 28
note: OK
tp: 18
tn: 10
fp: 0
fn: 0
assigned_count: 18
accuracy: 1.0
concordance: 1.0
sensitivity: 1.0
specificity: 1.0
PPV: 1.0
NPV: 1.0
f1_score: 1.0
mcc: 1.0
...

Metrics for CRO (First Line vs. Treatment):
Total Isolates Analyzed: 10
note: OK
...
```

### Interpreting clinical workflow metrics

In this mode, **each regimen block is evaluated only for the isolates assigned to it**. An isolate appears in only one regimen group — the first eligible one in the priority order.

- **PPV** for `CRO+AZM` — of isolates recommended for dual therapy, how many actually responded?
- **coverage_fraction** — fraction of isolates receiving any recommendation in this block
- **tn** — isolates correctly identified as not requiring this specific regimen (passed down the cascade)

> **Important:** Isolates assigned to `CRO` appear in the `CRO` block only because they were not assigned `CRO+AZM` first. Each block represents a sequential decision step.

---

## Step 3 — Confidence intervals

Confidence intervals are enabled with `--ci_flag`. The `hybrid` method (recommended) automatically applies:

- **Wilson score** for proportions: PPV, NPV, sensitivity, specificity, ME_rate, VME_rate
- **Parametric bootstrap** for composite metrics: F1, MCC, kappa, balanced accuracy, AUC

```bash
--ci_flag --ci_method hybrid --ci_level 0.95 --n_boot 2000 --seed 1
```

CI columns are appended to the output table as `{metric}_ci_lower` and `{metric}_ci_upper`.

---

## Step 4 — Sample size diagnostics (SSD)

SSD estimates the **minimum sample size** required to achieve a CI half-width of ±w for each metric.

```bash
--ssd_flag --ssd_width 0.05 --ssd_mode observed
```

| `--ssd_mode` | Description |
|-------------|-------------|
| `observed` | Uses the observed proportion as the true value |
| `worst` | Uses 0.5 (most conservative, largest sample size) |

SSD columns appear as `ssd_{metric}` in the output table. Useful for:

- Surveillance planning
- Power justification in manuscript methods sections
- Protocol design for prospective evaluations

---

## Step 5 — Radar plots

Radar plots provide a visual multi-metric summary per antibiotic or regimen.

```bash
--radar_flag --radar_metrics PPV,one_minus_FDR,coverage_fraction
```

Output is an **interactive Plotly HTML** file, suitable for:

- Supplementary material
- Internal dashboards
- Presentations

Recommended metric sets:

- **Model evaluation**: `sensitivity,specificity,PPV,NPV,coverage_fraction`
- **Clinical evaluation**: `PPV,one_minus_FDR,coverage_fraction,one_minus_ME`
- **Full panel**: `PPV,NPV,sensitivity,specificity,one_minus_ME,one_minus_VME,coverage_fraction`

---

## Output structure

```text
example_out/
├── metrics_summary.tab          # All metrics per antibiotic or regimen
├── CRO+AZM_id_extracted.tab    # Isolate IDs assigned to CRO+AZM block
├── CRO_id_extracted.tab        # Isolate IDs assigned to CRO block
├── AZM_id_extracted.tab        # Isolate IDs assigned to AZM block
└── radar_plot.html             # Interactive radar plot
```

The `metrics_summary.tab` is a tab-delimited table with one row per antibiotic/regimen and one column per metric (plus CI and SSD columns when enabled).

---

## Why both modes are required

| Aspect | `predicted_vs_treatment` | `first_line_vs_treatment` |
|------|------------------------|-------------------------|
| Evaluates model correctness | ✓ | ✗ |
| Evaluates clinical outcome | ✗ | ✓ |
| Per-antibiotic performance | ✓ | ✗ |
| Sequential decision logic | ✗ | ✓ |
| Reflects guideline practice | ✗ | ✓ |
| Suitable for methods papers | ✓ | ✗ |
| Suitable for policy / implementation | ✗ | ✓ |

Running only one would give a **biased or incomplete interpretation**.

---

## Complete example — both modes

```bash
# Mode 1: per-antibiotic model performance
python sensityping_metrics.py \
  -i sensiscript_with_phenotypes.tab \
  -o metrics_predicted.tab \
  -d metrics_predicted_dir \
  --analysis_type predicted_vs_treatment \
  --ci_flag --ci_method hybrid --ci_level 0.95 \
  --n_boot 2000 --seed 1 \
  --ssd_flag --ssd_width 0.05 --ssd_mode observed \
  --radar_flag \
  --radar_metrics PPV,sensitivity,specificity,coverage_fraction

# Mode 2: clinical cascade performance
python sensityping_metrics.py \
  -i sensiscript_with_phenotypes.tab \
  -o metrics_firstline.tab \
  -d metrics_firstline_dir \
  --analysis_type first_line_vs_treatment \
  --order 'CRO+AZM,CRO,AZM,AZM+SPC,SPC,ZOL' \
  --id_extraction CRO+AZM,CRO,AZM,AZM+SPC,SPC,ZOL \
  --ci_flag --ci_method hybrid --ci_level 0.95 \
  --n_boot 2000 --seed 1 \
  --ssd_flag --ssd_width 0.05 --ssd_mode observed \
  --radar_flag \
  --radar_metrics PPV,one_minus_FDR,coverage_fraction
```

---

## Intended use

This tool is designed for:

- Genomic AMR prediction pipelines
- Rule-based or ML-based susceptibility systems
- Surveillance frameworks (WHO / EUCAST-style)
- Manuscripts requiring **transparent, reproducible evaluation**

---

## Citation

If you use this framework in a manuscript, please cite the corresponding **Sensityping / SensiTyper** publication and reference this repository.
