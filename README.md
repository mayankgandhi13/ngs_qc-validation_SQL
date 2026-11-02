# clinical-ngs-lab-db

A small SQLite database modeled on a hospital pathology lab running clinical NGS panels. It covers samples, sequencing runs, QC metrics, variant calls and an assay validation truth set. A set of SQL queries answers the questions a lab bioinformatician gets asked every week: which samples failed QC and why, was that a bad run, are we hitting turnaround, what's reportable, and does the assay still meet its validated sensitivity.

All data is synthetic (`random.seed(42)`), so there's no PHI anywhere. MRNs start with `SYN`. The headline hotspots (KRAS G12D/V, BRAF V600E, EGFR L858R, PIK3CA H1047R, TP53 R175H/R273H, JAK2 V617F, IDH1 R132H, DNMT3A R882H and friends) sit at their real GRCh38 coordinates. A few secondary sites (some indels, MET splice sites, TET2) use approximate positions, since they only need to look like real panel data.

## Run it

Python 3.8+ and the standard library. Nothing to pip install.

```
python run.py
```

This will:

1. delete and rebuild `ngs_lab.db` from `schema.sql`
2. seed it with `seed_data.py` (60 patients, 12 weekly runs from 2026-07-06, about 120 samples)
3. run every query in `queries.sql` and print each result as a table in the terminal
4. write `report.html`, a single self-contained page with one section per query (open it in any browser)

To poke around afterwards: `sqlite3 ngs_lab.db`.

## What's in the data

- **3 tests:** ONCO-DNA (solid tumor hotspot, 10 day TAT target), HEME-DNA (myeloid, 14 days), FUSION-RNA (12 days)
- **FFPE vs blood:** FFPE samples get lower coverage and higher duplication, the way degraded input really behaves
- **One bad run:** run R260824 comes in at about 45% of normal coverage, so most of its samples fail QC and its reports go out late
- **Serial patients:** myeloid patients get follow-up tests with falling (or occasionally rising) VAFs, and a few solid tumor patients get a re-biopsy. EGFR-mutant re-biopsies pick up T790M.
- **Validation:** 3 replicates of a reference standard on ONCO-DNA, run on 3 different weeks. The truth set has 20 variants with expected VAFs from 1% to 35%. Detection drops off below 5%, and one replicate has a single low-VAF false positive.
- **QC cutoffs live in a table** (`qc_thresholds`), not in code: reads >= 2M, mapped >= 95%, duplicates <= 40%, mean coverage >= 500x, >= 90% of targets at 500x

## Schema

```
 patients                 tests                      runs
 ---------                ---------------            -----------
 patient_id PK            test_code PK  <----+       run_id PK
 mrn UNIQUE                test_name          |       instrument
 sex                      assay_type         |       run_date
 birth_year               panel_version      |       flowcell
     ^                    tat_target_days    |          ^
     |                        ^              |          |
     |                        |        panel_genes      |
     |                        |        (test_code, gene)|
     |                        |                         |
     |      samples           |                         |
     |      ---------------------------------           |
     +----- patient_id FK (NULL for ref standards)      |
            test_code FK ----+                          |
            run_id FK ------------------------------------+
            sample_id PK
            specimen_type, is_validation, tumor_purity
            received_date, reported_date
                 |
     +-----------+----------------+----------------------+
     |                            |                      |
 qc_metrics (1:1)              variants (1:many)      truth_variants
 ----------------              ------------------     ---------------
 sample_id PK/FK               variant_id PK          sample_id FK
 total_reads                   sample_id FK           chrom, pos, ref, alt
 pct_mapped                    chrom, pos, ref, alt   gene, var_type
 pct_duplicate                 gene, hgvs_p           expected_vaf
 mean_coverage                 var_type, vaf, depth
 pct_target_500x               tier (I-IV)

 qc_thresholds (config)        v_sample_qc (view)
 ----------------------        ------------------
 metric PK                     samples + qc_metrics + qc_pass flag
 min_value, max_value          built from qc_thresholds
```

## Queries

| id | Lab question | SQL concepts |
|----|--------------|--------------|
| `qc_status` | Which samples failed QC, and on which metrics? | CTE, pivot with `MAX(CASE ...)`, `CROSS JOIN` to a one-row config, string building with `\|\|` and `TRIM` |
| `run_summary` | Is a whole run bad, or just a few samples? | `GROUP BY` over a view, conditional aggregation, percent-of-total |
| `turnaround` | Are we meeting TAT targets per test? | `julianday()` date math, `SUM(boolean)` counting |
| `reportable` | What goes on patient reports? | Multi-table joins, filtering on tier, VAF, depth and QC pass |
| `hotspots` | What recurs most across patients? | `COUNT(DISTINCT)`, `RANK() OVER (PARTITION BY ...)`, top N per group |
| `concordance` | Do validation calls match the truth set? | Full outer join emulated with `LEFT JOIN` + `UNION ALL`, `NULLIF` for safe division, a hand-rolled ALL row |
| `lod` | At what VAF does detection fall off? | `CASE` binning, `LEFT JOIN` to keep misses, `COUNT(col)` vs `COUNT(*)` |
| `serial` | How did a patient's VAF move between tests? | `LAG()`, named `WINDOW`, window `COUNT` to keep only serial variants |

`concordance` and `lod` mirror the accuracy and LOD tables in a CAP/CLIA validation summary: TP/FP/FN per variant class, sensitivity and PPV, plus detection rate by VAF bin across replicate runs of a reference material.

## Files

- `schema.sql`: tables, foreign keys, indexes, and the `v_sample_qc` view
- `seed_data.py`: synthetic data generator
- `queries.sql`: the 8 analysis queries, each block starts with `-- name: <id> | <title>`
- `run.py`: rebuild, seed, run, print, write the report

## Ideas to extend

- Load a real GIAB HG002 VCF and its truth set into `variants` / `truth_variants` and rerun `concordance` against real benchmark data
- Port the schema to PostgreSQL (real `FULL OUTER JOIN`, `FILTER (WHERE ...)`, proper date types)
- A small Flask page for looking up a sample or MRN and seeing its QC and calls
