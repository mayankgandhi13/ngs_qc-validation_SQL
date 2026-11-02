-- Analysis queries for clinical-ngs-lab-db.
-- Each block starts with "-- name: <id> | <title>". Comment lines directly
-- under the name line are the lab question and get shown in the report.

-- name: qc_status | Samples failing QC
-- Which samples miss at least one QC cutoff, and which ones?
WITH t AS (
    -- pivot the threshold config table into a single row
    SELECT
        MAX(CASE WHEN metric = 'total_reads'     THEN min_value END) AS min_reads,
        MAX(CASE WHEN metric = 'pct_mapped'      THEN min_value END) AS min_mapped,
        MAX(CASE WHEN metric = 'pct_duplicate'   THEN max_value END) AS max_dup,
        MAX(CASE WHEN metric = 'mean_coverage'   THEN min_value END) AS min_cov,
        MAX(CASE WHEN metric = 'pct_target_500x' THEN min_value END) AS min_500x
    FROM qc_thresholds
),
flagged AS (
    SELECT
        s.sample_id, s.run_id, s.test_code, s.specimen_type,
        q.total_reads, q.pct_mapped, q.pct_duplicate, q.mean_coverage, q.pct_target_500x,
        TRIM(
            CASE WHEN q.total_reads     <  t.min_reads  THEN 'LOW_READS '    ELSE '' END ||
            CASE WHEN q.pct_mapped      <  t.min_mapped THEN 'LOW_MAPPING '  ELSE '' END ||
            CASE WHEN q.pct_duplicate   >  t.max_dup    THEN 'HIGH_DUP '     ELSE '' END ||
            CASE WHEN q.mean_coverage   <  t.min_cov    THEN 'LOW_COVERAGE ' ELSE '' END ||
            CASE WHEN q.pct_target_500x <  t.min_500x   THEN 'LOW_500X '     ELSE '' END
        ) AS fail_reasons
    FROM qc_metrics q
    JOIN samples s ON s.sample_id = q.sample_id
    CROSS JOIN t
)
SELECT
    sample_id, run_id, test_code, specimen_type,
    ROUND(total_reads / 1e6, 2) AS reads_m,
    pct_mapped, pct_duplicate,
    ROUND(mean_coverage, 0)     AS mean_cov,
    pct_target_500x             AS pct_500x,
    'FAIL'                      AS qc_status,
    fail_reasons
FROM flagged
WHERE fail_reasons <> ''
ORDER BY run_id, sample_id;

-- name: run_summary | Run-level QC summary
-- Is any sequencing run an outlier? A run where half the samples fail is a run problem, not a sample problem.
SELECT
    r.run_id,
    r.run_date,
    r.instrument,
    r.flowcell,
    COUNT(*)                                          AS n_samples,
    ROUND(AVG(v.mean_coverage), 0)                    AS avg_cov,
    ROUND(MIN(v.mean_coverage), 0)                    AS min_cov,
    SUM(1 - v.qc_pass)                                AS n_fail,
    ROUND(100.0 * SUM(1 - v.qc_pass) / COUNT(*), 1)   AS pct_fail,
    CASE WHEN 100.0 * SUM(1 - v.qc_pass) / COUNT(*) >= 50
         THEN 'FAIL' ELSE 'PASS' END                  AS run_status
FROM runs r
JOIN v_sample_qc v ON v.run_id = r.run_id
GROUP BY r.run_id, r.run_date, r.instrument, r.flowcell
ORDER BY r.run_date;

-- name: turnaround | Turnaround time vs target
-- Are we meeting the TAT target for each test? (received to reported, calendar days)
WITH tat AS (
    SELECT
        test_code,
        julianday(reported_date) - julianday(received_date) AS tat_days
    FROM samples
    WHERE is_validation = 0
      AND reported_date IS NOT NULL
)
SELECT
    t.test_code,
    t.test_name,
    t.tat_target_days                                             AS target_days,
    COUNT(*)                                                      AS n_reported,
    ROUND(AVG(tat.tat_days), 1)                                   AS mean_tat,
    CAST(MIN(tat.tat_days) AS INTEGER)                            AS min_tat,
    CAST(MAX(tat.tat_days) AS INTEGER)                            AS max_tat,
    SUM(tat.tat_days > t.tat_target_days)                         AS n_over_target,
    ROUND(100.0 * SUM(tat.tat_days > t.tat_target_days) / COUNT(*), 1) AS pct_over_target
FROM tests t
JOIN tat ON tat.test_code = t.test_code
GROUP BY t.test_code, t.test_name, t.tat_target_days
ORDER BY t.test_code;

-- name: reportable | Reportable variants (Tier I/II)
-- What goes on the patient reports? Tier I/II, VAF >= 5%, depth >= 250, QC-passing clinical samples only.
SELECT
    v.sample_id,
    p.mrn,
    s.test_code,
    v.gene,
    v.hgvs_p,
    v.var_type,
    ROUND(v.vaf * 100, 1) AS vaf_pct,
    v.depth,
    v.tier
FROM variants v
JOIN v_sample_qc s ON s.sample_id = v.sample_id
JOIN patients p    ON p.patient_id = s.patient_id
WHERE s.qc_pass = 1
  AND s.is_validation = 0
  AND v.tier IN ('I', 'II')
  AND v.vaf >= 0.05
  AND v.depth >= 250
ORDER BY s.test_code, v.tier, v.gene, v.sample_id;

-- name: hotspots | Recurrent variants per test
-- Which somatic variants recur most often across patients, per test? (Tier IV polymorphisms excluded)
WITH counts AS (
    SELECT
        s.test_code,
        v.gene,
        v.hgvs_p,
        v.chrom || ':' || v.pos || ' ' || v.ref || '>' || v.alt AS locus,
        COUNT(DISTINCT s.patient_id) AS n_patients,
        COUNT(*)                     AS n_calls
    FROM variants v
    JOIN samples s ON s.sample_id = v.sample_id
    WHERE s.is_validation = 0
      AND v.tier <> 'IV'
    GROUP BY s.test_code, v.gene, v.hgvs_p, v.chrom, v.pos, v.ref, v.alt
),
ranked AS (
    SELECT
        counts.*,
        RANK() OVER (PARTITION BY test_code ORDER BY n_patients DESC) AS rnk
    FROM counts
    WHERE n_patients >= 2
)
SELECT test_code, rnk, gene, hgvs_p, locus, n_patients, n_calls
FROM ranked
WHERE rnk <= 5
ORDER BY test_code, rnk, gene;

-- name: concordance | Validation concordance vs truth set
-- How well do the reference standard replicates match the known truth set? (TP/FP/FN, sensitivity, PPV)
WITH val_calls AS (
    SELECT v.*
    FROM variants v
    JOIN samples s ON s.sample_id = v.sample_id
    WHERE s.is_validation = 1
),
outcomes AS (
    -- truth side: every expected variant is either found (TP) or missed (FN)
    SELECT
        t.var_type,
        CASE WHEN c.variant_id IS NOT NULL THEN 'TP' ELSE 'FN' END AS outcome
    FROM truth_variants t
    LEFT JOIN val_calls c
           ON c.sample_id = t.sample_id AND c.chrom = t.chrom
          AND c.pos = t.pos AND c.ref = t.ref AND c.alt = t.alt
    UNION ALL
    -- call side: calls with no truth match are false positives
    SELECT c.var_type, 'FP'
    FROM val_calls c
    LEFT JOIN truth_variants t
           ON t.sample_id = c.sample_id AND t.chrom = c.chrom
          AND t.pos = c.pos AND t.ref = c.ref AND t.alt = c.alt
    WHERE t.sample_id IS NULL
),
tallies AS (
    SELECT var_type,
           SUM(outcome = 'TP') AS tp, SUM(outcome = 'FP') AS fp, SUM(outcome = 'FN') AS fn
    FROM outcomes
    GROUP BY var_type
    UNION ALL
    SELECT 'ALL',
           SUM(outcome = 'TP'), SUM(outcome = 'FP'), SUM(outcome = 'FN')
    FROM outcomes
)
SELECT
    var_type,
    tp, fp, fn,
    ROUND(100.0 * tp / NULLIF(tp + fn, 0), 1) AS sensitivity_pct,
    ROUND(100.0 * tp / NULLIF(tp + fp, 0), 1) AS ppv_pct
FROM tallies
ORDER BY var_type = 'ALL', var_type DESC;

-- name: lod | Limit of detection by VAF bin
-- Where does detection fall off? Detection rate and VAF bias (observed minus expected, detected calls only) per expected-VAF bin.
WITH val_calls AS (
    SELECT v.*
    FROM variants v
    JOIN samples s ON s.sample_id = v.sample_id
    WHERE s.is_validation = 1
),
binned AS (
    SELECT
        t.expected_vaf,
        c.vaf AS observed_vaf,
        CASE WHEN t.expected_vaf < 0.02 THEN 1
             WHEN t.expected_vaf < 0.05 THEN 2
             WHEN t.expected_vaf < 0.10 THEN 3
             ELSE 4 END AS bin_order
    FROM truth_variants t
    LEFT JOIN val_calls c
           ON c.sample_id = t.sample_id AND c.chrom = t.chrom
          AND c.pos = t.pos AND c.ref = t.ref AND c.alt = t.alt
)
SELECT
    CASE bin_order WHEN 1 THEN '<2%' WHEN 2 THEN '2-5%'
                   WHEN 3 THEN '5-10%' ELSE '>=10%' END          AS vaf_bin,
    COUNT(*)                                                     AS n_expected,
    COUNT(observed_vaf)                                          AS n_detected,
    ROUND(100.0 * COUNT(observed_vaf) / COUNT(*), 1)             AS detection_pct,
    ROUND(100 * AVG(expected_vaf), 2)                            AS mean_expected_pct,
    ROUND(100 * AVG(observed_vaf - expected_vaf), 2)             AS mean_bias_pct
FROM binned
GROUP BY bin_order
ORDER BY bin_order;

-- name: serial | Serial monitoring: VAF change over time
-- For patients tested more than once, how did each Tier I/II variant's VAF move between tests?
WITH calls AS (
    SELECT
        p.mrn, s.patient_id, s.sample_id, s.test_code, s.received_date,
        v.gene, v.hgvs_p, v.chrom, v.pos, v.ref, v.alt, v.vaf,
        COUNT(*) OVER (PARTITION BY s.patient_id, v.chrom, v.pos, v.ref, v.alt) AS n_timepoints
    FROM variants v
    JOIN samples s  ON s.sample_id = v.sample_id
    JOIN patients p ON p.patient_id = s.patient_id
    WHERE s.is_validation = 0
      AND v.tier IN ('I', 'II')
),
lagged AS (
    SELECT
        calls.*,
        LAG(vaf)           OVER w AS prev_vaf,
        LAG(received_date) OVER w AS prev_date
    FROM calls
    WHERE n_timepoints > 1
    WINDOW w AS (PARTITION BY patient_id, chrom, pos, ref, alt ORDER BY received_date)
)
SELECT
    mrn, gene, hgvs_p, test_code, sample_id, received_date,
    CAST(julianday(received_date) - julianday(prev_date) AS INTEGER) AS days_since_prev,
    ROUND(prev_vaf * 100, 1)             AS prev_vaf_pct,
    ROUND(vaf * 100, 1)                  AS vaf_pct,
    ROUND((vaf - prev_vaf) * 100, 1)     AS delta_pct,
    CASE WHEN vaf - prev_vaf <= -0.05 THEN 'falling'
         WHEN vaf - prev_vaf >=  0.05 THEN 'rising'
         ELSE 'stable' END               AS trend
FROM lagged
WHERE prev_vaf IS NOT NULL
ORDER BY mrn, gene, received_date;
