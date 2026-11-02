-- clinical-ngs-lab-db schema (SQLite)
-- All data is synthetic. VAF values are stored as fractions (0.05 = 5%).

PRAGMA foreign_keys = ON;

CREATE TABLE patients (
    patient_id  INTEGER PRIMARY KEY,
    mrn         TEXT    NOT NULL UNIQUE,
    sex         TEXT    NOT NULL CHECK (sex IN ('F', 'M')),
    birth_year  INTEGER NOT NULL CHECK (birth_year BETWEEN 1900 AND 2026)
);

CREATE TABLE tests (
    test_code        TEXT    PRIMARY KEY,
    test_name        TEXT    NOT NULL,
    assay_type       TEXT    NOT NULL CHECK (assay_type IN ('DNA', 'RNA')),
    panel_version    TEXT    NOT NULL,
    tat_target_days  INTEGER NOT NULL CHECK (tat_target_days > 0)
);

CREATE TABLE panel_genes (
    test_code  TEXT NOT NULL REFERENCES tests (test_code),
    gene       TEXT NOT NULL,
    PRIMARY KEY (test_code, gene)
);

CREATE TABLE runs (
    run_id      TEXT PRIMARY KEY,
    instrument  TEXT NOT NULL CHECK (instrument IN ('NovaSeq X', 'NextSeq 2000')),
    run_date    TEXT NOT NULL,              -- ISO date
    flowcell    TEXT NOT NULL UNIQUE
);

CREATE TABLE samples (
    sample_id      TEXT    PRIMARY KEY,
    patient_id     INTEGER REFERENCES patients (patient_id),   -- NULL for reference standards
    run_id         TEXT    NOT NULL REFERENCES runs (run_id),
    test_code      TEXT    NOT NULL REFERENCES tests (test_code),
    specimen_type  TEXT    NOT NULL CHECK (specimen_type IN ('FFPE', 'blood', 'reference_standard')),
    is_validation  INTEGER NOT NULL DEFAULT 0 CHECK (is_validation IN (0, 1)),
    tumor_purity   REAL    CHECK (tumor_purity BETWEEN 0 AND 1),
    received_date  TEXT    NOT NULL,
    reported_date  TEXT,
    CHECK (reported_date IS NULL OR reported_date >= received_date),
    CHECK ((is_validation = 1 AND patient_id IS NULL)
        OR (is_validation = 0 AND patient_id IS NOT NULL))
);

CREATE TABLE qc_metrics (
    sample_id        TEXT    PRIMARY KEY REFERENCES samples (sample_id),
    total_reads      INTEGER NOT NULL CHECK (total_reads >= 0),
    pct_mapped       REAL    NOT NULL CHECK (pct_mapped BETWEEN 0 AND 100),
    pct_duplicate    REAL    NOT NULL CHECK (pct_duplicate BETWEEN 0 AND 100),
    mean_coverage    REAL    NOT NULL CHECK (mean_coverage >= 0),
    pct_target_500x  REAL    NOT NULL CHECK (pct_target_500x BETWEEN 0 AND 100)
);

-- QC cutoffs live in data, not code. metric matches a qc_metrics column name.
CREATE TABLE qc_thresholds (
    metric     TEXT PRIMARY KEY,
    min_value  REAL,
    max_value  REAL,
    CHECK (min_value IS NOT NULL OR max_value IS NOT NULL)
);

CREATE TABLE variants (
    variant_id  INTEGER PRIMARY KEY,
    sample_id   TEXT    NOT NULL REFERENCES samples (sample_id),
    chrom       TEXT    NOT NULL,
    pos         INTEGER NOT NULL,
    ref         TEXT    NOT NULL,
    alt         TEXT    NOT NULL,
    gene        TEXT    NOT NULL,
    hgvs_p      TEXT,
    var_type    TEXT    NOT NULL CHECK (var_type IN ('SNV', 'INDEL')),
    vaf         REAL    NOT NULL CHECK (vaf BETWEEN 0 AND 1),
    depth       INTEGER NOT NULL CHECK (depth > 0),
    tier        TEXT    NOT NULL CHECK (tier IN ('I', 'II', 'III', 'IV')),   -- AMP/ASCO/CAP 2017
    UNIQUE (sample_id, chrom, pos, ref, alt)
);

-- Known variants in a reference standard, one row per validation sample.
CREATE TABLE truth_variants (
    sample_id     TEXT    NOT NULL REFERENCES samples (sample_id),
    chrom         TEXT    NOT NULL,
    pos           INTEGER NOT NULL,
    ref           TEXT    NOT NULL,
    alt           TEXT    NOT NULL,
    gene          TEXT    NOT NULL,
    var_type      TEXT    NOT NULL CHECK (var_type IN ('SNV', 'INDEL')),
    expected_vaf  REAL    NOT NULL CHECK (expected_vaf BETWEEN 0 AND 1),
    PRIMARY KEY (sample_id, chrom, pos, ref, alt)
);

CREATE INDEX idx_samples_run      ON samples (run_id);
CREATE INDEX idx_samples_patient  ON samples (patient_id);
CREATE INDEX idx_samples_test     ON samples (test_code);
CREATE INDEX idx_variants_sample  ON variants (sample_id);
CREATE INDEX idx_variants_locus   ON variants (chrom, pos, ref, alt);
CREATE INDEX idx_variants_gene    ON variants (gene);

-- One row per sequenced sample with its QC metrics and an overall pass flag.
-- Same threshold pivot as the qc_status query, so every query shares one
-- definition of "QC pass".
CREATE VIEW v_sample_qc AS
WITH t AS (
    SELECT
        MAX(CASE WHEN metric = 'total_reads'     THEN min_value END) AS min_reads,
        MAX(CASE WHEN metric = 'pct_mapped'      THEN min_value END) AS min_mapped,
        MAX(CASE WHEN metric = 'pct_duplicate'   THEN max_value END) AS max_dup,
        MAX(CASE WHEN metric = 'mean_coverage'   THEN min_value END) AS min_cov,
        MAX(CASE WHEN metric = 'pct_target_500x' THEN min_value END) AS min_500x
    FROM qc_thresholds
)
SELECT
    s.*,
    q.total_reads, q.pct_mapped, q.pct_duplicate, q.mean_coverage, q.pct_target_500x,
    CASE WHEN q.total_reads     >= t.min_reads
          AND q.pct_mapped      >= t.min_mapped
          AND q.pct_duplicate   <= t.max_dup
          AND q.mean_coverage   >= t.min_cov
          AND q.pct_target_500x >= t.min_500x
         THEN 1 ELSE 0 END AS qc_pass
FROM samples s
JOIN qc_metrics q ON q.sample_id = s.sample_id
CROSS JOIN t;
