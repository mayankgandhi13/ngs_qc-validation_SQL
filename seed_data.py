"""Synthetic data generator for clinical-ngs-lab-db.

Everything here is made up. MRNs use a SYN prefix and no record maps to a
real person. Headline hotspot coordinates are real GRCh38 positions; a few
secondary sites (some indels, MET splice sites, TET2) are approximate and
only there to make the data look like a real panel.

VAFs are fractions (0.05 = 5%).
"""

import math
import random
from datetime import date, timedelta

SEED = 42
FIRST_RUN = date(2026, 7, 6)
N_RUNS = 12
N_PATIENTS = 60
MIN_PER_RUN, MAX_PER_RUN = 8, 11
BAD_RUN = 7                       # index of the run with ~45% of normal coverage
VALIDATION_RUNS = (2, 3, 4)       # one reference standard replicate per run
FP_REPLICATE = 1                  # which replicate gets the low-VAF false positive
CALLER_MIN_VAF = 0.01
CALLER_MIN_ALT_READS = 5

TESTS = [
    ("ONCO-DNA", "Solid tumor hotspot panel", "DNA", "v3.2", 10),
    ("HEME-DNA", "Myeloid panel", "DNA", "v2.1", 14),
    ("FUSION-RNA", "Solid tumor fusion panel", "RNA", "v1.4", 12),
]

PANEL_GENES = {
    "ONCO-DNA": ["APC", "BRAF", "CTNNB1", "EGFR", "ERBB2", "FGFR3", "GNAS",
                 "IDH1", "KIT", "KRAS", "NRAS", "PIK3CA", "PTEN", "TP53"],
    "HEME-DNA": ["ASXL1", "CEBPA", "DNMT3A", "FLT3", "IDH1", "IDH2", "JAK2",
                 "KIT", "KRAS", "NPM1", "NRAS", "RUNX1", "SRSF2", "TET2", "TP53"],
    "FUSION-RNA": ["ALK", "BRAF", "FGFR2", "FGFR3", "MET", "NTRK1", "NTRK2",
                   "NTRK3", "RET", "ROS1"],
}

# Approximate GRCh38 gene spans, used only to place synthetic VUS calls.
GENE_SPAN = {
    "APC": ("chr5", 112707498, 112846239), "BRAF": ("chr7", 140719327, 140924929),
    "CTNNB1": ("chr3", 41194741, 41260096), "EGFR": ("chr7", 55019017, 55211628),
    "ERBB2": ("chr17", 39687914, 39730426), "FGFR3": ("chr4", 1793293, 1808872),
    "GNAS": ("chr20", 58839718, 58911192), "IDH1": ("chr2", 208236227, 208255071),
    "KIT": ("chr4", 54657918, 54740715), "KRAS": ("chr12", 25205246, 25250929),
    "NRAS": ("chr1", 114704469, 114716771), "PIK3CA": ("chr3", 179148114, 179240093),
    "PTEN": ("chr10", 87863113, 87971930), "TP53": ("chr17", 7661779, 7687538),
    "ASXL1": ("chr20", 32358062, 32439319), "CEBPA": ("chr19", 33299934, 33302534),
    "DNMT3A": ("chr2", 25227855, 25342590), "FLT3": ("chr13", 28003274, 28100592),
    "IDH2": ("chr15", 90083045, 90102504), "JAK2": ("chr9", 4985033, 5128183),
    "NPM1": ("chr5", 171387116, 171411810), "RUNX1": ("chr21", 34787801, 36004667),
    "SRSF2": ("chr17", 76734115, 76737374), "TET2": ("chr4", 105145875, 105279816),
}

# key: (gene, chrom, pos, ref, alt, hgvs_p, var_type, tier)
HOTSPOTS = {
    # solid tumor
    "KRAS_G12D":     ("KRAS", "chr12", 25245350, "C", "T", "p.G12D", "SNV", "I"),
    "KRAS_G12V":     ("KRAS", "chr12", 25245350, "C", "A", "p.G12V", "SNV", "I"),
    "KRAS_G13D":     ("KRAS", "chr12", 25245347, "C", "T", "p.G13D", "SNV", "I"),
    "NRAS_Q61K":     ("NRAS", "chr1", 114713909, "G", "T", "p.Q61K", "SNV", "I"),
    "BRAF_V600E":    ("BRAF", "chr7", 140753336, "A", "T", "p.V600E", "SNV", "I"),
    "EGFR_L858R":    ("EGFR", "chr7", 55191822, "T", "G", "p.L858R", "SNV", "I"),
    "EGFR_T790M":    ("EGFR", "chr7", 55181378, "C", "T", "p.T790M", "SNV", "I"),
    "EGFR_EX19DEL":  ("EGFR", "chr7", 55174771, "AGGAATTAAGAGAAGC", "A", "p.E746_A750del", "INDEL", "I"),
    "EGFR_EX20INS":  ("EGFR", "chr7", 55181319, "C", "CAGCGTGGAC", "p.A767_V769dup", "INDEL", "I"),
    "ERBB2_V777L":   ("ERBB2", "chr17", 39724747, "G", "T", "p.V777L", "SNV", "I"),
    "PIK3CA_H1047R": ("PIK3CA", "chr3", 179234297, "A", "G", "p.H1047R", "SNV", "I"),
    "PIK3CA_E545K":  ("PIK3CA", "chr3", 179218303, "G", "A", "p.E545K", "SNV", "I"),
    "PIK3CA_E542K":  ("PIK3CA", "chr3", 179218294, "G", "A", "p.E542K", "SNV", "I"),
    "TP53_R175H":    ("TP53", "chr17", 7675088, "C", "T", "p.R175H", "SNV", "II"),
    "TP53_R273H":    ("TP53", "chr17", 7673802, "C", "T", "p.R273H", "SNV", "II"),
    "TP53_Y220C":    ("TP53", "chr17", 7674872, "T", "C", "p.Y220C", "SNV", "II"),
    "TP53_G245S":    ("TP53", "chr17", 7674230, "C", "T", "p.G245S", "SNV", "II"),
    "TP53_R248Q":    ("TP53", "chr17", 7674220, "C", "T", "p.R248Q", "SNV", "II"),
    "TP53_R248W":    ("TP53", "chr17", 7674221, "G", "A", "p.R248W", "SNV", "II"),
    "TP53_R282W":    ("TP53", "chr17", 7673776, "G", "A", "p.R282W", "SNV", "II"),
    "APC_R1450X":    ("APC", "chr5", 112839942, "C", "T", "p.R1450*", "SNV", "II"),
    "APC_T1556FS":   ("APC", "chr5", 112840254, "T", "TA", "p.T1556Nfs*3", "INDEL", "II"),
    "FGFR3_S249C":   ("FGFR3", "chr4", 1801841, "C", "G", "p.S249C", "SNV", "I"),
    "IDH1_R132H":    ("IDH1", "chr2", 208248388, "C", "T", "p.R132H", "SNV", "I"),
    "KIT_D816V":     ("KIT", "chr4", 54733155, "A", "T", "p.D816V", "SNV", "I"),
    "GNAS_R201C":    ("GNAS", "chr20", 58909365, "C", "T", "p.R201C", "SNV", "II"),
    "CTNNB1_S33C":   ("CTNNB1", "chr3", 41224610, "C", "G", "p.S33C", "SNV", "II"),
    # myeloid
    "JAK2_V617F":    ("JAK2", "chr9", 5073770, "G", "T", "p.V617F", "SNV", "I"),
    "DNMT3A_R882H":  ("DNMT3A", "chr2", 25234373, "C", "T", "p.R882H", "SNV", "II"),
    "IDH2_R140Q":    ("IDH2", "chr15", 90088702, "C", "T", "p.R140Q", "SNV", "I"),
    "FLT3_ITD":      ("FLT3", "chr13", 28034105, "C", "CGAAGTACTCATTATCTGAGGA", "p.E598_Y599ins7", "INDEL", "I"),
    "NPM1_W288FS":   ("NPM1", "chr5", 171410539, "C", "CTCTG", "p.W288Cfs*12", "INDEL", "I"),
    "SRSF2_P95H":    ("SRSF2", "chr17", 76736877, "G", "T", "p.P95H", "SNV", "II"),
    "ASXL1_G646FS":  ("ASXL1", "chr20", 32434638, "A", "AG", "p.G646Wfs*12", "INDEL", "II"),
    "TET2_Q1034X":   ("TET2", "chr4", 105235829, "C", "T", "p.Q1034*", "SNV", "II"),
    # RNA panel (MET exon 14 skipping splice-site variants)
    "MET_EX14_D1":   ("MET", "chr7", 116771989, "G", "T", "c.3082+1G>T (exon 14 skipping)", "SNV", "I"),
    "MET_EX14_D2":   ("MET", "chr7", 116771990, "T", "C", "c.3082+2T>C (exon 14 skipping)", "SNV", "I"),
    # common germline polymorphism, benign
    "TP53_P72R":     ("TP53", "chr17", 7676154, "G", "C", "p.P72R", "SNV", "IV"),
}

# Reference standard truth set for ONCO-DNA validation: expected VAF 1% to 35%.
TRUTH_SET = [
    ("KRAS_G13D", 0.010), ("EGFR_T790M", 0.010), ("CTNNB1_S33C", 0.015), ("KIT_D816V", 0.018),
    ("EGFR_EX19DEL", 0.020), ("EGFR_L858R", 0.030), ("FGFR3_S249C", 0.040), ("GNAS_R201C", 0.045),
    ("KRAS_G12D", 0.050), ("ERBB2_V777L", 0.065), ("EGFR_EX20INS", 0.070), ("APC_R1450X", 0.080),
    ("PIK3CA_E545K", 0.090), ("BRAF_V600E", 0.105), ("NRAS_Q61K", 0.125), ("APC_T1556FS", 0.150),
    ("PIK3CA_H1047R", 0.175), ("TP53_R175H", 0.250), ("IDH1_R132H", 0.300), ("TP53_R273H", 0.350),
]
FALSE_POSITIVE = ("PIK3CA_E542K", 0.013)
SUBCLONES = ["TP53_Y220C", "TP53_G245S", "TP53_R248Q", "TP53_R248W", "TP53_R282W"]

SOLID_DRIVERS = [("KRAS_G12D", 3), ("KRAS_G12V", 2), ("BRAF_V600E", 2), ("EGFR_L858R", 2),
                 ("EGFR_EX19DEL", 2), ("PIK3CA_H1047R", 2), ("FGFR3_S249C", 1),
                 ("IDH1_R132H", 1), (None, 2)]
HEME_DRIVERS = [("JAK2_V617F", 3), ("DNMT3A_R882H", 3), ("NPM1_W288FS", 2), ("FLT3_ITD", 2),
                ("TET2_Q1034X", 2), ("IDH1_R132H", 1), ("IDH2_R140Q", 1), ("SRSF2_P95H", 1),
                ("ASXL1_G646FS", 1), ("TP53_R273H", 1)]

AMINO = "ACDEFGHIKLMNPQRSTVWY"


def weighted_choice(pairs):
    items, weights = zip(*pairs)
    return random.choices(items, weights=weights, k=1)[0]


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def random_vus(gene):
    chrom, start, end = GENE_SPAN[gene]
    ref = random.choice("ACGT")
    alt = random.choice([b for b in "ACGT" if b != ref])
    hgvs = f"p.{random.choice(AMINO)}{random.randint(20, 900)}{random.choice(AMINO)}"
    return (gene, chrom, random.randint(start, end), ref, alt, hgvs, "SNV", "III")


# --------------------------------------------------------------------------
# patients and their test orders
# --------------------------------------------------------------------------

def make_patients():
    patients, used_mrns = [], set()
    for pid in range(1, N_PATIENTS + 1):
        cohort = "heme" if random.random() < 0.4 else "solid"
        mrn = f"SYN{random.randint(1000000, 9999999)}"
        while mrn in used_mrns:
            mrn = f"SYN{random.randint(1000000, 9999999)}"
        used_mrns.add(mrn)
        birth_year = random.randint(1940, 1975) if cohort == "heme" else random.randint(1945, 1995)
        patients.append({"patient_id": pid, "mrn": mrn, "sex": random.choice("FM"),
                         "birth_year": birth_year, "cohort": cohort})
    return patients


def plan_orders(patient):
    """List of (week_offset, test_code, specimen_type) for one patient."""
    if patient["cohort"] == "solid":
        orders = [(0, "ONCO-DNA", "FFPE")]
        if random.random() < 0.45:
            orders.append((0, "FUSION-RNA", "FFPE"))        # DNA + RNA off the same block
        if random.random() < 0.25:
            orders.append((random.randint(5, 7), "ONCO-DNA", "FFPE"))   # re-biopsy at progression
    else:
        orders = [(0, "HEME-DNA", "blood")]
        r = random.random()
        if r < 0.65:
            first = random.randint(3, 5)
            orders.append((first, "HEME-DNA", "blood"))     # post-induction follow-up
            if r < 0.35:
                orders.append((first + random.randint(3, 5), "HEME-DNA", "blood"))
    return orders


def schedule(patients):
    """Pick a start week per patient so every run ends up with 8 to 11 samples."""
    for _ in range(200):
        plans = {p["patient_id"]: plan_orders(p) for p in patients}
        loads = [0] * N_RUNS
        for w in VALIDATION_RUNS:
            loads[w] += 1
        order = list(patients)
        random.shuffle(order)
        order.sort(key=lambda p: -max(o[0] for o in plans[p["patient_id"]]))
        starts, ok = {}, True
        for p in order:
            offsets = [o[0] for o in plans[p["patient_id"]]]
            best = None
            for w in range(N_RUNS - max(offsets)):
                trial = loads[:]
                for o in offsets:
                    trial[w + o] += 1
                if max(trial) > MAX_PER_RUN:
                    continue
                score = (max(trial[w + o] for o in offsets), sum(trial[w + o] for o in offsets), random.random())
                if best is None or score < best[0]:
                    best = (score, w, trial)
            if best is None:
                ok = False
                break
            starts[p["patient_id"]], loads = best[1], best[2]
        if ok and min(loads) >= MIN_PER_RUN:
            return plans, starts
    raise RuntimeError("could not build a balanced run schedule")


# --------------------------------------------------------------------------
# QC metrics
# --------------------------------------------------------------------------

def make_qc(specimen, bad_run):
    if specimen == "blood":
        reads, mapped, dup, cov = random.gauss(8.5e6, 1.2e6), random.uniform(98, 99.6), random.uniform(8, 20), random.gauss(1650, 200)
    elif specimen == "reference_standard":
        reads, mapped, dup, cov = random.gauss(8.0e6, 0.8e6), random.uniform(98.5, 99.6), random.uniform(8, 15), random.gauss(1550, 120)
    else:  # FFPE: degraded input, lower complexity
        reads, mapped, dup, cov = random.gauss(6.0e6, 1.3e6), random.uniform(95.5, 98.5), random.gauss(30, 5.5), random.gauss(1100, 200)
    if bad_run:  # under-clustered flowcell: roughly 45% of normal yield
        reads *= 0.45
        cov *= 0.45
        dup += random.uniform(3, 7)
    cov = max(cov, 50)
    pct_500x = clamp(100 / (1 + math.exp(-(cov - 500) / 130)) + random.gauss(0, 1), 0, 99.9)
    return (int(max(reads, 1e5)), round(mapped, 2), round(clamp(dup, 1, 95), 2),
            round(cov, 1), round(pct_500x, 2))


# --------------------------------------------------------------------------
# variants
# --------------------------------------------------------------------------

def make_genotype(patient):
    """Fixed per-patient somatic and germline variants, reused across serial samples."""
    g = {"drivers": [], "germline": [], "met14": None}
    p72r = random.random()
    if p72r < 0.45:
        g["germline"].append((HOTSPOTS["TP53_P72R"], 0.5))
    elif p72r < 0.55:
        g["germline"].append((HOTSPOTS["TP53_P72R"], 1.0))

    if patient["cohort"] == "solid":
        primary = weighted_choice(SOLID_DRIVERS)
        if primary:
            g["drivers"].append([primary, random.uniform(0.75, 1.0)])
        if random.random() < 0.4:
            g["drivers"].append([random.choice(["TP53_R175H", "TP53_R273H"]), 1.0])
        if random.random() < 0.15:
            g["drivers"].append(["APC_R1450X", random.uniform(0.6, 1.0)])
        genes = PANEL_GENES["ONCO-DNA"]
        if random.random() < 0.15:
            g["met14"] = random.choice(["MET_EX14_D1", "MET_EX14_D2"])
    else:
        pool = list(HEME_DRIVERS)
        for _ in range(random.choice([1, 1, 2, 2, 3])):
            key = weighted_choice(pool)
            pool = [x for x in pool if x[0] != key]
            g["drivers"].append([key, random.uniform(0.2, 0.9)])   # clonal fraction of blasts
        genes = PANEL_GENES["HEME-DNA"]
    for _ in range(random.choice([0, 0, 1, 1, 2])):
        g["germline"].append((random_vus(random.choice([x for x in genes if x in GENE_SPAN])), 0.5))
    return g


def observe(true_vaf, mean_cov):
    """Simulate sequencing a variant: returns (vaf, depth) or None if not called."""
    depth = max(30, int(mean_cov * random.gauss(1.0, 0.2)))
    alt = round(depth * clamp(true_vaf * random.gauss(1.0, 0.1), 0, 1))
    if alt < CALLER_MIN_ALT_READS or alt / depth < CALLER_MIN_VAF:
        return None
    return round(alt / depth, 4), depth


def clinical_calls(sample, genotype, mean_cov, heme_burden):
    """Variant calls for a clinical sample as a list of (hotspot_tuple, vaf, depth)."""
    calls = []
    test = sample["test_code"]
    if test == "FUSION-RNA":
        if genotype["met14"]:
            obs = observe(random.uniform(0.15, 0.6) * sample["tumor_purity"] * 1.5, mean_cov)
            if obs:
                calls.append((HOTSPOTS[genotype["met14"]], *obs))
        return calls

    for var, vaf in genotype["germline"]:
        obs = observe(vaf if vaf < 1 else 0.995, mean_cov)
        if obs:
            calls.append((var, *obs))

    for key, cf in genotype["drivers"]:
        if test == "ONCO-DNA":
            true_vaf = sample["tumor_purity"] * cf / 2
        else:
            true_vaf = clamp(cf / 2 * heme_burden, 0, 0.5)
        obs = observe(true_vaf, mean_cov)
        if obs:
            calls.append((HOTSPOTS[key], *obs))

    if random.random() < 0.2:   # low-level subclone, below the 5% reporting cutoff
        var = HOTSPOTS[random.choice(SUBCLONES)]
        obs = observe(random.uniform(0.02, 0.045), mean_cov)
        if obs and not any(c[0] == var for c in calls):
            calls.append((var, *obs))
    return calls


# --------------------------------------------------------------------------
# main entry point
# --------------------------------------------------------------------------

def flowcell_id(instrument):
    chars = "ABCDEFGHJKLMNPRSTVWXY23456789"
    body = "".join(random.choice(chars) for _ in range(6))
    return f"22{body}LT3" if instrument == "NovaSeq X" else f"AAF{body}M5"


def seed(conn):
    random.seed(SEED)
    cur = conn.cursor()

    cur.executemany("INSERT INTO tests VALUES (?, ?, ?, ?, ?)", TESTS)
    cur.executemany("INSERT INTO panel_genes VALUES (?, ?)",
                    [(t, g) for t, genes in PANEL_GENES.items() for g in genes])
    cur.executemany("INSERT INTO qc_thresholds VALUES (?, ?, ?)", [
        ("total_reads", 2_000_000, None),
        ("pct_mapped", 95.0, None),
        ("pct_duplicate", None, 40.0),
        ("mean_coverage", 500.0, None),
        ("pct_target_500x", 90.0, None),
    ])

    runs = []
    for i in range(N_RUNS):
        instrument = "NextSeq 2000" if i == BAD_RUN else random.choice(["NovaSeq X", "NovaSeq X", "NextSeq 2000"])
        run_date = FIRST_RUN + timedelta(weeks=i)
        runs.append({"run_id": f"R{run_date:%y%m%d}", "instrument": instrument,
                     "run_date": run_date, "flowcell": flowcell_id(instrument)})
    cur.executemany("INSERT INTO runs VALUES (?, ?, ?, ?)",
                    [(r["run_id"], r["instrument"], r["run_date"].isoformat(), r["flowcell"]) for r in runs])

    patients = make_patients()
    cur.executemany("INSERT INTO patients VALUES (?, ?, ?, ?)",
                    [(p["patient_id"], p["mrn"], p["sex"], p["birth_year"]) for p in patients])
    plans, starts = schedule(patients)
    genotypes = {p["patient_id"]: make_genotype(p) for p in patients}

    # Flatten orders into (week, patient_id, order_index, test, specimen).
    events = []
    for p in patients:
        pid = p["patient_id"]
        for idx, (offset, test, specimen) in enumerate(plans[pid]):
            events.append((starts[pid] + offset, pid, idx, test, specimen))
    for rep, w in enumerate(VALIDATION_RUNS):
        events.append((w, None, rep, "ONCO-DNA", "reference_standard"))
    events.sort(key=lambda e: (e[0], e[1] is None, e[1] or 0, e[2]))

    tat_model = {"ONCO-DNA": (8.5, 1.8), "HEME-DNA": (11.5, 2.2), "FUSION-RNA": (10.5, 2.0)}
    block_purity, received_by_week, heme_burden = {}, {}, {}
    sample_rows, qc_rows, variant_rows, truth_rows = [], [], [], []

    for n, (week, pid, idx, test, specimen) in enumerate(events, start=1):
        run = runs[week]
        sample_id = f"NGS26-{n:04d}"
        bad = week == BAD_RUN
        received = received_by_week.setdefault((pid, week), run["run_date"] - timedelta(days=random.randint(2, 6)))
        is_val = pid is None

        if specimen == "FFPE":
            purity = block_purity.setdefault((pid, week), round(random.uniform(0.2, 0.9), 2))
        else:
            purity = None

        if is_val:
            reported = None
        else:
            mean, sd = tat_model[test]
            floor = (run["run_date"] - received).days + 3
            tat = max(floor, round(random.gauss(mean, sd)))
            if bad:
                tat += random.randint(5, 8)       # failed run had to be re-sequenced
            reported = received + timedelta(days=tat)

        sample = {"sample_id": sample_id, "test_code": test, "tumor_purity": purity}
        sample_rows.append((sample_id, pid, run["run_id"], test, specimen, int(is_val), purity,
                            received.isoformat(), reported.isoformat() if reported else None))
        qc = make_qc(specimen, bad)
        qc_rows.append((sample_id, *qc))
        mean_cov = qc[3]

        if is_val:
            for key, exp_vaf in TRUTH_SET:
                gene, chrom, pos, ref, alt, _, vtype, _ = HOTSPOTS[key]
                truth_rows.append((sample_id, chrom, pos, ref, alt, gene, vtype, exp_vaf))
                p_detect = 1.0 if exp_vaf >= 0.05 else 0.8 if exp_vaf >= 0.02 else 0.45
                if random.random() < p_detect:
                    depth = max(200, int(random.gauss(1450, 200)))
                    alt_reads = max(CALLER_MIN_ALT_READS, round(depth * exp_vaf * random.gauss(1.0, 0.12)))
                    variant_rows.append((sample_id, chrom, pos, ref, alt, gene, HOTSPOTS[key][5],
                                         vtype, round(alt_reads / depth, 4), depth, HOTSPOTS[key][7]))
            if idx == FP_REPLICATE:
                key, vaf = FALSE_POSITIVE
                gene, chrom, pos, ref, alt, hgvs, vtype, tier = HOTSPOTS[key]
                variant_rows.append((sample_id, chrom, pos, ref, alt, gene, hgvs, vtype, vaf, 1388, tier))
            continue

        burden = 1.0
        if test == "HEME-DNA":
            if pid in heme_burden:        # follow-up: usually response, sometimes relapse
                factor = random.uniform(0.1, 0.6) if random.random() < 0.75 else random.uniform(1.1, 1.6)
                heme_burden[pid] = clamp(heme_burden[pid] * factor, 0.02, 1.0)
            else:
                heme_burden[pid] = 1.0
            burden = heme_burden[pid]
        elif test == "ONCO-DNA" and idx > 0:
            drivers = genotypes[pid]["drivers"]
            if any(k in ("EGFR_L858R", "EGFR_EX19DEL") for k, _ in drivers) and \
                    not any(k == "EGFR_T790M" for k, _ in drivers):
                drivers.append(["EGFR_T790M", random.uniform(0.3, 0.6)])   # acquired resistance

        for var, vaf, depth in clinical_calls(sample, genotypes[pid], mean_cov, burden):
            gene, chrom, pos, ref, alt, hgvs, vtype, tier = var
            variant_rows.append((sample_id, chrom, pos, ref, alt, gene, hgvs, vtype, vaf, depth, tier))

    cur.executemany("INSERT INTO samples VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", sample_rows)
    cur.executemany("INSERT INTO qc_metrics VALUES (?, ?, ?, ?, ?, ?)", qc_rows)
    cur.executemany("INSERT INTO variants (sample_id, chrom, pos, ref, alt, gene, hgvs_p, var_type, vaf, depth, tier) "
                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", variant_rows)
    cur.executemany("INSERT INTO truth_variants VALUES (?, ?, ?, ?, ?, ?, ?, ?)", truth_rows)
    conn.commit()
    return {"samples": len(sample_rows), "variants": len(variant_rows), "truth": len(truth_rows)}
