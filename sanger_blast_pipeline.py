#!/usr/bin/env python
import os
import re
import subprocess
import yaml
import pandas as pd
from Bio import SeqIO, SearchIO
from Bio.Seq import Seq

# ────────── Load Config ──────────
import argparse
p = argparse.ArgumentParser()
p.add_argument("-c", "--config", default="config.yaml",
               help="Path to YAML config file")
args = p.parse_args()
with open(args.config) as cf:
    cfg = yaml.safe_load(cf)

INPUT_DIR   = cfg["input_dir"]
OUTPUT_DIR  = cfg["output_dir"]
FIGS_DIR    = os.path.join(OUTPUT_DIR, "figs")
TMP_DIR     = os.path.join(OUTPUT_DIR, "tmp")

# Reference & BLAST DB
REF_CFG     = cfg["reference"]
BLAST_CFG   = cfg["blast"]
MUT_CFG     = cfg["mutations"]


# ────────── Step 0: Dependencies ──────────
# pip install biopython pyyaml pandas

def make_blast_db():
    """ Optionally build a local BLAST DB from the reference FASTA. """
    if not REF_CFG.get("make_blast_db", False):
        return os.path.join(REF_CFG["blast_db_dir"], REF_CFG["blast_db_name"])
    os.makedirs(REF_CFG["blast_db_dir"], exist_ok=True)
    prefix = os.path.join(REF_CFG["blast_db_dir"], REF_CFG["blast_db_name"])
    cmd = [
        "makeblastdb",
        "-in", REF_CFG["fasta"],
        "-dbtype", "nucl",
        "-out", prefix
    ]
    print("Building BLAST DB:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return prefix

BLAST_DB = make_blast_db()


# ────────── Step 1: Pair .ab1 reads ──────────
def parse_sample_key(fn):
    base = os.path.splitext(os.path.basename(fn))[0]
    m = re.search(r"(A\d+.*)", base, flags=re.IGNORECASE)
    if not m: return base, None
    sub = m.group(1)
    d = None
    if sub.upper().endswith(("-F","_F")):
        d, key = "F", sub[:-2]
    elif sub.upper().endswith(("-R","_R")):
        d, key = "R", sub[:-2]
    else:
        key = sub
    return key, d

def gather_pairs():
    pairs = {}
    for root, _, files in os.walk(INPUT_DIR):
        for f in files:
            if not f.lower().endswith(".ab1"): continue
            key, dirn = parse_sample_key(f)
            if dirn not in ("F","R"): continue
            pairs.setdefault(key, {})[dirn] = os.path.join(root, f)
    return {k:v for k,v in pairs.items() if "F" in v and "R" in v}


# ────────── Step 2: Build forward/reverse consensus ──────────
def load_abi_sequence(path):
    rec = SeqIO.read(path, "abi")
    return str(rec.seq).upper(), rec

def merge_with_shift(fwd, rev, shift):
    # see earlier example; for brevity we assume same code as above
    if shift >= 0:
        total = max(len(fwd), shift+len(rev))
        merged = ["N"]*total
        for i,b in enumerate(fwd): merged[i]=b
        for i,b in enumerate(rev):
            j = i+shift
            if merged[j]=="N": merged[j]=b
        ov_start, ov_end = shift, min(len(fwd),shift+len(rev))
        score = sum(1 for i in range(ov_start,ov_end)
                    if merged[i]==fwd[i]==rev[i-shift])
    else:
        total = max(len(rev), -shift+len(fwd))
        merged = ["N"]*total
        for i,b in enumerate(rev): merged[i]=b
        for i,b in enumerate(fwd):
            j = i-shift
            if merged[j]=="N": merged[j]=b
        ov_start, ov_end = -shift, min(len(rev),-shift+len(fwd))
        score = sum(1 for i in range(ov_start,ov_end)
                    if merged[i]==rev[i]==fwd[i+shift])
    return "".join(merged).strip("N"), score

def build_consensus(fwd_seq, rev_seq):
    revcomp = str(Seq(rev_seq).reverse_complement())
    best, best_score = None, -1
    for shift in range(-len(revcomp), len(fwd_seq)):
        m,s = merge_with_shift(fwd_seq, revcomp, shift)
        if s>best_score:
            best, best_score = m,s
    return best


# ────────── Step 3: BLAST the consensus ──────────
def run_blast(consensus, sample_key):
    os.makedirs(TMP_DIR, exist_ok=True)
    qfa  = os.path.join(TMP_DIR, f"{sample_key}.fa")
    qxml = os.path.join(TMP_DIR, f"{sample_key}.xml")
    with open(qfa,"w") as f:
        f.write(">qry\n"+consensus+"\n")
    cmd = [
        "blastn",
        "-task", BLAST_CFG["task"],
        "-query", qfa,
        "-db",    BLAST_DB,
        "-outfmt", BLAST_CFG["outfmt"],
        "-evalue", str(BLAST_CFG["evalue"]),
        "-out", qxml
    ]
    if cfg.get("debug", False):
        print("Running BLAST:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return qfa, qxml


# ────────── Step 4: Generic mutation check ──────────
def check_mutation(hsp, name, info):
    typ = info["type"]
    if typ=="SNV":
        p=info["pos"]
        if not (hsp.hit_start<=p<=hsp.hit_end): return False
        off, qi = p-hsp.hit_start, hsp.query_start + (p-hsp.hit_start)
        if qi<0 or qi>=len(hsp.query.seq): return False
        return hsp.query.seq[qi].upper()==info["mut"].upper()

    if typ=="deletion":
        a,b=info["start"],info["end"]
        if not(hsp.hit_start<=a and hsp.hit_end>=b): return False
        ref_len = b-a+1
        qs = hsp.query_start + (a-hsp.hit_start)
        qe = hsp.query_start + (b-hsp.hit_start)
        qr = hsp.query.seq[qs:qe+1]
        return len(qr) < ref_len * 0.5

    if typ=="duplication":
        a,b=info["start"],info["end"]
        if not(hsp.hit_start<=a and hsp.hit_end>=b): return False
        ref_len=b-a+1
        qs=hsp.query_start+(a-hsp.hit_start)
        qe=hsp.query_start+(b-hsp.hit_start)
        qr=hsp.query.seq[qs:qe+1]
        return len(qr) > ref_len + 2

    if typ=="compound":
        for sub in info["submutations"]:
            if not check_mutation(hsp, name, sub):
                return False
        return True

    return False


def parse_blast_for_all(qxml):
    qr = SearchIO.read(qxml, "blast-xml")
    if not qr.hits:
        return {k: False for k in MUT_CFG}
    hsp = qr.hits[0].hsps[0]
    out = {}
    for name,info in MUT_CFG.items():
        out[name] = check_mutation(hsp, name, info)
    return out


# ────────── Step 5: Worker ──────────
def process_sample(key, paths):
    seq_f, rec_f = load_abi_sequence(paths["F"])
    seq_r, rec_r = load_abi_sequence(paths["R"])
    cons = build_consensus(seq_f, seq_r)
    qfa, qxml = run_blast(cons, key)
    return {"Sample": key, **parse_blast_for_all(qxml)}


# ────────── Step 6: Main ──────────
if __name__=="__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(FIGS_DIR, exist_ok=True)
    pairs = gather_pairs()
    print(f"Found {len(pairs)} paired samples")

    results = []
    for key, paths in pairs.items():
        print(">> Processing", key)
        results.append(process_sample(key, paths))

    df = pd.DataFrame(results)
    out_csv = os.path.join(OUTPUT_DIR, "egfr_mutation_results.csv")
    df.to_csv(out_csv, index=False)
    print("✔️ Results written to", out_csv)