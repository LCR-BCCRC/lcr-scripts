#!/usr/bin/env python3
"""
Annotate acquired N-linked glycosylation sites (NxS/T, x != Pro) in IG/TCR
sequences with IMGT unique numbering (Lefranc et al.).

Reads sequence_alignment_aa and germline_alignment_aa from an AIRR-format TSV
(--source_tsv), keyed by sequence_id. Compatible with output from IgBLASTn
(igblast module) and IMGT V-QUEST (vquest module).

IMGT numbering is assigned by ANARCI (bioconda: conda install -c bioconda anarci).

Output TSV columns:
  sequence_id                       FASTA header ID
  aa_sequence                       Query amino acid sequence used for numbering
  num_glycosylation_sites           Count of NxS/T motifs in query (x != Pro)
  glycosylation_imgt_positions      IMGT positions of N residues in query, comma-separated
  glycosylation_motifs              Corresponding NxS/T triplets in query, comma-separated
  glycosylation_imgt_regions        IMGT FWR/CDR region of each query site, comma-separated
  num_acquired_glycosylation_sites  Sites in query absent from germline (SHM-acquired)
  acquired_glycosylation_imgt_positions  IMGT positions of acquired sites, comma-separated
  acquired_glycosylation_imgt_regions    IMGT FWR/CDR region of each acquired site, comma-separated
  num_acquired_glycosylation_sites_cdr   Acquired sites falling in CDR1/CDR2/CDR3
  manntype_ags                      POS/NEG: >=1 acquired site in a CDR (Tatterton et al. 2025
                                     AGS criterion). This is only the sequence-derived half of
                                     manntype classification; combine with a genomic FL signature
                                     (EZB LymphGen subtype or BCL2 translocation) downstream, as
                                     GAMBLR.results::collate_tatterton does for the published
                                     Tatterton cohort, to get the full manntype call.
  germline_aa_sequence              Germline amino acid sequence used for numbering
  num_germline_glycosylation_sites  Count of NxS/T motifs in germline
  germline_glycosylation_imgt_positions  IMGT positions of N residues in germline
  germline_glycosylation_motifs     Corresponding NxS/T triplets in germline
"""

import argparse
import csv
import sys

try:
    from anarci import anarci
except ImportError:
    sys.exit("ERROR: anarci is not installed. Install via: conda install -c bioconda anarci")


# ── TSV loading ───────────────────────────────────────────────────────────────

def load_sequences_by_id(tsv_path):
    """
    Build {sequence_id: {"seq": str, "germline": str}} from an AIRR TSV.
    Both sequences have alignment gap characters stripped.
    """
    data = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq_id = row["sequence_id"]
            aa_seq  = row.get("sequence_alignment_aa",  "").strip().replace("-", "").replace(".", "")
            gl_seq  = row.get("germline_alignment_aa",  "").strip().replace("-", "").replace(".", "")
            if aa_seq:
                data[seq_id] = {"seq": aa_seq, "germline": gl_seq}
    return data


# ── IMGT numbering ────────────────────────────────────────────────────────────

def _imgt_pos_str(pos_tuple):
    """Format ANARCI position (number, insertion_letter) as '27' or '111a'."""
    num, ins = pos_tuple
    return f"{num}{ins.strip()}" if ins.strip() else str(num)


# IMGT unique numbering V-domain region boundaries (Lefranc et al.), inclusive.
IMGT_REGION_BOUNDARIES = [
    (1,   26,  "FWR1"),
    (27,  38,  "CDR1"),
    (39,  55,  "FWR2"),
    (56,  65,  "CDR2"),
    (66,  104, "FWR3"),
    (105, 117, "CDR3"),
    (118, 128, "FWR4"),
]


def _imgt_region(pos_tuple):
    """Map an ANARCI IMGT position (number, insertion_letter) to its FWR/CDR region."""
    num, _ = pos_tuple
    for lo, hi, name in IMGT_REGION_BOUNDARIES:
        if lo <= num <= hi:
            return name
    return "NA"


def number_with_imgt(aa_seq):
    """
    Align aa_seq to IG/TCR germline HMMs with ANARCI using the IMGT scheme.
    Returns [(pos_tuple, aa), ...] with gap positions removed, or None on failure.
    """
    results, _, _ = anarci([("seq", aa_seq)], scheme="imgt", output=False)
    if results[0] is None:
        return None
    return [(pos, aa) for pos, aa in results[0][0][0] if aa != "-"]


def find_glycosylation_sites(numbered):
    """
    Scan IMGT-numbered residues for N-linked glycosylation motifs (NxS/T, x != Pro).

    Returns [(imgt_position_string, motif_string, region_string), ...] for each site.
    """
    sites = []
    for i in range(len(numbered) - 2):
        n_pos, n_aa  = numbered[i]
        _,     x_aa  = numbered[i + 1]
        _,     st_aa = numbered[i + 2]
        if n_aa == "N" and x_aa != "P" and st_aa in ("S", "T"):
            sites.append((_imgt_pos_str(n_pos), f"N{x_aa}{st_aa}", _imgt_region(n_pos)))
    return sites


def classify_sites(query_sites, germline_sites):
    """
    Return the subset of query_sites whose IMGT position is absent in germline_sites.
    These are SHM-acquired glycosylation sites.
    """
    germline_positions = {pos for pos, _, _ in germline_sites}
    return [site for site in query_sites if site[0] not in germline_positions]


def sites_in_cdr(sites):
    """Return the subset of sites (as returned by find_glycosylation_sites) in a CDR."""
    return [site for site in sites if site[2].startswith("CDR")]


def _make_lookup_fn(seq_lookup):
    """
    Return a lookup callable for seq_lookup that handles V-QUEST's fixed-length ID
    truncation. V-QUEST truncates input FASTA headers to 49 characters in AIRR output;
    if all TSV keys are the same length, fall back to truncating the query ID on a miss.
    """
    if not seq_lookup:
        return seq_lookup.get
    lengths = {len(k) for k in seq_lookup}
    if len(lengths) == 1:
        trunc = next(iter(lengths))
        def _lookup(seq_id, _d=seq_lookup, _t=trunc):
            return _d.get(seq_id) or _d.get(seq_id[:_t])
        return _lookup
    return seq_lookup.get


# ── I/O ───────────────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs = {}
    cur_id, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(parts)
                cur_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(parts)
    return seqs


FIELDS = [
    "sequence_id",
    "aa_sequence",
    "num_glycosylation_sites",
    "glycosylation_imgt_positions",
    "glycosylation_motifs",
    "glycosylation_imgt_regions",
    "num_acquired_glycosylation_sites",
    "acquired_glycosylation_imgt_positions",
    "acquired_glycosylation_imgt_regions",
    "num_acquired_glycosylation_sites_cdr",
    "manntype_ags",
    "germline_aa_sequence",
    "num_germline_glycosylation_sites",
    "germline_glycosylation_imgt_positions",
    "germline_glycosylation_motifs",
]


def _sites_str(sites):
    return ",".join(s[0] for s in sites) if sites else "NA"

def _motifs_str(sites):
    return ",".join(s[1] for s in sites) if sites else "NA"

def _regions_str(sites):
    return ",".join(s[2] for s in sites) if sites else "NA"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",      required=True,
                        help="Input FASTA. Used as the source of sequence IDs.")
    parser.add_argument("--source_tsv", required=True,
                        help="AIRR TSV with sequence_alignment_aa and "
                             "germline_alignment_aa columns (IgBLASTn or V-QUEST output).")
    parser.add_argument("--output",     required=True, help="Output TSV path")
    args = parser.parse_args()

    fasta_seqs  = read_fasta(args.fasta)
    seq_lookup  = load_sequences_by_id(args.source_tsv)
    lookup      = _make_lookup_fn(seq_lookup)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        for seq_id in fasta_seqs:
            entry = lookup(seq_id)

            if not entry:
                writer.writerow({
                    "sequence_id": seq_id, "aa_sequence": "NA",
                    "num_glycosylation_sites": 0,
                    "glycosylation_imgt_positions": "NA",
                    "glycosylation_motifs": "NA",
                    "glycosylation_imgt_regions": "NA",
                    "num_acquired_glycosylation_sites": 0,
                    "acquired_glycosylation_imgt_positions": "NA",
                    "acquired_glycosylation_imgt_regions": "NA",
                    "num_acquired_glycosylation_sites_cdr": 0,
                    "manntype_ags": "NA",
                    "germline_aa_sequence": "NA",
                    "num_germline_glycosylation_sites": 0,
                    "germline_glycosylation_imgt_positions": "NA",
                    "germline_glycosylation_motifs": "NA",
                })
                continue

            aa_seq = entry["seq"]
            gl_seq = entry["germline"]

            numbered    = number_with_imgt(aa_seq)
            gl_numbered = number_with_imgt(gl_seq) if gl_seq else None

            query_sites    = find_glycosylation_sites(numbered)    if numbered    else []
            germline_sites = find_glycosylation_sites(gl_numbered) if gl_numbered else []
            acquired_sites = classify_sites(query_sites, germline_sites)
            acquired_cdr_sites = sites_in_cdr(acquired_sites)

            writer.writerow({
                "sequence_id": seq_id,
                "aa_sequence": aa_seq,
                "num_glycosylation_sites": len(query_sites),
                "glycosylation_imgt_positions": _sites_str(query_sites),
                "glycosylation_motifs": _motifs_str(query_sites),
                "glycosylation_imgt_regions": _regions_str(query_sites),
                "num_acquired_glycosylation_sites": len(acquired_sites),
                "acquired_glycosylation_imgt_positions": _sites_str(acquired_sites),
                "acquired_glycosylation_imgt_regions": _regions_str(acquired_sites),
                "num_acquired_glycosylation_sites_cdr": len(acquired_cdr_sites),
                "manntype_ags": "POS" if acquired_cdr_sites else "NEG",
                "germline_aa_sequence": gl_seq if gl_seq else "NA",
                "num_germline_glycosylation_sites": len(germline_sites),
                "germline_glycosylation_imgt_positions": _sites_str(germline_sites),
                "germline_glycosylation_motifs": _motifs_str(germline_sites),
            })


if __name__ == "__main__":
    main()
