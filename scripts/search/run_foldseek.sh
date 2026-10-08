#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# Baculoviridae Foldseek structural similarity search
#
# Query:
# selected query structures
#
# Reference:
# all Baculoviridae proteins with usable structures
#
# Goal:
# Make output directly comparable with MMseqs2
# ============================================================


# ============================================================
# Check Foldseek
# ============================================================

if ! command -v foldseek >/dev/null 2>&1; then
    echo "ERROR: foldseek command not found."
    echo
    echo "Install Foldseek first, e.g."
    echo
    echo "conda install -c conda-forge -c bioconda foldseek"
    exit 1
fi


# ============================================================
# Paths
# ============================================================

BASE=/processed/Baculoviridae

OUT=/out

INPUT=${BASE}/foldseek/input

QUERY_STRUCTURES=${INPUT}/query_structures

REFERENCE_STRUCTURES=${INPUT}/reference_structures


FOLDSEEK_DIR=${OUT}/foldseek

DBDIR=${FOLDSEEK_DIR}/db

TMP=${FOLDSEEK_DIR}/tmp

RESULT=${FOLDSEEK_DIR}/results


THREADS=32


mkdir -p "${DBDIR}"
mkdir -p "${TMP}"
mkdir -p "${RESULT}"


echo "============================================================"
echo "Baculoviridae Foldseek structural similarity search"
echo "============================================================"

echo
echo "Query structures:"
echo "${QUERY_STRUCTURES}"

echo
echo "Reference structures:"
echo "${REFERENCE_STRUCTURES}"

echo


# ============================================================
# 1. Foldseek version
# ============================================================

echo "[1/6] Foldseek version"

foldseek version \
    | tee "${RESULT}/foldseek_version.txt"


# ============================================================
# 2. Count structures
# ============================================================

N_QUERY=$(
    find "${QUERY_STRUCTURES}" \
    -maxdepth 1 \
    -type l \
    -name "*.pdb" \
    | wc -l
)

N_REF=$(
    find "${REFERENCE_STRUCTURES}" \
    -maxdepth 1 \
    -type l \
    -name "*.pdb" \
    | wc -l
)

echo
echo "Query structures     : ${N_QUERY}"
echo "Reference structures : ${N_REF}"


# ============================================================
# 3. Build query database
# ============================================================

echo
echo "[2/6] Building query Foldseek database"

rm -f "${DBDIR}/queryDB"* 2>/dev/null || true


foldseek createdb \
    "${QUERY_STRUCTURES}" \
    "${DBDIR}/queryDB" \
    --threads "${THREADS}"


# ============================================================
# 4. Build reference database
# ============================================================

echo
echo "[3/6] Building reference Foldseek database"

rm -f "${DBDIR}/referenceDB"* 2>/dev/null || true


foldseek createdb \
    "${REFERENCE_STRUCTURES}" \
    "${DBDIR}/referenceDB" \
    --threads "${THREADS}"


# ============================================================
# 5. Structural search
# ============================================================

echo
echo "[4/6] Running Foldseek search"

rm -f "${DBDIR}/resultDB"* 2>/dev/null || true


foldseek search \
    "${DBDIR}/queryDB" \
    "${DBDIR}/referenceDB" \
    "${DBDIR}/resultDB" \
    "${TMP}" \
    -s 9.5 \
    -e 1e-3 \
    --max-seqs 7000 \
    -a \
    --threads "${THREADS}"


# ============================================================
# 6. Convert results
# ============================================================

echo
echo "[5/6] Exporting Foldseek TSV"

RAW=${RESULT}/foldseek_all_hits.no_header.tsv

ALL=${RESULT}/foldseek_all_hits.tsv


foldseek convertalis \
    "${DBDIR}/queryDB" \
    "${DBDIR}/referenceDB" \
    "${DBDIR}/resultDB" \
    "${RAW}" \
    --format-output \
"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov,alntmscore,qtmscore,ttmscore,lddt,prob"


# ------------------------------------------------------------
# Add header
# ------------------------------------------------------------

echo -e \
"query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqlen\ttlen\tqcov\ttcov\talntmscore\tqtmscore\tttmscore\tlddt\tprob" \
> "${ALL}"

cat "${RAW}" >> "${ALL}"

rm "${RAW}"


# ============================================================
# Normalize query/target identifiers
#
# Handles cases such as:
#
# AAN28026.1_12373.pdb
# AAN28026.1_12373
# AAN28026.1_12373_A
#
# We strip:
#   .pdb
#   optional final chain suffix
#
# IMPORTANT:
# because Viro3D_ID itself contains "_12373",
# only strip "_X" when X is a single chain-like token.
# ============================================================

NORMALIZED=${RESULT}/foldseek_all_hits.normalized.tsv


awk -F'\t' '
BEGIN {
    OFS="\t"
}

NR==1 {
    print
    next
}

{
    q=$1
    t=$2

    sub(/\.pdb$/, "", q)
    sub(/\.pdb$/, "", t)

    # If Foldseek adds a simple chain suffix such as "_A"
    if (q ~ /_[A-Za-z]$/) {
        sub(/_[A-Za-z]$/, "", q)
    }

    if (t ~ /_[A-Za-z]$/) {
        sub(/_[A-Za-z]$/, "", t)
    }

    $1=q
    $2=t

    print
}
' "${ALL}" \
> "${NORMALIZED}"


# ============================================================
# Remove self hits
# ============================================================

echo
echo "[6/6] Removing self hits"

NONSELF=${RESULT}/foldseek_nonself_hits.tsv


awk -F'\t' '
BEGIN {
    OFS="\t"
}

NR==1 {
    print
    next
}

$1 != $2 {
    print
}
' "${NORMALIZED}" \
> "${NONSELF}"


# ============================================================
# Statistics
# ============================================================

TOTAL_HITS=$(
    tail -n +2 "${NORMALIZED}" |
    wc -l
)

SELF_HITS=$(
    awk -F'\t' '
    NR>1 && $1==$2 {
        n++
    }
    END {
        print n+0
    }
    ' "${NORMALIZED}"
)

NONSELF_HITS=$(
    tail -n +2 "${NONSELF}" |
    wc -l
)

QUERIES_WITH_HITS=$(
    awk -F'\t' '
    NR>1 {
        print $1
    }
    ' "${NONSELF}" |
    sort -u |
    wc -l
)


# ============================================================
# Save parameters
# ============================================================

cat > "${RESULT}/foldseek_parameters.txt" << EOF
Dataset: Baculoviridae

Query structures:
${N_QUERY}

Reference structures:
${N_REF}

Search program:
Foldseek

Search parameters:
sensitivity = 9.5
E-value threshold = 1e-3
max-seqs = 7000
alignment backtrace = enabled (-a)
threads = ${THREADS}

Search-stage filters:
minimum structural coverage = none
minimum sequence identity = none

Self hits:
removed only after search

Query structures:
${QUERY_STRUCTURES}

Reference structures:
${REFERENCE_STRUCTURES}
EOF


# ============================================================
# Summary
# ============================================================

echo
echo "============================================================"
echo "Foldseek search completed"
echo "============================================================"

echo "Query structures          : ${N_QUERY}"
echo "Reference structures      : ${N_REF}"

echo "Total hits                : ${TOTAL_HITS}"
echo "Self hits                 : ${SELF_HITS}"
echo "Non-self hits             : ${NONSELF_HITS}"
echo "Queries with nonself hits : ${QUERIES_WITH_HITS}"

echo
echo "Main output:"
echo "${NONSELF}"

echo
echo "All normalized hits:"
echo "${NORMALIZED}"

echo
echo "Parameters:"
echo "${RESULT}/foldseek_parameters.txt"

echo "============================================================"