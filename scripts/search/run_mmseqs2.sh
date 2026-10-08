#运行相似度的计算mmseq2.sh


#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# Baculoviridae MMseqs2 sequence similarity search
#
# Query:
# selected_queries.tsv 中选择的代表性蛋白
#
# Reference:
# 全部 Baculoviridae proteins
#
# Important:
# No identity cutoff
# No coverage cutoff at search stage
# Keep weak sequence hits for later Foldseek comparison
# ============================================================


BASE=processed/Baculoviridae


# ============================================================
# Input
# ============================================================

REF_FASTA=${BASE}/mmseqs2/input/Baculoviridae_reference.fasta

QUERY_FASTA=${BASE}/mmseqs2/input/Baculoviridae_queries.fasta


# ============================================================
# Output
# ============================================================

out=/out

MMSEQS_DIR=${out}/mmseqs2

DBDIR=${MMSEQS_DIR}/db

TMP=${MMSEQS_DIR}/tmp

RESULT=${MMSEQS_DIR}/results


THREADS=32


mkdir -p "${DBDIR}"
mkdir -p "${TMP}"
mkdir -p "${RESULT}"


echo "============================================================"
echo "Baculoviridae MMseqs2 search"
echo "============================================================"

echo
echo "Reference:"
echo "${REF_FASTA}"

echo
echo "Queries:"
echo "${QUERY_FASTA}"

echo


# ============================================================
# 1. MMseqs2 version
# ============================================================

echo "[1/6] MMseqs2 version"

mmseqs version \
    | tee "${RESULT}/mmseqs_version.txt"


# ============================================================
# 2. Count FASTA sequences
# ============================================================

N_REF=$(grep -c "^>" "${REF_FASTA}")

N_QUERY=$(grep -c "^>" "${QUERY_FASTA}")

echo
echo "Reference proteins: ${N_REF}"
echo "Query proteins:     ${N_QUERY}"


# ============================================================
# 3. Create databases
# ============================================================

echo
echo "[2/6] Building query database"

rm -f \
    "${DBDIR}/queryDB"* \
    2>/dev/null || true

mmseqs createdb \
    "${QUERY_FASTA}" \
    "${DBDIR}/queryDB"


echo
echo "[3/6] Building reference database"

rm -f \
    "${DBDIR}/referenceDB"* \
    2>/dev/null || true

mmseqs createdb \
    "${REF_FASTA}" \
    "${DBDIR}/referenceDB"


# ============================================================
# 4. Search
# ============================================================

echo
echo "[4/6] Running MMseqs2 search"

rm -f \
    "${DBDIR}/resultDB"* \
    2>/dev/null || true

mmseqs search \
    "${DBDIR}/queryDB" \
    "${DBDIR}/referenceDB" \
    "${DBDIR}/resultDB" \
    "${TMP}" \
    -s 8.5 \
    -e 1e-3 \
    --max-seqs 7000 \
    -a \
    --threads "${THREADS}"


# ============================================================
# 5. Convert alignment results
# ============================================================

echo
echo "[5/6] Exporting TSV"

RAW_NO_HEADER=${RESULT}/mmseqs_all_hits.no_header.tsv

FINAL_ALL=${RESULT}/mmseqs_all_hits.tsv


mmseqs convertalis \
    "${DBDIR}/queryDB" \
    "${DBDIR}/referenceDB" \
    "${DBDIR}/resultDB" \
    "${RAW_NO_HEADER}" \
    --format-output \
"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov"


# Add our own header explicitly
echo -e \
"query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbits\tqlen\ttlen\tqcov\ttcov" \
> "${FINAL_ALL}"


cat "${RAW_NO_HEADER}" \
    >> "${FINAL_ALL}"


rm "${RAW_NO_HEADER}"


# ============================================================
# 6. Remove self hits
# ============================================================

echo
echo "[6/6] Removing self hits"

NONSELF=${RESULT}/mmseqs_nonself_hits.tsv


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
' \
"${FINAL_ALL}" \
> "${NONSELF}"


# ============================================================
# Statistics
# ============================================================

TOTAL_HITS=$(
    tail -n +2 "${FINAL_ALL}" |
    wc -l
)

NONSELF_HITS=$(
    tail -n +2 "${NONSELF}" |
    wc -l
)

SELF_HITS=$(
    awk -F'\t' '
    NR>1 && $1==$2 {n++}
    END {print n+0}
    ' "${FINAL_ALL}"
)

QUERIES_WITH_NONSELF_HITS=$(
    awk -F'\t' '
    NR>1 {print $1}
    ' "${NONSELF}" |
    sort -u |
    wc -l
)


# ============================================================
# Save run parameters
# ============================================================

cat > "${RESULT}/mmseqs_parameters.txt" << EOF
Dataset: Baculoviridae
Reference proteins: ${N_REF}
Query proteins: ${N_QUERY}

Search program:
MMseqs2

Search parameters:
sensitivity = 7.0
E-value threshold = 1e-3
max-seqs = 1000
alignment backtrace = enabled (-a)
threads = ${THREADS}

Search-stage filters:
minimum sequence identity = none
minimum query coverage = none
minimum target coverage = none

Self hits:
removed only after search

Reference FASTA:
${REF_FASTA}

Query FASTA:
${QUERY_FASTA}
EOF


echo
echo "============================================================"
echo "MMseqs2 search completed"
echo "============================================================"

echo "Reference proteins       : ${N_REF}"

echo "Query proteins           : ${N_QUERY}"

echo "Total hits               : ${TOTAL_HITS}"

echo "Self hits                : ${SELF_HITS}"

echo "Non-self hits            : ${NONSELF_HITS}"

echo "Queries with nonself hits: ${QUERIES_WITH_NONSELF_HITS}"

echo
echo "Main output:"
echo "${NONSELF}"

echo
echo "All hits:"
echo "${FINAL_ALL}"

echo
echo "Parameters:"
echo "${RESULT}/mmseqs_parameters.txt"

echo "============================================================"