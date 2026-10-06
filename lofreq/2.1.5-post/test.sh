#!/bin/bash
# Functional test for lofreq:2.1.5-post container.
# Run via:
#   docker run --rm \
#     -v /path/to/output:/data \
#     -v /path/to/this/dir:/scripts \
#     <image> bash /scripts/test.sh
#
# Output files are written to /data inside the container (host mount).

set -x
cd /data

# ── Reference ────────────────────────────────────────────────────────────────
python3 -c "
with open('ref.fa', 'w') as f:
    f.write('>chr1\n' + 'ACGT' * 50 + '\n')
    f.write('>chr2\n' + 'GCTA' * 50 + '\n')
" || { echo "FAIL: ref creation"; exit 1; }
samtools faidx ref.fa || { echo "FAIL: faidx"; exit 1; }

# ── Synthetic reads (50 per chromosome, 15/50 = 30% AF SNV) ─────────────────
python3 -c "
hdr = ['@HD\tVN:1.4\tSO:coordinate',
       '@SQ\tSN:chr1\tLN:200',
       '@SQ\tSN:chr2\tLN:200']
rows = []
for chrom, ref_unit, alt in [('chr1','ACGT','T'), ('chr2','GCTA','A')]:
    ref_seq = ref_unit * 50
    for i in range(50):
        seq = list(ref_seq[100:150])
        if i < 15:
            seq[24] = alt
        rows.append('\t'.join([
            'r{}_{}'.format(chrom, i), '0', chrom, '101', '40',
            '50M', '*', '0', '0', ''.join(seq), 'I' * 50
        ]))
print('\n'.join(hdr + rows))
" > reads.sam || { echo "FAIL: SAM creation"; exit 1; }

samtools sort -O BAM -o reads.bam reads.sam  || { echo "FAIL: samtools sort"; exit 1; }
samtools index reads.bam                     || { echo "FAIL: samtools index"; exit 1; }
lofreq indelqual --dindel -f ref.fa -o reads_iq.bam reads.bam || { echo "FAIL: indelqual"; exit 1; }
samtools index reads_iq.bam                  || { echo "FAIL: indelqual index"; exit 1; }

# ── Single-threaded call ─────────────────────────────────────────────────────
lofreq call -f ref.fa -o single.vcf reads_iq.bam || { echo "FAIL: lofreq call"; exit 1; }
n_single=$(grep -v '^#' single.vcf | wc -l)
echo "Single call: $n_single variants"

# ── Parallel call (exercises bcftools concat path from PR #109) ──────────────
lofreq call-parallel --pp-threads 2 -f ref.fa -o parallel.vcf.gz reads_iq.bam || { echo "FAIL: lofreq call-parallel"; exit 1; }
n_parallel=$(zgrep -v '^#' parallel.vcf.gz | wc -l)
echo "Parallel call: $n_parallel variants"

# ── Verify results agree ─────────────────────────────────────────────────────
if [ "$n_single" = "$n_parallel" ]; then
    echo "PASS: single and parallel agree ($n_single variants)"
else
    echo "FAIL: single=$n_single parallel=$n_parallel"
    exit 1
fi
