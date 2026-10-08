# TransHiC - Inter-chromosomal (trans) Hi-C contact analysis

A small, configurable pipeline that answers three questions from Hi-C `.cool` files:

1. **Which chromosome pairs contact each other more (or less) than their sizes predict?**
2. **Which chromosome pairs gain or lose contacts between conditions** (e.g. knockout vs wild type)?
3. **Which specific bins drive those changes, and what genes are in them?**

It was written for *Leishmania* (34–36 chromosomes, no CTCF, few TADs), where trans contacts mostly reflect
nuclear organisation, but nothing is organism-specific: any genome, any number of samples, any bin size.

Cestari Lab · McGill University · Lissa Cruz-Saavedra
```
.cool files ──► 1 extract ──► 2 enrichment ──► 3 differential ──► 4 candidate bins ──► 5 annotate ──► 6 plots
                chr×chr        obs vs exp        test vs control     bin-pair level        GFF overlap     heatmaps
```

## Installation

```bash
git clone https://github.com/<your-user>/hic-trans-interactions.git
cd hic-trans-interactions
python -m venv .venv && source .venv/bin/activate
pip install -e .            # or: pip install -r requirements.txt
```

Only standard Python packages are needed (numpy, pandas, scipy, statsmodels, h5py, pyyaml, matplotlib).
`cooler` and `bedtools` are **not** required.

On an HPC cluster, see [`slurm/run_transhic.sbatch`](slurm/run_transhic.sbatch).

## Try it in 30 seconds (synthetic data)

```bash
python tests/make_toy_data.py example_data
python -m transhic run --config example_data/config.yaml
```

The toy data has a knockout (`KO`) sequenced 3× deeper than WT with one planted gained interaction
(`toy_chr.02:100–105 kb` ↔ `toy_chr.05:200–205 kb`) and an addback that looks like WT. The pipeline should
report exactly that chromosome pair, rank the planted bin pair first, and show the addback reverses it.

## Running on your data

```bash
cp config/config.example.yaml config/config.yaml   # edit sample paths, comparisons, chromosomes, GFF
python -m transhic run --config config/config.yaml
```

Run individual steps (each step reads the previous step's output from `outdir`):

```bash
python -m transhic run --config config/config.yaml --steps 3-6
python -m transhic differential --config config/config.yaml
```

### Key config options

| Option | Meaning |
|---|---|
| `samples` | name → `.cool`/`.mcool` path; any number of samples |
| `reference` | default control |
| `comparisons` | `DAC3` (= DAC3 vs reference) or `{test: DAC3_ab, control: DAC3}` |
| `chromosomes` | list of names or `{regex: ...}` to exclude maxicircle / contigs |
| `resolution` | bin size to use from an `.mcool` |
| `normalize_by` | `total` (cis+trans, default) or `trans` — see *Statistics* |
| `fdr`, `min_abs_log2fc` | significance thresholds for step 3 |
| `top_bins_per_pair` | bin pairs reported per significant chromosome pair |
| `gff` | annotation for step 5 (`null` to skip) |

## Outputs (`outdir/`)

| File | Content |
|---|---|
| `qc_cis_trans_summary.tsv` | cis, trans totals and trans fraction per sample |
| `trans_contacts_per_pair.tsv` | summed trans contacts for every chromosome pair × sample |
| `trans_enrichment_stats.tsv` | observed, expected, log2(obs/exp), chi-square p, FDR per pair × sample |
| `differential_trans_interactions.tsv` | per comparison: counts, depth-normalised fractions, log2FC, p, FDR, gained/lost |
| `candidate_bin_pairs.tsv` | top bin pairs driving each significant chromosome pair, with bin-level log2FC and p |
| `candidate_anchors.bed` | unique anchor bins (BED, for IGV / bedtools) |
| `candidate_anchors_annotated.tsv` | anchors with overlapping gene IDs and descriptions |
| `candidate_bin_pairs_annotated.tsv` | **final ranked list**: bin pairs with genes on both anchors |
| `plots/enrichment_<sample>.png` | chr × chr heatmap of log2(obs/exp) |
| `plots/differential_<comparison>.png` | chr × chr heatmap of log2FC, dots = significant |

## Statistics

Full details and worked examples: [`docs/METHODS.md`](docs/METHODS.md).

**Step 2 — enrichment.** Under a size-only null, contacts between chromosomes A and B are proportional to
`size_A × size_B`:

```
Expected(A,B) = T × size_A × size_B / Σ_{i<j} size_i × size_j
```

`T` is that sample's own total trans contacts, so sequencing depth cancels and the expected values sum to `T`.
Observed vs expected is tested with a 1-df chi-square, then Benjamini–Hochberg FDR per sample.

**Step 3 — differential.** For each chromosome pair with `a` contacts in the test sample (library size `Na`) and
`b` in the control (`Nb`), the exact test for two Poisson rates with different exposures is used:
conditional on `n = a + b`, under H0 `a ~ Binomial(n, Na / (Na + Nb))`.
`log2FC = log2[((a+0.5)/Na) / ((b+0.5)/Nb)]`. BH-FDR per comparison; a pair is *significant* if
FDR < `fdr` **and** |log2FC| ≥ `min_abs_log2fc`.

`normalize_by: total` uses cis+trans contacts as library size, so a genome-wide gain of trans contacts
(e.g. loss of chromosome territories) appears as gains. `normalize_by: trans` asks a different question —
whether trans contacts are *redistributed* — and by construction gains in some pairs create relative losses in others.

**Step 4 — candidate bins.** Within each significant chromosome pair, every bin pair is tested the same way,
FDR-corrected within that chromosome pair, and the top bin pairs changing in the same direction are reported.

## Important caveats

- **Use raw counts for testing if you can.** Poisson/binomial tests assume integer counts. Balanced
  (KR/ICE) matrices are rescaled floats; they are rounded (`round_counts: true`), which is an approximation.
  Running steps 1–4 on the raw (unbalanced) `.cool` is statistically cleaner — `normalize_by` already
  handles depth differences.
- **Without biological replicates, p-values are optimistic.** Real Hi-C counts are over-dispersed relative
  to Poisson, and with millions of trans reads almost every pair becomes "significant". That is why
  `min_abs_log2fc` matters: rank by effect size and use p-values as a filter, not as proof. With replicates,
  feed `trans_contacts_per_pair.tsv` into DESeq2 / edgeR instead (pairs = features).
- **Chi-square assumes independent observations**, which Hi-C contacts are not; it is a standard approximation.
- **Check a control.** E.g. Cas9-only vs WT should give few or no hits; if it doesn't, the KO calls are suspect.
- **Copy number.** Aneuploidy (common in *Leishmania*) changes contact counts. Compare karyotypes between
  samples before interpreting a chromosome-wide change as reorganisation.

## Testing

```bash
pip install -e ".[test]"
python -m pytest -q
```

## Citation

If you use this code, please cite this repository (and the paper it accompanies, once published).

## License

MIT — see [LICENSE](LICENSE).
