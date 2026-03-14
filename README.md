# Conversion Notes: mad4hatter outputs to v1.0.0

This document describes how each file and column is handled when converting
mad4hatter results from older releases to v1.0.0 format.

Supported input versions: **v0.1.8, v0.1.9, v0.2.0, v0.2.1, v0.2.2**.
Releases prior to v0.1.8 are not supported.

Version is auto-detected from column headers in `allele_data.txt`. If detection
fails or is incorrect, it can be overridden by setting `release_version` manually
in `run_convert_mad4hatter.R`.

---

## Output structure

The converter produces a v1.0.0-compatible directory with the following structure:

```
output_dir/
├── allele_data.txt
├── allele_data_collapsed.txt
├── sample_coverage.txt
├── amplicon_coverage.txt
├── raw_dada2_output/
│   └── dada2.clusters.txt
├── resistance_marker_module/
│   ├── resmarker_table.txt
│   ├── resmarker_table_by_locus.txt        # v0.2.2 only
│   └── resmarker_microhaplotype_table.txt
├── panel_information/                       # generated from reference inputs
│   ├── amplicon_info.tsv
│   ├── reference.fasta
│   └── resmarker_info.tsv
├── run/                                     # copied as-is if present
└── legacy/
    └── quality_report/                      # copied as-is
```

---

## File-by-file mapping

### `allele_data.txt`

| v1.0.0 column        | v0.1.8                    | v0.1.9 / v0.2.0 / v0.2.1 / v0.2.2 | Notes                                         |
|----------------------|---------------------------|-------------------------------------|-----------------------------------------------|
| `sample_name`        | `sampleID`                | `SampleID`                          |                                               |
| `target_name`        | `locus`                   | `Locus`                             | Mapped via `target_id_conversion_table.tsv`   |
| `asv`                | `asv`                     | `ASV`                               |                                               |
| `pseudocigar_unmasked` | —                       | —                                   | **Set to NA** — requires realignment to reference; recompute with mad4hatter v1.0.0 |
| `asv_masked`         | —                         | —                                   | **Set to NA** — not available in these versions |
| `pseudocigar_masked` | `pseudo_cigar`            | `PseudoCIGAR`                       |                                               |
| `reads`              | `reads`                   | `Reads`                             |                                               |
| `pool`               | `locus`                   | `Locus`                             | Mapped via `target_id_conversion_table.tsv`   |
| `allele` / `Allele`  | dropped                   | dropped                             | Not present in v1.0.0                         |

### `allele_data_collapsed.txt`

Derived from the converted `allele_data.txt` by grouping on
`sample_name`, `target_name`, `asv_masked`, `pseudocigar_masked`, `pool`
and summing `reads`. Because `asv_masked` is NA for all pre-v1.0.0 outputs,
this file will contain NA values in those columns — this is expected.

### `sample_coverage.txt`

| v1.0.0 column | All old versions  |
|---------------|-------------------|
| `sample_name` | `SampleID`        |
| `stage`       | `Stage`           |
| `reads`       | `Reads`           |

### `amplicon_coverage.txt`

| v1.0.0 column          | All old versions       | Notes                                       |
|------------------------|------------------------|---------------------------------------------|
| `sample_name`          | `SampleID`             |                                             |
| `target_name`          | `Locus`                | Mapped via `target_id_conversion_table.tsv` |
| `reads`                | `Reads`                |                                             |
| `OutputDada2`          | `OutputDada2`          | Unchanged                                   |
| `OutputPostprocessing` | `OutputPostprocessing` | Unchanged                                   |

---

## `resistance_marker_module/`

### `resmarker_table.txt`

| v1.0.0 column  | v0.1.8 / v0.1.9       | v0.2.0 / v0.2.1        | v0.2.2          | Notes                                                                 |
|----------------|-----------------------|------------------------|-----------------|-----------------------------------------------------------------------|
| `sample_name`  | `SampleID`            | `SampleID`             | `SampleID`      |                                                                       |
| `gene_id`      | `GeneID`              | `GeneID`               | `GeneID`        |                                                                       |
| `gene`         | `Gene`                | `Gene`                 | `Gene`          |                                                                       |
| `aa_position`  | `CodonID`             | `CodonID`              | `CodonID`       |                                                                       |
| `ref_codon`    | `RefCodon`            | `RefCodon`             | `RefCodon`      |                                                                       |
| `codon`        | `Codon`               | `Codon`                | `Codon`         |                                                                       |
| `codon_ref_alt`| `CodonRefAlt`         | `CodonRefAlt`          | `CodonRefAlt`   |                                                                       |
| `ref_aa`       | `RefAA`               | `RefAA`                | `RefAA`         |                                                                       |
| `aa`           | `AA`                  | `AA`                   | `AA`            |                                                                       |
| `aa_ref_alt`   | `AARefAlt`            | `AARefAlt`             | `AARefAlt`      |                                                                       |
| `follows_indel`| —                     | —                      | `FollowsIndel`  | **Set to NA** for v0.1.8–v0.2.1                                       |
| `codon_masked` | —                     | —                      | `CodonMasked`   | **Set to NA** for v0.1.8–v0.2.1                                       |
| `multiple_loci`| —                     | —                      | `MultipleLoci`  | **Set to NA** for v0.1.8–v0.2.1                                       |
| `reads`        | `Reads` (summed)      | `Reads` (summed)       | `Reads`         | v0.1.8–v0.2.1 had a `CodonStart` column that created duplicate rows; reads are summed across `CodonStart` values |
| `Locus`        | —                     | present, dropped       | —               | Locus was present in v0.2.0/v0.2.1 resmarker_table but is not in v1.0.0 |

### `resmarker_table_by_locus.txt`

Only present in v0.2.2 input. Not produced for v0.1.8–v0.2.1.

| v1.0.0 column  | v0.2.2        | Notes                                       |
|----------------|---------------|---------------------------------------------|
| `sample_name`  | `SampleID`    |                                             |
| `gene_id`      | `GeneID`      |                                             |
| `gene`         | `Gene`        |                                             |
| `target_name`  | `Locus`       | Mapped via `target_id_conversion_table.tsv` |
| `aa_position`  | `CodonID`     |                                             |
| `ref_codon`    | `RefCodon`    |                                             |
| `codon`        | `Codon`       |                                             |
| `codon_ref_alt`| `CodonRefAlt` |                                             |
| `ref_aa`       | `RefAA`       | Note: source had a typo (`RefAAAA`) in some v0.2.2 outputs — both are handled |
| `aa`           | `AA`          |                                             |
| `aa_ref_alt`   | `AARefAlt`    |                                             |
| `follows_indel`| `FollowsIndel`|                                             |
| `codon_masked` | `CodonMasked` |                                             |
| `reads`        | `Reads`       |                                             |

### `resmarker_microhaplotype_table.txt`

| v1.0.0 column      | v0.1.8 / v0.1.9          | v0.2.0 / v0.2.1          | v0.2.2                    | Notes                                                              |
|--------------------|--------------------------|--------------------------|---------------------------|--------------------------------------------------------------------|
| `sample_name`      | `SampleID`               | `SampleID`               | `SampleID`                |                                                                    |
| `gene_id`          | `GeneID`                 | `GeneID`                 | `GeneID`                  |                                                                    |
| `gene`             | `Gene`                   | `Gene`                   | `Gene`                    |                                                                    |
| `target_name`      | —                        | `Locus`                  | `Locus`                   | **Set to NA** for v0.1.8/v0.1.9 — no Locus column in those versions |
| `mhap_aa_positions`| `MicrohapIndex`          | `MicrohapIndex`          | `MicrohaplotypeCodonIDs`  |                                                                    |
| `ref_mhap`         | `RefMicrohap`            | `RefMicrohap`            | `RefMicrohap`             |                                                                    |
| `mhap`             | `Microhaplotype`         | `Microhaplotype`         | `Microhaplotype`          |                                                                    |
| `mhap_ref_alt`     | `MicrohapRefAlt`         | `MicrohapRefAlt`         | `MicrohapRefAlt`          |                                                                    |
| `reads`            | `Reads`                  | `Reads`                  | `Reads`                   |                                                                    |

### Excluded files

The following files are intentionally **not carried over** to the v1.0.0 output,
as the data they contain was not correctly computed in older versions:

- `resmarker_new_mutations.txt` (present in v0.1.8–v0.2.1)
- `all_mutations_table.txt` (present in v0.2.2)

The equivalent correct output (`all_mutations_table.txt`) is only produced by
running mad4hatter v1.0.0 directly.

---

## `panel_information/`

This folder does not exist in any pre-v1.0.0 output and is generated by the
converter from the reference input files in `inputs/`.

| File                   | Source                                        | Filtering                                              |
|------------------------|-----------------------------------------------|--------------------------------------------------------|
| `amplicon_info.tsv`    | `inputs/amplicon_info.tsv`                    | Rows where `target_name` is present in the run data    |
| `reference.fasta`      | `inputs/references.fasta`                     | Sequences whose name matches a `target_name` in the run |
| `resmarker_info.tsv`   | `inputs/principal_resistance_marker_info_table.tsv` | Rows matching `gene_id`/`gene`/`aa_position` combinations present in `resmarker_table.txt` |

The `target_name` and `codon_start_in_target` columns in `resmarker_info.tsv`
are **set to NA** as this information is not available from the older inputs.

---

## `legacy/`

The `quality_report/` folder is copied as-is into `legacy/` since it no longer
exists in v1.0.0. Its contents (QC plots and reports) remain valid as a record
of the original run.

## Other copied folders

- `run/` — copied as-is to the output root
- `raw_dada2_output/` — copied as-is to the output root
- `dada2_analysis/` (v0.2.0 only) — copied to output root and **renamed to `raw_dada2_output/`** to match v1.0.0 structure

---

## Issues

If you encounter problems or have results from an unsupported version, please
[submit an issue](https://github.com/EPPIcenter/mad4hatter-tools/issues).
