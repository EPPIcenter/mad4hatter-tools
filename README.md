# mad4hatter-tools

Utility scripts for working with [mad4hatter](https://github.com/EPPIcenter/mad4hatter) pipeline outputs.

## Contents

### `convert_outputs/`

Scripts to convert mad4hatter results from older releases to v1.0.0 format. Supported input versions: v0.1.8, v0.1.9, v0.2.0, v0.2.1, v0.2.2. **Releases prior to v0.1.8 are not supported.**

See [`convert_outputs/CONVERSION_NOTES.md`](convert_outputs/CONVERSION_NOTES.md) for a full description of column mappings, files that are renamed or reorganized, and fields that cannot be recovered from older outputs.

#### Files

- **`run_convert_mad4hatter.R`** — main script to edit and run. Set your input/output paths here and source the functions file.
- **`run_convert_mad4hatter_local.R`** — same as above, for local path configurations.
- **`convert_mad4hatter_functions.R`** — all conversion logic. Do not need to edit this.
- **`inputs/`** — reference files required by the converter:
  - `target_id_conversion_table.tsv` — maps old locus names to new `target_name` and `pool`
  - `amplicon_info.tsv` — amplicon panel information for `panel_information/` output
  - `references.fasta` — reference sequences for `panel_information/` output
  - `principal_resistance_marker_info_table.tsv` — resistance marker reference for `panel_information/` output

#### Usage

Edit `run_convert_mad4hatter.R` with your paths and run it from RStudio or the R console. Requires the `tidyverse` package.

```r
# Minimal example
input_dir      = "path/to/old/results"
output_dir     = "path/to/new/results"
pipeline_path  = "path/to/mad4hatter-tools"
```

## Issues

If you encounter problems or results from an unsupported version, please [submit an issue](https://github.com/EPPIcenter/mad4hatter-tools/issues).
