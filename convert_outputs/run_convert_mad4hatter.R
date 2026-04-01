
# WRITE PATH TO DIRECTORY TO BE TRANSFORMED
input_dir   = "PATH/TO/DIR"

# WRITE PATH TO NEW DIRECTORY TO BE CREATED
output_dir  = "PATH/TO/NEW_DIR"

# SET PATH TO mad4hatter-tools REPO
repo_path = "PATH/TO/mad4hatter-tools"


# target_id_conversion_table.tsv WITH ALL TARGETS (PROVIDED)
locus_lookup = file.path(repo_path,"convert_outputs/inputs/target_id_conversion_table.tsv")
# references.fasta WITH ALL TARGETS (PROVIDED)
references = file.path(repo_path,"convert_outputs/inputs/references.fasta")
# amplicon_info.tsv WITH ALL TARGETS (PROVIDED)
amplicon_info = file.path(repo_path,"convert_outputs/inputs/amplicon_info.tsv")
# principal_resistance_marker_info_table.tsv (PROVIDED)
resmarker_info = file.path(repo_path,"convert_outputs/inputs/principal_resistance_marker_info_table.tsv")

source(file.path(repo_path, "convert_outputs/convert_mad4hatter_functions.R"))

# RUN
convert_mad4hatter(
  input_dir       = input_dir,
  output_dir      = output_dir,
  locus_lookup    = locus_lookup,
  amplicon_info   = amplicon_info,
  references      = references,
  resmarker_info  = resmarker_info,
  #release_version = OPTIONAL
  )

