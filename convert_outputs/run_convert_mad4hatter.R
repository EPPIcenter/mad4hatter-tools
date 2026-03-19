
# WRITE PATH TO DIRECTORY TO BE TRANSFORMED
input_dir   = "PATH/TO/DIR"

# WRITE PATH TO NEW DIRECTORY TO BE CREATED
output_dir  = "PATH/TO/NEW_DIR"

# SET PATH TO convert_outputs REPO
pipeline_path = "PATH/TO/convert_outputs"


# target_id_conversion_table.tsv WITH ALL TARGETS (PROVIDED)
locus_lookup = file.path(pipeline_path,"inputs/target_id_conversion_table.tsv")
# references.fasta WITH ALL TARGETS (PROVIDED)
references = file.path(pipeline_path,"inputs/references.fasta")
# amplicon_info.tsv WITH ALL TARGETS (PROVIDED)
amplicon_info = file.path(pipeline_path,"inputs/amplicon_info.tsv")
# principal_resistance_marker_info_table.tsv (PROVIDED)
resmarker_info = file.path(pipeline_path,"inputs/principal_resistance_marker_info_table.tsv")


source(file.path(pipeline_path,"aad/convert_mad4hatter_functions.R"))

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
