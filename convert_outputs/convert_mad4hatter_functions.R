library(tidyverse)

# =============================================================================
# convert_mad4hatter_functions.R
#
# Converts mad4hatter results to v1.0.0 format.
# Supported input versions: v0.1.8, v0.1.9, v0.2.0, v0.2.1, v0.2.2
#
# Usage: see run_convert.R
# =============================================================================

# -----------------------------------------------------------------------------
# VERSION DETECTION
# -----------------------------------------------------------------------------


detect_release_version <- function(input_dir) {
  
  # Read allele_data.txt headers to fingerprint release_version
  allele_file <- file.path(input_dir, "allele_data.txt")
  if (!file.exists(allele_file)) {
    stop("Cannot find allele_data.txt in input_dir. Is this a mad4hatter results folder?")
  }
  
  allele_cols <- names(read_tsv(allele_file, n_max = 0, show_col_types = FALSE))
  
  # v1.0.0: sample_name + target_name
  if ("sample_name" %in% allele_cols && "target_name" %in% allele_cols) {
    release_version <- "1.0.0"
    
    # v0.1.8: all-lowercase sampleID/locus/asv/allele/pseudo_cigar
  } else if ("sampleID" %in% allele_cols && "locus" %in% allele_cols && "allele" %in% allele_cols) {
    release_version <- "0.1.8"
    
    # v0.1.9 / v0.2.0 / v0.2.1 / v0.2.2: SampleID/Locus/ASV/Allele — distinguish by dada2 folder or allele_by_locus
  } else if ("SampleID" %in% allele_cols && "Locus" %in% allele_cols && "Allele" %in% allele_cols) {
    release_version <- detect_version_019_020_021_022(input_dir)
    
    } else {
    stop(paste(
      "Could not detect release_version from allele_data.txt column names.",
      "Found columns:", paste(allele_cols, collapse = ", ")
    ))
  }
  
  message(sprintf("Detected release_version: %s", release_version))
  return(release_version)
}


# Secondary fingerprinting for v0.1.9 vs v0.2.0 vs v0.2.1
# These three share identical allele_data.txt headers, distinguished by dada2 folder:
#   v0.1.9: neither dada2_analysis/ nor raw_dada2_output/ present
#   v0.2.0: has dada2_analysis/
#   v0.2.1: has raw_dada2_output/
#   v0.2.2: has resmarker_table_by_locus.txt
detect_version_019_020_021_022 <- function(input_dir) {
  has_dada2_analysis    <- dir.exists(file.path(input_dir, "dada2_analysis"))
  has_raw_dada2_output  <- dir.exists(file.path(input_dir, "raw_dada2_output"))
  has_resmarker_table_by_locus  <- file.exists(file.path(input_dir, "resistance_marker_module/resmarker_table_by_locus.txt"))
  
  if (has_resmarker_table_by_locus){
    return("0.2.2")
  }else if (has_dada2_analysis) {
    return("0.2.0")
  } else if (has_raw_dada2_output) {
    return("0.2.1")
  } else {
    return("0.1.9")
  }
}


# -----------------------------------------------------------------------------
# LOCUS LOOKUP HELPER
# -----------------------------------------------------------------------------

apply_locus_lookup <- function(df, locus_col, lookup) {
  # dplyr::rename the locus column to old_name for joining
  df <- df %>% dplyr::rename(old_name = !!sym(locus_col))
  
  # Check for unmatched loci
  unmatched <- setdiff(unique(df$old_name), lookup$old_name)
  if (length(unmatched) > 0) {
    warning(sprintf(
      "The following loci were NOT found in the lookup file and will be kept as-is with pool = NA:\n  %s",
      paste(unmatched, collapse = "\n  ")
    ))
  }
  
  df <- df %>%
    left_join(lookup, by = "old_name") %>%
    mutate(
      target_name = if_else(is.na(target_name), old_name, target_name),
      pool        = if_else(old_name %in% unmatched, NA_character_, pool)
    ) %>%
    select(-old_name)
  
  return(df)
}

# -----------------------------------------------------------------------------
# FILE CONVERTERS
# -----------------------------------------------------------------------------

convert_allele_data <- function(input_dir, output_dir, release_version, lookup) {
  f <- file.path(input_dir, "allele_data.txt")
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (release_version == "0.1.8") {
    df <- df %>%
      dplyr::rename(
        sample_name        = sampleID,
        pseudocigar_masked = pseudo_cigar
      ) %>%
      select(-allele) %>%
      mutate(asv_masked           = NA_character_,
             pseudocigar_unmasked = NA_character_)
    locus_col <- "locus"
    
  } else if (release_version %in% c("0.1.9", "0.2.0", "0.2.1","0.2.2")) {
    df <- df %>%
      dplyr::rename(
        sample_name        = SampleID,
        asv                = ASV,
        reads              = Reads,
        pseudocigar_masked = PseudoCIGAR
      ) %>%
      select(-Allele) %>%
      mutate(asv_masked           = NA_character_,
             pseudocigar_unmasked = NA_character_)
    locus_col <- "Locus"
    
  } 
  
  df <- apply_locus_lookup(df, locus_col, lookup)
  
  df <- df %>%
    select(sample_name, target_name, asv, pseudocigar_unmasked,
           asv_masked, pseudocigar_masked, reads, pool)
  
  write_tsv(df, file.path(output_dir, "allele_data.txt"))
  message("  [OK] allele_data.txt")
  return(df)
}

convert_allele_data_collapsed <- function(allele_data, output_dir) {
  df <- allele_data %>%
    group_by(sample_name, target_name, asv_masked, pseudocigar_masked, pool) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
    select(sample_name, target_name, asv_masked, pseudocigar_masked, reads, pool)
  
  write_tsv(df, file.path(output_dir, "allele_data_collapsed.txt"))
  message("  [OK] allele_data_collapsed.txt")
}

convert_sample_coverage <- function(input_dir, output_dir, release_version) {
  f  <- file.path(input_dir, "sample_coverage.txt")
  df <- read_tsv(f, show_col_types = FALSE) %>%
    dplyr::rename(sample_name = SampleID,
                  stage       = Stage,
                  reads       = Reads)
  
  write_tsv(df, file.path(output_dir, "sample_coverage.txt"))
  message("  [OK] sample_coverage.txt")
}

convert_amplicon_coverage <- function(input_dir, output_dir, release_version, lookup) {
  f  <- file.path(input_dir, "amplicon_coverage.txt")
  df <- read_tsv(f, show_col_types = FALSE) %>%
    dplyr::rename(sample_name = SampleID,
                  reads       = Reads)
  
  df <- apply_locus_lookup(df, "Locus", lookup) %>%
    select(sample_name, target_name, reads, OutputDada2, OutputPostprocessing)
  
  write_tsv(df, file.path(output_dir, "amplicon_coverage.txt"))
  message("  [OK] amplicon_coverage.txt")
}

convert_resmarker_table <- function(input_dir, output_dir, release_version, lookup) {
  f <- file.path(input_dir, "resistance_marker_module", "resmarker_table.txt")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (release_version == "0.1.8") {
    # No Locus column; has CodonStart — group_by to collapse it out
    df <- df %>%
      dplyr::rename(
        sample_name   = SampleID,
        gene_id       = GeneID,
        gene          = Gene,
        aa_position   = CodonID,
        ref_codon     = RefCodon,
        codon         = Codon,
        codon_ref_alt = CodonRefAlt,
        ref_aa        = RefAA,
        aa            = AA,
        aa_ref_alt    = AARefAlt,
        reads         = Reads
      ) %>%
      group_by(sample_name, gene_id, gene, aa_position, ref_codon, codon,
               codon_ref_alt, ref_aa, aa, aa_ref_alt) %>%
      summarize(reads = sum(reads, na.rm = FALSE), .groups = "drop") %>%
      mutate(follows_indel = NA,
             codon_masked  = NA,
             multiple_loci = NA)
    
  } else if (release_version == "0.1.9") {
    # No Locus column; has CodonStart — group_by to collapse it out
    df <- df %>%
      dplyr::rename(
        sample_name   = SampleID,
        gene_id       = GeneID,
        gene          = Gene,
        aa_position   = CodonID,
        ref_codon     = RefCodon,
        codon         = Codon,
        codon_ref_alt = CodonRefAlt,
        ref_aa        = RefAA,
        aa            = AA,
        aa_ref_alt    = AARefAlt,
        reads         = Reads
      ) %>%
      group_by(sample_name, gene_id, gene, aa_position, ref_codon, codon,
               codon_ref_alt, ref_aa, aa, aa_ref_alt) %>%
      summarize(reads = sum(reads, na.rm = FALSE), .groups = "drop") %>%
      mutate(follows_indel = NA,
             codon_masked  = NA,
             multiple_loci = NA)
    
  } else if (release_version %in% c("0.2.0", "0.2.1")) {
    # Has Locus column + CodonStart — apply lookup then group_by to collapse CodonStart
    df <- df %>%
      dplyr::rename(
        sample_name   = SampleID,
        gene_id       = GeneID,
        gene          = Gene,
        aa_position   = CodonID,
        ref_codon     = RefCodon,
        codon         = Codon,
        codon_ref_alt = CodonRefAlt,
        ref_aa        = RefAA,
        aa            = AA,
        aa_ref_alt    = AARefAlt,
        reads         = Reads
      )
    df <- apply_locus_lookup(df, "Locus", lookup)
    df <- df %>%
      group_by(sample_name, gene_id, gene, aa_position, ref_codon, codon,
               codon_ref_alt, ref_aa, aa, aa_ref_alt) %>%
      summarize(reads = sum(reads, na.rm = FALSE), .groups = "drop") %>%
      mutate(follows_indel = NA,
             codon_masked  = NA,
             multiple_loci = NA)
    
  } else if (release_version == "0.2.2") {
    # No Locus column; no CodonStart; has FollowsIndel/CodonMasked/MultipleLoci
    df <- df %>%
      dplyr::rename(
        sample_name   = SampleID,
        gene_id       = GeneID,
        gene          = Gene,
        aa_position   = CodonID,
        ref_codon     = RefCodon,
        codon         = Codon,
        codon_ref_alt = CodonRefAlt,
        ref_aa        = RefAA,
        aa            = AA,
        aa_ref_alt    = AARefAlt,
        follows_indel = FollowsIndel,
        codon_masked  = CodonMasked,
        multiple_loci = MultipleLoci,
        reads         = Reads
      )
  }
  
  df <- df %>%
    select(sample_name, gene_id, gene, aa_position, ref_codon, codon,
           codon_ref_alt, ref_aa, aa, aa_ref_alt, follows_indel,
           codon_masked, multiple_loci, reads)
  
  write_tsv(df, file.path(output_dir, "resistance_marker_module", "resmarker_table.txt"))
  message("  [OK] resistance_marker_module/resmarker_table.txt")
}

convert_resmarker_table_by_locus <- function(input_dir, output_dir, release_version, lookup) {
  f <- file.path(input_dir, "resistance_marker_module", "resmarker_table_by_locus.txt")
  
  # Only exists in v0.2.2
  if (!file.exists(f)) {
    message("  [SKIP] resistance_marker_module/resmarker_table_by_locus.txt — not present in this version")
    return(invisible(NULL))
  }
  
  df <- read_tsv(f, show_col_types = FALSE) %>%
    dplyr::rename(
      sample_name   = SampleID,
      gene_id       = GeneID,
      gene          = Gene,
      aa_position   = CodonID,
      ref_codon     = RefCodon,
      codon         = Codon,
      codon_ref_alt = CodonRefAlt,
      ref_aa        = RefAA,
      aa            = AA,
      aa_ref_alt    = AARefAlt,
      follows_indel = FollowsIndel,
      codon_masked  = CodonMasked,
      reads         = Reads
    )
  
  df <- apply_locus_lookup(df, "Locus", lookup)
  
  df <- df %>%
    select(sample_name, gene_id, gene, target_name, aa_position, ref_codon,
           codon, codon_ref_alt, ref_aa, aa, aa_ref_alt, follows_indel,
           codon_masked, reads)
  
  write_tsv(df, file.path(output_dir, "resistance_marker_module", "resmarker_table_by_locus.txt"))
  message("  [OK] resistance_marker_module/resmarker_table_by_locus.txt")
}


convert_resmarker_microhaplotype <- function(input_dir, output_dir, release_version, lookup) {
  
  if (release_version %in% c("0.1.8", "0.1.9", "0.2.0", "0.2.1")) {
    f <- file.path(input_dir, "resistance_marker_module", "resmarker_microhap_table.txt")
  } else {
    f <- file.path(input_dir, "resistance_marker_module", "resmarker_microhaplotype_table.txt")
  }
  if (!file.exists(f)) return(invisible(NULL))
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (release_version == "0.1.8") {
    # No Locus column — target_name will be NA
    df <- df %>%
      dplyr::rename(
        sample_name       = SampleID,
        gene_id           = GeneID,
        gene              = Gene,
        mhap_aa_positions = MicrohapIndex,
        ref_mhap          = RefMicrohap,
        mhap              = Microhaplotype,
        mhap_ref_alt      = MicrohapRefAlt,
        reads             = Reads
      ) %>%
      mutate(target_name = NA_character_)
    
  } else if (release_version == "0.1.9") {
    # No Locus column — target_name will be NA
    df <- df %>%
      dplyr::rename(
        sample_name       = SampleID,
        gene_id           = GeneID,
        gene              = Gene,
        mhap_aa_positions = MicrohapIndex,
        ref_mhap          = RefMicrohap,
        mhap              = Microhaplotype,
        mhap_ref_alt      = MicrohapRefAlt,
        reads             = Reads
      ) %>%
      mutate(target_name = NA_character_)
    
  } else if (release_version %in% c("0.2.0", "0.2.1")) {
    # Has Locus column
    df <- df %>%
      dplyr::rename(
        sample_name       = SampleID,
        gene_id           = GeneID,
        gene              = Gene,
        mhap_aa_positions = MicrohapIndex,
        ref_mhap          = RefMicrohap,
        mhap              = Microhaplotype,
        mhap_ref_alt      = MicrohapRefAlt,
        reads             = Reads
      )
    df <- apply_locus_lookup(df, "Locus", lookup)
    
  } else if (release_version == "0.2.2") {
    # Renamed file and columns
    df <- df %>%
      dplyr::rename(
        sample_name       = SampleID,
        gene_id           = GeneID,
        gene              = Gene,
        mhap_aa_positions = MicrohaplotypeCodonIDs,
        ref_mhap          = RefMicrohap,
        mhap              = Microhaplotype,
        mhap_ref_alt      = MicrohapRefAlt,
        reads             = Reads
      )
    df <- apply_locus_lookup(df, "Locus", lookup)
  }
  
  df <- df %>%
    select(sample_name, gene_id, gene, target_name, mhap_aa_positions,
           ref_mhap, mhap, mhap_ref_alt, reads)
  
  write_tsv(df, file.path(output_dir, "resistance_marker_module", "resmarker_microhaplotype_table.txt"))
  message("  [OK] resistance_marker_module/resmarker_microhaplotype_table.txt")
}

convert_all_mutations <- function(input_dir, output_dir, release_version, lookup) {
  # resmarker_new_mutations.txt (v0.1.8-v0.2.1) and all_mutations_table.txt (v0.2.2)
  # are intentionally excluded — prior versions were not correct
  old_files <- c(
    file.path(input_dir, "resistance_marker_module", "resmarker_new_mutations.txt"),
    file.path(input_dir, "resistance_marker_module", "all_mutations_table.txt")
  )
  for (f in old_files) {
    if (file.exists(f)) {
      message(sprintf("  [SKIP] %s — intentionally excluded as versions prior to v1.0.0 were not correct", basename(f)))
    }
  }
}

copy_passthrough_files <- function(input_dir, output_dir, release_version) {
  
  # quality_report goes into legacy/ as it no longer exists in v1.0.0
  legacy_dir <- file.path(output_dir, "legacy")
  src <- file.path(input_dir, "quality_report")
  dst <- file.path(legacy_dir, "quality_report")
  if (dir.exists(src)) {
    dir.create(dst, showWarnings = FALSE, recursive = TRUE)
    rel_paths <- list.files(src, recursive = TRUE)
    src_files <- file.path(src, rel_paths)
    for (i in seq_along(rel_paths)) {
      dest_file <- file.path(dst, rel_paths[i])
      dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
      file.copy(src_files[i], dest_file)
    }
    message("  [COPY -> legacy/] quality_report/")
  }
  
  # run and raw_dada2_output copied as-is to output root
  for (d in c("run", "raw_dada2_output")) {
    src <- file.path(input_dir, d)
    dst <- file.path(output_dir, d)
    if (dir.exists(src)) {
      dir.create(dst, showWarnings = FALSE, recursive = TRUE)
      rel_paths <- list.files(src, recursive = TRUE)
      src_files <- file.path(src, rel_paths)
      for (i in seq_along(rel_paths)) {
        dest_file <- file.path(dst, rel_paths[i])
        dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
        file.copy(src_files[i], dest_file)
      }
      message(sprintf("  [COPY] %s/", d))
    }
  }
  
  # dada2_analysis/ (v0.2.0 only) renamed to raw_dada2_output/ to match v1.0.0
  src <- file.path(input_dir, "dada2_analysis")
  dst <- file.path(output_dir, "raw_dada2_output")
  if (dir.exists(src)) {
    dir.create(dst, showWarnings = FALSE, recursive = TRUE)
    rel_paths <- list.files(src, recursive = TRUE)
    src_files <- file.path(src, rel_paths)
    for (i in seq_along(rel_paths)) {
      dest_file <- file.path(dst, rel_paths[i])
      dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
      file.copy(src_files[i], dest_file)
    }
    message("  [COPY] dada2_analysis/ -> raw_dada2_output/")
  }
}


# -----------------------------------------------------------------------------
# PANEL INFORMATION
# -----------------------------------------------------------------------------

create_panel_information <- function(output_dir, amplicon_info, references, resmarker_info) {
  
  panel_dir <- file.path(output_dir, "panel_information")
  dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get target_names present in this run from the already-converted amplicon_coverage.txt
  targets_in_data <- read_tsv(file.path(output_dir, "amplicon_coverage.txt"),
                              show_col_types = FALSE) %>%
    pull(target_name) %>%
    unique()
  
  # --- amplicon_info.tsv -----------------------------------------------------
  if (!file.exists(amplicon_info)) {
    warning("amplicon_info file not found, skipping: ", amplicon_info)
  } else {
    df_amplicon <- read_tsv(amplicon_info, show_col_types = FALSE)
    
    missing_targets <- setdiff(targets_in_data, df_amplicon$target_name)
    if (length(missing_targets) > 0) {
      warning(sprintf(
        "The following target_names from the data were NOT found in amplicon_info:\n  %s",
        paste(missing_targets, collapse = "\n  ")
      ))
    }
    
    df_amplicon %>%
      filter(target_name %in% targets_in_data) %>%
      write_tsv(file.path(panel_dir, "amplicon_info.tsv"))
    message("  [OK] panel_information/amplicon_info.tsv")
  }
  
  # --- reference.fasta -------------------------------------------------------
  if (!file.exists(references)) {
    warning("references fasta file not found, skipping: ", references)
  } else {
    fasta_lines <- readLines(references)
    
    header_idx <- which(startsWith(fasta_lines, ">"))
    seq_names  <- sub("^>([^ ]+).*", "\\1", fasta_lines[header_idx])
    
    entry_start <- header_idx
    entry_end   <- c(header_idx[-1] - 1, length(fasta_lines))
    
    missing_refs <- setdiff(targets_in_data, seq_names)
    if (length(missing_refs) > 0) {
      warning(sprintf(
        "The following target_names from the data were NOT found in the reference fasta:\n  %s",
        paste(missing_refs, collapse = "\n  ")
      ))
    }
    
    keep       <- seq_names %in% targets_in_data
    kept_lines <- unlist(mapply(
      function(s, e) fasta_lines[s:e],
      entry_start[keep], entry_end[keep],
      SIMPLIFY = FALSE
    ))
    
    writeLines(kept_lines, file.path(panel_dir, "reference.fasta"))
    message("  [OK] panel_information/reference.fasta")
  }
  
  # --- resmarker_info.tsv ----------------------------------------------------
  if (!file.exists(resmarker_info)) {
    warning("resmarker_info file not found, skipping: ", resmarker_info)
  } else {
    resmarker_table_path <- file.path(output_dir, "resistance_marker_module", "resmarker_table.txt")
    
    if (!file.exists(resmarker_table_path)) {
      message("  [SKIP] panel_information/resmarker_info.tsv — no resmarker_table.txt in output")
    } else {
      resmarker_keys <- read_tsv(resmarker_table_path, show_col_types = FALSE) %>%
        distinct(gene_id, gene, aa_position)
      
      df_resmarker <- read_tsv(resmarker_info, show_col_types = FALSE)
      
      missing_keys <- anti_join(resmarker_keys, df_resmarker,
                                by = c("gene_id", "gene", "aa_position"))
      if (nrow(missing_keys) > 0) {
        warning(sprintf(
          "The following gene_id/gene/aa_position combinations from resmarker_table.txt were NOT found in resmarker_info:\n  %s",
          paste(apply(missing_keys, 1, paste, collapse = " / "), collapse = "\n  ")
        ))
      }
      
      df_resmarker %>%
        semi_join(resmarker_keys, by = c("gene_id", "gene", "aa_position")) %>%
        mutate(target_name           = NA_character_,
               codon_start_in_target = NA_character_) %>%
        write_tsv(file.path(panel_dir, "resmarker_info.tsv"))
      message("  [OK] panel_information/resmarker_info.tsv")
    }
  }
}


# -----------------------------------------------------------------------------
# MAIN FUNCTION
# -----------------------------------------------------------------------------

convert_mad4hatter <- function(input_dir, output_dir, locus_lookup,
                               amplicon_info, references, resmarker_info,
                               release_version = NULL) {
  
  # --- Validate inputs -------------------------------------------------------
  if (!dir.exists(input_dir))     stop("input_dir does not exist: ", input_dir)
  if (!file.exists(locus_lookup)) stop("locus_lookup file not found: ", locus_lookup)
  
  # --- Validate outputs ------------------------------------------------------
  if (input_dir==output_dir){
    output_dir = paste0(output_dir,"_converted")
    message(paste("\033[1mOutput directory had the same name as the input directory, converted data will be saved into:\033[0m",
                  output_dir)
    )
  } 
  
  # --- Detect version --------------------------------------------------------
  if (is.null(release_version)) {
    release_version <- detect_release_version(input_dir)
  } else {
    message(sprintf("Using user-supplied release_version: %s", release_version))
  }
  
  if (release_version == "1.0.0") {
    message("Input already appears to be v1.0.0. No conversion needed.")
    return(invisible(NULL))
  }
  
  supported <- c("0.1.8", "0.1.9", "0.2.0", "0.2.1", "0.2.2")
  if (!release_version %in% supported) {
    stop(sprintf("Unsupported release_version '%s'. Supported: %s",
                 release_version, paste(supported, collapse = ", ")))
  }
  
  # --- Load lookup -----------------------------------------------------------
  lookup <- read_tsv(locus_lookup, show_col_types = FALSE)
  required_cols <- c("old_name", "new_name", "pool")
  missing_cols  <- setdiff(required_cols, names(lookup))
  if (length(missing_cols) > 0) {
    stop("locus_lookup is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  lookup <- lookup %>% dplyr::rename(target_name = new_name)
  
  # --- Create output directories ---------------------------------------------
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "resistance_marker_module"), showWarnings = FALSE)
  
  message(sprintf("\nConverting from v%s -> v1.0.0", release_version))
  message(sprintf("Input:  %s", input_dir))
  message(sprintf("Output: %s\n", output_dir))
  
  # --- Convert files ---------------------------------------------------------
  allele_data <- convert_allele_data(input_dir, output_dir, release_version, lookup)
  convert_allele_data_collapsed(allele_data, output_dir)
  convert_sample_coverage(input_dir, output_dir, release_version)
  convert_amplicon_coverage(input_dir, output_dir, release_version, lookup)
  convert_resmarker_table(input_dir, output_dir, release_version, lookup)
  convert_resmarker_table_by_locus(input_dir, output_dir, release_version, lookup)
  convert_resmarker_microhaplotype(input_dir, output_dir, release_version, lookup)
  convert_all_mutations(input_dir, output_dir, release_version, lookup)
  copy_passthrough_files(input_dir, output_dir, release_version)
  
  # --- Create panel_information ----------------------------------------------
  create_panel_information(output_dir, amplicon_info, references, resmarker_info)
  
  message(sprintf("\nDone. Output written to: %s", output_dir))
  invisible(output_dir)
}
