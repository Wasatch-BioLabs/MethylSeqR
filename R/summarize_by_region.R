summarize_by_region <- function(modseq_dat, 
                                annot_file)
{
  # Read annotation
  annotation <- 
    read_tsv(annot_file, col_names = c("chrom", "start", "end"))
  
  # Make sure annotation file is formatted correctly
  if (ncol(annotation) < 3 || ncol(annotation) > 4) {
    stop("Annotation must be in correct format. 
         ex. chrom, start, end, region_name OR chrom, start, end")
  }
  
  # Create regional dataframe
  modseq_dat <- 
    right_join(
      modseq_dat, annotation, 
      by = join_by(chrom, between(ref_position, start, end)),
      copy = TRUE)
  
  # Calculate the requested scores- summarize scores based on available columns
  modseq_dat |>
    summarize( 
      .by = c(sample_name, region_name),
      mean_mh_frac = if_else("mh_frac" %in% colnames(modseq_dat), 
                             sum(cov * mh_frac) / sum(cov), NA_real_),
      mean_m_frac = if_else("m_frac" %in% colnames(modseq_dat), 
                            sum(cov * m_frac) / sum(cov), NA_real_),
      mean_h_frac = if_else("h_frac" %in% colnames(modseq_dat), 
                            sum(cov * h_frac) / sum(cov), NA_real_),
      mean_cov = mean(cov, na.rm = TRUE)) |>
    select_if(
      ~ !all(is.na(.)))
}
