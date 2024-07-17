.fisher_pval <- function (q, m, n, k) {
  # derived from https://github.com/al2na/methylKit/issues/96
  
  mat <- cbind(q, m, n, k)
  
  apply(mat, 1, function(qmnk){
    dhyper_val <- 0.5 * dhyper(x = qmnk[1], m = qmnk[2], 
                               n = qmnk[3], k = qmnk[4])
    
    pval_right <- phyper(q = qmnk[1], m = qmnk[2], 
                         n = qmnk[3], k = qmnk[4], 
                         lower.tail = FALSE) + dhyper_val
    
    pval_left  <- phyper(q = qmnk[1] - 1, m = qmnk[2], 
                         n = qmnk[3], k = qmnk[4], 
                         lower.tail = TRUE) + dhyper_val
    
    return(ifelse(test = pval_right > pval_left, 
                  yes  = pval_left * 2, 
                  no   = pval_right * 2))
  })
}

.calc_diff_pos_fisher <- function(modseq_dat, 
                                  case, 
                                  control, 
                                  mod_type = "mh") 
{
  # Set stat to use
  mod_counts_col <- paste0(mod_type[1], "_counts")
  
  # Select necessary columns
  dat <- 
    modseq_dat |>
    select(
      sample_name, chrom, ref_position, c_counts, mod_counts = !!mod_counts_col) |>
    mutate(
      exp_group = case_when(
        sample_name %in% case ~ "case",
        sample_name %in% control ~ "control",
        TRUE ~ NA)) |>
    filter(
      !is.na(exp_group)) |>
    summarize(
      .by = c(exp_group, chrom, ref_position),
      c_counts = sum(c_counts),
      mod_counts = sum(mod_counts)) |>
    pivot_wider(
      id_cols = c(chrom, ref_position),
      names_from = exp_group, 
      values_from = c(c_counts, mod_counts),
      values_fill = 0) 
  
  # Extract Matrix and calculate pvals
  pvals <- 
    dat |>
    select(!chrom:ref_position) |>
    distinct() |>
    collect() |>
    mutate(
      mod_frac_case = mod_counts_case /
        sum(mod_counts_case, c_counts_case),
      mod_frac_control = mod_counts_control /
        sum(mod_counts_control, c_counts_control),
      meth_diff = mod_counts_case / sum(mod_counts_case, c_counts_case) -
        mod_counts_control / sum(mod_counts_control, c_counts_control),
      p_val = .fisher_pval(
        q = mod_counts_case,
        m = mod_counts_case + mod_counts_control,
        n = c_counts_case + c_counts_control,
        k = mod_counts_case + c_counts_case)) 
  
  dat |>
    inner_join(
      pvals,
      by = join_by(c_counts_case, c_counts_control, 
                   mod_counts_case, mod_counts_control),
      copy = TRUE) |>
    rename_with(
      ~ gsub("mod", mod_type[1], .x))
}

.calc_diff_pos_fisher2 <- function(modseq_dat, 
                                  case, 
                                  control, 
                                  mod_type = "mh") 
{
  # Set stat to use
  mod_counts_col <- paste0(mod_type[1], "_counts")
  
  # Select necessary columns
  modseq_dat |>
    select(
      sample_name, chrom, ref_position, c_counts, mod_counts = !!mod_counts_col) |>
    mutate(
      exp_group = case_when(
        sample_name %in% case ~ "case",
        sample_name %in% control ~ "control",
        TRUE ~ NA)) |>
    filter(
      !is.na(exp_group)) |>
    summarize(
      .by = c(exp_group, chrom, ref_position),
      c_counts = sum(c_counts),
      mod_counts = sum(mod_counts)) |>
    pivot_wider(
      id_cols = c(chrom, ref_position),
      names_from = exp_group, 
      values_from = c(c_counts, mod_counts),
      values_fill = 0) |>
    mutate(
      .by = c(chrom, ref_position),
      p_val = fisher.test(matrix(c_across(everything()), 2))$p.val,
      mod_frac_case = mod_counts_case /
        sum(mod_counts_case, c_counts_case),
      mod_frac_control = mod_counts_control /
        sum(mod_counts_control, c_counts_control),
      meth_diff = mod_counts_case / sum(mod_counts_case, c_counts_case) -
        mod_counts_control / sum(mod_counts_control, c_counts_control)) |>
    rename_with(
      ~ gsub("mod", mod_type[1], .x))
}
