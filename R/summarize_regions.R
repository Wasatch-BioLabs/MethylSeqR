#' Summarize Methylation Data by Regions
#'
#' This function summarizes methylation data from a DuckDB database based on specified regions 
#' defined in a BED, TSV, or CSV file. It performs a join operation between the methylation data and 
#' the regions specified in the annotation file, allowing for different types of joins.
#'
#' @param ch3_db A list containing the database file path. This should be a valid "ch3_db" class object.
#' @param region_file A string representing the path to the BED or CSV file that contains the region annotations.
#' @param join_type A string indicating the type of join to perform. Options are "inner", "left", 
#' "right", or "full". Default is "inner".
#'
#' @details
#' The function reads the region annotations from the regional annotation file and checks for its validity.
#' It connects to the DuckDB database, creates a summarized table of methylation data based on the specified 
#' regions, and performs the join operation according to the specified join type. A progress bar is displayed 
#' during the summarization process. The resulting data is stored in a table called "regions" within the database.
#'
#' @return The updated `ch3_db` object with the summarized regions data added to the DuckDB database.
#'
#' @import readr
#' @import dplyr
#' @import dbplyr
#' @import duckdb
#' @import duckplyr
#' @import progress
#'
#' @examples
#'  # Specify the path to the database
#'  ch3_db <- system.file("my_data.ch3.db", package = "MethylSeqR")
#'  
#'  region_bed = system.file("Islands_hg38_test.csv", package = "MethylSeqR")
#'  # Summarize Regions using annotation table
#'  summarize_regions(ch3_db, region_bed)
#'
#' @export
summarize_regions <- function(ch3_db,
                              region_file,
                              join_type = "inner")
{
  # Determine the file type (csv, tsv, or bed)
  file_ext <- tools::file_ext(region_file)
  
  if (file_ext == "csv") {
    annotation <- read_csv(region_file,
                           col_names = c("chrom", "start", "end", "region_name"),
                           show_col_types = FALSE)
  } else if (file_ext %in% c("bed", "tsv")) {
    annotation <- read_tsv(region_file,
                           col_names = c("chrom", "start", "end", "region_name"),
                           show_col_types = FALSE)
  } else {
    stop("Invalid file type. Only CSV, TSV, or BED files are supported.")
  }
  
  # Ensure proper column names and remove header row if needed
  if (annotation[1,1] %in% c("chr", "Chr", "chrom", "Chrom")) {
    annotation = annotation[-1, ]
  }
  
  # check format
  if (ncol(annotation) < 3 || ncol(annotation) > 4) {
    stop("Invalid annotation format. File must have 3 or 4 columns:
         chr, start, end, region_name (optional) annotation.")
  }
  
  if (ncol(annotation == 3)) {
    annotation <-
      annotation |>
      mutate(
        region_name = paste(chrom, start, end, sep = "_"))
  }
  
  # Open the database connection
  database <- .helper_connectDB(ch3_db)
  db_con <- database$db_con
    
    # Specify on exit what to do...
  on.exit(.helper_closeDB(database), add = TRUE)
  
  # Increase temp storage limit to avoid memory issues
  dbExecute(db_con, "PRAGMA max_temp_directory_size='100GiB';")
    
  # Drop the regions table if it already exists
  dbExecute(db_con, "DROP TABLE IF EXISTS regions;")
    
  cat("Building regions table...")
  
  db_tbl = tbl(db_con, "positions")
  
  # Upload annotation as a temporary table
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_annotation;")
  dbWriteTable(db_con, "temp_annotation", annotation, temporary = TRUE)
  
  # Upload positions (db_tbl) as a temporary table
  dbExecute(db_con, "DROP TABLE IF EXISTS temp_positions;")
  dbWriteTable(db_con, "temp_positions", collect(db_tbl), temporary = TRUE)
  
  # Perform the **inequality join** in DuckDB and create the `regions` table
  query <- "
  CREATE TABLE regions AS
  SELECT 
    p.sample_name, 
    a.region_name,
    COUNT(*) AS num_CpGs, 
    SUM(p.cov) AS cov, 
    SUM(p.c_counts) AS c_counts, 
    SUM(p.m_counts) AS m_counts, 
    SUM(p.h_counts) AS h_counts, 
    SUM(p.mh_counts) AS mh_counts, 
    SUM(p.m_counts * p.cov) / NULLIF(SUM(p.cov), 0) AS m_frac,
    SUM(p.h_counts * p.cov) / NULLIF(SUM(p.cov), 0) AS h_frac,
    SUM(p.mh_counts * p.cov) / NULLIF(SUM(p.cov), 0) AS mh_frac
  FROM temp_positions p
  JOIN temp_annotation a
    ON p.chrom = a.chrom 
    AND CAST(p.ref_position AS DOUBLE) BETWEEN CAST(a.start AS DOUBLE) AND CAST(a.end AS DOUBLE)
  GROUP BY p.sample_name, a.region_name;
"
  
  dbExecute(db_con, query)
            
    # Finish up: purge extra tables & update table list and close the connection
  keep_tables = c("positions", "regions", "windows", "meth_diff")
  .helper_purgeTables(db_con, keep_tables)
  
  # Finish Up
  database$last_table = "regions"
  invisible(database)
}
