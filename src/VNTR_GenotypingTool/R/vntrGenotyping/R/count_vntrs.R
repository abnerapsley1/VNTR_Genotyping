#' @importFrom reticulate import py_available py_module_available
NULL

# Module-level cache so we only import once per session
.vg_env <- new.env(parent = emptyenv())

.get_vg <- function() {
  if (is.null(.vg_env$vg)) {
    if (!py_available(initialize = TRUE)) {
      stop("Python is not available. Please configure reticulate (e.g. reticulate::use_python()).")
    }
    if (!py_module_available("vntr_genotyping")) {
      stop(
        "The vntr_genotyping Python package is not installed.\n",
        "Install it with:\n",
        "  pip install <path-to>/src/VNTR_GenotypingTool/Python\n",
        "or activate the conda/venv environment where it is installed."
      )
    }
    .vg_env$vg <- import("vntr_genotyping")
  }
  .vg_env$vg
}


#' Count reads overlapping VNTR regions
#'
#' Calls the \code{vntr_genotyping.count_vntrs()} Python function via
#' reticulate and returns the result as an R data frame.
#'
#' At least one region-selection argument (\code{default}, \code{gene},
#' \code{vntr}, \code{regions}) must be supplied.
#'
#' @param input_files Character vector of BAM/CRAM/SAM file paths.
#' @param default Logical. Use all VNTRs in the built-in BED file.
#' @param gene Character vector of gene names to subset from the built-in BED.
#' @param vntr Character vector of VNTR names to subset from the built-in BED.
#' @param regions Path to a custom BED file (combinable with other selectors).
#' @param gtf Path to a GTF or GTF.gz annotation file. When supplied, the
#'   function returns density ratios (VNTR / gene coverage); otherwise raw
#'   read counts are returned.
#' @param reference Path to a reference FASTA. Required for CRAM files
#'   without an embedded reference.
#' @param psl Path to the UCSC altSeqLiftOverPsl.txt(.gz) for alt-contig
#'   gene normalization. Defaults to the bundled GRCh38 file.
#' @param no_alt_contigs Logical. Disable alt-contig gene normalization.
#' @param norm_gene_types Character vector of gene biotypes eligible for
#'   nearest-gene assignment. Default: \code{"protein_coding"}. Use
#'   \code{"any"} to allow all biotypes.
#' @param output_csv Optional path to write results as a CSV file.
#'
#' @return A data frame with columns \code{sample}, \code{metric}, and one
#'   column per VNTR region. Values are density ratios (numeric) when
#'   \code{gtf} is supplied, or read counts (integer) otherwise. Cells where
#'   normalization could not be performed contain \code{NA}.
#'
#' @examples
#' \dontrun{
#' df <- count_vntrs(
#'   "sample.cram",
#'   default = TRUE,
#'   gtf     = "gencode.v38.annotation.gtf.gz",
#'   reference = "hg38.fa"
#' )
#'
#' # Subset to SLC6A3 VNTRs only
#' df <- count_vntrs(
#'   c("s1.cram", "s2.cram"),
#'   gene      = "SLC6A3",
#'   gtf       = "gencode.v38.annotation.gtf.gz",
#'   reference = "hg38.fa",
#'   output_csv = "slc6a3_results.csv"
#' )
#' }
#'
#' @export
count_vntrs <- function(
  input_files,
  default         = FALSE,
  gene            = NULL,
  vntr            = NULL,
  regions         = NULL,
  gtf             = NULL,
  reference       = NULL,
  psl             = NULL,
  no_alt_contigs  = FALSE,
  norm_gene_types = "protein_coding",
  output_csv      = NULL
) {
  vg <- .get_vg()

  # Build keyword argument list; only pass non-NULL optionals
  kwargs <- list(
    default         = default,
    no_alt_contigs  = no_alt_contigs,
    norm_gene_types = as.list(norm_gene_types)
  )
  if (!is.null(gene))        kwargs$gene        <- as.list(gene)
  if (!is.null(vntr))        kwargs$vntr        <- as.list(vntr)
  if (!is.null(regions))     kwargs$regions     <- regions
  if (!is.null(gtf))         kwargs$gtf         <- gtf
  if (!is.null(reference))   kwargs$reference   <- reference
  if (!is.null(psl))         kwargs$psl         <- psl
  if (!is.null(output_csv))  kwargs$output_csv  <- output_csv

  py_df <- do.call(
    vg$count_vntrs,
    c(list(as.list(input_files)), kwargs)
  )

  # Convert the pandas DataFrame to an R data frame
  reticulate::py_to_r(py_df)
}


#' Return the installed version of the vntr_genotyping Python package
#'
#' @return A character string with the version number, e.g. \code{"0.1.0"}.
#' @export
vntr_genotyping_version <- function() {
  .get_vg()$`__version__`
}
