#' Update FCS file references in a FlowJo workspace
#'
#' Replace FCS filenames or paths using an explicit mapping or matches to files
#' on disk, and save an updated workspace. Matching text is replaced without
#' reserializing the XML, preserving all other bytes, including line endings.
#' FCS files themselves are neither renamed nor modified.
#'
#' @param ws A single path to an existing FlowJo \code{.wsp} file.
#' @param replacements A nonempty named character vector, or \code{NULL}.
#'   Names identify old filenames, paths, or local file URIs in the workspace's
#'   \code{DataSet} references; values specify their replacements. Keys are
#'   literal strings, not regular expressions, and must be unique. A new
#'   filename without a directory preserves the old directory. A new path or
#'   file URI also changes the directory; relative replacement paths containing
#'   a slash are resolved against the current working directory. Explicit
#'   replacement targets need not exist on disk. Supply exactly one of
#'   \code{replacements} and \code{fcs_dir}.
#' @param fcs_dir A single existing directory to search for files with an
#'   \code{.fcs} extension, detected case-insensitively, or \code{NULL}.
#'   Successful matches replace both the filename and path with the matched
#'   file's location on disk.
#' @param out Output workspace path. Defaults to the input basename with
#'   \code{_updated.wsp} appended after removing its extension, in the input
#'   directory. The output directory must already exist. The output must differ
#'   from the input, including when \code{overwrite = TRUE}.
#' @param recursive Logical; search subdirectories of \code{fcs_dir}?
#'   Defaults to \code{TRUE}. Ignored for explicit replacements.
#' @param dry_run Logical; return a per-sample match report without writing a
#'   file? Defaults to \code{FALSE}. Ambiguous matches and collisions are
#'   reported rather than stopping a dry run. Input and output-path validation
#'   still apply. A dry run does not perform the final replacement, XML, or
#'   output-write checks.
#' @param overwrite Logical; allow replacing an existing output file?
#'   Defaults to \code{FALSE}. The input workspace is never overwritten.
#' @param exact_match Logical; require identical, case-sensitive filenames
#'   when searching \code{fcs_dir}? Defaults to \code{TRUE}. Directories may
#'   differ. If \code{FALSE}, select the unique candidate with the smallest
#'   finite string distance. Ignored for explicit replacements.
#' @param method String-distance method passed to
#'   \code{stringdist::stringdist()}. Defaults to \code{"osa"} (optimal string
#'   alignment). Other supported values are \code{"lv"}, \code{"dl"},
#'   \code{"hamming"}, \code{"lcs"}, \code{"qgram"}, \code{"cosine"},
#'   \code{"jaccard"}, \code{"jw"}, and \code{"soundex"}. Used only when
#'   \code{fcs_dir} is supplied and \code{exact_match = FALSE}.
#'
#' @details
#' In folder-search mode, matching compares complete basenames, including their
#' extensions, without stripping timestamps or applying a distance threshold.
#' Each workspace entry independently chooses its closest candidate. There is
#' no global assignment step or fallback to a second-best candidate. Tied best
#' candidates are ambiguous, and a disk file selected by multiple workspace
#' entries is a collision. Either condition prevents saving.
#'
#' Partial matches are allowed. Unmatched workspace entries retain their
#' original references, and unused disk files are ignored. A nearest match is
#' attempted whenever a finite distance is available, even if the filenames
#' are dissimilar. With an empty search directory, all entries are unmatched.
#'
#' Explicit mappings may update only a subset of samples. An old key absent
#' from the workspace, or multiple keys referring to the same sample, causes
#' an error. A filename key can refer to multiple samples. Assigning different
#' original file paths to a shared target also prevents saving.
#'
#' Replacements operate on matching filename and path strings throughout the
#' XML, including sample display names and other textual references. They are
#' simultaneous: an A-to-B replacement does not subsequently become C through
#' a B-to-C replacement. A shared string requiring conflicting replacements
#' causes an error. Existing \code{$FIL} values change only when their text
#' matches a replacement.
#'
#' Supports UTF-8 and ASCII workspaces, including CRLF line endings and a UTF-8
#' byte-order mark. Paths use forward slashes; remote file URI hosts are not
#' supported. The resulting XML and DataSet URIs are validated before writing
#' a verified temporary file and moving it to \code{out}. Unusual XML escaping
#' that prevents a required replacement causes validation to fail without
#' saving the output. The \pkg{xml2} package is required; nearest matching also
#' requires \pkg{stringdist}.
#'
#' @return With \code{dry_run = TRUE}, a data frame with one row per sample
#'   DataSet reference. Otherwise, invisibly returns a list containing
#'   \code{file}, the saved output path, and \code{changes}, the same report.
#'   The report has the following columns:
#'   \describe{
#'     \item{sampleID}{Sample identifier from the workspace.}
#'     \item{old_filename, old_path}{Original basename and decoded path.}
#'     \item{new_filename, new_path}{Proposed basename and path. Entries without
#'       a selected match retain their original values. Collision rows may
#'       contain proposed values that cannot be saved.}
#'     \item{method}{Matching route: \code{"explicit"}, \code{"exact"},
#'       \code{"nearest"}, or \code{"unchanged"}. This is not the name of the
#'       string-distance algorithm supplied through the method argument.}
#'     \item{status}{\code{"matched"} for a changed URI, \code{"unchanged"}
#'       for an unchanged URI, \code{"unmatched"} if no candidate is available,
#'       \code{"ambiguous"} for tied candidates, or \code{"collision"} for
#'       conflicting target assignments.}
#'     \item{distance}{Minimum finite string distance for nearest matching,
#'       zero for an exact match, or NA when no distance was calculated.}
#'     \item{candidates}{Newline-separated candidate paths in folder-search
#'       mode; an empty string when no candidate is available or when using
#'       explicit replacements.}
#'   }
#'
#' @seealso \code{stringdist::stringdist()}, \code{xml2::read_xml()}
#'
#' @examples
#' \dontrun{
#' ws <- "path/to/workspace.wsp"
#'
#' # Rename a file reference while retaining its original directory.
#' replacements <- c(
#'   "0037_-_40_-_2020.03.04-09.24.31_-_PBMCs_unstained.fcs" =
#'     "0037_-_PBMCs_unstained.fcs"
#' )
#' result <- wsx_relink_fcs(ws, replacements = replacements)
#'
#' # An absolute replacement path also changes the directory.
#' result <- wsx_relink_fcs(
#'   ws, replacements = c("old.fcs" = "/data/FCS/new.fcs"),
#'   out = "path/to/relocated.wsp"
#' )
#'
#' # Locate identical filenames in a different directory.
#' plan <- wsx_relink_fcs(
#'   ws, fcs_dir = "/data/FCS", exact_match = TRUE, dry_run = TRUE
#' )
#'
#' # Inspect nearest matches using optimal string alignment distance.
#' plan <- wsx_relink_fcs(
#'   ws, fcs_dir = "/data/FCS", exact_match = FALSE,
#'   method = "osa", dry_run = TRUE
#' )
#' subset(plan, status %in% c("ambiguous", "collision", "unmatched"))
#'
#' # Save after resolving any ambiguities or collisions.
#' result <- wsx_relink_fcs(
#'   ws, fcs_dir = "/data/FCS", exact_match = FALSE,
#'   out = "path/to/relinked.wsp"
#' )
#' result$file
#' result$changes
#' }
#' @export
wsx_relink_fcs <- function(ws,
                           replacements = NULL,
                           fcs_dir = NULL,
                           out = NULL,
                           recursive = TRUE,
                           dry_run = FALSE,
                           overwrite = FALSE,
                           exact_match = TRUE,
                           method = "osa") {
  if (!requireNamespace("xml2", quietly = TRUE))
    stop("Install xml2 first: install.packages('xml2')", call. = FALSE)
  if (length(ws) != 1L || is.na(ws) || !file.exists(ws))
    stop("ws must name one existing workspace.", call. = FALSE)
  if (is.null(replacements) == is.null(fcs_dir))
    stop("Supply exactly one of replacements or fcs_dir.", call. = FALSE)
  if (!is.logical(exact_match) || length(exact_match) != 1L || is.na(exact_match))
    stop("exact_match must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(fcs_dir) && !exact_match) {
    if (!requireNamespace("stringdist", quietly = TRUE))
      stop("Install stringdist first: install.packages('stringdist')", call. = FALSE)
    method <- match.arg(method, c("osa", "lv", "dl", "hamming", "lcs",
                                  "qgram", "cosine", "jaccard", "jw", "soundex"))
  }
  ws <- normalizePath(ws, mustWork = TRUE)
  if (is.null(out)) out <- paste0(tools::file_path_sans_ext(ws), "_updated.wsp")
  if (length(out) != 1L || is.na(out) || !nzchar(out))
    stop("out must be one output path.", call. = FALSE)
  out_dir <- normalizePath(dirname(path.expand(out)), mustWork = TRUE)
  out <- file.path(out_dir, basename(out))
  if (identical(ws, normalizePath(out, mustWork = FALSE)))
    stop("Choose an output different from the input workspace.", call. = FALSE)

  bytes <- readBin(ws, "raw", n = file.info(ws)$size)
  if (any(bytes == as.raw(0))) stop("Only UTF-8 workspaces are supported.")
  txt <- rawToChar(bytes)
  if (!validUTF8(txt)) stop("Only UTF-8 workspaces are supported.")
  doc <- xml2::read_xml(bytes, options = "NONET")
  # Do not allow a different declared encoding even if its current text is ASCII.
  encoding <- regmatches(txt, regexec("encoding=[\"']([^\"']+)[\"']", substr(txt, 1, 200)))[[1L]]
  if (length(encoding) && !toupper(encoding[2L]) %in% c("UTF-8", "UTF8", "US-ASCII"))
    stop("Only UTF-8/ASCII workspaces are supported.")
  datasets <- xml2::xml_find_all(doc, ".//Sample/DataSet")
  if (!length(datasets)) stop("No Sample/DataSet elements found.")
  old_uri <- xml2::xml_attr(datasets, "uri")
  if (anyNA(old_uri) || any(!nzchar(old_uri))) stop("Missing DataSet URI.")

  from_uri <- function(x) {
    if (grepl("^file:", x)) {
      x <- sub("^file://localhost/", "/", x)
      x <- sub("^file:", "", x)
      if (grepl("^//[^/]", x)) stop("Remote file URI hosts are unsupported: ", x)
      x <- sub("^/+", "/", x)
      x <- utils::URLdecode(x)
    }
    x
  }
  encode_path <- function(x) {
    # Encode path segments independently so / retains its separator meaning.
    paste(vapply(strsplit(x, "/", fixed = TRUE)[[1L]],
                 utils::URLencode, character(1), reserved = TRUE, repeated = TRUE),
          collapse = "/")
  }
  to_uri <- function(x) paste0("file:", encode_path(x))
  old_path <- vapply(old_uri, from_uri, character(1), USE.NAMES = FALSE)
  old_name <- basename(old_path)
  report <- data.frame(
    sampleID = xml2::xml_attr(datasets, "sampleID"),
    old_filename = old_name, old_path = old_path,
    new_filename = old_name, new_path = old_path,
    method = "unchanged", status = "unchanged", distance = NA_real_, candidates = "",
    stringsAsFactors = FALSE
  )
  new_uri <- old_uri
  targets <- rep(NA_character_, length(old_uri))

  if (!is.null(replacements)) {
    if (!is.character(replacements) || !length(replacements) ||
        is.null(names(replacements)) || anyNA(replacements) ||
        anyNA(names(replacements)) || any(!nzchar(replacements)) ||
        any(!nzchar(names(replacements))) || anyDuplicated(names(replacements)))
      stop("replacements must be a nonempty named character vector with unique old strings.")
    keys <- names(replacements)
    matched <- lapply(keys, function(k) {
      which(if (grepl("^file:", k)) old_uri == k | old_path == from_uri(k)
            else if (grepl("/", k, fixed = TRUE)) old_path == path.expand(k)
            else old_name == k)
    })
    if (any(lengths(matched) == 0L))
      stop("Mapping keys not found in DataSet references: ",
           paste(keys[lengths(matched) == 0L], collapse = ", "))
    if (anyDuplicated(unlist(matched)))
      stop("Multiple mapping keys refer to the same sample; use only one.")
    for (j in seq_along(keys)) {
      targets[matched[[j]]] <- replacements[[j]]
      report$method[matched[[j]]] <- "explicit"
    }
  } else {
    if (length(fcs_dir) != 1L || is.na(fcs_dir) || !dir.exists(fcs_dir))
      stop("fcs_dir must be an existing directory.")
    files <- list.files(fcs_dir, pattern = "[.]fcs$", ignore.case = TRUE,
                        full.names = TRUE, recursive = recursive)
    files <- files[!dir.exists(files)]
    files <- unique(normalizePath(files, mustWork = TRUE))
    disk_names <- basename(files)
    for (i in seq_along(old_name)) {
      report$method[i] <- if (exact_match) "exact" else "nearest"
      candidates <- integer()

      if (length(files)) {
        if (exact_match) {
          candidates <- which(disk_names == old_name[i])
          if (length(candidates)) report$distance[i] <- 0L
        } else {
          distances <- stringdist::stringdist(old_name[i], disk_names, method = method)
          finite <- which(is.finite(distances))
          if (length(finite)) {
            report$distance[i] <- min(distances[finite])
            candidates <- finite[distances[finite] == report$distance[i]]
          }
        }
      }

      report$candidates[i] <- paste(files[candidates], collapse = "\n")
      if (length(candidates) == 1L) {
        targets[i] <- files[candidates]
      } else {
        report$status[i] <- if (length(candidates)) "ambiguous" else "unmatched"
        # Preserve original references for entries without a selected match.
      }
    }
  }


  for (i in which(!is.na(targets))) {
    target <- from_uri(targets[i])
    if (grepl("/", target, fixed = TRUE)) {
      target <- path.expand(target)
      if (!startsWith(target, "/")) target <- file.path(getwd(), target)
      target <- normalizePath(target, mustWork = FALSE)
      new_uri[i] <- to_uri(target)
    } else {
      target <- file.path(dirname(old_path[i]), target)
      # Keep the directory's original URI spelling in filename-only mode.
      new_uri[i] <- paste0(sub("[^/]*$", "", old_uri[i]),
                           utils::URLencode(basename(target), reserved = TRUE, repeated = TRUE))
    }
    report$new_path[i] <- target
    report$new_filename[i] <- basename(target)
    report$status[i] <- if (old_uri[i] == new_uri[i]) "unchanged" else "matched"
  }
  # Folder mode: each disk file can be selected by at most one workspace entry,
  # even if the workspace contains repeated references to the same source path.
  if (!is.null(fcs_dir)) {
    reused <- !is.na(targets) &
      (duplicated(targets) | duplicated(targets, fromLast = TRUE))
    report$status[reused] <- "collision"
  }
  # Also detect a changed target that collides with an unchanged source file.
  for (target in unique(stats::na.omit(report$new_path))) {
    ii <- which(report$new_path == target)
    if (length(unique(report$old_path[ii])) > 1L) report$status[ii] <- "collision"
  }
  if (dry_run) return(report)
  bad <- report$status %in% c("ambiguous", "collision")
  if (any(bad)) {
    detail <- paste(report$old_filename[bad], report$status[bad], sep = ": ")
    stop("No file written. Resolve these matches (inspect dry_run = TRUE):\n",
         paste(detail, collapse = "\n"), call. = FALSE)
  }
  if (file.exists(out) && !overwrite) stop("Output exists; choose another out or use overwrite = TRUE.")

  escape_xml <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    gsub("'", "&apos;", x, fixed = TRUE)
  }
  pairs <- data.frame(old = character(), new = character())
  add_pair <- function(a, b) {
    pairs <<- rbind(pairs, data.frame(old = enc2utf8(a), new = enc2utf8(b)))
  }
  for (i in seq_along(old_uri)) {
    # Include unchanged references to detect conflicting basename replacements.
    add_pair(escape_xml(old_uri[i]), escape_xml(new_uri[i]))
    add_pair(escape_xml(old_path[i]), escape_xml(report$new_path[i]))
    add_pair(escape_xml(old_name[i]), escape_xml(report$new_filename[i]))
    add_pair(utils::URLencode(old_name[i], reserved = TRUE, repeated = TRUE),
             utils::URLencode(report$new_filename[i], reserved = TRUE, repeated = TRUE))
  }
  pairs <- unique(pairs)
  if (anyDuplicated(pairs$old))
    stop("A shared filename/path would need different replacements. No file written.")
  pairs <- pairs[pairs$old != pairs$new, , drop = FALSE]
  pairs <- pairs[order(nchar(pairs$old, type = "bytes"), decreasing = TRUE), , drop = FALSE]

  # Replace longest strings first through unique placeholders, preventing cascades
  # such as A -> B followed by B -> C. No XML serialization or newline conversion.
  token <- "__WSX_RELINK_PLACEHOLDER__"
  while (any(grepl(token, c(txt, pairs$old, pairs$new), fixed = TRUE))) token <- paste0(token, "X")
  tokens <- paste0(token, seq_len(nrow(pairs)), "__")
  edited <- txt
  for (i in seq_len(nrow(pairs)))
    edited <- gsub(pairs$old[i], tokens[i], edited, fixed = TRUE, useBytes = TRUE)
  for (i in seq_len(nrow(pairs)))
    edited <- gsub(tokens[i], pairs$new[i], edited, fixed = TRUE, useBytes = TRUE)
  edited_raw <- charToRaw(edited)
  check <- xml2::read_xml(edited_raw, options = "NONET")
  actual <- xml2::xml_attr(xml2::xml_find_all(check, ".//Sample/DataSet"), "uri")
  if (!identical(unname(actual), unname(new_uri)))
    stop("URI validation failed (possibly unusual XML escaping). No file written.")

  tmp <- tempfile(".wsx-relink-", tmpdir = out_dir)
  on.exit(unlink(tmp), add = TRUE)
  writeBin(edited_raw, tmp)
  if (!identical(readBin(tmp, "raw", n = file.info(tmp)$size), edited_raw))
    stop("Output verification failed.")
  if (file.exists(out) && !overwrite) stop("Output appeared during processing; no file replaced.")
  if (!file.rename(tmp, out)) stop("Could not move the verified output into place: ", out)
  message("Saved ", out, " (", sum(report$status == "matched"), " sample references changed).")
  invisible(list(file = out, changes = report))
}
