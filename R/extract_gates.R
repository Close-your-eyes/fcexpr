#' Extract sample-specific gate definitions from a FlowJo workspace
#'
#' Read a FlowJo workspace and return one row for each sample population with
#' a \code{Gate} element. Extracts gate metadata, dimension bounds, and polygon
#' vertices without modifying the workspace or reading FCS event data.
#'
#' @param ws A single path to a FlowJo \code{.wsp} file, or an XML document
#'   created by \code{xml2::read_xml()}.
#'
#' @details
#' Join the result to \code{fcexpr::wsx_get_popstats(ws)$counts} using
#' \code{FileName} and \code{PopulationFullPath}. Both tables must come from
#' the same workspace. \code{FileName} is derived from the URL-decoded
#' \code{DataSet} URI, rather than the sample display name. Duplicate join
#' keys cause an error to prevent a join from multiplying count rows.
#'
#' Only sample-specific \code{Population} elements containing a direct
#' \code{Gate} child are included. Sample totals, group gate templates, and
#' Boolean nodes represented only by \code{Dependents} are not returned.
#' Samples without such gates contribute no rows.
#'
#' Rectangle bounds and polygon vertices are extracted in their stored
#' coordinate system. Compensation, transformations, and event gating are
#' not performed. Other gate types retain their XML definition, but do not
#' necessarily have a complete parsed geometric representation.
#'
#' @return A data frame with the following columns. If no gates are found,
#'   a zero-row data frame with the same columns is returned.
#'   \describe{
#'     \item{FileName}{FCS basename from the sample's decoded DataSet URI.}
#'     \item{PopulationFullPath}{Slash-separated population hierarchy, excluding
#'       the sample name. Includes Boolean ancestors where present.}
#'     \item{sample_id}{Sample identifier stored in the workspace, as character.}
#'     \item{sample}{Sample display name from the SampleNode element.}
#'     \item{gate}{Population name.}
#'     \item{gate_path}{Alias of PopulationFullPath, retained for compatibility.}
#'     \item{parent_gate}{Name of the nearest ancestor Population, or NA.}
#'     \item{gate_id}{Stored gate identifier, or NA.}
#'     \item{parent_id}{Stored parent gate identifier, or NA.}
#'     \item{gate_type}{Gate element type, such as RectangleGate or PolygonGate.}
#'     \item{events_inside}{Stored eventsInside attribute, as character, or NA.}
#'     \item{definition_xml}{Serialized Gate element including its attributes
#'       and children. This is not a byte-preserving copy of the source text.}
#'     \item{dimensions}{List-column of data frames with character parameter
#'       names and numeric min and max bounds. Unspecified bounds are NA.}
#'     \item{vertices}{List-column of numeric matrices, with vertices in rows
#'       and parameter names as column names; NULL if there are no vertices.}
#'   }
#'
#' @examples
#' \dontrun{
#' ws <- "path/to/workspace.wsp"
#' gates <- extract_gates(ws)
#' gates$dimensions[[1]]
#' gates$vertices[[1]]
#'
#' # A parsed workspace can also be supplied.
#' gates <- extract_gates(xml2::read_xml(ws))
#'
#' counts <- fcexpr::wsx_get_popstats(ws)$counts
#' counts_with_gates <- dplyr::left_join(
#'   counts, gates,
#'   by = c("FileName", "PopulationFullPath"),
#'   relationship = "many-to-one"
#' )
#' }
#' @export
extract_gates <- function(ws) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Install xml2 first: install.packages('xml2')", call. = FALSE)
  }
  if (inherits(ws, "xml_document")) {
    doc <- ws
  } else if (is.character(ws) && length(ws) == 1L && !is.na(ws) &&
             file.exists(path.expand(ws)) && !dir.exists(path.expand(ws))) {
    doc <- xml2::read_xml(path.expand(ws), options = c("NONET", "NOBLANKS"))
  } else {
    stop("ws must be one existing workspace path or an xml2 XML document.",
         call. = FALSE)
  }
  ns <- xml2::xml_ns(doc)
  samples <- xml2::xml_find_all(doc, ".//Sample")
  rows <- list()
  for (s in samples) {
    sn <- xml2::xml_find_first(s, "./SampleNode")
    dataset <- xml2::xml_find_first(s, "./DataSet")
    filename <- basename(utils::URLdecode(sub("^file:", "", xml2::xml_attr(dataset, "uri"))))
    for (p in xml2::xml_find_all(sn, ".//Population[Gate]")) {
      gate <- xml2::xml_find_first(p, "./Gate")
      shape <- xml2::xml_child(gate)
      dims <- xml2::xml_find_all(shape, "./gating:dimension", ns)
      parameters <- xml2::xml_attr(
        xml2::xml_find_first(dims, "./data-type:fcs-dimension", ns),
        "data-type:name", ns
      )
      bounds <- data.frame(
        parameter = parameters,
        min = as.numeric(xml2::xml_attr(dims, "gating:min", ns)),
        max = as.numeric(xml2::xml_attr(dims, "gating:max", ns))
      )
      vertices <- xml2::xml_find_all(shape, "./gating:vertex", ns)
      coordinates <- lapply(vertices, function(v) {
        as.numeric(xml2::xml_attr(xml2::xml_children(v), "data-type:value", ns))
      })
      vertex_matrix <- if (length(coordinates)) {
        m <- do.call(rbind, coordinates)
        colnames(m) <- parameters
        m
      } else NULL
      ancestors <- xml2::xml_find_all(p, "ancestor-or-self::*[self::Population or self::AndNode or self::OrNode or self::NotNode]")
      parent <- xml2::xml_find_first(p, "ancestor::Population[1]")
      row <- data.frame(
        FileName = filename,
        PopulationFullPath = paste(xml2::xml_attr(ancestors, "name"), collapse = "/"),
        sample_id = xml2::xml_attr(dataset, "sampleID"),
        sample = xml2::xml_attr(sn, "name"),
        gate = xml2::xml_attr(p, "name"),
        gate_path = paste(xml2::xml_attr(ancestors, "name"), collapse = "/"),
        parent_gate = xml2::xml_attr(parent, "name"),
        gate_id = xml2::xml_attr(gate, "gating:id", ns),
        parent_id = xml2::xml_attr(gate, "gating:parent_id", ns),
        gate_type = xml2::xml_name(shape),
        events_inside = xml2::xml_attr(shape, "eventsInside"),
        definition_xml = as.character(gate),
        stringsAsFactors = FALSE
      )
      row$dimensions <- list(bounds)
      row$vertices <- list(vertex_matrix)
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (!length(rows)) {
    gates <- as.data.frame(stats::setNames(
      rep(list(character()), 12L),
      c("FileName", "PopulationFullPath", "sample_id", "sample", "gate",
        "gate_path", "parent_gate", "gate_id", "parent_id", "gate_type",
        "events_inside", "definition_xml")
    ), stringsAsFactors = FALSE)
    gates$dimensions <- vector("list", 0L)
    gates$vertices <- vector("list", 0L)
    return(gates)
  }
  gates <- do.call(rbind, rows)

  # Prevent joins from silently multiplying count rows.
  if (anyDuplicated(gates[c("FileName", "PopulationFullPath")])) {
    stop("Duplicate FileName + PopulationFullPath keys in this workspace.")
  }

  gates
}
