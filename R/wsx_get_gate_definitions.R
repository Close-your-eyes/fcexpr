#' Obtain geometric gate definitions (coordinates) from flowjo workspace
#'
#' @param ws
#' @param groups
#' @param invert_groups
#'
#' @return
#' @export
#'
#' @examples
wsx_get_gate_definitions <- function(ws,
                                     groups = NULL,
                                     invert_groups = F) {

  ws <- check_ws(ws)

  group_df <- wsx_get_groups(ws, collapse_groups = F)

  if (is.null(groups)) {
    groups <- unique(group_df[,"FlowJoGroup", drop=T])
  }

  if (invert_groups) {
    group_df <- group_df[which(!group_df[,"FlowJoGroup", drop=T] %in% groups),]
  } else {
    group_df <- group_df[which(group_df[,"FlowJoGroup", drop=T] %in% groups),]
  }

  if (nrow(group_df) == 0) {
    stop("Non of provided groups found.")
  }

  # each sample which may be in multiple groups only considered once
  ids <- unique(group_df[,"sampleID",drop=T])

  rel_nodes <- xml2::xml_children(xml2::xml_child(ws, "SampleList"))
  rel_nodes <- rel_nodes[which(purrr::map(rel_nodes, function(x) xml2::xml_attrs(xml2::xml_child(x, "DataSet"))[["sampleID"]]) %in% ids)]
  gg <- xml2::xml_find_all(rel_nodes, ".//Gate|.//Dependents")

  gates <- lapply(gg, function(n) xml2::as_list(n))
  names(gates) <- lapply(gates, function(x) attributes(x)[["id"]])

  gates_df <- stack(lapply(gates, names))
  gates_df$ind <- as.character(gates_df$ind)
  names(gates_df) <- c("gate_type", "gate_id")
  gates_df$parentgate_id

  # attach parent gate ids
  parentgates_df <- sapply(gates, function(x) attributes(x)[["parent_id"]])
  parentgates_df <- stack(parentgates_df[which(!sapply(parentgates_df, is.null))])
  names(parentgates_df) <- c("parentgate_id", "gate_id")
  gates_df <- dplyr::left_join(gates_df, parentgates_df, by = "gate_id")

  if (any(duplicated(gates_df$gate_id))) {
    warning("Duplicate gate found. Strange. Check!")
    gates_df <- gates_df[which(!duplicated(gates_df$gate_id)),]
    gates <- gates[which(!duplicated(names(gates)))]
  }


  # attach geometric gate definitions from flowjo
  gates_df$gate_def <- lapply(names(gates), function(x) {

    ## add routines for other gate types as well
    ## check how to save gate defs elegantly to have them easily passed to subsequent methods
    ## save numbers as numeric
    ## add channel desc if possible?

    if (names(gates[[x]]) == "RectangleGate") {
      list(dim1_name = attributes(gates[[x]]$RectangleGate[[1]]$`fcs-dimension`)$name,
           dim1_min = as.numeric(attributes(gates[[x]]$RectangleGate[[1]])$min),
           dim1_max = as.numeric(attributes(gates[[x]]$RectangleGate[[1]])$max),
           dim2_name = attributes(gates[[x]]$RectangleGate[[2]]$`fcs-dimension`)$name,
           dim2_min = as.numeric(attributes(gates[[x]]$RectangleGate[[2]])$min),
           dim2_max = as.numeric(attributes(gates[[x]]$RectangleGate[[2]])$max))
    } else if (names(gates[[x]]) == "PolygonGate") {

      inds <- 3:(sum(names(gates[[x]]$PolygonGate) == "vertex")+2)
      #dim1_vertices <- sapply(inds, function(y) attributes(gates[[x]]$PolygonGate[[y]][[1]])$value)
      #dim2_vertices <- sapply(inds, function(y) attributes(gates[[x]]$PolygonGate[[y]][[2]])$value)
      list(dim1_name = attributes(gates[[x]]$PolygonGate[[1]]$`fcs-dimension`)$name,
           dim1_vertices = as.numeric(sapply(inds, function(y) attributes(gates[[x]]$PolygonGate[[y]][[1]])$value)),
           dim2_name = attributes(gates[[x]]$PolygonGate[[2]]$`fcs-dimension`)$name,
           dim2_vertices = as.numeric(sapply(inds, function(y) attributes(gates[[x]]$PolygonGate[[y]][[2]])$value)))
    } else {
      message("Yet unknown gate type: ", names(gates[[x]]))
    }
  })

  return(list(list = gates, df = gates_df))
}




