
.model_ensure_is_tree <- function(d) {
  parent_ids <- unique(d$parent)
  child_ids <- unique(d$child)
  
  # make sure that the graph has exactly one root
  root = subset(d, !(parent %in% child))
  root_ids <- unique(root$parent)
  n_roots = length(root_ids)
  if(n_roots > 1)
    stop(sprintf("Model has %d root nodes: <%s>.", n_roots, paste(root_ids, collapse=",")))
}

.node_root_id <- function(d) {
  root = subset(d, !(parent %in% child))
  root_id = unique(root$parent)
  stopifnot(length(root_id) == 1)
  root_id
}

.node_is_terminal <- function(d, node_id) {
  stopifnot(node_id %in% d$parent || node_id %in% d$child)
  cur_d = subset(d, parent == node_id)
  if(nrow(cur_d) == 0)
    return(TRUE)
  else
    return(FALSE)
}


.harvest_terminals <- function(x) {
  if(typeof(x) == "character")
    x <- parse(text=x)
  
  if(typeof(x) == "expression")
    x <- x[[1]]
  
  if(length(x) == 1) {
    return(as.character(x))
    
  }
  if(typeof(x[[1]]) == "symbol")
    x[[1]] <- NULL
  
  
  unlist( lapply(x, .harvest_terminals) )
}
