
.generate_code_data <- function(data_vars, n_vars, data_length) {
  lines <- sprintf("%s %s[%s]", data_vars, names(data_vars), data_length)
  lines <- c(sprintf('int<lower=1> %s', n_vars), 
             lines)
  .format_block('data', lines)
}

.generate_code_par <- function(par_vars, ranef_vars) {
  ranef_lines <- lapply(names(ranef_vars), function(parameter_name) {
    ranef_vector <- sprintf("vector[n_%s] %s_%s", ranef_vars[[parameter_name]], parameter_name, ranef_vars[[parameter_name]])
    ranef_sd <- sprintf("real<lower=0> %s_%s_sd", parameter_name, ranef_vars[[parameter_name]])
    c(ranef_vector, ranef_sd)
  })
  ranef_lines <- unlist(ranef_lines)
  fixedeff_lines <- sprintf("%s %s", par_vars, names(par_vars))
  lines <- c(fixedeff_lines, ranef_lines)
  .format_block('parameters', lines)
}


.generate_code_likHyperpar <- function(ranef_vars) {
  res <- lapply(names(ranef_vars), function(parameter_name) {
    sprintf("%s_%s ~ normal(0, %s_%s_sd)", parameter_name, ranef_vars[[parameter_name]], parameter_name, ranef_vars[[parameter_name]])
  })
  unlist(res)
}


.model_harvest_code <- function(d, root_id=NULL) {
  if(is.null(root_id))
    root_id <- .node_root_id(d)
  
  is_tree_root = (.node_root_id(d) == root_id)
  cur_children <- subset(d, parent == root_id)
  
  # return edge info if node is terminal
  if( nrow(cur_children) == 0 ) {
    return(list(list(path=root_id, prob=NULL, code=NULL, condition=NULL)))
  }
  
  # append edge info unless node is the tree root
  res <- list()
  for(i in 1:nrow(cur_children)) {
    cur_child <- cur_children[i,]
    cur_res <- .model_harvest_code(d, root_id=cur_child$child)
    for(j in 1:length(cur_res)) {
      cur_res[[j]] <- with(cur_res[[j]], 
                           if(cur_child$is_logical)
                              list(path=c(root_id, path), code=c(cur_child$code, code),
                                   prob=prob, condition=c(cur_child$prob, condition))
                           else
                             list(path=c(root_id, path), code=c(cur_child$code, code),
                                  prob=c(cur_child$prob, prob), condition=condition)
      )
    }
    res <- append(res, cur_res)
  }
  res
}

.model_reshape_code <- function(yield)
{
  ldply(yield, function(y) {
    path <- paste(y$path, collapse='-')
    condition <- paste(y$condition, collapse=' && ')
    prob <- setdiff(y$prob, "1") # y$prob
    prob <- paste(sprintf("(%s)", prob), collapse='*')
    y$code <- y$code[!(y$code %in% c('', " "))]
    code <- paste(y$code, collapse=';')
    c(id=path, prob=prob, condition=condition, code=code)
  })
}

.format_lines <- function(lines, indentation=0) {
  indentation <- paste0(rep('\t', indentation), collapse="")
  lines <- gsub('\n', paste('\n', indentation, sep=""), lines)
  lines <- paste(indentation, lines, ifelse(grepl(";$", lines), "", ";"), sep='')
}

.format_block <- function(section, lines, indentation=0) {
  lines <- .format_lines(lines, indentation+1)
  code_block <- paste(lines, collapse="\n")
  paste('',section, '{', code_block, '}', sep='\n')
}

.append_to_symbol <- function(lst, symbols, suffix) {
  .circumfix_symbol(lst=lst, symbols=symbols, prefix='', suffix=suffix)
}

.circumfix_symbol <- function(lst, symbols, prefix, suffix)
{
  if(typeof(lst) == "character")
  {
    ret <- lapply(lst, function(cur_lst) as.list(parse(text=cur_lst)) %>% 
                    .circumfix_symbol(symbols=symbols, prefix=prefix, suffix=suffix) %>% 
                    paste(., collapse="; ") )
    return(unlist(ret))
  }
  
  if(is.null(lst) || length(lst) == 0)
    return(lst)
  
  if(typeof(lst) == "symbol")
  {
    if(lst == "=")
      return( as.symbol('=') )
    
    if(as.character(lst) %in% symbols) {
      res <- parse(text=paste0(prefix,lst,suffix))[[1]]
      return( res )  
    }
    return(lst)
  }
  
  if( is.numeric(lst) ) {
    return(lst)
  }
  
  for(i in 1:length(lst)) {
    lst[[i]] <- .circumfix_symbol(lst[[i]], symbols, prefix, suffix)
  }
  lst
}


.check_variables <- function(branch_conditions, branch_cprob, branch_code, par_vars, iv_vars, dv_vars, logLik, ranef, ignore_vars = NULL)
{
  lhs_expr = unique( sapply(unlist(branch_code), function(eq) eq[[2]]) )
  lhs_vars = as.character(lhs_expr)
  
  all_symbols = unique(.serialize_expression( unlist(c(branch_code,branch_cprob,branch_conditions)) ))
  all_vars =  setdiff(all_symbols, all_operators) %>% setdiff(ignore_vars)
#  all_vars = c(names(par_vars), names(iv_vars), names(dv_vars)) 
#  if(!all(all_expr_vars %in%  all_vars)) {
#    stop(sprintf("Type of variables in the model was not specified for: %s.", paste(all_expr_vars[!(all_expr_vars %in% all_vars)], collapse=", ")))
#  }
  
  rhs_vars = setdiff(all_vars, lhs_vars)
  
  if( !setequal(c(names(par_vars),names(iv_vars),names(dv_vars)), rhs_vars) ) {
    extra_in_model <- setdiff( rhs_vars, c(names(par_vars),names(iv_vars),names(dv_vars)) )
    extra_in_par_vars <- setdiff( names(par_vars), rhs_vars )
    extra_in_iv_vars <- setdiff( names(iv_vars), rhs_vars )
    extra_in_dv_vars <- setdiff( names(dv_vars), rhs_vars )
    
    if(length(extra_in_model) > 0) {
      stop(sprintf("Type of variable in the model was not specified for: %s.", paste(extra_in_model, collapse=',') ))
    } 
    #else if(length(extra_in_iv_vars) > 0) {
    #  stop(sprintf("Variables in iv_vars do not exist in the model: %s.", paste(extra_in_iv_vars, collapse=',') ))
    #}
    #
    #else if(length(extra_in_par_vars) > 0) {
    #  stop(sprintf("Variables in par_vars do not exist in the model: %s.", paste(extra_in_par_vars, collapse=',') ))
    #}
    #
    #else if(length(extra_in_dv_vars) > 0) {
    #  stop(sprintf("Variables in dv_vars do not exist in the model: %s.", paste(extra_in_dv_vars, collapse=',') ))
    #}
  }
  
  all_vars = unique(c(all_vars, names(par_vars), names(iv_vars), names(dv_vars)))
#print(all_vars)

  # find out what the random effects level variables are, and add them to the rest
  ran_eff_vars = lapply(ranef, function(expr) { 
    #print( setdiff(.harvest_terminals(expr), all_vars) )
    setdiff(.harvest_terminals(expr), all_vars)
  })
#print(ran_eff_vars)  

  list(all=all_vars, lhs=lhs_vars, rhs=rhs_vars, ran_eff=ran_eff_vars,
       par=names(par_vars), data=c(names(iv_vars),names(dv_vars)) )
}


.modify_symbol <- function(lst, symbols_from, symbols_to) {
  if(typeof(lst) == "character") {
    ret <- lapply(lst, function(cur_lst) as.list(parse(text=cur_lst)) %>% 
                    .modify_symbol(symbols_from=symbols_from, symbols_to=symbols_to) %>% 
                    paste(., collapse="; ") )
    return(unlist(ret))
  }
  
  if(is.null(lst) || length(lst) == 0)
    return(lst)
  
  if(typeof(lst) == "symbol") {
    if(lst == "=")
      return( as.symbol('=') )
    if(as.character(lst) %in% symbols_from) {
      new_symbol <- symbols_to[which(as.character(lst) == symbols_from)]
      res <- parse(text=new_symbol)[[1]]
      return( res )  
    }
    return(lst)
  }
  
  if( is.numeric(lst) ) {
    return(lst)
  }
  
  for(i in 1:length(lst)) {
    lst[[i]] <- .modify_symbol(lst[[i]],  symbols_from, symbols_to)
  }
  lst
}

.serialize_expression <- function(lst) {
  if(typeof(lst) == "symbol") {
    return(as.character(lst))
  }
  if(is.numeric(lst)) {
    return(NULL)
  }
  res <- c()
  for(i in 1:length(lst))
    res <- c(res, .serialize_expression(lst[[i]]))
  res
}
