#Nonredundant Helper Functions
#Brian Chivers, DARC team

#' Nonredundant
#'
#' \code{nonredundant} flags for redundant dummy levels
#'
#' @param iid A vector of group dummy indicators
#' @param tid A vector of time dummy indicators
#' @param w A vector of non-negative weights
#' @return A list will be returned with the following named values:
#'         flag - Are there redundant dummy levels?
#'         nr - a listing of 
#' @export

nonredundant <- function(iid, tid, w) {
  obs <- length(iid)
  iid <- iid[w != 0]
  tid <- tid[w != 0]
  
  for (k in 1:max(length(iid), length(tid))) {
    kid <- where_id_with_single_obs(iid)
    iid <- iid[!kid]
    tid <- tid[!kid]
    kid <- where_id_with_single_obs(tid)
    
    if(!any(kid)) {
      break
    }
    iid <- iid[!kid]
    tid <- tid[!kid]
  }
  return_list <- list()
  flag <- obs != length(iid)
  nr <- list()
  nr$iid <- unique(iid)
  nr$tid <- unique(tid)
  
  return_list$flag <- flag
  return_list$nr <- nr
  return(return_list)
}

#' \code{where_id_with_single_obs} returns a true/false vector, 
#' Does this dummy occur multiple times?
#'
#' @param id A vector of dummy indicators
#' @return A true/false vector
#' @export
where_id_with_single_obs <- function(id) {
  return(!(id %in% ids_with_multiple_obs(id)))
}

#' \code{ids_with_multiple_obs} returns a vector of values that occur multiple times 
#'
#' @param id A vector of dummy indicators
#' @return A vector of repeated values
#' @export
ids_with_multiple_obs <- function(id) {
  return(unique(id[duplicated(id)]))
}
