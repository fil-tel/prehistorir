## ---------------------------
##
##
## Purpose of script: This script contains functions used to obtain gene flows
##
## Author: Filippo Tell
##
## Date Created: 2026-01-29
##
##
## ---------------------------
##
## Notes: many will be removed/modify
##
##
## ---------------------------


######################################
######################################
# These functions might be useful ####
######################################
######################################

# function to check whether two populations are alive at the sam time
check_if_alive <- function(model, pop1, pop2){
  # extract time pop1 and pop2
  pop1_t <- model %>% dplyr::filter(pop==pop1) %>% dplyr::select(t_start, t_end)
  pop2_t <- model %>% dplyr::filter(pop==pop2) %>% dplyr::select(t_start, t_end)
  # set to 0 if t_end is NA, i.e. still alive
  # pop1_t[is.na(pop1_t)] <- 0
  # pop2_t[is.na(pop2_t)] <- 0
  if((pop1_t["t_end"]-pop2_t["t_start"]<0) & (pop2_t["t_end"]-pop1_t["t_start"]<0)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# function to compute distance between populations
compute_dist <- function(model, pop1, pop2){
  # extract coordi pop1 and pop2
  pop1_coord <- model %>% dplyr::filter(pop==pop1) %>% dplyr::select(x, y)
  pop2_coord <- model %>% dplyr::filter(pop==pop2) %>% dplyr::select(x, y)
  # return euclidean distance
  return(norm(as.matrix(pop1_coord - pop2_coord), "F"))
}

# sigmoid function
sigmoid <- function(x, a=1){
  1/(1+exp(x = (-a*x)))
}
######################################
######################################
######################################
######################################


##################################################
####################EXTRA#########################
##################################################

# function to get the starting time of a gene flow
# between two populations
get_t_start <- function(model, pop1, pop2){
  # we need to add a small offset
  # because for slendr a geneflow cannot start at the same time
  # as the splitting of the population
  offset <- 1
  t_start1 <- model %>% dplyr::filter(pop==pop1) %>% .[["t_start"]] - offset
  t_start2 <- model %>% dplyr::filter(pop==pop2) %>% .[["t_start"]] - offset
  # get the min between the two cause bigger means before
  min(t_start1, t_start2)
}

# function to get the ending time of
# a gene flow
get_t_end <- function(model, pop1, pop2){
  # if both the populations are alive then the gene flow terminates at time 0
  if(all(model %>% dplyr::filter(pop==pop1) %>% .[["alive"]], model %>% dplyr::filter(pop==pop2) %>% .[["alive"]])){
    t_end1 <- t_end2 <- 0
  }else{
    # note that the offset is not needed for the
    # end of the gene flow
    t_end1 <- model %>% dplyr::filter(pop==pop1) %>% .[["t_end"]]
    t_end2 <- model %>% dplyr::filter(pop==pop2) %>% .[["t_end"]]
  }
  max(t_end1, t_end2)
}

##################################################
##################################################
##################################################
##################################################


# function that given a certain distance_df
# created by create_distances_df returns a set of geneflows
# given the location of the populations
create_gf_df <- function(dist_df, t_cutoff, gf_fun, ...) {
  if(!inherits(dist_df, "prehistorik_dist_df")) stop("dist_df do not belong to the class prehistorik_dist_df\n Are you sure you use the function create_distances_df to generate it?")
  model <- get_model(dist_df)
  generation_time <- get_generation_time(model)
  # a gene flow cannot be shorter than the generation time
  # so we need to check that the time cutoff, that filter geneflows shorter than t_cutoff
  # is greater than the generation time
  if(t_cutoff<generation_time) stop("The t_cutoff cannot be smaller than the generation time.")
  if(!is.function(gf_fun)) stop("gf_fun has to be a function")
  # add gf rate column
  dist_df$gf_rate <- gf_fun(dist_df$dist_m, ...)
  # now we add the t_start and t_end of the gf
  # it needs to happen in the overlapping period of the two populations
  # that we already know are existing simultaneously
  # add t_start and t_end for the pops
  # by this means that the model has only
  pop1_ind <- match(dist_df$pop1, model$pop)
  dist_df$t_start_1 <- model$t_start[pop1_ind]
  dist_df$t_end_1 <- model$t_end[pop1_ind]
  pop2_ind <- match(dist_df$pop2, model$pop)
  dist_df$t_start_2 <- model$t_start[pop2_ind]
  dist_df$t_end_2 <- model$t_end[pop2_ind]
  # first define offset
  # because for slendr a geneflow cannot start at the same time
  # as the splitting of the population
  offset <- 1
  gf_df <- dist_df %>% mutate(
    gf_start = ifelse(t_start_1-t_start_2<=0, t_start_1-offset, t_start_2-offset),
    gf_end = ifelse(t_end_1-t_end_2<=0, t_end_2, t_end_1)
  )
  gf_df <- gf_df[(gf_df$gf_start-gf_df$gf_end)>=t_cutoff,]
  attr(gf_df, "model") <- model
  attr(gf_df, "t_cutoff") <- t_cutoff
  attr(gf_df, "gf_fun") <- gf_fun
  class(gf_df) <- set_class(gf_df, "geneflows_df")
  gf_df
}


# function to get the dist df at a specific time point
get_dist_df_at_t <- function(model, time, threshold, resolution = NULL, as_df = TRUE){
  if(!inherits(model, "prehistorik_model")) stop("model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?")
  if(is.null(resolution)) resolution <- get_resolution(model)
  alive_df <- get_alive_pop(model, time)
  # get distance in meters
  dist_vec <- dist(alive_df[, c("x", "y")])*resolution
  comb_mat <- RcppAlgos::comboGeneral(1:nrow(alive_df), 2)
  gf_vec <- ifelse(dist_vec<=threshold, dist_vec, NA)
  # populations are ordered so I can use the second column
  # and named the combinations after the population names
  gf_mat <- matrix(c(alive_df[,2][comb_mat], gf_vec, rep(time, length(gf_vec))) , ncol = 4)
  gf_df <- gf_mat[!is.na(gf_vec), ]
  if(as_df){
    gf_df <- as.data.frame(gf_df)
    names(gf_df) <- c("pop1", "pop2", "dist_m", "time")
    attr(gf_df, "model") <- model
    attr(gf_df, "threshold") <- threshold
    class(gf_df) <- set_class(gf_df, "dist_df")
  }
  gf_df
}


# functions to create a df of distances
# between tribes living at the same time
# threshold gives the threshold in meters after the one
# gene flows are not possible
create_dist_df <- function(model, threshold) {
  if(!inherits(model, "prehistorik_model")) stop("model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?")
  # check that the pops have coordinates
  if (any(is.na(model$x)) ||
      any(is.na(model$y)))
    stop("The model is not spatial.")
  # extract resolution of the model to get distances in meters
  resolution <- get_resolution(model)
  # check whether the threshold is meaningful
  if (threshold<resolution) stop("The threshold is less than the resolution.")
  # at each generation extract population that are alive
  param <- get_parameters(model)
  ncores <- parallel::detectCores()-1
  time_points <- seq(param$t_start, param$t_stop, -param$generation_time)
  gf_df <- do.call(
    rbind,
    parallel::mclapply(time_points, function(t)
      get_dist_df_at_t(
        model,
        time = t,
        threshold = threshold,
        resolution = resolution,
        as_df = FALSE
      )
      , mc.cores = ncores)
  )
  # keep only uniques
  gf_df <- unique(gf_df)
  gf_df <- as.data.frame(gf_df)
  names(gf_df) <- c("pop1", "pop2", "dist_m", "time")
  attr(gf_df, "model") <- model
  attr(gf_df, "threshold") <- threshold
  class(gf_df) <- set_class(gf_df, "dist_df")
  gf_df
}

# function to extract pops alive at time t
# extract pops that are alive at time t
get_alive_pop <- function(model, time){
  if(!inherits(model, "prehistorik_model")) stop("model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?")
  params <- get_parameters(model)
  generations <- seq(params$t_start, params$t_stop, -params$generation_time)
  # extract pops that are alive at time t
  if (time==0) {
    model[model$alive, ]
  }else{
    model[model$t_end<time & time<=model$t_start, ]
  }
}

plot_gf_as_graph <- function(model, time, threshold, lat_range=NULL, lon_range=NULL){
  if(!inherits(model, "prehistorik_model")) stop("Model do not belong to the class prehistorik_model\n Are you sure you use the function create_distances_df to generate it?")
  # remove African populations
  my_model_f <- model[!model$pop==0,]
  dist_df <- get_dist_df_at_t(my_model_f, time = time, threshold = 50e3)
  # get pop in this time period
  my_pop_df <- get_alive_pop(model = my_model_f, time)
  edges <- dist_df[, c("pop1", "pop2")]
  # switch first two column so to have first the pop
  # otherwise this creates problems
  my_pop_df <- my_pop_df[,c(2,1, 3:ncol(my_pop_df))]
  # the conditions in this way it is just because I do not wnat to modify
  # other functions
  if(is.null(lat_range) && is.null(lon_range)){
    g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = my_pop_df)
  }else if(!(is.null(lat_range) && is.null(lat_range))){
    map_mat <- get_map_mat(model)
    xy_min <- get_xy_from_latlon(lat = lat_range[1], lon = lon_range[1], map_mat)
    xy_max <- get_xy_from_latlon(lat = lat_range[2], lon = lon_range[2], map_mat)
    # choose only pop tha live in that area
    my_pop_df <- my_pop_df[between(my_pop_df$x, xy_min$x, xy_max$x) & between(my_pop_df$y, xy_min$y, xy_max$y),]
    # clean edges df
    edges <- edges[edges$pop1 %in% my_pop_df$pop & edges$pop2 %in% my_pop_df$pop, ]
    g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = my_pop_df)
  }else stop("Either lat-lon are NULL, or both have to be present.")
  location <- matrix(c(igraph::V(g)$x, igraph::V(g)$y), ncol = 2)
  igraph::V(g)$label.cex <- 0.5
  igraph::plot.igraph(g, vertex.size = 80,
       layout = location, rescale=FALSE,
       xlim = range(V(g)$x), ylim = range(V(g)$y))
}



