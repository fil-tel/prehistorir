## ---------------------------
##
##
## Purpose of script: This script contains functions
## that creates/help slendr functions to work in our setting of
## data structures/models
##
## Author: Filippo Tell
##
## Date Created: 2026-01-29
##
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

suppressWarnings({
  library(magrittr)
  library(sf)
  library(stars)
  library(dplyr)
  library(slendr)
  library(ggplot2)
})

# function to create the populations schedule for the slendr_model
# given a dataframe as returned by the function get_samples_from_location
get_schedule_sampling_from_df <- function(model, samples_df, ...) {
  # extract populations list from model
  populations <- model$populations
  schedules <- lapply(
    1:nrow(samples_df),
    FUN = function(row_id) {
      row <- samples_df[row_id,]
      schedule_sampling(model = model,
                        times = row$times,
                        list(populations[[paste0("POP", row$sampling_pop)]], row$n),
                        ...)
    }
  )
  # bind as a final df
  do.call(rbind, schedules)
}

# THIS FUNCTION MIGHT BE PROBLEMATIC
# function to sample individuals from tribes in a specific
# location of the map at a specific time
get_samples_from_location <- function(model, n, time, lat, lon){
  if(!inherits(model, "prehistorik_model"))
    stop(
      "model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?"
    )
  # first extract the layer corresponding to the time
  # closest to the sampling time
  map_mat <- get_map_mat(model)
  pop_df <- get_alive_pop(model, time = time)
  resolution <- get_resolution(model)
  # convert the location from latlon to xy
  xy_tbl <- get_xy_from_latlon(lat = lat, lon = lon, map_mat = map_mat)
  # now sample the closest exisiting tribe in space to that location
  # (if two or more are exisitng at the same distance just pick one)
  # place on top the point
  dist_m <- dist(rbind(xy_tbl, pop_df[, c("x", "y")]))
  # remove rows corresponding to
  # extract the pop closest
  pop <- pop_df[which.min(dist_m[1:nrow(pop_df)]),]
  # might have to do it better, for example a cutoff might be useful otherwise it could pick a population that is super far just because there is nothing nearby
  # or select a neighbourhood
  cat(sprintf("The selected population is %s meters distant from the given sampling point.\n", min(dist_m[1:nrow(pop_df)])*resolution), "\n")
  # get lon and lat from the exact sampled cell
  latlon <- as.numeric(get_latlon_from_xy(x = pop$x, y = pop$y, map_mat = map_mat))
  samples_df <- data.frame(sampling_pop = as.character(pop$pop),  times = time, n = n, x=pop$x, y=pop$y, lat = latlon[1], lon = latlon[2])
  samples_df
}


# helper function that given the df indicating the gf and the list of pop creates a gf object in slendr
create_gf_list <- function(gf_df, pop_list) {
  # now I create a list of slendr gf
  # do invert pops so to ave symmetric geneflows
  n_cores <- parallel::detectCores() - 5
  gf_list <- parallel::mclapply(
    1:nrow(gf_df),
    FUN = function(i) {
      row <- gf_df[i, ]
      list(
        gene_flow(
          from = pop_list[[as.character(row["pop1"])]],
          to = pop_list[[as.character(row["pop2"])]],
          proportion = as.numeric(row["gf_rate"]),
          start = as.numeric(row["gf_start"]),
          end = as.numeric(row["gf_end"]),
          overlap = FALSE
        ),
        gene_flow(
          from = pop_list[[as.character(row["pop2"])]],
          to = pop_list[[as.character(row["pop1"])]],
          proportion = as.numeric(row["gf_rate"]),
          start = as.numeric(row["gf_start"]),
          end = as.numeric(row["gf_end"]),
          overlap = FALSE
        )
      )
    },
    mc.cores = n_cores
  )
  return(unlist(gf_list, recursive = FALSE))
}

# function to create a list of slendr populations
# given a data frame
# that is parento of all the others, this makes sense in case there are more than one starting populations
create_poplist <- function(model, ...){
  if(!inherits(model, "prehistorik_model")) stop("model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?")
  pop_list <- vector(mode = "list", length = nrow(model))
  names(pop_list) <- as.character(model$pop)
  # iterate over the row and create each population
  for (row_id in 1:nrow(model)) {
    # extract the population details
    row <- model[row_id,]
    # if row$parent==row$pop this is the first population
    if(row$parent==row$pop){
      # check if the pop is alive
      # do this way cause ifelse gives problems
      if(row$alive){
        t_end <- NULL
      }else{
        t_end <- row$t_end
      }
      pop_tmp <- population(name = paste0("POP", row$pop), time = row$t_start, remove = t_end, N = row$N, ...)
    }else{
      # check if the pop is alive
      # do this way cause ifelse gives problems
      if(row$alive){
        t_end <- NULL
      }else{
        t_end <- row$t_end
      }
      pop_tmp <- population(parent = pop_list[[as.character(row$parent)]], name = paste0("POP", row$pop), time = row$t_start, remove = t_end, N = row$N, ...)
    }
    pop_list[[as.character(row$pop)]] <- pop_tmp
  }
  return(pop_list)
}
