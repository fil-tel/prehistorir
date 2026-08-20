## ---------------------------
##
##
## Purpose of script: This script contains all the functions
## used to generate and work on map_mat object
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

# Function that given a map obtain a raster of it
# and a matrix of the location as well
# it returns it as a map_mat object
# set_free if TRUE set whatever is not NA to -1, hence free habitable cells
# in the map
#' Convert an \code{sf} object into a matrix
#'
#' Given an \code{sf} object, it rasterizes it and converts it into a \code{prehistorik_map}.
#'
#' @param map Object of the class \code{sf}.
#' @param resolution Resolution in meters of a side of a matrix cell.
#' @param attr Which attribute of the \code{sf} object should be kept to define the value of a matrix cell?
#' @param set_free Should not NA values be set to -1, hence be defined as free habitable cells of the map?
#'
#' @returns An object of the class \code{prehistorik_map}, corresponding to the a raster version of the map.
#' @export
#'
get_map_as_matrix <- function(map, resolution, attr, set_free=TRUE){
  # get bounding box to convert into template
  bbox <- sf::st_bbox(map)
  template <- stars::st_as_stars(bbox, dx = resolution, dy = resolution, values = NA_real_)
  # extract raster as stars
  raster_st <- stars::st_rasterize(map, template)
  # convert it to sf object but not remove NA values
  # those can always be removed afterwarss if needed
  map_mat <- raster_st[[attr]]
  # flip the columns to get the right orientation
  map_mat <- map_mat[, ncol(map_mat):1]
  # if flag on set -1 whatever is not NA
  # otherwise keep the value of the attrribute
  if(set_free)   map_mat[!is.na(map_mat)] <- -1
  colnames(map_mat) <- 1:ncol(map_mat)
  rownames(map_mat) <- 1:nrow(map_mat)
  attr(map_mat, "resolution") <- resolution
  attr(map_mat, "crs") <- sf::st_crs(map)
  attr(map_mat, "bbox") <- sf::st_bbox(map)
  class(map_mat) <- set_class(map_mat, "map")
  return(map_mat)
}

#' Plot prehistorik map
#'
#' Plot a prehistorik map
#'
#' @param map_mat Object of the class \code{prehistorik_map} to be plotted.
#' @param plotly Should the plot be interactive?
#'
#' @returns Plot.
#' @export
#'
plot_map_matrix <- function(map_mat, plotly=FALSE){
  if(!inherits(map_mat, "prehistorik_map")) stop("map_mat do not belong to the class prehistorik_map\n Are you sure you use the function get_map_as_matrix to generate it?")
  # melt matrix into df
  map_df <- reshape2::melt(map_mat, c("x", "y"), value.name = "value")
  # x are the columns (j) and y the rows(i)
  p <- ggplot2::ggplot(data=map_df,aes(x=x,y=y,fill=value))+
    ggplot2::geom_tile()
  if(plotly) p <- plotly::ggplotly(p = p)
  p
}

#' Extract a temporal slice from a prehistorik model
#'
#' Obtain a temporal slice of the map of a prehistorik model.
#'
#' @param model Object of the class \code{prehistorik_model}.
#' @param time Time point at which to extract the map.
#'
#' @returns An object of the class \code{prehistorik_map}.
#' @export
#'
get_map_slice <- function(model, time){
  if(!inherits(model, "prehistorik_model")) stop("model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?")
  # extract map
  slice <- get_map_mat(model)
  # get generations
  params <- get_parameters(model)
  generations <- seq(params$t_start, params$t_stop, -params$generation_time)
  time_i <- which.min(abs(generations-time))
  time <- generations[time_i]
  # make it all -1, cause then we repopulate it
  slice[!is.na(slice)] <- -1
  # extract slice
  pops <- model %>% dplyr::filter(t_end<time, time<=t_start)
  slice[as.matrix(pops[c("x", "y")])] <- pops[["pop"]]
  # assign class
  attr(slice, "time") <- time
  slice
}


# function to extract the pops living at a certain time points
# place them on a map and plotting them
# plot_map_slice <- function(model, time, plotly = FALSE){
#   if(!inherits(model, "prehistorik_model")) stop("model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?")
#   slice <- get_map_slice(model, time)
#   p <- plot_map_matrix(slice)+ggtitle(sprintf("Years before present: %s", attr(slice, "time")))
#   if(plotly) p <- plotly::ggplotly(p)
#   p
# }
#' Plot a temporal slice of a prehistorik_model
#'
#' @param model Object of the class \code{prehistorik_model}.
#' @param time Time point at which to extract the map.
#' @param plotly Should the plot be interactive?
#'
#' @returns Plot
#' @export
#'
plot_map_slice <- function(model, time, plotly = FALSE){
  if(!inherits(model, "prehistorik_model")) stop("model do not belong to the class prehistorik_model\n Are you sure you use the function create_prehistorik_model to generate it?")
  # extract slice
  slice <- get_map_slice(model, time)
  map_df <- reshape2::melt(slice, c("x", "y"), value.name = "value")
  map_df_t <- map_df %>% filter(!is.na(value))%>% mutate(bin = as.factor(ifelse(value>0, 1, 0)))
  p <- ggplot2::ggplot(data = map_df_t, aes(x = x, y = y, fill = bin)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c("0" = "#EFEFEF", "1" = "black")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")+
    ggplot2::ggtitle(sprintf("Years before present: %s", attr(slice, "time")))
  if(plotly) p <- plotly::ggplotly(p)
  p
}


# Function to set the value of some pixels in the matrix map
#' Set entries of a \code{prehistorik_map}
#'
#' Function to modify entries of a \code{prehistorik_map}.
#'
#' @param map_mat Object of the class \code{prehistorik_map}.
#' @param value Value to insert at the given coordinates.
#' @param x x coordinate(s) to modify.
#' @param y y coordinate(s) to modify.
#' @param subset Should the (x,y) coordinates be treated as a list of points or as an area?
#' @param na If the existing value is NA, should we modify it?
#'
#' @returns \code{map} given as input, but with the entries modified.
#' @export
#'
set_map_mat_values <- function(map_mat,
                               value,
                               x,
                               y,
                               subset = FALSE,
                               na = FALSE) {
  if(!inherits(map_mat, "prehistorik_map")) stop("map_mat do not belong to the class prehistorik_map\n Are you sure you use the function get_map_as_matrix to generate it?")
  if (!na) {
    # if subset is true apply subset
    # otherwise treat them as a list of points
    if (!subset) {
      map_mat[cbind(x, y)][!is.na(map_mat[cbind(x, y)])] <- value
    } else{
      map_mat[x, y][!is.na(map_mat[x, y])] <- value
    }
  } else{
    # if subset is true apply subset
    # otherwise treat them as a list of points
    if (!subset) {
      map_mat[cbind(x, y)] <- value
    } else{
      map_mat[x, y] <- value
    }
  }
  map_mat
}



# function that given a matrix, a value, and some coordinate in the matrix
# applies the flood fill algorithm (https://en.wikipedia.org/wiki/Flood_fill)
# to populate all the connected cells with the same value
# with the given value
#
#' Flood fill algorithm
#'
#' Implementation of the flood fill algorithm to fill a section of a matrix with a given value.
#'
#' @param map_mat Object of the class \code{prehistorik_map}.
#' @param new_value Value to be set in the area.
#' @param x x coordinate at which to start the filling.
#' @param y y coordinate at which to start the filling.
#'
#' @returns Object of the class \code{prehistorik_map} where the region has been filled by the given value.
#' @export
#'
flood_fill <- function(map_mat, new_value, x, y){
  if(!inherits(map_mat, "prehistorik_map")) stop("map_mat do not belong to the class prehistorik_map\n Are you sure you use the function get_map_as_matrix to generate it?")
  # stop if we are in a NA place (no map)
  if(is.na(map_mat[x,y])) stop(sprintf("The location %s, %s is NA.", x, y))
  # if already that value return the original
  if(map_mat[x,y]==new_value) return(map_mat)
  # directions
  dir <-  list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))

  q <- list(c(x, y))
  old_value <- map_mat[x, y]
  map_mat[x, y] <- new_value

  while (length(q) > 0) {
    pos <- q[[1]]
    q <- q[-1]
    # move around
    for (dxy in dir) {
      new_x <- pos[1] + dxy[1]
      new_y <- pos[2] + dxy[2]

      if (new_y >= 1 &
          new_y <= dim(map_mat)[2] &
          new_x >= 1 &
          new_x <= dim(map_mat)[1]) {
        if (!is.na(map_mat[new_x, new_y])) {
          if (map_mat[new_x, new_y] == old_value) {
            map_mat[new_x, new_y] <- new_value
            q <- append(q, list(c(new_x, new_y)))
          }
        }
      }
    }
  }
  map_mat
  # browser()
}

# OPTIMIZED VERSION
# function to convolce a matrix with a given kernel
# (apparently the built-in convolve is different)
# note that it is not a super general convolve function but it
# slightly tuned for our scenario, in particular
# we are dealing with NA and we always binarized the matrix
# according if there are population or not
# convolve_mat_cpp <- function(mat, kernel) {
#   # transpose the kernel, cause since we are dealing
#   # with a map and not a picture its orientation
#   # is different
#   kernel <- t(kernel)
#   mat_dim <- dim(mat)
#   k_dim <- dim(kernel)
#   # check that the kernel has odd dimension
#   if(!all(k_dim%%2==1)){
#     stop("The kernel has to have odd dimensions.")
#   }
#   # pad the matrix
#   p_dim <- (k_dim-1)/2
#   # save na positions
#   is_na <- is.na(mat)
#   # convert them to -1 to facilitate the math
#   mat[is_na] <- -1
#   # make matrix binary
#   # -1 means no populations
#   mat[mat!=-1] <- 1
#   mat[mat==-1] <- 0
#   res <- conv_2d_cpp(mat, kernel)
#   attr(res, "resolution") <- get_resolution(mat)
#   class(res) <- set_class(res, "map")
#   # set as NA also cell that are already occupied
#   # as no populations can spawn there
#   res[mat==1] <- NA
#   # this to keep the positions that are NA
#   res[is_na] <- NA
#   dimnames(res) <- dimnames(mat)
#   res
# }

# function to update only the affected area
# of the convolved matrix
# update_conv_map_mat <- function(conv_map_mat, map_mat, x_new, y_new, radius, kernel){
#   # update the part of the convolved matrix that will chnage
#   # first extract the submatrix that will be affected
#   map_mat_neighb <- extract_neighbourhood(map_mat, x_new, y_new, radius = radius*2, meters = FALSE)
#   # convolve only this
#   conv_map_mat_neighb <- convolve_mat_cpp(map_mat_neighb, kernel)
#   # update that part of the convolved mat with this new one
#   # browser()
#   x_ys <- extract_names_aff_area(map_mat, x_new, y_new, radius)
#   xs <- as.character(x_ys$xs)
#   ys <- as.character(x_ys$ys)
#   conv_map_mat[xs, ys] <- conv_map_mat_neighb[xs, ys]
#   conv_map_mat
# }

# function to generate of a given radius in meters
# populate by all 1s, except the center that it is 0
#' Generate a 1s kernel
#'
#' @param map_mat An object of the class \code{prehistorik_map}.
#' @param radius Radius in meters for the kernel.
#'
#' @returns A matrix
#' @export
#'
generate_kernel <- function(map_mat, radius) {
  # to define the dimension of the kernel we need
  # to start from the resolution
  resolution <- get_resolution(map_mat)
  # the kernel have a radius in number of cells
  # as radius(m)/resolution
  radius_cells <- ceiling(radius / resolution)
  if (radius_cells == 1)
    stop(
      "The radius provided is to small.\n Either you increase the radius or decrase the resolution of the map."
    )
  # the dimensions of the kernel
  # is twice the radius_cell -1.
  # -1 because otherwise the central cell is counted twice
  # HAVE TO THINK BETTER HERE IF -1 OR +1
  dim <- 2 * radius_cells - 1
  kernel <- matrix(1, ncol = dim, nrow = dim)
  # set the center as 0
  kernel[radius_cells, radius_cells] <- 0
  kernel
}

# generate_ring_kernel <- function(map_mat, radius) {
#   # to define the dimension of the kernel we need
#   # to start from the resolution
#   resolution <- get_resolution(map_mat)
#   # the kernel have a radius in number of cells
#   # as radius(m)/resolution
#   radius_cells <- ceiling(radius / resolution)
#   if (radius_cells == 1)
#     stop(
#       "The radius provided is to small.\n Either you increase the radius or decrase the resolution of the map."
#     )
#   # the dimensions of the kernel
#   # is twice the radius_cell -1.
#   # -1 because otherwise the central cell is counted twice
#   # HAVE TO THINK BETTER HERE IF -1 OR +1
#   dim <- 2 * radius_cells - 1
#   kernel <- diag(dim)
#   kernel <- radius_cells - pmax(abs(row(kernel) - radius_cells), abs(col(kernel) - radius_cells))
#   # set the center as 0
#   kernel[radius_cells, radius_cells] <- 0
#
#   kernel
# }



# gauss_kernel <- function(map_mat, radius, sd=1){
#   resolution <- get_resolution(map_mat)
#   radius_cells <- ceiling(radius / resolution)
#   if (radius_cells == 1)
#     stop(
#       "The radius provided is to small.\n Either you increase the radius or decrase the resolution of the map."
#     )
#   # define the values of the kernel as if it
#   # was on the x-axis
#   axis <- -(radius_cells-1):(radius_cells-1)
#   # get the gaussian values for the axis
#   gauss <- dnorm(axis, sd=sd)
#   # get the kernel with the outer product
#   gauss <- gauss%o%gauss
#   gauss[radius_cells, radius_cells] <- 0
#   gauss
# }

# normalize matrix between 0.01 and 0.99
# 0-0.95 because I do not want the 1 to be there,
# as later on it would be interpreted as 0 probability
# normalize <- function(mat){
#   # check special case when min and mat have the same value
#   if(min(mat, na.rm = TRUE)==max(mat, na.rm = TRUE)){
#     mat/min(mat, na.rm = TRUE)*0.99
#   }else{
#     (mat-min(mat, na.rm = TRUE))/(max(mat, na.rm = TRUE)-min(mat, na.rm = TRUE))*(0.99-0.01)+0.01
#   }
# }

# function that samples the cell in the map
# where the next population will spawn
# sample_cell <- function(mat){
#   # melt the matrix
#   mat_melted <- reshape2::melt(mat, na.rm = TRUE, value.name = "weight")
#   mat_melted %>% dplyr::sample_n(size = 1, weight = weight) %>% dplyr::select(Var1, Var2) %>% dplyr::rename(x = Var1, y = Var2) %>% unlist()
# }

# helper function to extract the names of xs that are affected by the convolution
# radius is in cells
# extract_names_aff_area <- function(map_mat, x, y, radius){
#   # the area affected is the (radius+1)
#   xl <- x-radius
#   xr <- x+radius
#   yl <- y-radius
#   yr <- y+radius
#   # check whether we go out of the border
#   if(xl<1) xl <- 1
#   if(yl<1) yl <- 1
#   if(xr>dim(map_mat)[1]) xr <- dim(map_mat)[1]
#   if(yr>dim(map_mat)[2]) yr <- dim(map_mat)[2]
#   xs <- xl:xr
#   ys <- yl:yr
#   return(list(xs=xs, ys=ys))
# }

# function that given the coordinate of a population
# and a neighbourhood size it extracts the neighbour cells
# extract_neighbourhood <- function(map_mat, x, y, radius, meters=TRUE){
#   resolution <- get_resolution(map_mat)
#   # if the radius is in meters convert it
#   # otherwise it is the number of cells
#   # convert radius from meters to cells of kernel
#   # note I do not need -1 cause I am using floor
#   if(meters)  radius_cells <- floor(radius / resolution)
#   else radius_cells <- radius
#   if (radius_cells == 0)
#     stop(
#       "The radius provided is to small.\n Either you increase the radius or decrase the resolution of the map."
#     )
#   # set the limits of the neighbouroods
#   xl <- x-radius_cells
#   xr <- x+radius_cells
#   yl <- y-radius_cells
#   yr <- y+radius_cells
#   # check whether we go out of the border
#   if(xl<1) xl <- 1
#   if(yl<1) yl <- 1
#   if(xr>dim(map_mat)[1]) xr <- dim(map_mat)[1]
#   if(yr>dim(map_mat)[2]) yr <- dim(map_mat)[2]
#   # extract neigbourhood
#   neighbourhood <- map_mat[xl:xr, yl:yr]
#   # name it so to have the righ coordinates
#   rownames(neighbourhood) <- xl:xr
#   colnames(neighbourhood) <- yl:yr
#   attr(neighbourhood, "resolution") <- resolution
#   class(neighbourhood) <- set_class(neighbourhood, "map")
#   neighbourhood
# }


# coordinates conversion --------------------------------------------------

# function to convert the coordinate in the xy form to
# coordinates in map units
# note that the name is a bit weird cause I have another
# function called convert_coord
# xy_to_coords <- function(coord, resolution){
#   (coord-1)*resolution+resolution/2
# }
#
# # functions to convert x and y to latitude and longitude
# x_to_long <- function(map, map_mat, x){
#   # get xrange of slendr map
#   xrange <- attr(map, "xrange")
#   # extract dimensin matrix (rows are the longitude)
#   x_len <- dim(map_mat)[1]
#   # dimension of a cell in long
#   cell_len <- diff(xrange)/x_len
#   # long of the cell
#   xrange[1]+x*cell_len-cell_len/2
# }
#
# y_to_lat <- function(map, map_mat, y){
#   # get yrange of slendr map
#   yrange <- attr(map, "yrange")
#   # extract dimensin matrix (rows are the longitude)
#   y_len <- dim(map_mat)[2]
#   # dimension of a cell in long
#   cell_len <- diff(yrange)/y_len
#   # long of the cell
#   yrange[1]+y*cell_len-cell_len/2
# }


# given a point in latitude and longitude it gives back the corresponding cell
#' Title
#'
#' @param lat aa
#' @param lon aa
#' @param map_mat aa
#'
#' @returns aa
#' @export
#'
get_xy_from_latlon <- function(lat, lon, map_mat){
  if(!inherits(map_mat, "prehistorik_map")) stop("Not a prehisotrik_map.")
  # map_mat has to be prehistorik_map
  crs <- get_crs(map_mat)
  resolution <- get_resolution(map_mat)
  bbox <- get_bbox(map_mat)
  # transform lat lon into the crs
  coords <- data.frame(x = lon, y = lat)
  # authority_compliant = TRUE for having lat and lon -> y and x
  coords_proj <- ggplot2::sf_transform_xy(data = coords, source_crs = "EPSG:4326", target_crs = crs, authority_compliant = FALSE)
  # now we need to get the closest centroid to that point
  # subtract the xmin and ymin of the bbox and divide for the resolution
  xy <- unlist((coords_proj-bbox[c("xmin", "ymin")])/resolution)
  # get available position in the map
  cells <- data.frame(which(!is.na(map_mat), arr.ind = T))
  names(cells) <- c("x", "y")
  # and pick the closes to xy
  cell_min <- cells %>% dplyr::rowwise() %>% dplyr::mutate(dist_m=sqrt((x-xy[1])**2+(y-xy[2])**2)*resolution) %>% dplyr::ungroup() %>% dplyr::slice(which.min(dist_m))
  cell_min %>% select(x,y)
}


#' Title
#'
#' @param x aa
#' @param y aa
#' @param map_mat aa
#'
#' @returns aa
#' @export
#'
get_latlon_from_xy <- function(x, y, map_mat){
  if(!inherits(map_mat, "prehistorik_map")) stop("Not a prehisotrik_map.")
  # map_mat has to be prehistorik_map
  crs <- get_crs(map_mat)
  resolution <- get_resolution(map_mat)
  bbox <- get_bbox(map_mat)
  # transfor x y into meters
  coords <- c(x,y)*resolution+bbox[c("xmin", "ymin")]-resolution/2
  coords_df <- data.frame(x = coords[1], y = coords[2])
  # authority_compliant = TRUE for having lat and lon -> y and x
  coords_lonlat <- ggplot2::sf_transform_xy(data = coords_df, source_crs = crs, target_crs = "EPSG:4326", authority_compliant = FALSE)
  # coords_proj is in lon, lat
  names(coords_lonlat) <- c("lon", "lat")
  # and pick the closes to xy
  unlist(coords_lonlat)
}
