#' @title Parallel - Generalised geographically weighted regression
#' @name get
#' 
#' @description This is a parallel version adpted from spgwr::ggwr. The function runs a 
#' generalised geographically weighted regression approach to exploring spatial non-stationarity 
#' for given global bandwidth and chosen weighting scheme. 
#' 
#' @inheritParams twdtwMatches-class
#' @param n_cores an integer. Default is all available cores minus one 
#' @param max_dist maximum distance from the sample point in the metric of the points if longlat = FALSE, 
#' or in kilometers if longlat = TRUE. Default is NULL.
#' @param min_weight minimum weight for a sample point be considered in the regression. Default is NULL.
#' @praam ... for other arguments see ?spgwr::ggwr 
#' 
par_ggwr <- function (formula, data = list(), coords, bandwidth, gweight = gwr.Gauss, 
                   adapt = NULL, fit.points, family = gaussian, longlat = NULL,
                   type = c("working", "deviance", "pearson", "response"), n_cores = parallel::detectCores() - 1, max_dist = NULL, min_weight = NULL){
  
  type <- match.arg(type)
  resid_name <- paste(type, "resids", sep = "_")
  this.call <- match.call()
  p4s <- as.character(NA)
  Polys <- NULL
  if (is(data, "SpatialPolygonsDataFrame")) 
    Polys <- as(data, "SpatialPolygons")
  if (is(data, "Spatial")) {
    if (!missing(coords)) 
      warning("data is Spatial* object, ignoring coords argument")
    coords <- coordinates(data)
    p4s <- proj4string(data)
    if (is.null(longlat) || !is.logical(longlat)) {
      if (!is.na(is.projected(data)) && !is.projected(data)) {
        longlat <- TRUE
      }
      else {
        longlat <- FALSE
      }
    }
    data <- as(data, "data.frame")
  }
  if (is.null(longlat) || !is.logical(longlat)) 
    longlat <- FALSE
  if (missing(coords)) 
    stop("Observation coordinates have to be given")
  if (is.null(colnames(coords))) 
    colnames(coords) <- c("coord.x", "coord.y")
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  offset <- model.offset(mf)
  if (!is.null(offset) && length(offset) != NROW(y)) 
    stop("number of offsets should equal number of observations")
  if (is.null(offset)) 
    offset <- rep(0, length(c(y)))
  
  glm_fit <- glm.fit(x = x, y = y, offset = offset, family = family)
  if (missing(fit.points)) {
    fp.given <- FALSE
    fit.points <- coords
    colnames(fit.points) <- colnames(coords)
  }
  else fp.given <- TRUE
  griddedObj <- FALSE

  if (is(fit.points, "Spatial")) {
    Polys <- NULL
    if (is(fit.points, "SpatialPolygonsDataFrame")) {
      Polys <- Polygons(fit.points)
      fit.points <- coordinates(fit.points)
    }
    else {
      griddedObj <- gridded(fit.points)
      fit.points <- coordinates(fit.points)
    }
  }
  n <- NROW(fit.points)
  rownames(fit.points) <- NULL
  if (is.null(colnames(fit.points))) 
    colnames(fit.points) <- c("x", "y")
  m <- NCOL(x)
  if (NROW(x) != NROW(coords)) 
    stop("Input data and coordinates have different dimensions")
  
  if(Sys.info()["sysname"] %in% c("Windows", "windows")){
    cl <- parallel::makePSOCKcluster(names = getOption("mc.cores", n_cores))
  } else {
    cl <- parallel::makeForkCluster(nnodes = getOption("mc.cores", n_cores))
  }
  
  # Set cluster env
  parallel::clusterEvalQ(cl, library(stats))
  parallel::clusterEvalQ(cl, library(sp))
  parallel::clusterEvalQ(cl, library(spgwr))
  varlist <- list("coords", "x", "y", "family", "offset", "type", "fp.given", "longlat", "adapt", "n", "max_dist")
  env <- new.env()
  assign("n", n, envir = env)
  assign("m", m, envir = env)
  assign("coords", coords, envir = env)
  assign("y", y, envir = env)
  assign("x", x, envir = env)
  assign("family", family, envir = env)
  assign("offset", offset, envir = env)
  assign("type", type, envir = env)
  assign("fp.given", fp.given, envir = env)
  assign("adapt", adapt, envir = env)
  assign("max_dist", max_dist, envir = env)
  
  if (is.null(adapt)) {
    if (!missing(bandwidth)) {
      bw <- bandwidth
      fit.points <- cbind(fit.points, rep(bandwidth, n))
    }
    else stop("Bandwidth must be given for non-adaptive weights")
  }
  else {

    cat("\nComputing bandwidth using",n_cores,"cores...")
    start_time <- Sys.time()
    bw <- parallel::parRapply(cl, fit.points, function(fp) spgwr::gw.adapt(dp = coords, fp = cbind(fp[1], fp[2]), quant = adapt, longlat = longlat))
    cat(" Done\n")
    fit.points <- cbind(fit.points, bw)
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  
  if (any(bw < 0) | any(is.na(bw))) 
    stop("Invalid bandwidth NA or lower than 0 (zero). It must be greater than 0 (zero)")

  fun <- function(pto){
    
    dxs <- sp::spDistsN1(coords, cbind(pto[1], pto[2]), longlat = longlat)
    
    if (any(!is.finite(dxs)))
      dxs[which(!is.finite(dxs))] <- .Machine$double.xmax/2

    I <- seq_along(dxs)
    
    if(!is.null(max_dist))
      I <- which(dxs <= max_dist)

    xx <- matrix(x[I,], ncol = 2)
    yy <- y[I]
    of <- offset[I]
    w.i <- gweight(dxs[I]^2, pto[3])
    
    J <- seq_along(w.i)
    if(!is.null(min_weight))
      J <- which(w.i >= min_weight)
    
    lm.i <- glm.fit(y = yy[J], x = matrix(xx[J,], ncol = 2), weights = w.i[J], offset = of[J], family = family)
    
    sum.w <- sum(w.i)
    gwr.b <- matrix(coefficients(lm.i), nrow = 1, ncol = m)
    colnames(gwr.b) <- colnames(x)
    
    v_resids <- 0
    if (!fp.given)
      v_resids <- residuals.glm(lm.i, type = type)
    else is.na(v_resids) <- TRUE
    
    df.r <- lm.i$df.residual
    
    if (lm.i$family$family %in% c("poisson", "binomial"))
      dispersion <- 1
    
    else {
      if (df.r > 0) {
        dispersion <- sum((lm.i$weights * lm.i$residuals^2)[lm.i$weights > 0])/df.r
      }
      else {
        dispersion <- NaN
      }
    }
    df <- data.frame(x = pto[1], y = pto[2], bw = pto[3], sum.w = sum.w, gwr.b, dispersion = dispersion)
    df[[resid_name]] <- v_resids
    return(df)
  }

  cat("\nProcessing using",n_cores,"cores...")
  start_time <- Sys.time()
  df <- parallel::parRapply(cl, fit.points, FUN = fun)
  # df <- apply(fit.points, 1, FUN = fun)
  cat(" Done\n")
  end_time <- Sys.time()
  print(end_time - start_time)

  cat("\nCreating output data frame...")
  start_time <- Sys.time()
  df <- dplyr::bind_rows(df)
  cat(" Done\n")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # Clean cluster env
  parallel::clusterEvalQ(cl, rm(varlist))
  rm(env)

  # Stop cluster 
  stopCluster(cl)

  xy <- df %>% 
    dplyr::select(x, y)
  
  bw <- df %>% 
    dplyr::select(bw)
  
  df <- df %>% 
    dplyr::select(-c(x, y, bw))

  cat("\nBuilding Spatial object...")  
  start_time <- Sys.time()
  SDF <- SpatialPointsDataFrame(coords = xy, data = df, proj4string = CRS(p4s))
  cat(" Done\n")
  end_time <- Sys.time()
  print(end_time - start_time)
  
  if (griddedObj) {
    gridded(SDF) <- TRUE
  }
  else {
    if (!is.null(Polys)) {
      cat("\nBuilding polygons...")
      start_time <- Sys.time()
      df <- data.frame(SDF@data)
      rownames(df) <- sapply(slot(Polys, "polygons"), function(i) slot(i,"ID"))
      SDF <- SpatialPolygonsDataFrame(Sr = Polys, data = df)
      cat(" Done\n")
      end_time <- Sys.time()
      print(end_time - start_time)
    }
  }

  z <- list(SDF = SDF, lhat = NA, lm = glm_fit, results = NULL,
            bandwidth = bw, adapt = adapt, hatmatrix = FALSE, gweight = deparse(substitute(gweight)),
            fp.given = fp.given, this.call = this.call)

  class(z) <- "gwr"

  invisible(z)
  
}

