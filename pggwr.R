#' @title Parallel - Generalised geographically weighted regression
#' @name get
#' @author Victor Maus, \email{vwmaus1@@gmail.com}
#' 
#' @description The function implements generalised geographically weighted regression 
#' approach to exploring spatial non-stationarity for given global bandwidth and chosen 
#' weighting scheme.
#' 
#' @inheritParams twdtwMatches-class
#' @param n_cores an integer. Default all available colors minus one 
#' @praam ... other arguments see ?spgwr::ggwr 
#' 
pggwr <- function (formula, data = list(), coords, bandwidth, gweight = gwr.Gauss, 
                   adapt = NULL, fit.points, family = gaussian, longlat = NULL,
                   type = c("working", "deviance", "pearson", "response"), n_cores = parallel::detectCores() - 1){
  
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
  
  if (is.null(adapt)) {
    if (!missing(bandwidth)) {
      bw <- bandwidth
      bandwidth <- rep(bandwidth, n)
    }
    else stop("Bandwidth must be given for non-adaptive weights")
  }
  else {
    # Too slow ... NOT SURE IT WORKS 
    # bandwidth <- gw.adapt(dp = coords, fp = fit.points, quant = adapt, longlat = longlat)
    
    warning(cat("\nStarting parallel bandwidth processing using",n_cores,"cores\n"))
    cl <- parallel::makeForkCluster(nnodes = getOption("mc.cores", n_cores))
    
    # Set cluster env
    parallel::clusterEvalQ(cl, library(spgwr))
    varlist <- list("coords", "fit.points", "adapt", "longlat")
    env <- new.env()
    assign("coords", coords, envir = env)
    assign("fit.points", fit.points, envir = env)
    assign("adapt", adapt, envir = env)
    assign("longlat", longlat, envir = env)
    parallel::clusterExport(cl, varlist, env)
    
    # Split fit points 
    l_fp <- lapply(parallel::splitIndices(nrow(fit.points), length(cl)), function(i) fit.points[i, , drop = FALSE])
    
    # Run parallel processing 
    bandwidth <- unlist(parallel::parLapply(cl, l_fp, function(fp) gw.adapt(dp = coords, fp = fp, quant = adapt, longlat = longlat)))
    
    # Clean cluster env
    parallel::clusterEvalQ(cl, rm(varlist))
    rm(env)
    
    # Stop cluster 
    stopCluster(cl)
    
    bw <- bandwidth
  }
  
  if (any(bandwidth < 0)) 
    stop("Invalid bandwidth")
  
  fun <- function(i){
    dxs <- spDistsN1(coords, fit.points[i, ], longlat = longlat)
    if (any(!is.finite(dxs)))
      dxs[which(!is.finite(dxs))] <- .Machine$double.xmax/2
    
    w.i <- gweight(dxs^2, bandwidth[i])
    
    if (any(w.i < 0 | is.na(w.i)))
      stop(paste("Invalid weights for i:", i))
    
    lm.i <- glm.fit(y = y, x = x, weights = w.i, offset = offset, family = family)
    sum.w <- sum(w.i)
    gwr.b <- coefficients(lm.i)
    
    v_resids <- 0
    if (!fp.given)
      v_resids <- residuals.glm(lm.i, type = type)[i]
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
    return(data.frame(sum.w = sum.w, gwr.b, dispersion = dispersion, v_resids = v_resids))
  }
  
  # Start cluster 
  warning(cat("\nStarting parallel processing using",n_cores,"cores\n"))
  cl <- parallel::makeForkCluster(nnodes = getOption("mc.cores", n_cores))
  
  # Set cluster env
  parallel::clusterEvalQ(cl, library(spgwr))
  varlist <- list("coords", "fit.points", "bandwidth", "x", "y", "family", "offset", "type", "fp.given", "longlat")
  env <- new.env()
  assign("coords", coords, envir = env)
  assign("fit.points", fit.points, envir = env)
  assign("y", y, envir = env)
  assign("x", x, envir = env)
  assign("bandwidth", weights, envir = env)
  assign("family", family, envir = env)
  assign("offset", offset, envir = env)
  assign("type", type, envir = env)
  assign("fp.given", fp.given, envir = env)
  assign("longlat", longlat, envir = env)
  parallel::clusterExport(cl, varlist, env)
  
  # Run parallel processing 
  df <- dplyr::bind_rows(parallel::parLapply(cl, 1:n, fun = fun))
  
  # Clean cluster env
  parallel::clusterEvalQ(cl, rm(varlist))
  rm(env)

  # Stop cluster 
  stopCluster(cl)
  
  SDF <- SpatialPointsDataFrame(coords = fit.points, data = df, proj4string = CRS(p4s))
  if (griddedObj) {
    gridded(SDF) <- TRUE
  }
  else {
    if (!is.null(Polys)) {
      df <- data.frame(SDF@data)
      rownames(df) <- sapply(slot(Polys, "polygons"), function(i) slot(i,"ID"))
      SDF <- SpatialPolygonsDataFrame(Sr = Polys, data = df)
    }
  }

  z <- list(SDF = SDF, lhat = NA, lm = glm_fit, results = NULL,
            bandwidth = bw, adapt = adapt, hatmatrix = FALSE, gweight = deparse(substitute(gweight)),
            fp.given = fp.given, this.call = this.call)

  class(z) <- "gwr"

  invisible(z)
  
}

