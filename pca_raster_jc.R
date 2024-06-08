M_simulationR_jc=function (data, simulation_variables, starting_proportion = 0.5, 
                           sampling_rule = "random", barriers = NULL, scale = TRUE, 
                           center = TRUE, project = FALSE, projection_variables = NULL, 
                           dispersal_kernel = "normal", kernel_spread = 1, max_dispersers = 4, 
                           suitability_threshold = 5, replicates = 10, dispersal_events = 25, 
                           access_threshold = 5, simulation_period = 50, stable_lgm = 7, 
                           transition_to_lgm = 100, lgm_to_current = 7, stable_current = 13, 
                           scenario_span = 1, out_format = "GTiff", set_seed = 1, write_all_scenarios = FALSE, 
                           output_directory) 
{
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  }
  if (missing(simulation_variables)) {
    stop("Argument 'simulation_variables' must be defined")
  }
  if (missing(output_directory)) {
    stop("Argument 'output_directory' must be defined")
  }
  if (!sampling_rule %in% c("random", "suitability")) {
    stop("Argument 'sampling_rule' is not valid")
  }
  n <- raster::nlayers(simulation_variables)
  if (n < 2) {
    stop("At least 2 variables are needed in 'simulation_variables'")
  }
  if (project == TRUE) {
    if (missing(projection_variables)) {
      stop("If 'project' = TRUE, argument 'projection_variables' must be defined")
    }
    if (!all(simulation_variables@extent == projection_variables@extent)) {
      stop("'projection_variables' and 'simulation_variables' must have the same extent")
    }
    if (!all(names(simulation_variables) == names(projection_variables))) {
      stop("Variable names in 'projection_variables' and 'simulation_variables' must bee the same")
    }
  }
  if (!is.null(barriers)) {
    if (simulation_variables@extent != barriers@extent) {
      stop("'barriers' and 'simulation_variables' must have the same extent")
    }
  }
  message("Preparing data to run simulation...")
  dir.create(output_directory)
  ftype <- rformat_type(out_format)
  sp_nam <- as.character(data[1, 1])
  data <- data[,2:3]
  suit_fol <- paste0(output_directory, "/Suitability_results")
  dir.create(suit_fol)
  npcs <- ifelse(n > 3, 3, n)
  opca_fol <- paste0(output_directory, "/PCA_results")
  if (project == FALSE) {
    simulation_variables <- pca_raster_jc(pca_variables = simulation_variables, 
                                       scale = scale, center = center, n_pcs = npcs, project = project, 
                                       write_to_directory = TRUE, output_directory = opca_fol)[[2]]
    suit_mod <- ellipsoid_suitability_jc(data, simulation_variables, 
                                      suitability_threshold)
    suit_layer <- suit_mod[[5]]
    if (!is.null(barriers)) {
      message("\nSuitability layer will be corrected using barriers")
      barr <- is.na(barriers)
      suit_layer <- suit_layer * barr
    }
    emodfile <- paste0(suit_fol, "/Ellipsoid_metadata.rds")
    write_ellmeta(suit_mod, emodfile)
    s_name <- paste0(suit_fol, "/suitability", ftype)
    raster::writeRaster(suit_layer, filename = s_name, format = out_format)
    suit_name <- normalizePath(s_name)
  } else {
    print("Not yet formatted for projection.")
    break
    pcs <- pca_raster(variables = simulation_variables, scale = scale, 
                      center = center, n_pcs = npcs, project = project, 
                      projection_variables = projection_variables, return_projection = TRUE, 
                      write_to_directory = TRUE, out_format = out_format, 
                      output_directory = opca_fol)[c(2, 3)]
    simulation_variables <- pcs[[1]]
    projection_variables <- pcs[[2]][[1]]
    suit_mod <- ellipsoid_suitability(data, 
                                      simulation_variables, 
                                      suitability_threshold, 
                                      project = project, 
                                      projection_variables = projection_variables)
    suit_layer <- suit_mod[[5]][[1]]
    suit_lgm <- suit_mod[[5]][[2]]
    if (!is.null(barriers)) {
      message("\nSuitability layers will be corrected using barriers")
      barr <- is.na(barriers)
      suit_layer <- suit_layer * barr
      suit_lgm <- suit_lgm * barr
    }
    else {
      barr <- barriers
    }
    emodfile <- paste0(suit_fol, "/Ellipsoid_metadata")
    write_ellmeta(suit_mod, emodfile)
    s_name <- paste0(suit_fol, "/suitability_current", ftype)
    raster::writeRaster(suit_layer, filename = s_name, format = out_format)
    l_name <- paste0(suit_fol, "/suitability_lgm", ftype)
    raster::writeRaster(suit_lgm, filename = l_name, format = out_format)
    sp_name <- normalizePath(s_name)
    lp_name <- normalizePath(l_name)
    message("\nPreparing interpolations:")
    int_vals <- interpolation_values(simulation_period, transition_to_lgm, 
                                     stable_lgm, lgm_to_current, stable_current, scenario_span)
    suit_name <- interpolation(suit_mod, suitability_threshold, 
                               int_vals, barr, simulation_variables, projection_variables, 
                               sp_name, lp_name, out_format, suit_fol)
  }
  occ_suit <- suit_mod[[1]][, 1:2]
  suit_lay <- raster::raster(suit_name[1])
  oca <- data.frame(Species = sp_nam, occ_suit)
  oca <- suitable_cells(suit_lay, data = oca)
  oca_nam <- paste0(output_directory, "/occ_simulation.csv")
  write.csv(oca[, 1:3], oca_nam, row.names = FALSE)
  save_nicheplot(suit_mod, suitability_threshold, simulation_variables, 
                 size_proportion = 0.55, suit_fol)
  out_dir <- normalizePath(output_directory)
  message("")
  res <- dispersal_simulationR(data = oca, suit_layers = suit_name, 
                               starting_proportion = starting_proportion, proportion_to_disperse = 1, 
                               sampling_rule = sampling_rule, dispersal_kernel = dispersal_kernel, 
                               kernel_spread = kernel_spread, max_dispersers = max_dispersers, 
                               dispersal_events = dispersal_events, replicates = replicates, 
                               threshold = access_threshold, results_by = "scenario", 
                               set_seed = set_seed, return = "accessed", write_to_directory = TRUE, 
                               write_all = write_all_scenarios, raster_format = out_format, 
                               output_directory = out_dir)
  message("\nPreparing M as a spatial polygon...")
  m <- M_preparation(output_directory, res$A, raster_format = out_format)
  m_poly <- m[[2]]
  m <- m[[1]]
  save_Mplot(suit_mod, suit_layer, m_poly, size_proportion = 0.55, 
             output_directory)
  message("\nM simulation finished\nCheck results in:  ", out_dir, 
          "\n")
  return(list(Simulation_occurrences = oca, Simulation_scenarios = suit_name, 
              Summary = res$Summary, A_raster = m, A_polygon = m_poly, 
              A_mean = res$A_mean, A_var = res$A_var, Barriers = barriers))
}

pca_raster_jc=function (pca_variables, in_format = NULL, scale = TRUE, center = TRUE, 
                        n_pcs = NULL, project = FALSE, projection_pca_variables, return_projection = FALSE, 
                        write_to_directory = FALSE, out_format = "GTiff", output_directory) 
{
  if (missing(pca_variables)) {
    stop("Argument 'pca_variables' must be defined")
  }
  clsv <- class(pca_variables)[1]
  if (!clsv %in% c("character", "RasterStack", "RasterBrick")) {
    stop("'pca_variables' must be of class 'character' or 'RasterStack'")
  }
  if (clsv == "character") {
    if (is.null(in_format)) {
      stop("If class of 'pca_variables' = 'character', 'in_format' must be defined")
    }
  }
  if (project == TRUE) {
    if (missing(projection_pca_variables)) {
      stop("If 'project' = TRUE, argument 'projection_pca_variables' must be defined")
    }
    if (return_projection == FALSE & write_to_directory == 
        FALSE) {
      message("Setting 'return_projection' = TRUE to keep results from projections")
      return_projection <- TRUE
    }
  }
  if (write_to_directory & missing(output_directory)) {
    stop("If 'write_to_directory' = TRUE, 'output_directory' must be defined")
  }
  patt1 <- rformat_type(out_format)
  if (write_to_directory == TRUE) {
    dir.create(output_directory)
    pca_fol <- paste0(output_directory, "/Initial")
    dir.create(pca_fol)
  }
  if (clsv == "character") {
    patt <- paste0(rformat_type(in_format), "$")
    var <- list.files(pca_variables, pattern = patt, full.names = TRUE)
    pca_variables <- raster::stack(var)
  }
  var_points <- pca_variables[]
  to_fill <- which(complete.cases(var_points))
  var_points <- na.omit(var_points)
  pcra <- pca_variables[[1]]
  pcra[] <- NA
  if (is.null(n_pcs)) {
    n_pcs <- ncol(var_points)
  }
  pca <- prcomp(var_points, center = center, scale = scale)
  pcras <- list()
  message("Preparing raster PCs...")
  for (i in 1:n_pcs) {
    pcra2 <- pcra
    pcra2[to_fill] <- pca$x[, i]
    pcras[[i]] <- pcra2
    if (write_to_directory == TRUE) {
      filenam <- paste0(pca_fol, "/PC", i, patt1)
      raster::writeRaster(pcra2, filenam, format = out_format)
    }
  }
  pcras <- do.call(raster::stack, pcras)
  names(pcras) <- paste0("PC", 1:(raster::nlayers(pcras)))
  SumPCAMat <- summary(pca)$importance
  if (write_to_directory == TRUE) {
    txtfile <- paste0(pca_fol, "/PCA_summary.txt")
    cat("Principal component analysis:\n\nPCA summary\n", 
        file = txtfile)
    suppressWarnings(write.table(SumPCAMat, sep = "\t", file = txtfile, 
                                 append = TRUE, quote = FALSE))
    cat("\n\nPCA loadings\n", file = txtfile, append = TRUE)
    suppressWarnings(write.table(pca$rotation, sep = "\t", 
                                 file = txtfile, append = TRUE, quote = FALSE))
    rd_file <- paste0(pca_fol, "/PCA_results.RData")
    save(pca, file = rd_file)
  }
  if (project == TRUE) {
    clspr <- class(projection_pca_variables)[1]
    if (return_projection == TRUE) {
      ppcrass <- list()
    }
    message("Projecting PCs to distinct scenarios...")
    if (clspr == "character") {
      proj_dirs <- list.dirs(projection_pca_variables, recursive = FALSE)
      proj_names <- list.dirs(projection_pca_variables, recursive = FALSE, 
                              full.names = FALSE)
      fol_names <- paste(output_directory, proj_names, 
                         sep = "/")
    }
    else {
      proj_dirs <- "projection"
      proj_names <- "Projected_PCs"
      fol_names <- paste(output_directory, proj_names, 
                         sep = "/")
    }
    for (h in 1:length(proj_dirs)) {
      if (clspr == "character") {
        pvar <- list.files(proj_dirs[h], pattern = patt, 
                           full.names = TRUE)
        projection_pca_variables <- raster::stack(pvar)
      }
      pcra <- projection_pca_variables[[1]]
      to_fillp <- !is.na(pcra[])
      if (write_to_directory == TRUE) {
        dir.create(fol_names[h])
      }
      if (return_projection == TRUE) {
        ppcras <- list()
      }
      p_stackp <- na.omit(projection_pca_variables[])
      p_pcs <- predict(pca, newdata = p_stackp)
      for (i in 1:n_pcs) {
        pcra[to_fillp] <- p_pcs[, i]
        if (return_projection == TRUE) {
          ppcras[[i]] <- pcra
        }
        if (write_to_directory == TRUE) {
          filenam <- paste0(fol_names[h], "/PC", i, patt1)
          raster::writeRaster(pcra, filenam, format = out_format)
        }
      }
      if (return_projection == TRUE) {
        ppcrass[[h]] <- do.call(raster::stack, ppcras)
        names(ppcrass[[h]]) <- paste0("PC", 1:(raster::nlayers(ppcrass[[h]])))
      }
    }
    if (return_projection == TRUE) {
      names(ppcrass) <- paste0("PCRaster_", proj_names)
    }
  }
  if (project == TRUE & return_projection == TRUE) {
    results <- list(PCA_results = pca, PCRaster_initial = pcras, 
                    PCRaster_projection = ppcrass)
  }
  else {
    results <- list(PCA_results = pca, PCRaster_initial = pcras, 
                    PCRaster_projection = NULL)
  }
  message("Raster PCA finished")
  if (write_to_directory == TRUE) {
    message("Check results in:  ", normalizePath(output_directory))
  }
  return(results)
}

ellipsoid_suitability_jc=function(data, ell_variables, suitability_threshold = 5, project = FALSE, 
                                   projection_ell_variables, tolerance = 1e-60){
  raster_formats <- c("RasterStack", "RasterBrick")
  if (missing(data)) {
    stop("Argument 'data' must be defined")
  } else {
    if (ncol(data) != 2) {
      stop("Columns in 'data' must be: longitude and latitude; in that order")
    }
  }
  if (missing(ell_variables)) {
    stop("Argument 'ell_variables' must be defined")
  } else {
    if (!class(ell_variables)[1] %in% raster_formats) {
      stop("Argument 'ell_variables' must be of class RasterStack")
    }
  }
  if (project == TRUE) {
    if (missing(projection_ell_variables)) {
      stop("If projections are needed, argument 'projection_ell_variables' must be defined")
    }
    if (!class(projection_ell_variables)[1] %in% c(raster_formats, 
                                               "list")) {
      stop("Argument 'projection_ell_variables' must be of class RasterStack or list")
    }
  }
  occ_data <- raster::extract(ell_variables, 
                              data, 
                              na.rm=F) # keep to index
  full.index=which(complete.cases(occ_data))
  occ_data=occ_data[full.index,] # amended to remove NAs
  data2=data[full.index,] # creates dataframe with same row number
  centroid <- colMeans(occ_data)
  covariance <- cov(occ_data)
  if (project == FALSE) {
    suit_model <- predict_esuitability_jc(centroid = centroid, 
                                       covariance_matrix = covariance, 
                                       suit_variables = ell_variables, 
                                       suitability_threshold = suitability_threshold, 
                                       tolerance = tolerance)
    suit <- na.omit(raster::extract(suit_model[[4]], data))
    occ_comp <- cbind(data2, occ_data, suit)
    colnames(occ_comp) <- c("Longitude", "Latitude", names(ell_variables), 
                            "Suitability")
    results <- c(list(occurrences = occ_comp), suit_model)
  } else {
    not_suitable <- list()
    suit_layer <- list()
    if (class(projection_ell_variables)[1] == "RasterStack") {
      projection_ell_variables <- list(ell_variables, projection_ell_variables)
      proj_names <- "projection"
    } else {
      if (is.null(names(projection_ell_variables))) {
        proj_names <- paste0("projection_", 1:length(projection_ell_variables))
      } else {
        proj_names <- paste0("projection_", names(projection_ell_variables))
      }
      projection_ell_variables <- c(ell_variables, projection_ell_variables)
    }
    suit_names <- c("suitability_layer", proj_names)
    for (i in 1:length(projection_ell_variables)) {
      suit_model <- predict_esuitability(centroid = centroid, 
                                         covariance_matrix = covariance, suit_variables = projection_ell_variables[[i]], 
                                         suitability_threshold = suitability_threshold, 
                                         tolerance = tolerance)
      not_suitable[[i]] <- suit_model[[3]]
      suit_layer[[i]] <- suit_model[[4]]
    }
    scenarios <- rep(c("Initial", paste("Transfer area", 
                                        1:(length(projection_ell_variables) - 1))), each = 2)
    not_suitable <- data.frame(Scenario = scenarios, do.call(rbind, 
                                                             not_suitable))
    names(suit_layer) <- suit_names
    suit <- na.omit(raster::extract(suit_layer[[1]], data))
    occ_comp <- cbind(data, occ_data, suit)
    colnames(occ_comp) <- c("Longitude", "Latitude", names(ell_variables), 
                            "Suitability")
    results <- list(occurrences = occ_comp, centroid = suit_model[[1]], 
                    covariance_matrix = suit_model[[2]], suitable_area_prop = not_suitable, 
                    suitability_layer = suit_layer)
  }
  return(results)
}

predict_esuitability_jc=function (ellipsoid_model = NULL, centroid = NULL, covariance_matrix = NULL, 
                                  suit_variables, suitability_threshold = 5, tolerance = 1e-60) 
{
  raster_formats <- c("RasterStack", "RasterBrick")
  if (!missing(ellipsoid_model)) {
    occ <- ellipsoid_model[[1]]
    centroid <- ellipsoid_model[[2]]
    covariance_matrix <- ellipsoid_model[[3]]
  } else {
    if (missing(centroid) | missing(covariance_matrix)) {
      stop("Argument 'ellipsoid_model' missing, 'centroid' and 'covariance_matrix' must be defined")
    }
  }
  if (missing(suit_variables)) {
    stop("Argument 'suit_variables' must be defined")
  } else {
    if (!class(suit_variables)[1] %in% raster_formats) {
      stop("Argument 'suit_variables' must be of class RasterStack")
    }
  }
  back <- na.omit(suit_variables[])
  # back <- suit_variables[]
  # back <- back[which(complete.cases(back)),]
  # back <- as.matrix(back)
  if(ncol(back)>3){back <- back[,1:3]}
  maha <- mahalanobis(x = back, center = centroid, cov = covariance_matrix, 
                      tol = tolerance)
  alpha <- (100 - suitability_threshold)/100
  chi_sq <- qchisq(alpha, ncol(back))
  suitability <- exp(-0.5 * maha)
  suitability <- ifelse(maha/chi_sq <= 1, suitability, 0)
  p_no_suit_g <- sum(suitability != 0)/length(suitability)
  u_suit <- suitability[!duplicated(back)]
  p_no_suit_e <- sum(u_suit != 0)/length(u_suit)
  suitable <- data.frame(c("Geographic_space", "Environmental_space"), 
                         c(p_no_suit_g, p_no_suit_e))
  colnames(suitable) <- c("Space", "Proportion_suitable")
  suit_layer <- suit_variables[[1]]
  suit_layer[!is.na(suit_layer[])] <- suitability
  results <- list(centroid = centroid, covariance_matrix = covariance_matrix, 
                  suitable_area_prop = suitable, suitability_layer = suit_layer)
  return(results)
}

suitable_cells <- function(suit_layer, data = NULL) {
  if (is.null(data)) {
    noz <- which((suit_layer[] > 0))
    suit_bar <- suit_layer[noz]
    noz <- raster::xyFromCell(suit_layer, noz)
  } else {
    suit_bar <- raster::extract(suit_layer, data[, 2:3])
    tokeep <- suit_bar > 0 & !is.na(suit_bar)
    noz <- data[tokeep, 2:3]
    suit_bar <- suit_bar[tokeep]
  }
  sp <- ifelse(is.null(data), "Species", as.character(data[1, 1]))
  
  
  return(data.frame(species = sp, longitude = noz[, 1], latitude = noz[, 2],
                    suitability = suit_bar))
}
