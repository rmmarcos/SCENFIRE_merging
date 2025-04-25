# Function to merge all simulated fire perimeters into a single layer. It calcultes area size in hectares (epsg:25830)
read_perimeters <- function(n, shapefiles){
  
  foo <- read_sf(shapefiles[n])
  scenario <- str_split(shapefiles[n],pattern = '/')[[1]][7]
  foo$WD <- str_split(scenario,pattern = '_')[[1]][3]
  foo$Duration <- str_split(scenario,pattern = '_')[[1]][4]
  
  return(foo)
  
}

# Function to calculate the relative distance between the target distribution and the selected perimeters
calculate_discrepancy <- function(selected_surfaces, target_hist, bins, logaritmic=T) {
  # Calcular la densidad del histograma del conjunto seleccionado
  if(logaritmic==T){
    hist_selected <- hist(log(selected_surfaces + 1e-6), breaks = bins, plot = FALSE)$density
  }else{
    hist_selected <- hist((selected_surfaces), breaks = bins, plot = FALSE)$density
  }
  
  
  # Asegurarse de que las densidades tengan la misma longitud que el histograma objetivo
  if (length(hist_selected) != length(target_hist)) {
    stop("Las densidades del histograma seleccionado no coinciden en longitud con el histograma objetivo.")
  }
  
  # Calcular la diferencia absoluta entre densidades
  abs_diff <- abs(hist_selected - target_hist)
  
  # Calcular la discrepancia como el área relativa
  total_area_target <- sum(target_hist)  # Área total del histograma objetivo
  if (total_area_target == 0) {
    stop("El área total del histograma objetivo es cero.")
  }
  
  relative_discrepancy <- sum(abs_diff) / total_area_target  # Discrepancia relativa
  
  return(relative_discrepancy)
}

# Function to select simulated perimeters
#' Select events based on surface distribution matching
#' 
#' @param event_surfaces Numeric vector of surface values
#' @param event_probabilities Numeric vector of probabilities for each event
#' @param target_hist Numeric vector representing target histogram distribution
#' @param bins Numeric vector of bin breakpoints
#' @param reference_surface Numeric value for reference surface
#' @param surface_threshold Numeric value between 0 and 1 for surface threshold
#' @param tolerance Numeric value for acceptable tolerance
#' @param max_it Integer for maximum iterations (default: 100)
#' @param logaritmic Logical for logarithmic binning (default: TRUE)
#' @return List containing selected surfaces, indices, total surface, and final discrepancy
#' @export
select_events <- function(event_surfaces, event_probabilities, target_hist, bins,
                                          reference_surface, surface_threshold, tolerance, max_it = 100, logaritmic = T) {
  max_iterations <- max_it
  
  # Clasificar eventos en bins
  if(logaritmic == T) {
    bin_indices <- cut(log(event_surfaces + 1e-6), breaks = bins, include.lowest = TRUE, labels = FALSE)
  } else {
    bin_indices <- cut(event_surfaces, breaks = bins, include.lowest = TRUE, labels = FALSE)  
  }
  
  results <- ({
    foreach(block = 1:max_iterations,
            .combine = list,
            .multicombine = TRUE,
            .packages = c("stats"),
            .export = "calculate_discrepancy") %dopar% {
              
              best_selection <- NULL
              best_discrepancy <- Inf
              
              for (i in 1:max_it) {
                # Usar tryCatch para manejar posibles errores
                result <- tryCatch({
                  selected_surfaces <- c()
                  surface_index <- c()
                  total_surface <- 0
                  current_hist <- rep(0, length(target_hist))  # Inicializar histograma de la selección actual
                  
                  # Calcular la frecuencia de los bins con la tabla de contingencia
                  bin_counts <- table(bin_indices)
                  bin_priorities <- rep(0, length(bins) - 1)  # Inicializar prioridades para cada bin
                  
                  for (bin in 1:(length(bins) - 1)) {
                    if (bin %in% names(bin_counts)) {
                      bin_discrepancy <- abs(target_hist[bin] - current_hist[bin])
                      bin_priorities[bin] <- bin_discrepancy
                    }
                  }
                  
                  # Mantener un registro de índices seleccionados
                  selected_indices_global <- c()
                  
                  while ((total_surface < reference_surface * surface_threshold)) {
                    # Selección probabilística basada en la prioridad de los bins
                    if (any(bin_priorities > 0)) {  # Solo seleccionar si hay algún bin con prioridad positiva
                      selected_bin <- sample(1:length(bin_priorities), 1, prob = bin_priorities + 1e-6)
                      
                      # Filtrar índices elegibles excluyendo los seleccionados previamente
                      eligible_indices <- unique(which(bin_indices == selected_bin))  # Asegurar unicidad
                      eligible_indices <- setdiff(eligible_indices, selected_indices_global)  # Excluir ya seleccionados
                      
                      while (length(eligible_indices) > 0) {
                        if (length(eligible_indices) == 1) {
                          selected_index <- eligible_indices[1]  # Selección directa si solo hay un índice elegible
                        } else if (length(eligible_indices) > 1) {
                          selected_index <- sample(eligible_indices, 1, prob = event_probabilities[eligible_indices])
                        } else {
                          break  # Salir si no hay índices elegibles
                        }
                        
                        # Verificar si el índice ya ha sido seleccionado
                        if (!(selected_index %in% selected_indices_global)) {
                          # Si no está repetido, procesar la selección
                          selected_surface <- event_surfaces[selected_index]
                          
                          # Actualizar selección
                          selected_surfaces <- c(selected_surfaces, selected_surface)
                          surface_index <- c(surface_index, selected_index)
                          total_surface <- sum(selected_surfaces)
                          
                          # Registrar el índice seleccionado globalmente
                          selected_indices_global <- c(selected_indices_global, selected_index)
                          
                          # Actualizar histograma de la selección actual sin transformación logarítmica
                          if(logaritmic == T) {
                            current_hist <- hist(log(selected_surfaces + 1e-6), breaks = bins, plot = FALSE)$density
                          } else {
                            current_hist <- hist(selected_surfaces, breaks = bins, plot = FALSE)$density
                          }
                          
                          break
                        } else {
                          # Si está repetido, eliminarlo de los índices elegibles y continuar buscando
                          eligible_indices <- setdiff(eligible_indices, selected_indices_global)
                          
                          if (length(eligible_indices) == 0) {
                            break
                          }
                        }
                      }
                    } else {
                      break  # Salir si no hay bins con prioridad
                    }
                    
                    temp_discrepancy <- calculate_discrepancy(selected_surfaces, target_hist, bins, logaritmic)
                  }
                  
                  # Calcular la discrepancia
                  final_discrepancy <- calculate_discrepancy(selected_surfaces, target_hist, bins, logaritmic)
                  
                  # Verificar si es un resultado válido
                  if (final_discrepancy < best_discrepancy) {
                    best_selection <- list(
                      selected_surfaces = selected_surfaces,
                      surface_index = surface_index,
                      total_surface = total_surface,
                      final_discrepancy = final_discrepancy
                    )
                    best_discrepancy <- final_discrepancy
                  }
                  
                  best_selection
                  
                }, error = function(e) {
                  cat("Error en la iteración", i, ":", conditionMessage(e), "\n")
                  next
                })
                
                # Almacenar el resultado de la iteración
                if (!is.null(result)) {
                  return(result)
                }
              }
              
              # Retornar la mejor selección de este bloque
              best_selection
            }
  })
  
  # Filtrar resultados nulos y seleccionar el mejor
  valid_results <- Filter(Negate(is.null), results)
  
  if (length(valid_results) == 0) {
    return(NULL)  # Si no hay resultados válidos
  }
  
  # Seleccionar el mejor resultado
  best_result <- valid_results[[which.min(sapply(valid_results, function(x) x$final_discrepancy))]]
  
  # Verificar y eliminar duplicados
  if (!is.null(best_result)) {
    # Encontrar índices únicos manteniendo el orden original
    unique_indices <- !duplicated(best_result$surface_index)
    
    # Actualizar todos los componentes del resultado
    best_result$surface_index <- best_result$surface_index[unique_indices]
    best_result$selected_surfaces <- best_result$selected_surfaces[unique_indices]
    best_result$total_surface <- sum(best_result$selected_surfaces)
    
    # Recalcular la discrepancia final con los eventos únicos
    best_result$final_discrepancy <- calculate_discrepancy(
      best_result$selected_surfaces, 
      target_hist, 
      bins,
      logaritmic
    )
  }
  
  return(best_result)
}

# Function to remove duplicated perimeters
cleanse_duplicates <- function(candidates){
  
  candidate_surfaces <- candidates
  
  duplicate_indices <- list(-1)
  
  while(length(duplicate_indices)>0){
    # Comprobar si las geometrías están repetidas
    duplicates <- st_equals(candidate_surfaces)
    
    # Mostrar índices de geometrías duplicadas
    duplicate_indices <- which(sapply(duplicates, length) > 1)
    
    # Imprimir los resultados
    if (length(duplicate_indices) > 0) {
      print("Se encontraron polígonos duplicados en los siguientes índices:\n")
      print(duplicate_indices)
      candidate_surfaces <- candidate_surfaces[-duplicate_indices,]
    } else {
      print("No se encontraron polígonos duplicados.\n")
    }
    
  }
  
  return(candidate_surfaces)
}