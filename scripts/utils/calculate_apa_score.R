#' Calculate APA score from a normalized APA matrix
#' 
#' @param matrix A numeric matrix representing normalized APA data
#' @param foreground_size Size of foreground region (default 1 for single pixel, 3 for 3x3)
#' @return APA score as numeric value
calculate_apa_score <- function(matrix, foreground_size = 1) {
  
  # Get matrix dimensions
  n_row <- nrow(matrix)
  n_col <- ncol(matrix)
  
  # Calculate center position
  center_row <- ceiling(n_row / 2)
  center_col <- ceiling(n_col / 2)
  
  # Extract foreground (center pixel or 3x3 region)
  if (foreground_size == 1) {
    foreground <- matrix[center_row, center_col]
  } else if (foreground_size == 3) {
    # 3x3 region centered on loop pixel
    row_range <- (center_row - 1):(center_row + 1)
    col_range <- (center_col - 1):(center_col + 1)
    foreground <- matrix[row_range, col_range]
  } else {
    stop("foreground_size must be 1 or 3")
  }
  
  # Extract background (corners)
  # Use 5x5 corners as in mariner and Figure 2 script
  corner_size <- 5
  
  # Top-left corner
  top_left <- matrix[1:corner_size, 1:corner_size]
  
  # Bottom-right corner
  bottom_right <- matrix[(n_row - corner_size + 1):n_row, 
                         (n_col - corner_size + 1):n_col]
  
  # Combine background values
  background <- c(as.vector(top_left), as.vector(bottom_right))
  
  # Remove NA values
  foreground <- foreground[!is.na(foreground)]
  background <- background[!is.na(background)]
  
  # Calculate APA score
  if (length(foreground) == 0 || length(background) == 0) {
    return(NA)
  }
  
  apa_score <- median(foreground + 1) / median(background + 1)
  
  return(apa_score)
}