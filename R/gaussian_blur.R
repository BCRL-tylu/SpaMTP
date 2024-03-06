#' @gaussian_blur 
#' @description This is an application that apply guassian blur to a given image matrix
#' @param matrix is a 2-D matrix like object
#' @param sigma is a numerical value indicates the weight of gussian, the higher the sigma, the higher weight will be applied to kernel edges
#' @param kernel_size is an odd integer value, indicate the side length of the kernel
#' @param return_vector is a boolean vector indicates that whether a matrix if returned (False) or a vector is returned (True)

#' @return Returns the gaussian blurred matrix

gaussian_blur <- function(matrix, sigma = 1,
                          kernel_size = 3,
                          return_vector = F) {
  # Create a kernel for Gaussian blur
  if((kernel_size %% 2) == 0){
    stop("Please use odd kernel size, even number can leads to asymmetric blurring")
  }
  size <- floor(kernel_size/2)
  kernel <- outer(
    seq(-size, size),
    seq(-size, size),
    function(x, y) {
      exp(-(x^2 + y^2) / (2 * as.numeric(sigma)^2)) / (2 * pi * as.numeric(sigma)^2)
    }
  )
  
  # Normalize the kernel
  kernel <- kernel / sum(kernel)
  
  # Apply convolution with the kernel
  blurred_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  for (i in 1:(nrow(matrix))) {
    # Define the boundaries for convolution
    row_index = (i-size):(i+size)
    if(i<=size){
      row_index[which(row_index<=0)]= abs(row_index[which(row_index<=0)]-1-i)
    }
    if(i>(nrow(matrix)-size)){
      row_index[which(row_index>nrow(matrix))]= 2*nrow(matrix) - row_index[which(row_index>nrow(matrix))] 
    }
    for (j in 1:(ncol(matrix))) {
      col_index = (j-size):(j+size)
      if(j<=size){
        col_index[which(col_index<=0)]= abs(col_index[which(col_index<=0)]-1-j) 
      }
      if(j>(ncol(matrix)-size)){
        col_index[which(col_index>ncol(matrix))]= 2*ncol(matrix) - col_index[which(col_index>ncol(matrix))] 
      }
      # Get the unprocessed kernel
      matrix_kernel = matrix[row_index,
                             col_index]
      # Convolution operation within boundaries
      
      blurred_matrix[i, j] <- sum(matrix_kernel* kernel)
    }
  }
  if(return_vector == F){
    return(blurred_matrix) 
  }else{
    return(as.vector(blurred_matrix))
  }
}
