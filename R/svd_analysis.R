#' @svd_analysis
#' @description  This is an application that use single value decomposition to evaluate the potential importance of specific m/z values by analyzing the singular values of the original features. It integrate multiple image processing techniques to filter the noise and improve effiency of SVD calculation.
#' @param mass_matrix is a matrix like object with each row corresponding to a pixel, and each column corresponding to a m/z value (Colnames = m/z values), support sparse/dense matrix,m each cell is the intensity of at the pixel; Note that the mass matrix should be cleaned against backgrounds and potential artifacts
#' @param width is a integer value indicates the width of single imaging data in each column of mass matrix,
#' @param height is a integer value indicates the height of single imaging data in each column of mass matrix,
#' @param plot is a boolean variable indicates whether a plot is included in the output.
#' @param number_to_plot an integer value only needed when plot = True, the number of m/z values, which explain the highest proportion of variances, to be included in plots.
#' @param randomised_svd a boolean variable indicates whether use randomized SVD(True) or conventional deterministic SVD(False), if matrix decomposition take too long consider set this to True
#' @param k is a integer value goes to rsvd::rsvd, specifies the target rank of the low-rank decomposition. Default is full rank which equal to the number of input features which can be really slow.
#' @param gaussian_blur is a boolean variable whether Gaussian blur will be applied before downsampling. A typical 1000 pixel *1000 pixel image take roughly 4 seconds to process, if you want to speed up, set 'use_parallel' as TRUE, and enter a numeric value 'num_cores' to indicate number of your aviliable cores
#' @param kernel_size is a integer value indicates the kernel size of Gaussian blur.
#' @param sigma is the weight of gussian blur,the Standard deviation of isotropic Gaussian smoothing kernel.
#' @param use_paralle is a boolean factor indicates whether to use parallel to boost performance
#' @param num_cores is a number indicates tha cores you want to use for the job, if not stated, it will use 'parallelly::availableCores()' as default
#' @param resampling_factor A numerical value bigger than 0 to indicate the extend you want to resample you images, for example if you want to resize a 200x200 image into 100x100, you should input 2 in this case, and 0.5 for the reverse case
#' @param variance_explained Is a numeric value between 0 and 1 which indicates the cut-off variance explained by singular vectors from svd analysis

#' @return characteristic_mz_matices: is a dataframe contains the most influential m/z values on the singular values 
#' @return svd: is single value decomposition analysis result, referred to base::svd 
#' @return new_width,new_height: the new dimension of data after transformation
#' @return kernel_size,resampling_factor,use_gaussian_blur,sigma: are the parameters that input by user



svd_analysis = function(mass_matrix,
                                  width,
                                  height,
                                  variance_explained = 0.99,
                                  plot = F,
                                  number_to_plot = 8,
                                  randomised_svd = F,
                                  k = NULL,
                                  kernel_size = 3,
                                  resampling_factor = NULL,
                                  sigma = 1,
                                  use_gaussian_blur = F,
                                  use_paralle = F,
                                  num_cores =NULL){
  require(Matrix)
  mass_matrix = Matrix(mass_matrix,
                       sparse = T)

  
  ##########################################################
  ###################Gaussian blur##########################
  ##########################################################
  
  if(use_gaussian_blur  == T & use_paralle == F){
    print("Running gaussian blur")
    require(utils)
    pb = txtProgressBar(min = 0, max = ncol(mass_matrix), initial = 0, style = 3) 
    blur_matrix = matrix(nrow = width*height)
    source("./network_based_pathway_ananlysis/gaussian_blur.R")
    for(i in 1:ncol(mass_matrix)){
      temp_mz_matrix = matrix(mass_matrix[,i],
                              ncol = width,
                              nrow = height, byrow = T)
      blured_temp = gaussian_blur(temp_mz_matrix, 
                                  sigma = sigma,
                                  kernel_size = kernel_size,
                                  return_vector = T)
      blur_matrix = cbind(blur_matrix,blured_temp)
      setTxtProgressBar(pb,i)
    }
    blur_matrix = blur_matrix[,-1]
    colnames(blur_matrix) =colnames(mass_matrix)
    close(pb)
    print("Guassian Blur finished")
  }
  
  if(use_gaussian_blur == T & use_paralle == T){
    print("Running gaussian blur on parallel cores")
    require(parallel)
    require(parallelly)
    print("Making cluster on your aviliable cores")
    if(is.null(num_cores)){
      n.cores = parallelly::availableCores()
    }else if(is.numeric(num_cores) == T){
      n.cores = as.integer(num_cores)
    }else{
      stop("Please enter correct number of cores 'num_cores'")
    }
    source("./network_based_pathway_ananlysis/gaussian_blur.R")
    clust <- makeCluster(n.cores)
    clusterExport(clust, varlist = c("mass_matrix",
                                     "gaussian_blur",
                                     "width",
                                     "height",
                                     "kernel_size",
                                     "sigma"), envir = environment())
    print("Running gaussian blur (Takes roughly 2 min on 16 cores processing 500 frames of 1000pix*1000pix matrices)")
    blur_list <- parLapply(clust, split(t(mass_matrix),seq_len(nrow(t(mass_matrix)))), function(x){
      temp_mz_matrix = matrix(x,
                              ncol = width,
                              nrow = height, byrow = T)
      blured_temp = gaussian_blur(temp_mz_matrix, 
                                  sigma = sigma,
                                  kernel_size = kernel_size,
                                  return_vector = T)
      return(blured_temp)
    })
    stopCluster(clust)
    gc()
    blur_matrix =do.call(cbind,blur_list)
    colnames(blur_matrix) =colnames(mass_matrix)
    print("Guassian Blur finished")
  }
  ##########################################################
  #####################Resampling###########################
  ##########################################################
  rescale <- function(x, newrange=range(x)){
    xrange <- range(x)
    mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
    newrange[1]+(x-xrange[1])*mfac
  }
  
  ResizeMat <- function(mat, ndim=dim(mat)){
    if(!require(fields)) stop("`fields` required.")
    
    # input object
    odim <- dim(mat)
    obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
    
    # output object
    ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
    ndim <- dim(ans)
    
    # rescaling
    ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
    loc <- ncord
    loc[,1] = rescale(ncord[,1], c(1,odim[1]))
    loc[,2] = rescale(ncord[,2], c(1,odim[2]))
    
    # interpolation
    ans[ncord] <- interp.surface(obj, loc)
    ans
  }
  
  if(use_gaussian_blur == T){
    blur_matrix = Matrix::Matrix(blur_matrix,
                                 sparse = T)
  }else{
    blur_matrix = mass_matrix
  }
  
  if(!is.null(resampling_factor)){
    print("Running matrix resampling")
    pb = txtProgressBar(min = 0, max = ncol(mass_matrix), initial = 0, style = 3) 
    if(!is.numeric(resampling_factor)){
      stop("Please enter correct resampling_factor")
    }
    new.width = as.integer(width/resampling_factor)
    new.height = as.integer(height/resampling_factor)
    
    resampled_mat = matrix(nrow =  new.height*new.width)
    for(i in 1:ncol(blur_matrix)){
      temp_mz_matrix = matrix(blur_matrix[,i],
                              ncol = width,
                              nrow = height, byrow = F)
      resampled_temp = ResizeMat(temp_mz_matrix, c(new.width,
                                                   new.height))
      resampled_mat = cbind(resampled_mat,as.vector(resampled_temp))
      setTxtProgressBar(pb,i)
    }
    close(pb)
    resampled_mat = resampled_mat[,-1]
    colnames(resampled_mat) = colnames(blur_matrix)
    print("Resampling finished!")
    gc()
  }
  
  # if(!is.null(resampling_factor) & use_paralle == T){
  #   if(!is.numeric(resampling_factor)){
  #     stop("Please enter correct resampling_factor")
  #   }
  #   new.width = as.integer(width/resampling_factor)
  #   new.height = as.integer(height/resampling_factor)
  #   print("Running matrix resampling")
  #   print("Making cluster on your aviliable cores")
  #   clust <- makeCluster(n.cores)
  #   clusterExport(clust, varlist = c("blur_matrix",
  #                                    "ResizeMat",
  #                                    "new.width",
  #                                    "new.height",
  #                                    "rescale",
  #                                    "width",
  #                                    "height"), envir = environment())
  #   print("Running matrix resampling on multiple cores")
  #   resampled_list <- parLapply(clust, split(t(blur_matrix),seq_len(nrow(t(blur_matrix)))), function(x){
  #     temp_mz_matrix = matrix(x,
  #                             ncol = width,
  #                             nrow = height, byrow = F)
  #     resampled_temp = ResizeMat(temp_mz_matrix, c(new.width,
  #                                                  new.height))
  #     return(as.vector(resampled_temp))
  #   })
  #   stopCluster(clust)
  #   gc()
  #   resampled_mat =do.call(cbind,resampled_list)
  #   colnames(resampled_mat) =colnames(blur_matrix)
  #   print("Resampling finished!")
  # }
  
  ###################################################
  #####################SVD###########################
  ###################################################
  if(is.null(resampling_factor)){
    resampled_mat = blur_matrix
    new.width  = width
    new.height = height
    sigma = NULL
  }
  
  if(randomised_svd == F){
    print("Running exact matrix single value decomposition")
    svd = svd(resampled_mat)
    print("SVD done")
  }else{
    require(rsvd)
    print("Running random matrix decomposition")
    if(is.null(k)){
      k = ncol(resampled_mat)
    }
    svd = rsvd(resampled_mat, k)
    print("SVD done")
  }
  print("Searching out most characteristic m/z which contribute most to the variance explanation")
  cumsum = cumsum(svd$d^2/sum(svd$d^2))
  # The retained part of the data
  retain = length(which(cumsum<= variance_explained))
  # column vectors V stores the composition of variance explained by each components
  retained_component = svd$v[,1:retain]
  # Get the distribution matrix of the characteristic m/z values
  order_matrix = data.frame()
  value_matrix = c()
  print("Evaluating the relative proportion of variance explained by individual m/z values")
  for(i in 1:retain){
    if(i ==1L){
      percentage_explain =  cumsum[i]
      linear_factor = abs(retained_component[,i])
      order_matrix = order(linear_factor*percentage_explain, decreasing = T)
      value_matrix = linear_factor*percentage_explain
    }else{
      percentage_explain =  cumsum[i]-cumsum[i-1]
      linear_factor = abs(retained_component[,i])
      order_matrix = cbind(order_matrix ,
                           order(linear_factor * percentage_explain, decreasing = T))
      value_matrix = value_matrix + linear_factor*percentage_explain
    }
  }
  if(plot == T){
    plot.new()
    print("Plotting figures")
    require(fields)
    require(spam)
    if (number_to_plot %% 2 == 0) {
      result <- 1:number_to_plot
    } else {
      result <- c(1:number_to_plot, 0)
    }
    # layout.matrix <- matrix(result, ncol = 2)
    # layout(mat = layout.matrix,
    #        heights = rep(1.5, times = nrow(layout.matrix)), # Heights of the two rows
    #        widths = c(2, 2))
    par(mfrow = c(ceiling(number_to_plot/2), 2))
    par(mar = c(2, 2, 1, 1))
    for(i in 1:number_to_plot){
      plot_matrix = matrix(resampled_mat[,order(value_matrix,decreasing = T)[i]],
                           nrow = new.height, ncol = new.width, byrow = T)
      image = fields::image.plot(plot_matrix, useRaster = T,
                                 main = paste0("m/z:",colnames(resampled_mat[,order(value_matrix,decreasing = T)])[i]),axes=FALSE, xlab = "", ylab = "")
      image
    }
  }
  
  returned_item = list(characteristic_mz_matices = resampled_mat[,order(value_matrix,decreasing = T)],
                       svd = svd,
                       randomised_svd = randomised_svd,
                       k = k,
                       kernel_size = kernel_size,
                       resampling_factor = resampling_factor,
                       new_width = new.width,
                       new_height = new.height,
                       sigma = sigma,
                       use_gaussian_blur = use_gaussian_blur)
  return(returned_item)
}
