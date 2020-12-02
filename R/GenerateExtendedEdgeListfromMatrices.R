#' Generate Extended Edges List that is Compatible with MuxViz
#'
#' @param Matrix_list List of matrices to be included in the multilayer network
#'
#' @param Layers Number of Layers in Network
#' @param Nodes Number of Nodes in Network
#' @param Threshold_percentile Percentile at which to threshold each matrix
#' @param doOMST Perform OMST thresholding
#' @param ... If stucturally constrained, Bin_struct = Bin_struct (case sensitive)
#'
#' @return Returns Extended Edges List to be used with muxViz
#' @export
#'
#' @importFrom gtools combinations
#' @importFrom pracma repmat
#'
#' @examples
#' \dontrun{
#' Ext_edges <- GenerateExtendedEdgeListfromMatrices <- function(Matrix_list, Layers, Nodes, Threshold_percentile = 0,  doOMST = FALSE, Bin_struct = Bin_struct)
#' }
#'
GenerateExtendedEdgeListfromMatrices <- function(Matrix_list, Layers, Nodes, Threshold_percentile = 0,  doOMST = FALSE, ...){
  combs <- combinations(Layers, 2, repeats.allowed = T)
  A <- numeric(0)
  for(i in 1:length(Matrix_list)){
    matrix <- Matrix_list[[i]]

    start_layer <- combs[i, 1]
    end_layer <- combs[i, 2]
    matrix[is.nan(matrix)] <- 0


    args <- list(...)
    exist <- "Bin_struct" %in% names(args)
    if(exist == TRUE){
      diag(Bin_struct) <- 1
      matrix <- matrix*Bin_struct
      print("Constrained Matrix")
    }

    if(doOMST == TRUE){
      omst_matrix <- calc_omst_matlab(abs(matrix))
      omst_matrix[omst_matrix > 0] <- 1
      matrix <- omst_matrix * matrix
    }
    else{
      matrix <- matrix
    }

    input <- matrix(0, length(matrix), 5)
    input[,1] <- rep(c(1:Nodes), each = Nodes)
    input[,2] <- start_layer
    input[,3] <- repmat(matrix(1:Nodes, Nodes, 1), Nodes, 1)
    input[,4] <- end_layer
    input[,5] <- reshapeR(matrix, length(matrix), 1)

    Threshold <- quantile(abs(input[,5][abs(input[,5]) > 0]), probs = Threshold_percentile)
    input_thresh <- input[abs(input[,5]) > Threshold,]

    A <- rbind(A, input_thresh)
    print(paste("Finished Processing Matrix", i, sep=''))
  }
  A = A[complete.cases(A), ]
  return(A)
}
