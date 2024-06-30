# FUNCTION MY CODE
normalize_tan <- function(X) {
  
  #' For each sample, values of every variable are divided by the total sum of
  #' all feature intensity values measured in that sample (\code{NA} values
  #' ignored by default), before multiplication by 100; the unit is \%.

  # X is your data matrix, where samples are in rows and variables are in columns.
  # X.norm is the data matrix after normalization; the order of the samples is the same as in the original data X.
  
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  
  # Calculate the mean of each column
  sums <- as.numeric(apply(X,1,sum))
  
  for(i in 1:nrow(X)) {
    X.norm[i, ] <- as.numeric((X[i, ] / sums[i]) * 100)
  }
  return(X.norm)
}



pqn <- function(X) {
  # X is your data matrix, where samples are in rows and variables are in columns.
  # X.norm is the data matrix after normalization; the order of the samples is the same as in the original data X.
  X.norm <- matrix(nrow = nrow(X), ncol = ncol(X))
  colnames(X.norm) <- colnames(X)
  rownames(X.norm) <- rownames(X)
  
  # Calculate the mean of each column
  mX <- as.numeric(apply(X,2,mean))
  
  # Choose non zero variables
  #ii <- which(mX > 0.001 * max(mX))
  
  # Normalize the data matrix
  for(i in 1:nrow(X)) {
    # Exception:
    if(median(as.numeric(X[i, ] / mX)) == 0) {
      tmp <- sort(as.numeric(X[i, ] / mX))
      X.norm[i, ] <- as.numeric(X[i, ] / tmp[which(tmp != 0)][1])
    }
    else {
      X.norm[i, ] <- as.numeric(X[i, ] / median(as.numeric(X[i, ] / mX)))
    }
  }
  
  return(X.norm)
}
# PQN CODE FROM AGI:
# pqn_agi <- function(X) {
#   # Get the dimensions of the data matrix
#   dims <- dim(X)
#   m <- dims[1]
#   n <- dims[2]
#   
#   # Calculate the mean of each column
#   mx <- colMeans(X)
#   
#   # Identify columns with mean greater than 0.001 times the maximum mean
#   ii <- which(mx > 0.001 * max(mx)) # to reduce the effect of zeros
#   
#   # Normalize the data matrix
#   rr <- X[, ii] / matrix(rep(mx[ii], each = m), nrow = m, byrow = FALSE)
#   
#   # Calculate the median for each row
#   w <- apply(rr, 1, median)
#   
#   ###... do sth with which(w == 0)
#   
#   # Normalize the original data matrix
#   Xnew <- X / matrix(rep(w, each = n), nrow = m, byrow = FALSE)
#   
#   # Return the normalized data matrix and the row-wise medians
#   return(list(Xnew = Xnew, w = w))
# }

### URF FUNCTIONS FROM AGI (from Matlab to R)
script_URF <- function(nr_itteration, real_data, nr_trees, nr_samples) {
  #Input arguments:
  # nr_itteration --> number of times the procedure should be repeated; e.g.50, 100
  # real_data --> data set (samples in row and variables in column
  # nr_trees --> number of trees in the RF; 1500 or more
  # nr_samples --> number of samples in a  terminal nodes, default is 1 but better to put higher number 5 or 8
  
  #Outputs:
  #pc --> pca score plot
  #pr--> % of variance per PC
  #mean proximity matrix
  require(randomForest)
  require(factoextra)
  
  proximity_all <- array(NA, dim = c(nrow(real_data), nrow(real_data), nr_itteration))
  for(i in 1:nr_itteration) {
    cat(i, "\n")
    result <- unsupervisedRF_proximityadd(real_data, nr_trees, nr_samples) #[b_tree_3_classes2,proximity_realdata]=unsupervisedRF_proximityadd(real_data,nr_trees,nr_samples);
    proximity_all[, ,i] <- result$proximity_realdata  #proximity_all(:,:,i)=proximity_realdata;
  }
  
  # PCA on the mean value of proximity
  mean_proximity <- apply(proximity_all, c(1,2), mean)
  
  pc <- prcomp(doublecentering(1-mean_proximity),
               center = FALSE,
               scale = FALSE)
  
  pc.eigval <- get_eigenvalue(pc)
  
  return(list(pc = pc,
              pr = pc.eigval,
              mean_proximity = mean_proximity))
}

unsupervisedRF_proximityadd <- function(real_data, nr_trees, nr_samples) {
  #the input of the function is data matrix and original class vector,  nr_trees (I usually do 1000) and nr_samples (i.e. number of samples in tree leaf;note that the default is 1 in Matlab for classification)
  #output is Random Forest tree, pca scores (pc) and PCA loadings (d), pr (% variance per pc), proximity matrix of real data
  #unsupervised RF
  
  # permute the columns of the real data
  artificial_data <- apply(real_data, 2, function(x) sample(x))
  
  Data_for_RF <- rbind(real_data, artificial_data)
  class_data_RF <- rep(1:2, each = nrow(real_data))
  
  # RF
  set.seed(123)  # set seed for reproducibility
  b_tree_3_classes2 <- randomForest(x = Data_for_RF, y = as.factor(class_data_RF),
                                    ntree = nr_trees, proximity = TRUE, 
                                    nodesize = nr_samples, importance = TRUE,
                                    mtry = round(sqrt(ncol(real_data))))
  
  proximity_realdata <- b_tree_3_classes2$proximity[1:nrow(real_data), 1:nrow(real_data)]
  
  # PCA analysis
  #pca_result <- PCA(1 - proximity_realdata)
  #pc <- pca_result$li
  
  return(list(b_tree_3_classes2 = b_tree_3_classes2, 
              proximity_realdata = proximity_realdata))
}

# Double centering function 
# the dissimilarity (proximity matrix)  is double mean-centered, 
# i.e. the equivalent of the usual mean-centering for a symmetric matrix: 
# the mean of rows and columns are subtracted, and the overall mean added. 
doublecentering <- function(data) {
  nr_objects <- nrow(data)
  ones_mat <- matrix(1, nrow = nr_objects, ncol = nr_objects)
  
  column_mean <- t(ones_mat * colMeans(data))
  row_mean <- ones_mat * rowMeans(data)
  overall_mean <- mean(data) * ones_mat
  
  data_center <- data - column_mean - row_mean + overall_mean
  
  return(data_center)
}