IKPCA_UL <- function(separation_method = c("SVM", "LDA"),
                     data_train,
                     data_test,
                     label_train,
                     label_test,
                     kernel_agg_param=0.5,
                     kernel_type = c("gaussian", "linear", "polynomial"),
                     sigma_method = c("user","best"),
                     a,b,c,
                     sigma = NULL,                                                # sigma corresponds to sigma^2 (variance term)
                     poly_degree = 4,
                     coef0 = 1) {
  
  #-------------------- INITIAL SETUP --------------------
  sigma_method <- match.arg(sigma_method)
  kernel_type <- match.arg(kernel_type)
  separation_method <- match.arg(separation_method)
  
  if (!inherits(data_train, "intData")) stop("Data must be in class 'intData'.")
  if (!inherits(data_test, "intData")) stop("Data must be in class 'intData'.")
  
  C <- as.matrix(as.data.frame(data_train@Centers))
  R <- as.matrix(as.data.frame(data_train@Ranges))/2
  n <- nrow(C)
  U <- C+R
  L <- C-R
  agg <- kernel_agg_param
  
  #-----------------------Best------------------------
  if (sigma_method == "best" && kernel_type == "gaussian") {
    ac <- 0
    best_sigma <- sigma
    best_inds <- 0
    for (inds in seq(0, 1, length.out = 5)) {
      for (i in seq(a, b, length.out = c)) {
        auxi <- IKPCA_UL_aux(separation_method = separation_method,
                             data_train = data_train,
                             data_test = data_test,
                             label_train = label_train,
                             label_test = label_test,
                             kernel_agg_param=inds,
                             kernel_type = "gaussian",
                             sigma_method = "user",
                             sigma = i)
      
        if (auxi>ac) {
          ac <- auxi
          best_sigma <- i
          best_inds <- inds
        }
      }
    }
    sigma <- best_sigma
    agg <- best_inds
  }
  if (kernel_type == "polynomial") {
    ac <- 0
    best_inds <- 0
    for (inds in seq(0, 1, length.out = 5)) {
       auxi <- IKPCA_UL_aux(separation_method = separation_method,
                           data_train = data_train,
                           data_test = data_test,
                           label_train = label_train,
                           label_test = label_test,
                           kernel_agg_param=inds,
                           kernel_type = "polynomial",
                           poly_degree = poly_degree,
                           coef0 = coef0)
        
      if (auxi>ac) {
        ac <- auxi
        best_inds <- inds
      }
      
    }
    agg <- best_inds
  }
  
  if (kernel_type == "linear") {
    ac <- 0
    best_inds <- 0
    for (inds in seq(0, 1, length.out = 5)) {
      auxi <- IKPCA_UL_aux(separation_method = separation_method,
                           data_train = data_train,
                           data_test = data_test,
                           label_train = label_train,
                           label_test = label_test,
                           kernel_agg_param=inds,
                           kernel_type = "linear")
      
      if (auxi>ac) {
        ac <- auxi
        best_inds <- inds
      }
      
    }
    agg <- best_inds
  }
  
  
  #-------------------- KERNEL MATRIX U --------------------
  K_u <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- U[i, ]
    for (j in seq_len(n)) {
      xj <- U[j, ]
      if (kernel_type == "gaussian") {
        K_u[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_u[i, j] <- sum(xi * xj)
      } else {
        K_u[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }

  #-------------------- KERNEL MATRIX L --------------------
  K_l <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- L[i, ]
    for (j in seq_len(n)) {
      xj <- L[j, ]
      if (kernel_type == "gaussian") {
        K_l[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_l[i, j] <- sum(xi * xj)
      } else {
        K_l[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }

  #-------------------- KERNEL MATRIX UL --------------------
  K_ul <- agg*K_u + (1-agg)*K_l
  
  K_ul_full_mean <- mean(K_ul)
  k_ul_mean <- colMeans(K_ul)
  H <- diag(n) - matrix(1, n, n) / n
  K_ul <- H %*% K_ul %*% H
  #-------------------- EIGEN DECOMPOSITION --------------------
  eig_ul <- eigen(K_ul, symmetric = TRUE)
  lambda_ul <- eig_ul$values
  V_ul <- eig_ul$vectors
  
  tol <- 1e-10
  keep <- which(lambda_ul > tol)
  lambda_ul <- lambda_ul[keep]
  V_ul <- V_ul[, keep, drop = FALSE]
  
  projected_data <- K_ul %*% V_ul %*% diag(1 / sqrt(lambda_ul))
  
  #-------------------- VARIANCE EXPLAINED --------------------
  expl <- cumsum(lambda_ul) / sum(lambda_ul)
  pcexp <- which(expl >= 0.8)[1]
  
  #-------------------- TEST DATA PROJECTION --------------------
  C_test <- as.matrix(as.data.frame(data_test@Centers))
  R_test <- as.matrix(as.data.frame(data_test@Ranges))/2
  n_test <- nrow(C_test)
  U_test <- C_test+R_test
  L_test <- C_test-R_test
  
  projected_data_test <- matrix(0, nrow = n_test, ncol = length(lambda_ul))
 
   #-------------------- TEST DATA PROJECTION--------------------
  for (d in 1:n_test) {
    k_u_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- U[i,]
      xd <- U_test[d,]
      if (kernel_type == "gaussian") {
        k_u_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_u_new[i] <- sum(xi * xd)
      } else {
        k_u_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    k_l_new <- numeric(n)
      
    for (i in 1:n) {
      xi <- L[i,]
      xd <- L_test[d,]
      if (kernel_type == "gaussian") {
        k_l_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_l_new[i] <- sum(xi * xd)
      } else {
        k_l_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    k_ul_new <- agg*k_u_new + (1-agg)*k_l_new
      
    k_ul_new_centered <- k_ul_new - k_ul_mean - mean(k_ul_new) + K_ul_full_mean
    
    projected_data_test[d, ] <- t(k_ul_new_centered) %*% V_ul %*% diag(1 / sqrt(lambda_ul))
  } 
  
  #-------------------- CLASSIFICATION --------------------
  if (separation_method == "SVM") {
    svm_models <- e1071::svm(x = projected_data[, 1:pcexp, drop = FALSE],
                             y = as.factor(label_train), kernel = "linear")
    svm_predi <- predict(svm_models, projected_data_test[, 1:pcexp, drop = FALSE])
    
    all_labels <- sort(unique(c(label_test, svm_predi)))
    label_test <- factor(label_test, levels = all_labels)
    label_predi <- factor(svm_predi, levels = all_labels)
    
    conf_pos <- caret::confusionMatrix(label_predi, label_test)
    metrics_pos <- conf_pos$byClass
    class_pos <- svm_predi
    
  } else if (separation_method == "LDA") {
    lda_models <- MASS::lda(x = projected_data[, 1:pcexp, drop = FALSE],
                            grouping = as.factor(label_train))
    lda_predi <- predict(lda_models, projected_data_test[, 1:pcexp, drop = FALSE])
    
    all_labels <- sort(unique(c(label_test, lda_predi$class)))
    label_test <- factor(label_test, levels = all_labels)
    label_predi <- factor(lda_predi$class, levels = all_labels)
    
    conf_pos <- caret::confusionMatrix(label_predi, label_test)
    metrics_pos <- conf_pos$byClass
    class_pos <- lda_predi$class
  }
  
  #-------------------- PLOT --------------------
  classif <- as.numeric(class_pos)
  classif <- as.factor(classif)
  palette_colors <- rainbow(length(unique(classif)))
  
  plot( projected_data_test[,1],
        projected_data_test[,2],
        main = "Interval Kernel Principal Components Analysis",
        xlab = "Kernel Principal Component 1",
        ylab = "Kernel Principal Component 2",
        col = palette_colors[classif], pch = 19 )
  #-------------------- RETURN --------------------
  return(list(  
    sigma = sigma,
    kernel_agg_para=agg,
    kernel = kernel_type,
    eigenvalues = lambda_ul,
    var_exp = expl,
    pcexp = pcexp,
    principal_components = V_ul,
    projected_data = projected_data,
    projected_data_test = projected_data_test,
    pred_labels_pos = class_pos,
    conf_pos = conf_pos,
    metrics_pos = metrics_pos
  ))
}

IKPCA_UL_aux <- function(separation_method = c("SVM", "LDA"),
                     data_train,
                     data_test,
                     label_train,
                     label_test,
                     kernel_agg_param=0.5,
                     kernel_type = c("gaussian","linear","polynomial"),
                     sigma_method = c("user"),
                     sigma = NULL,                                                # sigma corresponds to sigma^2 (variance term)
                     poly_degree = 4,
                     coef0 = 1) {
  
  #-------------------- INITIAL SETUP --------------------
  sigma_method <- match.arg(sigma_method)
  kernel_type <- match.arg(kernel_type)
  separation_method <- match.arg(separation_method)
  
  if (!inherits(data_train, "intData")) stop("Data must be in class 'intData'.")
  if (!inherits(data_test, "intData")) stop("Data must be in class 'intData'.")
  
  C <- as.matrix(as.data.frame(data_train@Centers))
  R <- as.matrix(as.data.frame(data_train@Ranges))/2
  n <- nrow(C)
  U <- C+R
  L <- C-R
  agg <- kernel_agg_param
  
  #-------------------- KERNEL MATRIX U --------------------
  K_u <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- U[i, ]
    for (j in seq_len(n)) {
      xj <- U[j, ]
      if (kernel_type == "gaussian") {
        K_u[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_u[i, j] <- sum(xi * xj)
      } else {
        K_u[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  
  #-------------------- KERNEL MATRIX L --------------------
  K_l <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- L[i, ]
    for (j in seq_len(n)) {
      xj <- L[j, ]
      if (kernel_type == "gaussian") {
        K_l[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_l[i, j] <- sum(xi * xj)
      } else {
        K_l[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  
  #-------------------- KERNEL MATRIX UL --------------------
  K_ul <- agg*K_u + (1-agg)*K_l
  
  K_ul_full_mean <- mean(K_ul)
  k_ul_mean <- colMeans(K_ul)
  H <- diag(n) - matrix(1, n, n) / n
  K_ul <- H %*% K_ul %*% H
  #-------------------- EIGEN DECOMPOSITION --------------------
  eig_ul <- eigen(K_ul, symmetric = TRUE)
  lambda_ul <- eig_ul$values
  V_ul <- eig_ul$vectors
  
  tol <- 1e-10
  keep <- which(lambda_ul > tol)
  lambda_ul <- lambda_ul[keep]
  V_ul <- V_ul[, keep, drop = FALSE]
  
  projected_data <- K_ul %*% V_ul %*% diag(1 / sqrt(lambda_ul))
  
  #-------------------- VARIANCE EXPLAINED --------------------
  expl <- cumsum(lambda_ul) / sum(lambda_ul)
  pcexp <- which(expl >= 0.8)[1]
  
  #-------------------- TEST DATA PROJECTION --------------------
  C_test <- as.matrix(as.data.frame(data_test@Centers))
  R_test <- as.matrix(as.data.frame(data_test@Ranges))/2
  n_test <- nrow(C_test)
  U_test <- C_test+R_test
  L_test <- C_test-R_test
  
  projected_data_test <- matrix(0, nrow = n_test, ncol = length(lambda_ul))
  
  #-------------------- TEST DATA PROJECTION--------------------
  for (d in 1:n_test) {
    k_u_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- U[i,]
      xd <- U_test[d,]
      if (kernel_type == "gaussian") {
        k_u_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_u_new[i] <- sum(xi * xd)
      } else {
        k_u_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    k_l_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- L[i,]
      xd <- L_test[d,]
      if (kernel_type == "gaussian") {
        k_l_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_l_new[i] <- sum(xi * xd)
      } else {
        k_l_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    k_ul_new <- agg*k_u_new + (1-agg)*k_l_new
    
    k_ul_new_centered <- k_ul_new - k_ul_mean - mean(k_ul_new) + K_ul_full_mean
    
    projected_data_test[d, ] <- t(k_ul_new_centered) %*% V_ul %*% diag(1 / sqrt(lambda_ul))
  } 
  
  #-------------------- CLASSIFICATION --------------------
  if (separation_method == "SVM") {
    svm_models <- e1071::svm(x = projected_data[, 1:pcexp, drop = FALSE],
                             y = as.factor(label_train), kernel = "linear")
    svm_predi <- predict(svm_models, projected_data_test[, 1:pcexp, drop = FALSE])
    
    all_labels <- sort(unique(c(label_test, svm_predi)))
    label_test <- factor(label_test, levels = all_labels)
    label_predi <- factor(svm_predi, levels = all_labels)
    
    conf_pos <- caret::confusionMatrix(label_predi, label_test)
    metrics_pos <- conf_pos$byClass
    class_pos <- svm_predi
    
  } else if (separation_method == "LDA") {
    lda_models <- MASS::lda(x = projected_data[, 1:pcexp, drop = FALSE],
                            grouping = as.factor(label_train))
    lda_predi <- predict(lda_models, projected_data_test[, 1:pcexp, drop = FALSE])
    
    all_labels <- sort(unique(c(label_test, lda_predi$class)))
    label_test <- factor(label_test, levels = all_labels)
    label_predi <- factor(lda_predi$class, levels = all_labels)
    
    conf_pos <- caret::confusionMatrix(label_predi, label_test)
    metrics_pos <- conf_pos$byClass
    class_pos <- lda_predi$class
  }
  
  #-------------------- RETURN --------------------
  return(conf_pos$overall[1])
}
