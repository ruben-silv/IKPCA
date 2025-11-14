MR_KPCA <- function(separation_method = c("SVM", "LDA"),
                   data_train,
                   data_test,
                   label_train,
                   label_test,
                   kernel_type = c("gaussian", "linear", "polynomial"),
                   sigma_method = c("user","best"),
                   a,b,c,
                   sigma = NULL,                                                 # sigma corresponds to sigma^2 (variance term)
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
  
  #-----------------------Best------------------------
  if (sigma_method == "best" && kernel_type == "gaussian") {
    ac <- 0
    best_sigma <- sigma
    for (i in seq(a, b, length.out = c)) {
      auxi <- MR_KPCA_aux(separation_method = separation_method,
                           data_train = data_train,
                           data_test = data_test,
                           label_train = label_train,
                           label_test = label_test,
                           kernel_type = "gaussian",
                           sigma_method = "user",
                           sigma = i)
      
      if (auxi>ac) {
        ac <- auxi
        best_sigma <- i
      }
    }
    sigma <- best_sigma
  }
  
  #-------------------- KERNEL MATRIX C --------------------
  K_c <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- C[i, ]
    for (j in seq_len(n)) {
      xj <- C[j, ]
      if (kernel_type == "gaussian") {
        K_c[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_c[i, j] <- sum(xi * xj)
      } else {
        K_c[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  K_c_full_mean <- mean(K_c)
  k_c_mean <- colMeans(K_c)
  H <- diag(n) - matrix(1, n, n) / n
  K_c <- H %*% K_c %*% H
  
  #-------------------- KERNEL MATRIX R --------------------
  K_r <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- R[i, ]
    for (j in seq_len(n)) {
      xj <- R[j, ]
      if (kernel_type == "gaussian") {
        K_r[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_r[i, j] <- sum(xi * xj)
      } else {
        K_r[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  K_r_full_mean <- mean(K_r)
  k_r_mean <- colMeans(K_r)
  K_r <- H %*% K_r %*% H
  
  #-------------------- KERNEL MATRIX CR --------------------
  K_cr <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- C[i, ]
    for (j in seq_len(n)) {
      xj <- R[j, ]
      if (kernel_type == "gaussian") {
        K_cr[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_cr[i, j] <- sum(xi * xj)
      } else {
        K_cr[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  K_cr_full_mean <- mean(K_cr)
  k_cr_mean <- colMeans(K_cr)
  K_cr <- H %*% K_cr %*% H
  
  #-------------------- SVD & A_phi --------------------
  svd_K_cr <- svd(K_cr)
  P_cr <- svd_K_cr$u    
  Q_cr <- svd_K_cr$v     
  A_cr <- Q_cr %*% t(P_cr)
  
  #-------------------- GLOBAL KERNEL MATRIX --------------------
  K_g <- K_c + K_r - 2 * K_cr %*% A_cr
  K_g_full_mean <- mean(K_g)
  k_g_mean <- colMeans(K_g)
  K_g <- H %*% K_g %*% H
  K_g <- 0.5 * (K_g + t(K_g))
  
  #-------------------- EIGEN DECOMPOSITION --------------------
  eig_g <- eigen(K_g, symmetric = TRUE)
  lambda_g <- eig_g$values
  V_g <- eig_g$vectors
  
  tol <- 1e-10
  keep <- which(lambda_g > tol)
  lambda_g <- lambda_g[keep]
  V_g <- V_g[, keep, drop = FALSE]
  
  projected_data <- K_g %*% V_g %*% diag(1 / sqrt(lambda_g))
  
  #-------------------- VARIANCE EXPLAINED --------------------
  expl <- cumsum(lambda_g) / sum(lambda_g)
  pcexp <- which(expl >= 0.8)[1]
  
  #-------------------- TEST DATA PROJECTION --------------------
  C_test <- as.matrix(as.data.frame(data_test@Centers))
  R_test <- as.matrix(as.data.frame(data_test@Ranges))/2
  n_test <- nrow(C_test)
  
  projected_data_test <- matrix(0, nrow = n_test, ncol = length(lambda_g))
  
  for (d in 1:n_test) {
    k_c_new <- numeric(n)
    k_r_new <- numeric(n)
    k_cr_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- C[i,]; xd <- C_test[d,]
      if (kernel_type == "gaussian") {
        k_c_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_c_new[i] <- sum(xi * xd)
      } else {
        k_c_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    for (i in 1:n) {
      xi <- R[i,]; xd <- R_test[d,]
      if (kernel_type == "gaussian") {
        k_r_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_r_new[i] <- sum(xi * xd)
      } else {
        k_r_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    for (i in 1:n) {
      xi <- R[i,]; xd <- C_test[d,]
      if (kernel_type == "gaussian") {
        k_cr_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_cr_new[i] <- sum(xi * xd)
      } else {
        k_cr_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    # Kernel global (centrar apenas no fim)
    k_g_new <- k_c_new + k_r_new - 2 * as.numeric(k_cr_new %*% A_cr)
    k_g_new_centered <- k_g_new - k_g_mean - mean(k_g_new) + K_g_full_mean
    
    projected_data_test[d, ] <- t(k_g_new_centered) %*% V_g %*% diag(1 / sqrt(lambda_g))
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
    kernel = kernel_type,
    eigenvalues = lambda_g,
    var_exp = expl,
    pcexp = pcexp,
    principal_components = V_g,
    projected_data = projected_data,
    projected_data_test = projected_data_test,
    pred_labels_pos = class_pos,
    conf_pos = conf_pos,
    metrics_pos = metrics_pos
  ))
}

MR_KPCA_aux <- function(separation_method = c("SVM", "LDA"),
                    data_train,
                    data_test,
                    label_train,
                    label_test,
                    kernel_type = c("gaussian"),
                    sigma_method = c("user"),
                    sigma = NULL) {
  
  #-------------------- INITIAL SETUP --------------------
  sigma_method <- match.arg(sigma_method)
  kernel_type <- match.arg(kernel_type)
  separation_method <- match.arg(separation_method)
  
  if (!inherits(data_train, "intData")) stop("Data must be in class 'intData'.")
  if (!inherits(data_test, "intData")) stop("Data must be in class 'intData'.")
  
  C <- as.matrix(as.data.frame(data_train@Centers))
  R <- as.matrix(as.data.frame(data_train@Ranges))/2
  n <- nrow(C)
  
  
  #-------------------- KERNEL MATRIX C --------------------
  K_c <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- C[i, ]
    for (j in seq_len(n)) {
      xj <- C[j, ]
      if (kernel_type == "gaussian") {
        K_c[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_c[i, j] <- sum(xi * xj)
      } else {
        K_c[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  K_c_full_mean <- mean(K_c)
  k_c_mean <- colMeans(K_c)
  H <- diag(n) - matrix(1, n, n) / n
  K_c <- H %*% K_c %*% H
  
  #-------------------- KERNEL MATRIX R --------------------
  K_r <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- R[i, ]
    for (j in seq_len(n)) {
      xj <- R[j, ]
      if (kernel_type == "gaussian") {
        K_r[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_r[i, j] <- sum(xi * xj)
      } else {
        K_r[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  K_r_full_mean <- mean(K_r)
  k_r_mean <- colMeans(K_r)
  K_r <- H %*% K_r %*% H
  
  #-------------------- KERNEL MATRIX CR --------------------
  K_cr <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- C[i, ]
    for (j in seq_len(n)) {
      xj <- R[j, ]
      if (kernel_type == "gaussian") {
        K_cr[i, j] <- exp(-sum((xi - xj)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        K_cr[i, j] <- sum(xi * xj)
      } else {
        K_cr[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  K_cr_full_mean <- mean(K_cr)
  k_cr_mean <- colMeans(K_cr)
  K_cr <- H %*% K_cr %*% H
  
  #-------------------- SVD & A_phi --------------------
  svd_K_cr <- svd(K_cr)
  P_cr <- svd_K_cr$u    
  Q_cr <- svd_K_cr$v     
  A_cr <- Q_cr %*% t(P_cr)
  
  #-------------------- GLOBAL KERNEL MATRIX --------------------
  K_g <- K_c + K_r - 2 * K_cr %*% A_cr
  K_g_full_mean <- mean(K_g)
  k_g_mean <- colMeans(K_g)
  K_g <- H %*% K_g %*% H
  K_g <- 0.5 * (K_g + t(K_g))
  
  #-------------------- EIGEN DECOMPOSITION --------------------
  eig_g <- eigen(K_g, symmetric = TRUE)
  lambda_g <- eig_g$values
  V_g <- eig_g$vectors
  
  tol <- 1e-10
  keep <- which(lambda_g > tol)
  lambda_g <- lambda_g[keep]
  V_g <- V_g[, keep, drop = FALSE]
  
  projected_data <- K_g %*% V_g %*% diag(1 / sqrt(lambda_g))
  
  #-------------------- VARIANCE EXPLAINED --------------------
  expl <- cumsum(lambda_g) / sum(lambda_g)
  pcexp <- which(expl >= 0.8)[1]
  
  #-------------------- TEST DATA PROJECTION --------------------
  C_test <- as.matrix(as.data.frame(data_test@Centers))
  R_test <- as.matrix(as.data.frame(data_test@Ranges))/2
  n_test <- nrow(C_test)
  
  projected_data_test <- matrix(0, nrow = n_test, ncol = length(lambda_g))
  
  for (d in 1:n_test) {
    k_c_new <- numeric(n)
    k_r_new <- numeric(n)
    k_cr_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- C[i,]; xd <- C_test[d,]
      if (kernel_type == "gaussian") {
        k_c_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_c_new[i] <- sum(xi * xd)
      } else {
        k_c_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    for (i in 1:n) {
      xi <- R[i,]; xd <- R_test[d,]
      if (kernel_type == "gaussian") {
        k_r_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_r_new[i] <- sum(xi * xd)
      } else {
        k_r_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    for (i in 1:n) {
      xi <- R[i,]; xd <- C_test[d,]
      if (kernel_type == "gaussian") {
        k_cr_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_cr_new[i] <- sum(xi * xd)
      } else {
        k_cr_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    # Kernel global (centrar apenas no fim)
    k_g_new <- k_c_new + k_r_new - 2 * as.numeric(k_cr_new %*% A_cr)
    k_g_new_centered <- k_g_new - k_g_mean - mean(k_g_new) + K_g_full_mean
    
    projected_data_test[d, ] <- t(k_g_new_centered) %*% V_g %*% diag(1 / sqrt(lambda_g))
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

