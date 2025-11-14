IKPCA_CR <- function(separation_method = c("SVM", "LDA"),
                    data_train,
                    data_test,
                    label_train,
                    label_test,
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
  CR <- cbind(C,R)
  
  #-----------------------Best------------------------
  if (sigma_method == "best" && kernel_type == "gaussian") {
    ac <- 0
    best_sigma <- sigma
    for (i in seq(a, b, length.out = c)) {
      auxi <- IKPCA_CR_aux(separation_method = separation_method,
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
  #-------------------- KERNEL MATRIX CR --------------------
  K_cr <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- CR[i, ]
    for (j in seq_len(n)) {
      xj <- CR[j, ]
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
  H <- diag(n) - matrix(1, n, n) / n
  K_cr <- H %*% K_cr %*% H
  
  #-------------------- EIGEN DECOMPOSITION --------------------
  eig_cr <- eigen(K_cr, symmetric = TRUE)
  lambda_cr <- eig_cr$values
  V_cr <- eig_cr$vectors
  
  tol <- 1e-10
  keep <- which(lambda_cr > tol)
  lambda_cr <- lambda_cr[keep]
  V_cr <- V_cr[, keep, drop = FALSE]
  
  projected_data <- K_cr %*% V_cr %*% diag(1 / sqrt(lambda_cr))
  
  #-------------------- VARIANCE EXPLAINED --------------------
  expl <- cumsum(lambda_cr) / sum(lambda_cr)
  pcexp <- which(expl >= 0.8)[1]
  
  #-------------------- TEST DATA PROJECTION --------------------
  C_test <- as.matrix(as.data.frame(data_test@Centers))
  R_test <- as.matrix(as.data.frame(data_test@Ranges))/2
  n_test <- nrow(C_test)
  CR_test <- cbind(C_test,R_test)
  
  projected_data_test <- matrix(0, nrow = n_test, ncol = length(lambda_cr))
  
  for (d in 1:n_test) {
    k_cr_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- CR[i,]
      xd <- CR_test[d,]
      if (kernel_type == "gaussian") {
        k_cr_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_cr_new[i] <- sum(xi * xd)
      } else {
        k_cr_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    k_cr_new_centered <- k_cr_new - k_cr_mean - mean(k_cr_new) + K_cr_full_mean
    
    projected_data_test[d, ] <- t(k_cr_new_centered) %*% V_cr %*% diag(1 / sqrt(lambda_cr))
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
    eigenvalues = lambda_cr,
    var_exp = expl,
    pcexp = pcexp,
    principal_components = V_cr,
    projected_data = projected_data,
    projected_data_test = projected_data_test,
    pred_labels_pos = class_pos,
    conf_pos = conf_pos,
    metrics_pos = metrics_pos
  ))
}

IKPCA_CR_aux <- function(separation_method = c("SVM", "LDA"),
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
  CR <- cbind(C,R)
  
  #-------------------- KERNEL MATRIX CR --------------------
  K_cr <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- CR[i, ]
    for (j in seq_len(n)) {
      xj <- CR[j, ]
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
  H <- diag(n) - matrix(1, n, n) / n
  K_cr <- H %*% K_cr %*% H
  
  #-------------------- EIGEN DECOMPOSITION --------------------
  eig_cr <- eigen(K_cr, symmetric = TRUE)
  lambda_cr <- eig_cr$values
  V_cr <- eig_cr$vectors
  
  tol <- 1e-10
  keep <- which(lambda_cr > tol)
  lambda_cr <- lambda_cr[keep]
  V_cr <- V_cr[, keep, drop = FALSE]
  
  projected_data <- K_cr %*% V_cr %*% diag(1 / sqrt(lambda_cr))
  
  #-------------------- VARIANCE EXPLAINED --------------------
  expl <- cumsum(lambda_cr) / sum(lambda_cr)
  pcexp <- which(expl >= 0.8)[1]
  
  #-------------------- TEST DATA PROJECTION --------------------
  C_test <- as.matrix(as.data.frame(data_test@Centers))
  R_test <- as.matrix(as.data.frame(data_test@Ranges))/2
  n_test <- nrow(C_test)
  CR_test <- cbind(C_test,R_test)
  
  projected_data_test <- matrix(0, nrow = n_test, ncol = length(lambda_cr))
  
  for (d in 1:n_test) {
    k_cr_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- CR[i,]
      xd <- CR_test[d,]
      if (kernel_type == "gaussian") {
        k_cr_new[i] <- exp(-sum((xi - xd)^2) / (2 * sigma))
      } else if (kernel_type == "linear") {
        k_cr_new[i] <- sum(xi * xd)
      } else {
        k_cr_new[i] <- (sum(xi * xd) + coef0)^poly_degree
      }
    }
    
    k_cr_new_centered <- k_cr_new - k_cr_mean - mean(k_cr_new) + K_cr_full_mean
    
    projected_data_test[d, ] <- t(k_cr_new_centered) %*% V_cr %*% diag(1 / sqrt(lambda_cr))
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
