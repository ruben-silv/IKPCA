interval_kpca <- function(separation_method = c("SVM", "LDA"),
                                    data_train,
                                    data_test,
                                    label_train,
                                    label_test,
                                    kernel_type = c("gaussian", "linear", "polynomial"),
                                    sigma_method = c("user", "medist", "best_medist"),
                                    grid_length=40,
                                    sigma = NULL,                                #sigma, in this code, corresponds to what is normally called sigma^2
                                    poly_degree = 4,
                                    coef0 = 1) {
  
  #Initial Setup---------------
  sigma_method <- match.arg(sigma_method)
  kernel_type <- match.arg(kernel_type)
  separation_method <- match.arg(separation_method)

  if (!inherits(data_train, "intData")) {
    stop("Data must be in class 'intData'.")
  }
  
  if (!inherits(data_test, "intData")) {
    stop("Data must be in class 'intData'.")
  }
  
  C <- as.matrix(as.data.frame(data_train@Centers))
  R <- as.matrix(as.data.frame(data_train@Ranges))
  n <- nrow(C)
  p <- ncol(C)
  
  E_U <- as.matrix(data_train@LatentParam[[2]])
  E_UU <- 4*as.matrix(data_train@LatentParam[[1]])
  SD_U <- sqrt(E_UU - E_U^2)
  latent_means <- diag(as.numeric(E_U / 2), p, p)
  latent_sd <- diag(as.numeric(SD_U / 2), p, p)
  
  #Data Transformation (G Matrix)-------------
  transformed_data <- matrix(0, n, 2*p)
  Gt <- rbind(
    cbind(diag(p), latent_means),
    cbind(matrix(0, p, p), latent_sd)
  )
  for (i in 1:n) {
    macro_vector <- matrix(c(C[i,], R[i,]), ncol = 1)
    transformed_data[i, ] <- t(Gt %*% macro_vector)
  }
  #Sigma Definition, if Gaussian------------
  if (kernel_type == "gaussian") {
    if (sigma_method == "medist" || sigma_method == "best_medist") {
      dist_matrix <- as.matrix(dist(transformed_data))
      sigma <- mean(dist_matrix[lower.tri(dist_matrix)])
      sigma_sd <- sd(dist_matrix[lower.tri(dist_matrix)])
    } else if (is.null(sigma)) {
      stop("Sigma must be given if kernel_type = 'gaussian' and sigma_method = 'user'.")
    }
  }
  
  #Best_Medist
  if (sigma_method == "best_medist" && kernel_type == "gaussian") {
    ac <- 0
    best_sigma <- sigma
    for (i in seq(max(0.01,sigma-2*sigma_sd), sigma+2*sigma_sd, length.out = grid_length)) {
      auxi <- ikpca_aux(separation_method = separation_method,
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
  
  # Find & Center Kernel Matrix-----------
  K <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- transformed_data[i, ]
    for (j in seq_len(n)) {
      xj <- transformed_data[j, ]
      if (kernel_type == "gaussian") {
        K[i, j] <- exp(- sum((xi - xj)^2 / (2 * sigma)))
      } else if (kernel_type == "linear") {
        K[i, j] <- sum(xi * xj)
      } else if (kernel_type == "polynomial") {
        K[i, j] <- (sum(xi * xj) + coef0)^poly_degree
      }
    }
  }
  
  K_full_mean <- mean(K)
  k_mean <- colMeans(K)
  ones_n <- matrix(1, n, n) / n
  K <- (diag(n) - ones_n) %*% K %*% (diag(n) - ones_n)
  
  #Data projection--------------
  eig <- eigen(K, symmetric = TRUE)
  lambda <- eig$values
  V <- eig$vectors
  projected_data <- K %*% V %*% diag(1/sqrt(lambda))
  
  #80% definition--------------
  expl <- cumsum(lambda) / sum(lambda)
  pcexp <- which(expl >= 0.8)[1]
  
  #Data Test Projection---------------
    #Initial Setup
    C_test <- as.matrix(as.data.frame(data_test@Centers))
    R_test <- as.matrix(as.data.frame(data_test@Ranges))
    n_test <- nrow(C_test)
    p_test <- ncol(C_test)

  
    #Data Test Transformation (G Matrix)
    transformed_data_test <- matrix(0, n_test, 2*p_test)
    Gt <- rbind(
      cbind(diag(p), latent_means),
      cbind(matrix(0, p, p), latent_sd)
    )
    for (i in 1:n_test) {
      macro_vector <- matrix(c(C_test[i,], R_test[i,]), ncol = 1)
      transformed_data_test[i, ] <- t(Gt %*% macro_vector)
    }
  
    #Definição K_new
    projected_data_test <- matrix(0, nrow = n_test, ncol = n)
    for (d in 1:n_test) {
      
      k_new <- numeric(n)
      
      for (i in 1:n) {
        xi <- transformed_data[i,]
        xd <- transformed_data_test[d,]
        if (kernel_type == "gaussian") {
          norm_diff <- sum((xi - xd)^2)
          k_new[i] <- exp(-norm_diff / (2 * sigma))
        } else if (kernel_type == "linear") {
          k_new[i] <- sum(xi * xd)
        } else if (kernel_type == "polynomial") {
          k_new[i] <- (sum(xi * xd) + coef0)^poly_degree
        }
      }
      
      
      k_new_centered <- k_new - k_mean - mean(k_new) + K_full_mean
      projected_data_test[d, ] <- (t(k_new_centered) %*% V %*% diag(1 / sqrt(lambda)))
    }
  
  #SVM vs LDA------------------
    if (separation_method == "SVM") {
      
      svm_model <- svm(x = transformed_data, y = as.factor(label_train), kernel = "linear")
      svm_pred <- predict(svm_model, transformed_data_test)
      
      all_labels <- sort(unique(c(label_test, svm_pred)))
      label_test <- factor(label_test, levels = all_labels)
      label_pred <- factor(svm_pred, levels = all_labels)
      
      conf_pre <- confusionMatrix(label_pred, label_test)
      metrics_pre <- conf_pre$byClass
      
      svm_models <- svm(x = projected_data[,1:pcexp], y = as.factor(label_train), kernel = "linear")
      svm_predi <- predict(svm_models, projected_data_test[,1:pcexp])
      
      all_labels <- sort(unique(c(label_test, svm_predi)))
      label_test <- factor(label_test, levels = all_labels)
      label_predi <- factor(svm_predi, levels = all_labels)
      
      conf_pos <- confusionMatrix(label_predi, label_test)
      metrics_pos <- conf_pos$byClass
      
      class_pre <- svm_pred 
      class_pos <- svm_predi
      } else if (separation_method == "LDA") {
        
        lda_model <- lda(x = transformed_data, grouping = as.factor(label_train))
        lda_pred <- predict(lda_model, transformed_data_test)
        
        all_labels <- sort(unique(c(label_test, lda_pred$class)))
        label_test <- factor(label_test, levels = all_labels)
        label_pred <- factor(lda_pred$class, levels = all_labels)
        
        conf_pre <- confusionMatrix(label_pred, label_test)
        metrics_pre <- conf_pre$byClass
        
        lda_models <- lda(x = projected_data[,1:pcexp], grouping = as.factor(label_train))
        lda_predi <- predict(lda_models,projected_data_test[,1:pcexp])
        
        all_labels <- sort(unique(c(label_test, lda_predi$class)))
        label_test <- factor(label_test, levels = all_labels)
        label_predi <- factor(lda_predi$class, levels = all_labels)
        
        conf_pos <- confusionMatrix(label_predi, label_test)
        metrics_pos <- conf_pos$byClass
      
        class_pre <- lda_pred$class 
        class_pos <- lda_predi$class
        }

  #Projection Plot------------
  classif <- as.numeric(class_pos)
  classif <- as.factor(classif)
  
  palette_colors <- rainbow(length(unique(classif)))
  if (any(classif == 0)) palette_colors[classif == 0] <- "black"
  
  plot(
    projected_data_test[,1], projected_data_test[,2],
    main = paste("Interval Kernel Principal Components Analysis -", kernel_type),
    xlab = "Principal Component 1",
    ylab = "Principal Component 2",
    col = palette_colors[classif],
    pch = 19
  )
  
  
  return(list(  
    sigma=sigma,
    transformed_data = transformed_data,
    transformed_data_test = transformed_data_test,
    kernel = kernel_type,
    eigenvalues = lambda,
    var_exp = expl ,
    pcexp=pcexp,
    principal_components = V,
    projected_data = projected_data,
    projected_data_test = projected_data_test,
    pred_labels_pre = class_pre,
    pred_labels_pos = class_pos,
    conf_pre = conf_pre,
    metrics_pre = metrics_pre,
    conf_pos = conf_pos,
    metrics_pos = metrics_pos
  ))
}

ikpca_aux <- function(separation_method = c("SVM", "LDA"),
                      data_train,
                      data_test,
                      label_train,
                      label_test,
                      kernel_type = c("gaussian"),
                      sigma_method = c("user"),
                      sigma = NULL) {
  
  sigma_method <- match.arg(sigma_method)
  kernel_type <- match.arg(kernel_type)
  separation_method <- match.arg(separation_method)
  
  C <- as.matrix(data_train@Centers)
  R <- as.matrix(data_train@Ranges)
  n <- nrow(C)
  p <- ncol(C)
  
  E_U <- as.numeric(data_train@LatentParam[[2]])
  E_UU <- 4*as.numeric(data_train@LatentParam[[1]])
  SD_U <- sqrt(E_UU - E_U^2)
  latent_means <- diag(as.numeric(E_U / 2), p, p)
  latent_sd <- diag(as.numeric(SD_U / 2), p, p)
  
  transformed_data <- matrix(0, n, 2*p)
  Gt <- rbind(
    cbind(diag(p), latent_means),
    cbind(matrix(0, p, p), latent_sd)
  )
  for (i in 1:n) {
    macro_vector <- matrix(c(C[i,], R[i,]), ncol = 1)
    transformed_data[i, ] <- t(Gt %*% macro_vector)
  }

  K <- matrix(0, n, n)
  for (i in seq_len(n)) {
    xi <- transformed_data[i, ]
    for (j in seq_len(n)) {
      xj <- transformed_data[j, ]
      if (kernel_type == "gaussian") {
        K[i, j] <- exp(- sum((xi - xj)^2 / (2 * sigma)))
      }
    }
  }
  
  K_full_mean <- mean(K)
  k_mean <- colMeans(K)
  ones_n <- matrix(1, n, n) / n
  K <- (diag(n) - ones_n) %*% K %*% (diag(n) - ones_n)

  eig <- eigen(K, symmetric = TRUE)
  lambda <- eig$values
  V <- eig$vectors
  projected_data <- K %*% V %*% diag(1/sqrt(lambda))

  expl <- cumsum(lambda) / sum(lambda)
  pcexp <- which(expl >= 0.8)[1]

  C_test <- as.matrix(data_test@Centers)
  R_test <- as.matrix(data_test@Ranges)
  n_test <- nrow(C_test)
  p_test <- ncol(C_test)
  
  transformed_data_test <- matrix(0, n_test, 2*p_test)
  for (i in 1:n_test) {
    Gt <- rbind(
      cbind(diag(p), latent_means),
      cbind(matrix(0, p, p), latent_sd)
    )
    macro_vector <- matrix(c(C_test[i,], R_test[i,]), ncol = 1)
    transformed_data_test[i, ] <- t(Gt %*% macro_vector)
  }
  
  projected_data_test <- matrix(0, nrow = n_test, ncol = n)
  for (d in 1:n_test) {
    
    k_new <- numeric(n)
    
    for (i in 1:n) {
      xi <- transformed_data[i,]
      xd <- transformed_data_test[d,]
      if (kernel_type == "gaussian") {
        norm_diff <- sum((xi - xd)^2)
        k_new[i] <- exp(-norm_diff / (2 * sigma))
      }
    }
    
    k_new_centered <- k_new - k_mean - mean(k_new) + K_full_mean
    projected_data_test[d, ] <- (t(k_new_centered) %*% V %*% diag(1 / sqrt(lambda)))
  }

  if (separation_method == "SVM") {
    
    svm_model <- svm(x = transformed_data, y = as.factor(label_train), kernel = "linear")
    svm_pred <- predict(svm_model, transformed_data_test)
    
    all_labels <- sort(unique(c(label_test, svm_pred)))
    label_test <- factor(label_test, levels = all_labels)
    label_pred <- factor(svm_pred, levels = all_labels)
    
    conf_pre <- confusionMatrix(label_pred, label_test)
    metrics_pre <- conf_pre$byClass
    
    svm_models <- svm(x = projected_data[,1:pcexp], y = as.factor(label_train), kernel = "linear")
    svm_predi <- predict(svm_models, projected_data_test[,1:pcexp])
    
    all_labels <- sort(unique(c(label_test, svm_predi)))
    label_test <- factor(label_test, levels = all_labels)
    label_predi <- factor(svm_predi, levels = all_labels)
    
    conf_pos <- confusionMatrix(label_predi, label_test)
    metrics_pos <- conf_pos$byClass
    
    class_pre <- svm_pred 
    class_pos <- svm_predi
  } else if (separation_method == "LDA") {
    
    lda_model <- lda(x = transformed_data, grouping = as.factor(label_train))
    lda_pred <- predict(lda_model, transformed_data_test)
    
    all_labels <- sort(unique(c(label_test, lda_pred$class)))
    label_test <- factor(label_test, levels = all_labels)
    label_pred <- factor(lda_pred$class, levels = all_labels)
    
    conf_pre <- confusionMatrix(label_pred, label_test)
    metrics_pre <- conf_pre$byClass
    
    lda_models <- lda(x = projected_data[,1:pcexp], grouping = as.factor(label_train))
    lda_predi <- predict(lda_models,projected_data_test[,1:pcexp])
    
    all_labels <- sort(unique(c(label_test, lda_predi$class)))
    label_test <- factor(label_test, levels = all_labels)
    label_predi <- factor(lda_predi$class, levels = all_labels)
    
    conf_pos <- confusionMatrix(label_predi, label_test)
    metrics_pos <- conf_pos$byClass
    
    class_pre <- lda_pred$class 
    class_pos <- lda_predi$class
  }
  
  return(conf_pos$overall[1])
}