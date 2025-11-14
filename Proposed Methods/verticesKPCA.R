generate_vertices <- function(centers, ranges) {
  
  centers <- as.matrix(centers)
  ranges  <- as.matrix(ranges)
  
  n <- nrow(centers)  
  p <- ncol(centers)  
  
  half_ranges <- ranges / 2
  
  mins <- centers - half_ranges
  maxs <- centers + half_ranges
  
  
  combos <- expand.grid(rep(list(c(0, 1)), p))
  combos_mat <- as.matrix(combos)
  num_verts <- 2^p
  
  all_vertices <- matrix(numeric(0), nrow = 0, ncol = p)
  
  for (i in seq_len(n)) {
    
    verts <- matrix(0, nrow = num_verts, ncol = p)
    
    for (k in seq_len(num_verts)) {
      for (j in seq_len(p)) {
        if (combos_mat[k, j] == 0) {
          verts[k, j] <- mins[i, j]
        } else {
          verts[k, j] <- maxs[i, j]
        }
      }
    }
    
    all_vertices <- rbind(all_vertices,verts)
  }
  
  return(all_vertices)
}

reconstruct_vertices_v <- function(M, p) {
  M <- as.matrix(M)
  block_size <- 2^p
  total_linhas <- nrow(M)
  
  if (total_linhas %% block_size != 0) {
    stop("Number of rows in M  is not a multiple of 2^p.")
  }
  
  nblocks <- total_linhas / block_size
  result_list <- matrix(numeric(0), nrow = 0, ncol = p)
  
  combos <- as.matrix(expand.grid(rep(list(c(0, 1)), p)))
  nverts <- nrow(combos)  
  
  for (i in seq_len(nblocks)) {
    idx_bloc <- ((i - 1) * block_size + 1):(i * block_size)
    bloco   <- M[idx_bloc, , drop = FALSE] 
    
    mins <- apply(bloco, 2, min)
    maxs <- apply(bloco, 2, max)
    
    verts <- matrix(0, nrow = nverts, ncol = p)
    for (j in seq_len(p)) {
      verts[, j] <- ifelse(combos[, j] == 0, mins[j], maxs[j])
    }
    
    result_list <- rbind(result_list,verts)
  }
  
  return(result_list)
}

vertices_kpca <- function(separation_method = c("SVM", "LDA"),
                                data_train,
                                data_test,
                                label_train,
                                label_test,
                                kernel_type = c("gaussian", "linear", "polynomial"),
                                sigma_method = c("user", "medist"),
                                grid_length=40,
                                sigma = NULL,
                                poly_degree = 4,
                                coef0 = 1) {
  
  #Setup Inicial---------------
  sigma_method <- match.arg(sigma_method)
  kernel_type <- match.arg(kernel_type)
  separation_method <- match.arg(separation_method)
  
  if (!inherits(data_train, "intData")) {
    stop("O conjunto de dados deve estar na classe 'intData'.")
  }
  
  if (!inherits(data_test, "intData")) {
    stop("O conjunto de dados deve estar na classe 'intData'.")
  }
  
  C <- as.matrix(data_train@Centers)
  R <- as.matrix(data_train@Ranges)
  n <- nrow(C)
  p <- ncol(C)


  #Gera os microdados---------
  
  micro <- generate_vertices(C,R)
  
  #Definição do Sigma------------
  if (kernel_type == "gaussian") {
    if (sigma_method == "medist") {
      dist_matrix <- as.matrix(dist(micro))
      sigma <- mean(dist_matrix[lower.tri(dist_matrix)])
    }else if (is.null(sigma)) {
      stop("O parâmetro sigma deve ser fornecido se kernel_type = 'gaussian' e sigma_method = 'user'.")
    }
  }
  
  # Definir e centrar Kernel Matrix------------
  m <- n*(2^p)
  K <- matrix(0, m, m)
  for (i in seq_len(m)) {
    xi <- micro[i, ]
    for (j in seq_len(m)) {
      xj <- micro[j, ]
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
  ones_m <- matrix(1, m, m) / m
  K <- (diag(m) - ones_m) %*% K %*% (diag(m) - ones_m)
  
  
  #Projeção dos dados--------------
  eig <- eigen(K, symmetric = TRUE)
  lambda <- eig$values
  V <- eig$vectors
  projected_data <- K %*% V %*% diag(1 / sqrt(abs(lambda)))
  
  
  
  #Definição dos 80%--------------
  expl <- cumsum(lambda) / sum(lambda)
  pcexp <- which(expl >= 0.8)[1]
  
 
  #Projeção data_test---------------
    #Setup Inicial
    C_test <- as.matrix(data_test@Centers)
    R_test <- as.matrix(data_test@Ranges)
    n_test <- nrow(C_test)
    p_test <- ncol(C_test)
    
    #Gera os microdados---------
    micro_test <- generate_vertices(C_test,R_test)
    
    #Definição K_new
    m_test <- n_test*(2^p)
    projected_data_test <- matrix(0, nrow = m_test, ncol = m)
    for (d in 1:m_test) {
      
      k_new <- numeric(m)
      
      for (i in 1:m) {
        xi <- micro[i,]
        xd <- micro_test[d,]
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
    
    #Construção do modelo SVM vs LDA------------------
    if (separation_method == "SVM") {
      
      svm_model <- svm(x = micro, y = as.factor(label_train), kernel = "linear")
      svm_pred <- predict(svm_model, micro_test)
      
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
      
      lda_model <- lda(x = micro, grouping = as.factor(label_train))
      lda_pred <- predict(lda_model, micro_test)
      
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
    

  #Reconstruir os retângulos  
  mcar <- reconstruct_vertices_v(projected_data_test[,1:pcexp],p_test)
  classif <- class_pos

  block_size <- 2^p
  nblocks <- nrow(mcar) / block_size
  if (nblocks != floor(nblocks)) {
    stop("nrow(mcar) não é múltiplo inteiro de 2^p_test.")
  }
  
  rect_classif <- integer(nblocks)
  for (b in seq_len(nblocks)) {
    idx <- ((b - 1) * block_size + 1):(b * block_size)
    bloc_cls <- classif[idx]
    freqs <- table(bloc_cls)
    rect_classif[b] <- as.integer(names(which.max(freqs)))
  }
  
  # --- Plot
  if (ncol(mcar) >= 2) {
    xlim <- range(mcar[, 1]); ylim <- range(mcar[, 2])
    
    rect_fac <- factor(rect_classif)
    pal <- rainbow(nlevels(rect_fac))
    border_cols <- pal[as.integer(rect_fac)]
    
    plot(
      NA, xlim = xlim, ylim = ylim,
      xlab = "PC1", ylab = "PC2",
      main = paste0("Plot of ", nblocks, " rectangles")
    )
    
    for (b in seq_len(nblocks)) {
      idx <- ((b - 1) * block_size + 1):(b * block_size)
      bloco <- mcar[idx, , drop = FALSE]
      ord <- chull(bloco); ord <- c(ord, ord[1])
      polygon(
        x = bloco[ord, 1],
        y = bloco[ord, 2],
        border = border_cols[b],
        col = NA,
        lwd = 2
      )
    }
  }
  
  #Definição label do macrodado
  label_rect <- label_test[seq(1, length(label_test), by = 2^p)]
  
  all_labels <- sort(unique(c(label_rect, rect_classif)))
  label_rect <- factor(label_rect, levels = all_labels)
  rect_classif <- factor(rect_classif, levels = all_labels)
  
  conf_rect <- confusionMatrix(label_rect, rect_classif)
  metrics_rect <- conf_rect$byClass

  
  return(list(
    sigma=sigma,
    micro = micro,
    micro_test = micro_test,
    mcar = mcar,
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
    metrics_pos = metrics_pos,
    conf_rect=conf_rect,
    metrics_rect=metrics_rect
  ))
}