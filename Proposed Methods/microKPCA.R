library(triangle)
#' Gera microdados para cada observação como matrizes simples sem nomes de linhas/colunas
#'
#' @param C          Matriz ou data.frame de centros, dimensões n × p
#' @param R          Matriz ou data.frame de ranges, dimensões n × p
#' @param latentDist Caractere ou vetor (comprimento p) com distribuição em [-1,1]:
#'                   "Unif", "Triang", "TNorm", "InvTri", "Beta", "KDE", "Degenerated"
#' @param nMicro     Número inteiro de micro-pontos por observação
#' @param TriangParam Vetor (ou único) com modo da triangular (em [-1,1])
#' @param BetaParam.a Vetor (ou único) ?? da Beta (em [0,1] antes de mapear para [-1,1])
#' @param BetaParam.b Vetor (ou único) ?? da Beta
#' @return Lista de comprimento n (observações). Cada elemento é uma matriz nMicro × p, sem nomes.
generateMicroData <- function(
    C, R,
    latentDist,
    nMicro,
    TriangParam = 0,
    BetaParam.a  = 1,
    BetaParam.b  = 1
) {
  C <- as.matrix(C)
  R <- as.matrix(R)
  if (!all(dim(C) == dim(R))) {
    stop("As matrizes C e R devem ter as mesmas dimensões.")
  }
  n <- nrow(C)
  p <- ncol(C)
  if (!is.numeric(nMicro) || length(nMicro) != 1 || nMicro < 1) {
    stop("nMicro deve ser um inteiro ??? 1.")
  }
  
  if (length(latentDist) == 1) {
    latentDist <- rep(latentDist, p)
  }
  if (length(latentDist) != p) {
    stop("latentDist deve ser caractere único ou vetor de comprimento p.")
  }
  if (length(TriangParam) == 1) TriangParam <- rep(TriangParam, p)
  if (length(TriangParam) != p) {
    stop("TriangParam deve ser escalar ou vetor de comprimento p.")
  }
  if (length(BetaParam.a) == 1) BetaParam.a <- rep(BetaParam.a, p)
  if (length(BetaParam.b) == 1) BetaParam.b <- rep(BetaParam.b, p)
  if (!all(length(BetaParam.a) == p, length(BetaParam.b) == p)) {
    stop("BetaParam.a e BetaParam.b devem ser escalares ou vetores de comprimento p.")
  }
  
  sampleU <- function(dist, n, mode_tri, a_beta, b_beta) {
    dist <- match.arg(dist,
                      c("Unif","Triang","TNorm","InvTri","Beta","KDE","Degenerated")
    )
    if (dist == "Unif") {
      return(runif(n, -1, 1))
    }
    if (dist == "Degenerated") {
      return(rep(0, n))
    }
    if (dist == "Triang") {
      m <- mode_tri  
      x <- rtriangle(n, -1, 1, m)
      return(x)
    }
    if (dist == "InvTri") {
      u0 <- runif(n)
      x <- numeric(n)
      idx <- which(u0 < 0.5)
      if (length(idx) > 0) {
        x[idx] <- -sqrt(1 - 2 * u0[idx])
      }
      idx2 <- which(u0 >= 0.5)
      if (length(idx2) > 0) {
        x[idx2] <-  sqrt(2 * u0[idx2] - 1)
      }
      return(x)
    }
    if (dist == "TNorm") {
      out <- numeric(0)
      while (length(out) < n) {
        candidate <- rnorm(n * 2)
        keep <- candidate[candidate >= -1 & candidate <= 1]
        if (length(keep) > 0) {
          needed <- min(length(keep), n - length(out))
          out <- c(out, keep[1:needed])
        }
      }
      return(out[1:n])
    }
    if (dist == "Beta") {
      X <- rbeta(n, shape1 = a_beta, shape2 = b_beta)
      U <- 2 * X - 1
      return(U)
    }
    if (dist == "KDE") {
      stop("Distribuição 'KDE' requer microdados para estimar densidade. Não implementada.")
    }
    stop("Distribuição desconhecida em sampleU().")
  }
  
  microdata_list <-  matrix(numeric(0), nrow = 0, ncol = p)
  for (i in seq_len(n)) {
    mat_i <- matrix(NA_real_, nrow = nMicro, ncol = p)
    for (j in seq_len(p)) {
      U_ij <- sampleU(
        dist      = latentDist[j],
        n         = nMicro,
        mode_tri  = TriangParam[j],
        a_beta    = BetaParam.a[j],
        b_beta    = BetaParam.b[j]
      )
      centro_ij <- C[i, j]
      meio_raio <- R[i, j] / 2
      mat_i[, j] <- centro_ij + U_ij * meio_raio
    }

    dimnames(mat_i) <- NULL
    microdata_list <- rbind(microdata_list,mat_i)
  }
  
  return(microdata_list)
}


reconstruct_vertices <- function(M, p, nmic) {
  M <- as.matrix(M)
  block_size <- nmic
  total_linhas <- nrow(M)
  
  if (total_linhas %% block_size != 0) {
    stop("Número de linhas em M não é múltiplo exato de nmic.")
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


micro_kpca <- function(      microdata=NULL,
                             microdata_test=NULL,
                             nmic=10,
                             separation_method = c("SVM", "LDA"),
                             data_train,
                             data_test,
                             label_train, #label_train e label_test devem ter a dimensão de micro e micro_teste
                             label_test,
                             kernel_type = c("gaussian", "linear", "polynomial"),
                             sigma_method = c("user", "medist"),
                             sigma = NULL,
                             poly_degree = 3,
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
  latentdist <- data_train@LatentDist
  TriangParam <- 3*as.numeric(data_train@LatentParam[[2]])
  
  #Definição do microdata
  if (is.null(microdata)) {
    micro <- generateMicroData(
      C = C,
      R = R,
      latentDist = latentdist,
      TriangParam = TriangParam,
      nMicro = nmic
    )
  } else {
    micro <- as.matrix(microdata)
    dimnames(micro) <- NULL
  }
  
  #Definição do Sigma
  if (kernel_type == "gaussian") {
    if (sigma_method == "medist") {
      dist_matrix <- as.matrix(dist(micro))
      sigma <- mean(dist_matrix[lower.tri(dist_matrix)])
    } else if (is.null(sigma)) {
      stop("O parâmetro sigma deve ser fornecido se kernel_type = 'gaussian' e sigma_method = 'user'.")
    }
  }
  
  # Definir e centrar Kernel Matrix
  m <- nrow(micro)
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
    latentdist_test <- data_test@LatentDist
    
    #Definição do microdata_test
    if (is.null(microdata_test)) {
      micro_test <- generateMicroData(
        C = C_test,
        R = R_test,
        latentDist = latentdist_test,
        nMicro = nmic
      )
    } else {
      micro_test <- as.matrix(microdata_test)
      dimnames(micro_test) <- NULL
    }
    
    #Definição K_new
    projected_data_test <- matrix(0, nrow = nrow(micro_test), ncol = nrow(micro))
    for (d in 1:nrow(micro_test)) {
      
      k_new <- numeric(nrow(micro))
      
      for (i in 1:nrow(micro)) {
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
  mcar <- reconstruct_vertices(projected_data_test[,1:pcexp],p,nmic)
  classif <- class_pos
  
  block_size_micro <- nmic
  nrect <- length(classif) / block_size_micro
  if (nrect != floor(nrect)) {
    stop("length(classif) não é múltiplo de nmic.")
  }
  nrect <- as.integer(nrect)
  
  rect_classif <- integer(nrect)
  for (i in seq_len(nrect)) {
    idx_micro <- ((i - 1) * block_size_micro + 1):(i * block_size_micro)
    freqs <- table(classif[idx_micro])
    rect_classif[i] <- as.integer(names(which.max(freqs)))
  }
  
  block_size_vertices <- 2^p
  nblocks <- nrow(mcar) / block_size_vertices
  if (nblocks != floor(nblocks)) {
    stop("nrow(mcar) não é múltiplo de 2^p (número de vértices por retângulo).")
  }
  nblocks <- as.integer(nblocks)
  
  if (nblocks != nrect) {
    stop("Número de retângulos inconsistente entre micro (nrect) e vértices (nblocks).")
  }
  
  if (ncol(mcar) >= 2) {
    xlim <- range(mcar[, 1]); ylim <- range(mcar[, 2])
    rect_fac <- factor(rect_classif)
    pal <- rainbow(nlevels(rect_fac))
    border_cols <- pal[as.integer(rect_fac)]
    
    plot(NA, xlim = xlim, ylim = ylim, xlab = "PC1", ylab = "PC2",
         main = paste0("Plot de ", nblocks, " retângulos"))
    
    for (i in seq_len(nblocks)) {
      idx_vert <- ((i - 1) * block_size_vertices + 1):(i * block_size_vertices)
      bloco <- mcar[idx_vert, , drop = FALSE]
      ord <- chull(bloco); ord <- c(ord, ord[1])
      polygon(x = bloco[ord, 1], y = bloco[ord, 2],
              border = border_cols[i], col = NA, lwd = 2)
    }
  }
  
  
  #Definição label do macrodado
  label_rect <- label_test[seq(1, length(label_test), by = nmic)]
  
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
