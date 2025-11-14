library(class)

#' @param train_mat   Matriz de treino (n × p)
#' @param train_labels Fator ou vetor de labels de comprimento n
#' @param test_mat    Matriz de teste (m × p)
#' @param k           Número de vizinhos (padrão = 5)
#' @return Fator com as labels previstas para cada linha de test_mat
#' @examples

i_knn_intervalar <- function(data_train, data_test, label_train, label_test, k = 5) {

  
  if (!inherits(data_train, "intData")) {
    stop("O conjunto de dados deve estar na classe 'intData'.")
  }
  
  if (!inherits(data_test, "intData")) {
    stop("O conjunto de dados deve estar na classe 'intData'.")
  }
  
  if (!is.numeric(k) || k < 1 || k != as.integer(k)) {
    stop("k deve ser um inteiro ??? 1.")
  }
  
  C <- as.matrix(data_train@Centers)
  R <- as.matrix(data_train@Ranges)
  n <- nrow(C)
  p <- ncol(C)
  
  C_test <- as.matrix(data_test@Centers)
  R_test <- as.matrix(data_test@Ranges)
  n_test <- nrow(C_test)
  p_test <- ncol(C_test)
  
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
  
  preds <- knn(
    train = transformed_data,
    test  = transformed_data_test,
    cl    = label_train,
    k     = k
  )
  
  conf_mat <- confusionMatrix(as.factor(preds), as.factor(label_test))
  metrics_cm <- conf_mat$byClass
  
  return(list(
    preds = preds,
    cm = conf_mat,
    metrics_cm = metrics_cm
  ))
}
