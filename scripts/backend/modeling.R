loocv_validation <- function(df, logistic_model) {
  n_samples <- nrow(df)
  predicted_probabilities <- numeric(n_samples)
  
  for (i in 1:n_samples) {
    test_set <- df[i,]
    test_prob <- predict(logistic_model, newdata = test_set[-1], type = "response")
    predicted_probabilities[i] <- test_prob
  }
  
  return(predicted_probabilities)
  
}

roc_curve <- function(df, report, P.Val = 0.05, stepwise = FALSE, validation = "none") {
  lrep_sigs <- report[report$adj.P.Val < P.Val,]
  analytes <- row.names(lrep_sigs)
  df_selected <- df[, c("group", analytes)]
  df_selected$group = as.numeric(df_selected$group)
  logistic_model <- glm(group ~ ., data = df_selected, family = "binomial")
  
  if (stepwise) {
    logistic_model <- step(logistic_model, direction = "backward")
  }
  
  model_summary <- summary(logistic_model)
  accuracy <- NULL
  if (validation == "LOOCV") {
    predicted_probabilities <- loocv_validation(df, logistic_model)
    true_labels <- as.numeric(df$group)
    threshold <- 0.5
    predicted_labels <- ifelse(predicted_probabilities >= threshold, 1, 0)
    correct_predictions <- predicted_labels == true_labels
    accuracy <- mean(correct_predictions)
    roc_obj <- roc(true_labels, predicted_probabilities)
  } else {
    predicted_probabilities <- predict(logistic_model, type = "response")
    roc_obj <- roc(df_selected$group, predicted_probabilities)
  }
  
  auc <- auc(roc_obj)
  plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc, 3), ")"),asp=1)
  roc_plot <- recordPlot()
  
  logistic_regression <- list(model_summary = model_summary, logistic_model = logistic_model)
  validation_results <- list(validation_method = validation, predicted_probabilities = predicted_probabilities, accuracy = accuracy)
  roc_results <- list(roc_obj = roc_obj, auc = auc, plot = roc_plot)
  
  if (validation == "none") {
    results <- list(logistic_regression = logistic_regression, roc = roc_results)
  } else {
    results <- list(logistic_regression = logistic_regression, validation = validation_results, roc = roc_results)
  }
  
  return(results)
}

fancy_roc <- function(df, report, P.Val = 0.05, validation = "none") {
  library(pROC)
  lrep_sigs <- report[report$adj.P.Val < P.Val,]
  analytes <- row.names(lrep_sigs)
  df_selected <- df[, c("group", analytes)]
  df_selected$group = as.numeric(df_selected$group)
  logistic_model <- glm(group ~ ., data = df_selected, family = "binomial")
  
  logistic_model_stepwise <- NULL
  auc_stepwise <- NULL
  predicted_probabilities_stepwise <- NULL
  if (length(analytes) > 1) {
    logistic_model_stepwise <- step(logistic_model, direction = "backward")
    predicted_probabilities_stepwise <- predict(logistic_model_stepwise, type = "response")
    roc_obj_stepwise <- roc(df_selected$group, predicted_probabilities_stepwise)
    auc_stepwise <- auc(roc_obj_stepwise)
  }
  
  predicted_probabilities <- predict(logistic_model, type = "response")
  roc_obj <- roc(df_selected$group, predicted_probabilities)
  auc <- auc(roc_obj)
  
  par(cex.axis = 1.5)
  plot(roc_obj, col="red", main = "", xlim=c(1,0), ylim=c(0,1))

  if (!is.null(logistic_model_stepwise)) {
    predicted_probabilities_stepwise <- predict(logistic_model_stepwise, type = "response")
    roc_obj_stepwise <- roc(df_selected$group, predicted_probabilities_stepwise)
    auc_stepwise <- auc(roc_obj_stepwise)
    lines(roc_obj_stepwise, col="blue")
    legend("bottomright", legend=c("Logistic Regression", "Logistic Regression (Backwards Step)"), col=c("red", "blue"), lty=1, cex=0.8)
  }
  else {
    legend("bottomright", legend="Logistic Regression", col="red", lty=1, cex=0.8)
  }
  text(x=0.4, y=0.4, labels=paste("AUC: ", paste0(round(auc, 3))), pos=1, cex=1, col="red")
  text(x=0.4, y=0.35, labels=paste("AUC: ", paste0(round(auc_stepwise, 3))), pos=1, cex=1, col="blue")}

