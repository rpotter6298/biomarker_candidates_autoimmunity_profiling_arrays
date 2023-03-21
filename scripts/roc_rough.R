roc_curve <- function(df, report, P.Val=0.05, stepwise=TRUE){
  lrep_sigs <- report[report$adj.P.Val<P.Val,]
  analytes <- row.names(lrep_sigs)
  df_selected <- df[,c("group", analytes)]
  df_selected$group = as.numeric(df_selected$group)
  logistic_model <- glm(group ~ ., data = df_selected, family = "binomial")
  
  if (stepwise) {
    logistic_model <- step(logistic_model, direction = "backward")
  }
  
  predicted_probabilities <- predict(logistic_model, type = "response")
  roc_obj <- roc(df_selected$group, predicted_probabilities)
  auc <- auc(roc_obj)
  return(plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc, 3), ")")))
}


