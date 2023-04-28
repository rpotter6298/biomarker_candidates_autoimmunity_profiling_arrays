fancy_roc <- function(df, report, P.Val = 0.05, validation = "none") {
  library(pROC)
  library(ggplot2)
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
  
  roc_data <- data.frame(
    FPR = roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Model = "Logistic Regression"
  )
  
  if (!is.null(logistic_model_stepwise)) {
    roc_data_stepwise <- data.frame(
      FPR = roc_obj_stepwise$specificities,
      TPR = roc_obj_stepwise$sensitivities,
      Model = "Logistic Regression (Backwards Step)"
    )
    roc_data <- rbind(roc_data, roc_data_stepwise)
  }
  
    plot <- ggplot(data = roc_data, aes(x = FPR, y = TPR, color = Model)) +
      geom_line(size = 1) +
      labs(
        x = "1 - Specificity",
        y = "Sensitivity",
        title = ""
      ) 
    
    +
      theme(
        panel.background = element_rect("white"),
        legend.position = c(0.75,0.15), # Change legend location
        #legend.justification = "bottom", 
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=14), # Change legend text size
        legend.background = element_rect(color="black"),
        axis.text = element_text(size = 14) # Change axis tick text size
      ) +
      coord_cartesian(xlim = c(1, 0), ylim = c(0, 1)) +
      scale_color_manual(values = c("red", "blue")) +
      annotate(
        "text",
        x = 0.4,
        y = 0.3,
        label = paste("AUC: ", paste0(round(auc, 3))),
        color = "red",
        size = 6 #AUC text size
      )
    
    if (!is.null(logistic_model_stepwise)) {
      plot <- plot +
        annotate(
          "text",
          x = 0.4,
          y = 0.25,
          label = paste("AUC: ", paste0(round(auc_stepwise, 3))),
          color = "blue",
          size = 6 #AUC text size
        )
    }
    
    print(plot)
}
