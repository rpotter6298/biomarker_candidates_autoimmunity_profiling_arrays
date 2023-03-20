df = P03_Transformed$Set_3

limma_report = reports_P03_Transformed_bcrsn$Set_3_limma
lrep = limma_report[limma_report$adj.P.Val<0.05,]
analytes = row.names(reports_P03_Transformed_bcrsn$Set_3_limma[reports_P03_Transformed_bcrsn$Set_3_limma$adj.P.Val<0.05,])
df_selected <-df[,c("group", analytes)]
df_selected$group = as.numeric(df_selected$group)
logistic_model <- glm(group ~ ., data = df_selected, family = "binomial")
predicted_probabilities <- predict(logistic_model, type = "response")
roc_obj <- roc(df$group, predicted_probabilities)
auc <- auc(roc_obj)
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc, 3), ")"))
