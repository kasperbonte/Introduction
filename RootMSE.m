function RMSE = RootMSE(H_true,H_model)
RMSE = sqrt(mean((H_true-H_model).^2))
end
