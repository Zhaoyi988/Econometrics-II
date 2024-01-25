% Function to calculate MAE and RMSE
% Input: predictions and true values
function [mae, rmse] = calcAccuracy(Ypred, Ytrue)
    Yerror = abs(Ypred - Ytrue);
    mae = nanmean(Yerror, 2);
    rmse = sqrt(nanmean(Yerror.^2, 2));
end