% Function to predict using AR(1)
% Here we input a regression model, and one newx, then we will do rolling
% prediction for pred_periods using the predicted y
function pred = predAR1(model, newX, pred_period)
    ypred = zeros(pred_period, 1);
    for t = 1:pred_period
        if t == 1 % use the provided new x
            newx = newX;
        else % use the predicted y
            newx = ypred(t-1, 1);
        end
        ypred(t, 1) = predict(model, newx);
    end
    pred = ypred;
end