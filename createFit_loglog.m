function [fitresult, gof] = createFit1(X, Y)
%CREATEFIT1(X,Y)
%  Create a fit.
%
%  Data for 'angular_l' fit:
%      X Input : X
%      Y Output: Y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 07-Nov-2018 20:40:37


%% Fit: 'angular_l'.
[xData, yData] = prepareCurveData( X, Y );

% Set up fittype and options.
ft = fittype( 'a*x + b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.132420604790519 0.765735252974052];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'angular_l' );
h = plot( fitresult, xData, yData );
legend( h, 'Y vs. X', 'angular_l', 'Location', 'NorthEast' );
% Label axes
xlabel X
ylabel Y
grid on


