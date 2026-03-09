%% Demonstrate the effects of measurement error (both classical and Berkson, for controlled and uncontrolled exp)
rng(0);
width = 900; height = 700; left = 50; bottom = 80;
fs = 34;

% CI appearance settings 
ciLevel   = 0.95;   % 95% confidence interval
fillAlpha = 0.20;  

%% Case 1 - Non-controlled Exp, classical error

% X - true data,
% W - observed data
% Y - response
%
% Error model is W = X + U where U is independent of X with mean zero and
% variance 0.25
% X has mean zero and variance 1
% Y = a + bX + eps, with a = 0, b = 1 and epsilon has mean zero and
% variance 0.25

numData = 100; % total number of data points
% X parameters 
muX = 0; sigmaX = 1;
% epsilon parameters
muE = 0; sigmaE = sqrt(0.25);
% W (error U) parameters
muW = 0; sigmaW = sqrt(0.25);

% linear model to calculate Y
model = @(a,b,X,epsilon) a + b*X + epsilon;

% Generate true X data
X = normrnd(muX,sigmaX,numData,1);
% Generate observed W data
U = normrnd(muW,sigmaW,numData,1);
W = X + U;
% Generate epsilon error
epsilon = normrnd(muE,sigmaE,numData,1);
% Generate Y data
aTrue = 0; bTrue = 1;
Y = model(aTrue,bTrue,X,epsilon);

% Fit a linear model to true data
fo = fitoptions('Method','LinearLeastSquares');
[trueFit, trueGOF] = fit(X,Y,'poly1',fo); 

% Fit with observed data
[obsFit, obsGOF]   = fit(W,Y,'poly1',fo); 

% Plot the fits
fig = figure(1);
fig.Position = [left, bottom, width, height];

scatter(X,Y,100,'bo','LineWidth',2)
hold on

% True fit (blue solid)
pTrue = plot(trueFit,'b'); % preserves your style
set(pTrue,'LineWidth',3);
xg_true  = get(pTrue,'XData'); xg_true = xg_true(:); % ensure column
yci_true = predint(trueFit, xg_true, ciLevel, 'functional'); % [low high]
ciPatch(gca, xg_true, yci_true(:,1), yci_true(:,2), 'b', fillAlpha);

% Observed fit (red dashed)
scatter(W,Y,100,'r^','LineWidth',2)
pObs = plot(obsFit,'r--');
set(pObs,'LineWidth',3);
xg_obs  = get(pObs,'XData'); xg_obs = xg_obs(:);
yci_obs = predint(obsFit, xg_obs, ciLevel, 'functional');
ciPatch(gca, xg_obs, yci_obs(:,1), yci_obs(:,2), 'r', fillAlpha);

legend off
xlabel('Independent Variable','FontSize',fs)
ylabel('Dependent Variable','FontSize',fs)
hold off
set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);
xlim([-3,5])
exportgraphics(gcf,'Figures/linRegression_addError_ci_a.png')

%% Case 2 - Non-controlled exp, Berkson error

clear X U W Y
% X - true data,
% W - observed data
% Y - response
%
% Error model is X = W + U 
% X has mean zero and variance 1
% Y = a + bX + eps, with a = 0, b = 1 and epsilon has mean zero and
% variance 0.25

numData = 50; % total number of data points
% X parameters 
muW = 0; sigmaX = 1; 
% epsilon parameters
muE = 0; sigmaE = sqrt(0.25);
% U (Berkson) parameters
muU = 0; sigmaU = sqrt(0.25);
% ensure sigmaW exists (self-contained)
sigmaW = sqrt(0.25);

% linear model to calculate Y
model = @(a,b,X,epsilon) a + b*X + epsilon;

% Generate W data
W = normrnd(muW,sigmaW,numData,1);
% Generate true X data
U = normrnd(muU,sigmaU,numData,1);
X = W + U;
% Generate epsilon error
epsilon = normrnd(muE,sigmaE,numData,1);
% Generate Y data
aTrue = 0; bTrue = 1;
Y = model(aTrue,bTrue,X,epsilon);

% Fit a linear model to true data
fo = fitoptions('Method','LinearLeastSquares');
[trueFit, trueGOF] = fit(X,Y,'poly1',fo); 

% Fit with observed data
[obsFit, obsGOF]   = fit(W,Y,'poly1',fo); 

% Plot the fits
fig = figure(2);
fig.Position = [left, bottom, width, height];

scatter(X,Y,100,'bo','LineWidth',2)
hold on

% True fit (blue)
pTrue = plot(trueFit,'b');
set(pTrue,'lineWidth',3);
xg_true  = get(pTrue,'XData'); xg_true = xg_true(:);
yci_true = predint(trueFit, xg_true, ciLevel, 'functional');
ciPatch(gca, xg_true, yci_true(:,1), yci_true(:,2), 'b', fillAlpha);

% Observed fit (red dashed)
scatter(W,Y,100,'r^','LineWidth',2)
pObs = plot(obsFit,'r--');
set(pObs,'lineWidth',3);
xg_obs  = get(pObs,'XData'); xg_obs = xg_obs(:);
yci_obs = predint(obsFit, xg_obs, ciLevel, 'functional');
ciPatch(gca, xg_obs, yci_obs(:,1), yci_obs(:,2), 'r', fillAlpha);

legend off
xlabel('Independent Variable','FontSize',fs)
ylabel('Dependent Variable','FontSize',fs)
hold off
set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);

exportgraphics(gcf,'Figures/linRegression_addError_ci_b.png')

%% Case 3 - controlled exp - Berkson error

clear X U W Y
% X - true data,
% W - observed data
% Y - response
%
% Error model is X = W + U 
% X has mean zero and variance 1
% Y = a + bX + eps, with a = 0, b = 1 and epsilon has mean zero and
% variance 0.25

numData = 50; % total number of data points
% X parameters 
muW = 0; sigmaX = 1; 
% epsilon parameters
muE = 0; sigmaE = sqrt(0.25);
% U (Berkson) parameters
muU = 0; sigmaU = sqrt(0.25);

% linear model to calculate Y
model = @(a,b,X,epsilon) a + b*X + epsilon;

% Generate observed W data
% Data is taken at 5 equally spaced points with 10 observations at each point
datPerPoint = numData / 5;
W = repmat(linspace(-2,2,5)',datPerPoint,1);
% Generate true X data
U = normrnd(muU,sigmaU,numData,1);
X = W + U;
% Generate epsilon error
epsilon = normrnd(muE,sigmaE,numData,1);
% Generate Y data
aTrue = 0; bTrue = 1;
Y = model(aTrue,bTrue,X,epsilon);

% Fit a linear model to true data
fo = fitoptions('Method','LinearLeastSquares');
[trueFit, trueGOF] = fit(X,Y,'poly1',fo); 

% Fit with observed data
[obsFit, obsGOF]   = fit(W,Y,'poly1',fo); 

% Plot the fits
fig = figure(3);
fig.Position = [left, bottom, width, height];

scatter(X,Y,100,'bo','LineWidth',2)
hold on

% True fit (blue)
pTrue = plot(trueFit,'b');
set(pTrue,'LineWidth',3);
xg_true  = get(pTrue,'XData'); xg_true = xg_true(:);
yci_true = predint(trueFit, xg_true, ciLevel, 'functional');
ciPatch(gca, xg_true, yci_true(:,1), yci_true(:,2), 'b', fillAlpha);

% Observed fit (red dashed)
scatter(W,Y,100,'r^','LineWidth',2)
pObs = plot(obsFit,'r--');
set(pObs,'LineWidth',3);
xg_obs  = get(pObs,'XData'); xg_obs = xg_obs(:);
yci_obs = predint(obsFit, xg_obs, ciLevel, 'functional');
ciPatch(gca, xg_obs, yci_obs(:,1), yci_obs(:,2), 'r', fillAlpha);

legend off
xlabel('Independent Variable','FontSize',fs)
ylabel('Dependent Variable','FontSize',fs)
hold off

% Preserve your tiledlayout label code safely if t exists
if exist('t','var') && isa(t,'matlab.graphics.layout.TiledChartLayout')
    t.XLabel.String = 'Predictor';
    t.YLabel.String = 'Response';
    t.XLabel.FontSize = 32;
    t.YLabel.FontSize = 32;
end

lgd = legend('True Data','','Observed Data','','Location','southeast');
lgd.FontSize = 32;

set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);

exportgraphics(gcf,'Figures/linRegression_addError_ci_c.png')
end

%% --------- Local helper (safe CI patch) ----------
function ciPatch(ax, x, ylo, yhi, colorCharOrRGB, alphaVal)
% Ensures shapes are column vectors and makes a shaded CI polygon behind plots.
    x   = x(:);
    ylo = ylo(:);
    yhi = yhi(:);
    if ~(numel(x)==numel(ylo) && numel(x)==numel(yhi))
        error('ciPatch:SizeMismatch', 'x, ylo, yhi must have the same number of elements.');
    end
    xpoly = [x; flipud(x)];
    ypoly = [ylo; flipud(yhi)];
    h = patch('Parent', ax, ...
              'XData', xpoly, 'YData', ypoly, ...
              'FaceColor', colorCharOrRGB, ...
              'FaceAlpha', alphaVal, ...
              'EdgeColor', 'none', ...
              'HandleVisibility','off');
    % Keep the shading behind markers/lines
    try
        uistack(h,'bottom');
    catch
        % uistack might not be available in very old releases; ignore if so.
    end
end
