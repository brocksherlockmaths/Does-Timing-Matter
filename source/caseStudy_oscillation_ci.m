%% A case study to show the effects of measurement error on oscillating systems
rng(2)
width = 900; height = 700; left = 50; bottom = 80;

% number of data points
numData = 50;

% Always plot between [-3,3]
extendedX = linspace(-2.5,2.5,1001); 
extendedX_col = extendedX(:);        % column vector for predint

fs = 26;

% CI appearance settings
ciLevel   = 0.95;   % 95% confidence interval for the mean response
fillAlpha = 0.20;   % transparent shading

%% Case 1 - Non-controlled Exp, classical error

% epsilon parameters
muE = 0; sigmaE = sqrt(0.05);
% W parameters
muW = 0; sigmaW = sqrt(4/9);
% U parameters
muU = 0; sigmaU = sqrt(0.05);

% nonlinear model to calculate Y
model = @(a,b,X,epsilon) a.*cos(b.*X) + epsilon;

% Generate true X data
X = unifrnd(-2,2,[numData,1]);
% Generate observed W data
U = normrnd(muU,sigmaU,numData,1);
W = X + U;
% Generate epsilon error
epsilon = normrnd(muE,sigmaE,numData,1);
% Generate Y data
aTrue = 3; bTrue = 4;
Y = model(aTrue,bTrue,X,epsilon);

% Fit a model to true data (fit parameter 'a', keep b fixed)
ft = fittype(@(a,x) model(a,bTrue,x,0), 'independent','x');
fo = fitoptions('Method','NonlinearLeastSquares');
[trueFit, trueGOF] = fit(X,Y,ft,fo); 

% Fit with observed data
[obsFit, obsGOF] = fit(W,Y,ft,fo);   

% Plot the fits
fig = figure(1);
fig.Position = [left, bottom, width, height];
scatter(X,Y,100,'bo','LineWidth',2)
hold on

% True fit (blue)
extendedY = feval(trueFit, extendedX);
p = plot(extendedX,extendedY,'b');
set(p,'LineWidth',3);

% CI band for mean response (seaborn-style)
yci_true = predint(trueFit, extendedX_col, ciLevel, 'functional'); % [low high]
ciPatch(gca, extendedX_col, yci_true(:,1), yci_true(:,2), 'b', fillAlpha);

% Observed fit (red dashed)
scatter(W,Y,100,'r^','LineWidth',2)
extendedY = feval(obsFit, extendedX);
p = plot(extendedX,extendedY,'r--');
set(p,'LineWidth',3);

% CI band (observed fit)
yci_obs = predint(obsFit, extendedX_col, ciLevel, 'functional');
ciPatch(gca, extendedX_col, yci_obs(:,1), yci_obs(:,2), 'r', fillAlpha);

legend off
xlabel('Independent Variable')
ylabel('Dependent Variable')
xlim([-2.5,2.5])
hold off
set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);

exportgraphics(gcf,'Figures/oscillating_ci_a.png')

%% Case 2 - Non-controlled exp, Berkson error

clear X U W Y
% X - true data,
% W - observed data
% Y - response

% Generate true X data
W = unifrnd(-2,2,numData,1);
% Generate observed W data
U = normrnd(muU,sigmaU,numData,1);
X = W + U;
% Generate epsilon error
epsilon = normrnd(muE,sigmaE,numData,1);
% Generate Y data
Y = model(aTrue,bTrue,X,epsilon);

% Fit a model to true data
[trueFit, trueGOF] = fit(X,Y,ft,fo); 

% Fit with observed data
[obsFit, obsGOF] = fit(W,Y,ft,fo);   

% Plot the fits
fig = figure(2);
fig.Position = [left, bottom, width, height];
scatter(X,Y,100,'bo','LineWidth',2)
hold on

% True fit (blue)
extendedY = feval(trueFit, extendedX);
p = plot(extendedX,extendedY,'b');
set(p,'LineWidth',3);
% CI band
yci_true = predint(trueFit, extendedX_col, ciLevel, 'functional');
ciPatch(gca, extendedX_col, yci_true(:,1), yci_true(:,2), 'b', fillAlpha);

% Observed fit (red dashed)
scatter(W,Y,100,'r^','LineWidth',2)
extendedY = feval(obsFit, extendedX);
p = plot(extendedX,extendedY,'r--');
set(p,'LineWidth',3);
% CI band (observed)
yci_obs = predint(obsFit, extendedX_col, ciLevel, 'functional');
ciPatch(gca, extendedX_col, yci_obs(:,1), yci_obs(:,2), 'r', fillAlpha);

legend off
xlabel('Independent Variable')
ylabel('Dependent Variable')
xlim([-2.5,2.5])
hold off
set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);

exportgraphics(gcf,'Figures/oscillating_ci_b.png')

%% Case 3 - controlled exp - Berkson error

clear X U W Y

% Generate observed W data
% Data is taken at 5 equally spaced points with 10 observations at each point
datPerPoint = numData / 5;
U = normrnd(muU,sigmaU,numData,1);
for l = 1:2

    if l == 1 % misaligned sample times
        W = repmat(linspace(-2,2,5)',datPerPoint,1);
    else % align sample times to peaks
        W = repmat((-pi/2:pi/4:pi/2)',datPerPoint,1);
    end

    % Generate true X data
    X = W + U;
    % Generate epsilon error
    epsilon = normrnd(muE,sigmaE,numData,1);
    % Generate Y data
    Y = model(aTrue,bTrue,X,epsilon);

    % Fit a model to true data
    [trueFit, trueGOF] = fit(X,Y,ft,fo); 
    % Fit with observed data
    [obsFit, obsGOF] = fit(W,Y,ft,fo);   

    % Plot the fits
    fig = figure(2+l);
    fig.Position = [left, bottom, width, height];
    nexttile
    scatter(X,Y,100,'bo','LineWidth',2)
    hold on

    % True fit (blue)
    extendedY = feval(trueFit, extendedX);
    p = plot(extendedX,extendedY,'b');
    set(p,'LineWidth',3);
    % CI band
    yci_true = predint(trueFit, extendedX_col, ciLevel, 'functional');
    ciPatch(gca, extendedX_col, yci_true(:,1), yci_true(:,2), 'b', fillAlpha);

    % Observed fit (red dashed)
    scatter(W,Y,100,'r^','LineWidth',2)
    extendedY = feval(obsFit, extendedX);
    p = plot(extendedX,extendedY,'r--');
    set(p,'LineWidth',3);
    % CI band (observed)
    yci_obs = predint(obsFit, extendedX_col, ciLevel, 'functional');
    ciPatch(gca, extendedX_col, yci_obs(:,1), yci_obs(:,2), 'r', fillAlpha);

    legend off
    xlabel('Independent Variable')
    ylabel('Dependent Variable')
    xlim([-2.5,2.5])

    hold off

    set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);
    if l == 1
        exportgraphics(gcf,'Figures/oscillating_ci_c.png')
    else
        lgd = legend('True Data','','Observed Data','');
        lgd.Location = 'southeast';
        lgd.Orientation = 'horizontal';
        lgd.FontSize = 22;
        exportgraphics(gcf,'Figures/oscillating_ci_d.png')
    end
end

lgd = legend('True Data','','Observed Data','');
lgd.Location = 'southeast';
lgd.Orientation = 'horizontal';
lgd.FontSize = 22;

set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);

%% ---------- Local helper ----------
function ciPatch(ax, x, ylo, yhi, colorCharOrRGB, alphaVal)
% ciPatch draws a shaded confidence interval polygon behind other graphics.
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
        % ignore if uistack unavailable
    end
end
