%% A case study to show the effects of measurement error on oscilatting systems

width = 1200;
height = 1000;
left = 50;
bottom = 80;

% number of data points
numData = 50;

% Always plot between [-3,3]
extendedX = linspace(-2.5,2.5,1001);

numRepeats = 10000;
trueFits = nan(numRepeats,4);
observedFits = trueFits;

%% Case 1 - Non-controlled Exp, classical error

% epsilon parameters
muE = 0; sigmaE = sqrt(0.05);
% W parameters
muW = 0; sigmaW = sqrt(4/9);
% U parameters
muU = 0; sigmaU = sqrt(0.05);

% linear model to calculate Y
model = @(a,b,X,epsilon) a.*cos(b.*X) + epsilon;

for rep=1:numRepeats
    % Generate true X data
    X = unifrnd(-2,2,[numData,1]);% Generate observed W data
    U = normrnd(muU,sigmaU,numData,1);
    W = X + U;
    % Generate epsilon error
    epsilon = normrnd(muE,sigmaE,numData,1);
    % Generate Y data
    aTrue = 3; bTrue = 4;
    Y = model(aTrue,bTrue,X,epsilon);
    
    % Fit a model to true data
    ft = fittype(@(a,x) model(a,bTrue,x,0), independent='x');
    fo = fitoptions('Method','NonlinearLeastSquares');
    trueFit = fit(X,Y,ft,fo);
    
    % Fit with observed data
    obsFit = fit(W,Y,ft,fo);
    
    trueFits(rep,1) = coeffvalues(trueFit);
    observedFits(rep,1) = coeffvalues(obsFit);
    
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
    % % X parameters 
    % muW = 0; sigmaX = 1;
    % % epsilon parameters
    % muE = 0; sigmaE = sqrt(0.25);
    % % W parameters
    % muU = 0; sigmaU = sqrt(0.25);
    
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
    trueFit = fit(X,Y,ft,fo);
    
    % Fit with observed data
    obsFit = fit(W,Y,ft,fo);
    
    trueFits(rep,2) = coeffvalues(trueFit);
    observedFits(rep,2) = coeffvalues(obsFit);
        
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
    
    % numData = 50; % total number of data points
    % % X parameters 
    % muW = 0; sigmaX = 1;
    % % epsilon parameters
    % muE = 0; sigmaE = sqrt(0.25);
    % % W parameters
    % muU = 0; sigmaU = sqrt(0.25);
    
    
    % Generate observed W data
    % Data is taken at 5 equally spaced points with 10 observations at each
    % point
    datPerPoint = numData / 5;
    U = normrnd(muU,sigmaU,numData,1);
    for l = 1:2
        
        if l == 1 % misaligned sample times
            W = repmat(linspace(-2,2,5)',datPerPoint,1);
        else % align sample times to peaks
            W = repmat([-pi/2:pi/4:pi/2]',datPerPoint,1);
        end
    
        % Generate true X data
        X = W + U;
        % Generate epsilon error
        epsilon = normrnd(muE,sigmaE,numData,1);
        % Generate Y data
        Y = model(aTrue,bTrue,X,epsilon);
        
        % Fit a model to true data
        trueFit = fit(X,Y,ft,fo);
        
        % Fit with observed data
        obsFit = fit(W,Y,ft,fo);
        
        trueFits(rep,2+l) = coeffvalues(trueFit);
        observedFits(rep,2+l) = coeffvalues(obsFit);
    end

end

%% Make plots
fig = figure(1);
fig.Position = [left, bottom, width, height];
t = tiledlayout(2,2);
capLabels = {'(a)','(b)','(c)','(d)'};
for i = 1:4
    nexttile
    hold on
    
    [f,xi] = ksdensity(trueFits(:,i));
    plot(xi,f,'b')
    [f,xi] = ksdensity(observedFits(:,i));
    plot(xi,f,'r--')
    xline(aTrue,'LineWidth',3)
    hold off
    xlim([0.5,3.5])
    xlabel(capLabels{i});
end



t.XLabel.String = 'Parameter Estimate';
t.YLabel.String = 'Density';
t.XLabel.FontSize = 32;
t.YLabel.FontSize = 32;

lgd = legend('True Data','Observed Data');

lgd.Layout.Tile = 'north';
lgd.Orientation = 'horizontal';
lgd.FontSize = 22;
fs = 26;
set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);
set(findobj(gcf,'type','line'),'LineWidth',3);

% exportgraphics(gcf,'Figures/oscillating_density.png')

%% Individual plots
capLabels = {'a','b','c','d'};
fs = 26;  % tick label font size
width = 900; height = 700; left = 50; bottom = 80;

for i = 1:4
    fig = figure(i);
    fig.Position = [left, bottom, width, height];
    hold on
    
    % KDE for true fits
    [f,xi] = ksdensity(trueFits(:,i));
    plot(xi, f, 'b', 'LineWidth', 3)
    
    % KDE for observed fits
    [f,xi] = ksdensity(observedFits(:,i));
    plot(xi, f, 'r--', 'LineWidth', 3)
    
    % Vertical line at aTrue
    xline(aTrue, 'LineWidth', 3)
    
    hold off
    xlim([0.5, 3.5])
    xlabel('Parameter Estimate')
    ylabel('Density')
    % title(capLabels{i}, 'FontWeight', 'normal')  % panel label
    
    % Legend
    if i == 1
        lgd = legend('True Data', 'Observed Data');
        lgd.Location = 'northwest';
        lgd.Orientation = 'horizontal';
        lgd.FontSize = 22;
    end
    
    % Axes styling
    set(findobj(gcf, 'type', 'axes'), ...
        'FontSize', fs, ...
        'LineWidth', 2, ...
        'LabelFontSizeMultiplier', 1.4);

    % Export each figure
    exportgraphics(gcf, fullfile('Figures', sprintf('oscillating_density_%s.png', capLabels{i})))
end
