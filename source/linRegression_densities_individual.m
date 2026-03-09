%% Reproducibility (optional)
rng(1);

%% Simulation parameters
numSims = 1000;

% Storage table template
makeResultTable = @() table( ...
    'Size',[numSims 4], ...
    'VariableTypes',{'double','double','double','double'}, ...
    'VariableNames',{'a_hat','b_hat','CIwidth_a','CIwidth_b'});

% Results for each case and predictor type (True X vs Observed W)
results_case1_true = makeResultTable();
results_case1_obs  = makeResultTable();
results_case2_true = makeResultTable();
results_case2_obs  = makeResultTable();
results_case3_true = makeResultTable();
results_case3_obs  = makeResultTable();

% Helper to extract from a single fit object
extractFromFit = @(fobj) deal(coeffvalues(fobj), diff(confint(fobj)));

%% Case 1 - Non-controlled experiment, Classical Error
for s = 1:numSims
    % Parameters
    numData = 100;
    muX = 0; sigmaX = 1;
    muE = 0; sigmaE = sqrt(0.25);
    muU = 0; sigmaU = sqrt(0.25);

    model = @(a,b,X,e) a + b*X + e;

    % Generate data
    X = normrnd(muX, sigmaX, numData, 1);
    U = normrnd(muU, sigmaU, numData, 1);
    W = X + U;
    epsilon = normrnd(muE, sigmaE, numData, 1);
    Y = model(0, 1, X, epsilon);

    % Fits (same iteration; separate single objects for X and W)
    f_true = fit(X, Y, 'poly1');
    f_obs  = fit(W, Y, 'poly1');

    [betaT, ciwT] = extractFromFit(f_true);
    [betaO, ciwO] = extractFromFit(f_obs);

    % Store (poly1 returns [slope, intercept])
    results_case1_true.a_hat(s)     = betaT(2);
    results_case1_true.b_hat(s)     = betaT(1);
    results_case1_true.CIwidth_a(s) = ciwT(2);
    results_case1_true.CIwidth_b(s) = ciwT(1);

    results_case1_obs.a_hat(s)      = betaO(2);
    results_case1_obs.b_hat(s)      = betaO(1);
    results_case1_obs.CIwidth_a(s)  = ciwO(2);
    results_case1_obs.CIwidth_b(s)  = ciwO(1);
end

%% Case 2 - Non-controlled experiment, Berkson Error
for s = 1:numSims
    numData = 50;
    muW = 0; sigmaW = 1;
    muE = 0; sigmaE = sqrt(0.25);
    muU = 0; sigmaU = sqrt(0.25);

    model = @(a,b,X,e) a + b*X + e;

    W = normrnd(muW, sigmaW, numData, 1);
    U = normrnd(muU, sigmaU, numData, 1);
    X = W + U;
    epsilon = normrnd(muE, sigmaE, numData, 1);
    Y = model(0, 1, X, epsilon);

    f_true = fit(X, Y, 'poly1');
    f_obs  = fit(W, Y, 'poly1');

    [betaT, ciwT] = extractFromFit(f_true);
    [betaO, ciwO] = extractFromFit(f_obs);

    results_case2_true.a_hat(s)     = betaT(2);
    results_case2_true.b_hat(s)     = betaT(1);
    results_case2_true.CIwidth_a(s) = ciwT(2);
    results_case2_true.CIwidth_b(s) = ciwT(1);

    results_case2_obs.a_hat(s)      = betaO(2);
    results_case2_obs.b_hat(s)      = betaO(1);
    results_case2_obs.CIwidth_a(s)  = ciwO(2);
    results_case2_obs.CIwidth_b(s)  = ciwO(1);
end

%% Case 3 - Controlled experiment, Berkson Error
for s = 1:numSims
    numData = 50;
    muU = 0; sigmaU = sqrt(0.25);
    muE = 0; sigmaE = sqrt(0.25);

    model = @(a,b,X,e) a + b*X + e;

    datPerPoint = numData / 5;
    W = repmat(linspace(-2, 2, 5)', datPerPoint, 1);
    U = normrnd(muU, sigmaU, numData, 1);
    X = W + U;
    epsilon = normrnd(muE, sigmaE, numData, 1);
    Y = model(0, 1, X, epsilon);

    f_true = fit(X, Y, 'poly1');
    f_obs  = fit(W, Y, 'poly1');

    [betaT, ciwT] = extractFromFit(f_true);
    [betaO, ciwO] = extractFromFit(f_obs);

    results_case3_true.a_hat(s)     = betaT(2);
    results_case3_true.b_hat(s)     = betaT(1);
    results_case3_true.CIwidth_a(s) = ciwT(2);
    results_case3_true.CIwidth_b(s) = ciwT(1);

    results_case3_obs.a_hat(s)      = betaO(2);
    results_case3_obs.b_hat(s)      = betaO(1);
    results_case3_obs.CIwidth_a(s)  = ciwO(2);
    results_case3_obs.CIwidth_b(s)  = ciwO(1);
end

aTrue = 1;
width = 900; height = 700; left = 50; bottom = 80;

fig = figure(1);
fig.Position = [left, bottom, width, height];
t = tiledlayout(1,1);
%capLabels = {'(a)','(b)','(c)','(d)'};
    nexttile
    hold on
    
    [f,xi] = ksdensity(results_case3_true.b_hat);
    plot(xi,f,'b')
    [f,xi] = ksdensity(results_case3_obs.b_hat);
    plot(xi,f,'r--')
    xline(aTrue,'LineWidth',3)
    hold off
    %xlim([0.5,3.5])
    xlabel(capLabels{i});



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

exportgraphics(gcf,'Figures/linRegression_addError_a_density_individual.png')

%% Save for plotting
save('measurement_error_sim_results.mat', ...
    'results_case1_true','results_case1_obs', ...
    'results_case2_true','results_case2_obs', ...
    'results_case3_true','results_case3_obs');

%% Load simulation results
load('measurement_error_sim_results.mat', ...
    'results_case1_true','results_case1_obs', ...
    'results_case2_true','results_case2_obs', ...
    'results_case3_true','results_case3_obs');

cases_true = {results_case1_true, results_case2_true, results_case3_true};
cases_obs  = {results_case1_obs,  results_case2_obs,  results_case3_obs};
caseTitles = { ...
    'Case 1: Non-controlled, Classical', ...
    'Case 2: Non-controlled, Berkson', ...
    'Case 3: Controlled, Berkson'};

colors = struct('X',[0 0.4470 0.7410], 'W',[0.8500 0.3250 0.0980]); % MATLAB blue/red

%% Helper: plot parameter panels (2 rows x 3 cols)
function plot_param_panels(param_name, ci_varname, true_value, ...
    cases_true, cases_obs, caseTitles, colors, figTitle)

width = 900; height = 700; left = 50; bottom = 80;
fig = figure('Position',[0.05 0.08 0.5 0.9]);
%t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
%t.TileIndexing = "columnmajor";

% Precompute x-limits for row 1 and row 2 so they're consistent across columns
allVals_top = []; allVals_bottom = [];
for c = 1:3
    allVals_top    = [allVals_top;    cases_true{c}.(param_name); cases_obs{c}.(param_name)];
    allVals_bottom = [allVals_bottom; cases_true{c}.(ci_varname); cases_obs{c}.(ci_varname)];
end
xlim_top    = [min(allVals_top)    max(allVals_top)];
xlim_bottom = [min(allVals_bottom) max(allVals_bottom)];

% Keep handles to top-row axes for later alignment
ax_top = gobjects(1,3);

for c = 1:3
    T = cases_true{c};
    O = cases_obs{c};

    % --- Top row: parameter estimates ---
    ax_top(c) = nexttile(c);
    [fT,xiT] = ksdensity(T.(param_name));
    [fO,xiO] = ksdensity(O.(param_name));
    plot(xiT, fT, 'Color', colors.X, 'LineWidth', 3); hold on
    plot(xiO, fO, '--', 'Color', colors.W, 'LineWidth', 3);
    grid on
    ylabel(caseTitles{c}, 'FontSize', 32)
    set(gca,'FontSize',22,'LineWidth',2)
    % Reference line at true value
    yL = ylim; plot([true_value true_value], yL, 'k--', 'LineWidth', 2);
    ylim(yL); xlim(xlim_top);

    % --- Bottom row: CI widths ---
    nexttile(c+3);
    [fTw,xiTw] = ksdensity(T.(ci_varname));
    [fOw,xiOw] = ksdensity(O.(ci_varname));
    plot(xiTw, fTw, 'Color', colors.X, 'LineWidth', 3); hold on
    plot(xiOw, fOw, '--', 'Color', colors.W, 'LineWidth', 3);
    grid on
    set(gca,'FontSize',22,'LineWidth',2)
    xlim(xlim_bottom);
end

% Shared legend at top
nexttile(4)
lg = legend({'True Data','Observed Data'}, 'Orientation','vertical',Location='northeast');
lg.FontSize = 22;
nexttile(3)
xlabel('Parameter Estimate','FontSize',22);
nexttile(6)
xlabel('Confidence Interval Width','FontSize',22);


t.YLabel.String = 'Density';  % hide built-in shared label
t.YLabel.FontSize = 34;

% % --- Align all y‑labels in the top row ---
% drawnow  % ensure positions are computed
% minX = max(arrayfun(@(ax) ax.YLabel.Position(1), ax_top)); % leftmost position
% for k = 1:numel(ax_top)
%     pos = ax_top(k).YLabel.Position;
%     pos(1) = minX;
%     ax_top(k).YLabel.Position = pos;
% end

% exportgraphics(gcf,'Figures/linRegression_Densities.png')
% saveas(gcf,'Figures/linRegression_Densities.eps','epsc')


end


% Plot slope (b_hat)87
plot_param_panels('b_hat', 'CIwidth_b', 1, ...
    cases_true, cases_obs, caseTitles, colors, ...
    'Slope estimates (top) and CI widths (bottom)');

%% Plot intercept (a_hat)
plot_param_panels('a_hat', 'CIwidth_a', 0, ...
    cases_true, cases_obs, caseTitles, colors, ...
    'Intercept estimates (top) and CI widths (bottom)');
