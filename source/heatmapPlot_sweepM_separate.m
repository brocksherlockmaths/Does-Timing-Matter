model = @(theta,t) (theta(1) - theta(2)) .* exp(-theta(3).*t) + theta(2);
time = linspace(0,60,101);

% Get teh param space
a_range = [0.05, 0.1, 0.15, 0.2, 0.25];
M_range = 1:7:35;
k_range = [0.05, 0.1, 0.15, 0.2, 0.25];

paramSpace = combvec(a_range,M_range,k_range);

% Get the indices where a and M are fixed and k is altered
[r,c] = find(paramSpace(1,:) == 0.1 & paramSpace(3,:) == 0.1);

indices = c;

% Set all interpreters to LaTeX by default
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, 'defaultTextInterpreter',          'latex');

%% size of figures to display
width = 1500;
height = 900;
left = 50;
bottom = 80;

load batlow.mat

filePrefix = 'results_heatmap/results_timeError_lognormalData_heatMapCalc_paramSweep';

load(strcat(filePrefix,'1.mat'),'sigma','epsilon','num_params','params');

numGridpoints = 31;

%% Collect interval data for each param set
intervalWidth = nan(numGridpoints^2, 4, length(indices));
mlEstimate = intervalWidth;
intervalContained_cont1 = intervalWidth;
intervalContained_cont2 = intervalWidth;
allTrueVals = intervalWidth;


% Loop through each param set
for i = 1:length(indices)
    load(strcat(filePrefix,num2str(indices(i)),'.mat'),'lowerBound','upperBound','MLE','trueVals','SIGMA_SET');
    thisTrueVals = [repmat(trueVals,numGridpoints^2,1),SIGMA_SET];
    allTrueVals(:,:,i) = thisTrueVals;
    intervalWidth(:,:,i) = upperBound - lowerBound;
    truthContained = isbetween(thisTrueVals,lowerBound,upperBound);
    mlEstimate(:,:,i) = MLE;
    intervalMidpoint = (lowerBound + upperBound) / 2;
    intervalContained_cont1(:,:,i) = abs(thisTrueVals - intervalMidpoint) ./ intervalWidth(:,:,i);
    tmpIntervalContained = min(abs(thisTrueVals - lowerBound), abs(thisTrueVals-upperBound)); % distance to closest interval bound
    tmpIntervalContained(truthContained == 1) = 0;
    intervalContained_cont2(:,:,i) = tmpIntervalContained;
end

%% Plot the interval width
%% === SETTINGS for square, matching tiles ===
tileSize_cm = 4;   % side length of each tile in cm
nCols = 5;
nRows_top = 3;     % top figure rows
nRows_bottom = 1;  % bottom figure rows

figWidth_cm  = tileSize_cm * nCols;
figHeight_top_cm    = tileSize_cm * nRows_top;
figHeight_bottom_cm = tileSize_cm * nRows_bottom;

%% === TOP FIGURE ===
fig1 = figure(1);
fig1.Units = 'centimeters';
fig1.Position = [2, 2, figWidth_cm, figHeight_top_cm];

t1 = tiledlayout(fig1, nRows_top, nCols);
t1.TileIndexing = "columnMajor";
t1.Padding = 'loose';
t1.TileSpacing = 'loose';
colormap(batlow)

clims = [min(intervalWidth,[],'all'), max(intervalWidth,[],'all')];
levels = linspace(clims(1), clims(2), 21);

% Plot only the top 3 rows (i = 1:3)
for j = 1:length(indices)
    for i = 1:(num_params-1)
        hm = reshape(intervalWidth(:,i,j), numGridpoints, numGridpoints);
        nexttile
        if i ~= 4
            imagesc(epsilon, sigma, hm, clims);
            yticks(0:0.25:0.5)
        end
        set(gca,'YDir','normal');
        if j == 1
            if i ~= 4
                ylabel(strcat('$',params{i},'$'));
            else
                ylabel('PM GLUT4')
            end
        end
        if i == 1
            title(strcat("$M = ", num2str(M_range(j)),'$'));
        end
    end
end

% Tick label suppression
ax = flipud(findall(t1, 'Type', 'axes'));
leftColIdx   = 1:nRows_top;
bottomRowIdx = nRows_top:nRows_top:(nRows_top*nCols);
for k = 1:numel(ax)
    if ~ismember(k, leftColIdx)
        ax(k).YTickLabel = [];
    end
    if ~ismember(k, bottomRowIdx)
        ax(k).XTickLabel = [];
    end
end
set(ax, 'TickDir', 'out');

% Shared labels
%t1.Title.String = 'Interval Width';
t1.XLabel.String = '\Delta';
t1.YLabel.String = '\sigma_Y';
t1.Title.FontSize = 24;
t1.XLabel.FontSize = 24;
t1.YLabel.FontSize = 24;

% Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';

% Axes font/line
set(findobj(fig1,'type','axes'), 'FontSize', 16, 'LineWidth', 1);

% Export
set(fig1, 'PaperPositionMode', 'auto');
print(fig1, 'Figures_withDynamics/heatmap_sweepM_intervalWidth_top.png', ...
      '-dpng', '-r300');

% === BOTTOM FIGURE ===
fig2 = figure(8);
fig2.Units = 'centimeters';
fig2.Position = [2, 2, fig1.Position(3), figHeight_bottom_cm+2];

t2 = tiledlayout(fig2, nRows_bottom, nCols, 'TileIndexing', 'columnMajor');
t2.Padding = t1.Padding;
t2.TileSpacing = t1.TileSpacing;

% Plot only the bottom row (i = 4)
for j = 1:length(indices)
    nexttile
    plot(time, model(allTrueVals(1,1:3,j), time), 'LineWidth', 2);
    title(strcat("$M = ", num2str(M_range(j)),'$'));
    xticks(0:20:60)
    xtickangle(0)
    yticks(0:7:29)
    ylim([0,29])
    grid on
end

% Tick label suppression
ax2 = flipud(findall(t2, 'Type', 'axes'));
leftColIdx   = 1:nRows_bottom;
bottomRowIdx = nRows_bottom:nRows_bottom:(nRows_bottom*nCols);
for k = 1:numel(ax2)
    if ~ismember(k, leftColIdx)
        ax2(k).YTickLabel = [];
    end
    if ~ismember(k, bottomRowIdx)
        ax2(k).XTickLabel = [];
    end
end

% Dummy colorbar for alignment
cb2 = colorbar;
cb2.Layout.Tile = 'east';
cb2.Visible = 'off';

% Match horizontal alignment
t2.InnerPosition(1) = t1.InnerPosition(1);
t2.InnerPosition(3) = t1.InnerPosition(3);

% Adjust outer position for shared X label
outerPos = t2.OuterPosition;
outerPos(2) = outerPos(2) - 0.05;
outerPos(4) = outerPos(4) + 0.05;
t2.OuterPosition = outerPos;

set(findobj(fig2,'type','axes'), 'FontSize', 24, 'LineWidth', 1);


% Tick direction out
set(ax2, 'TickDir', 'out');

% Shared labels
t2.XLabel.String = '$t$';
t2.XLabel.Interpreter = 'latex';
t2.YLabel.String = '$P(t)$';
t2.YLabel.Interpreter = 'latex';
t2.Title.FontSize = 24;
t2.XLabel.FontSize = 24;
t2.YLabel.FontSize = 24;

% Export
set(fig2, 'PaperPositionMode', 'auto');
print(fig2, 'Figures_withDynamics/heatmap_sweepM_intervalWidth_bottom.png', ...
      '-dpng', '-r300');

% Set region where true value not in interval to nan so it appears blank

%% Plot the MLE estimate relative error
fig = figure(2);
fig.Position = [left, bottom, width, height];
t = tiledlayout(4,5);
t.TileIndexing = "columnMajor";
colormap(batlow)
mlLims_bottom = min((mlEstimate - allTrueVals) ./ allTrueVals,[],'all');
mlLims_top = (mlEstimate - allTrueVals) ./ allTrueVals;
mlLims_top(mlLims_top == Inf) = 0;
clims = [mlLims_bottom,...
    max(mlLims_top,[],'all','omit')]
% clims = [-0.1,5]
for j = 1:length(indices)
    for i = 1:num_params
        % Reshape the data to grid
        hm = reshape(mlEstimate(:,i,j),numGridpoints,numGridpoints);
        % get relative error
        if i <= 3 % param not sigma
            hm = (hm - allTrueVals(1,i,j)) ./ allTrueVals(1,i,j);
        else % param is sigma
            hm = (hm - sigma') ./ sigma';
        end
        nexttile
        if i ~= 4
            imagesc(epsilon,sigma,hm,clims);
            yticks(0:0.25:0.5)
        else % plot the dynamics
            plot(time,model(allTrueVals(1,1:3,j),time),'LineWidth',2);
            xticks(0:10:60)
            yticks(0:7:29)
            ylim([0,29])
            grid on
        end
        set(gca,'YDir','normal');
        if j == 1
            if i ~= 4
                ylabel(strcat('$',params{i},'$'));
            else
                ylabel('PM GLUT4')
            end
        end
        if i == 1
            title(strcat("M = ",num2str(M_range(j))));
        end
    end
end
t.Title.String = 'MLE Estimate Relative Error';
t.XLabel.String = 'Potential Delay (\Delta) - Bottom Row Time';
t.YLabel.String = '\sigma';
t.Title.FontSize = 24;
t.XLabel.FontSize = 24;
t.YLabel.FontSize = 24;
cb = colorbar;
cb.Layout.Tile = 'east';

%% Plot the MLE estimate absolute error
%% === SETTINGS for square, matching tiles ===
tileSize_cm = 4;   % side length of each tile in cm
nCols = 5;
nRows_top = 3;     % top figure rows
nRows_bottom = 1;  % bottom figure rows

figWidth_cm  = tileSize_cm * nCols;
figHeight_top_cm    = tileSize_cm * nRows_top;
figHeight_bottom_cm = tileSize_cm * nRows_bottom;

%% === TOP FIGURE ===
fig1 = figure(3);
fig1.Units = 'centimeters';
fig1.Position = [2, 2, figWidth_cm, figHeight_top_cm];

t1 = tiledlayout(fig1, nRows_top, nCols);
t1.TileIndexing = "columnMajor";
t1.Padding = 'loose';
t1.TileSpacing = 'loose';
colormap(batlow)

% Colour limits
mlLims_bottom = min(abs(mlEstimate - allTrueVals),[],'all');
mlLims_top = abs(mlEstimate - allTrueVals);
mlLims_top(mlLims_top == Inf) = 0;
clims = [mlLims_bottom, max(mlLims_top,[],'all','omit')];

% Plot only the top 3 rows (i = 1:3)
for j = 1:length(indices)
    for i = 1:(num_params-1)
        hm = reshape(mlEstimate(:,i,j), numGridpoints, numGridpoints);
        if i <= 3
            hm = (hm - allTrueVals(1,i,j));
        else
            hm = (hm - sigma');
        end
        nexttile
        if i ~= 4
            imagesc(epsilon, sigma, hm, clims);
            yticks(0:0.25:0.5)
        end
        set(gca,'YDir','normal');
        if j == 1
            if i ~= 4
                ylabel(strcat('$',params{i},'$'));
            else
                ylabel('PM GLUT4')
            end
        end
        if i == 1
            title(strcat("$M = ", num2str(M_range(j)),'$'));
        end
    end
end

% Tick label suppression
ax = flipud(findall(t1, 'Type', 'axes'));
leftColIdx   = 1:nRows_top;
bottomRowIdx = nRows_top:nRows_top:(nRows_top*nCols);
for k = 1:numel(ax)
    if ~ismember(k, leftColIdx)
        ax(k).YTickLabel = [];
    end
    if ~ismember(k, bottomRowIdx)
        ax(k).XTickLabel = [];
    end
end
set(ax, 'TickDir', 'out');

% Shared labels
%t1.Title.String = 'MAP Estimate Absolute Error';
t1.XLabel.String = '\Delta';
t1.YLabel.String = '\sigma_Y';
t1.Title.FontSize = 24;
t1.XLabel.FontSize = 24;
t1.YLabel.FontSize = 24;

% Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';

% Axes font/line
set(findobj(fig1,'type','axes'), 'FontSize', 16, 'LineWidth', 1);

% Export
set(fig1, 'PaperPositionMode', 'auto');
print(fig1, 'Figures_withDynamics/heatmap_sweepM_mleAbsolute_top.png', ...
      '-dpng', '-r300');

% === BOTTOM FIGURE ===
fig2 = figure(9);
fig2.Units = 'centimeters';
fig2.Position = [2, 2, fig1.Position(3), figHeight_bottom_cm];

t2 = tiledlayout(fig2, nRows_bottom, nCols, 'TileIndexing', 'columnMajor');
t2.Padding = t1.Padding;
t2.TileSpacing = t1.TileSpacing;

% Plot only the bottom row (i = 4)
for j = 1:length(indices)
    nexttile
    plot(time, model(allTrueVals(1,1:3,j), time), 'LineWidth', 2);
    xticks(0:20:60)
    yticks(0:7:29)
    ylim([0,29])
    grid on
end

% Tick label suppression
ax2 = flipud(findall(t2, 'Type', 'axes'));
leftColIdx   = 1:nRows_bottom;
bottomRowIdx = nRows_bottom:nRows_bottom:(nRows_bottom*nCols);
for k = 1:numel(ax2)
    if ~ismember(k, leftColIdx)
        ax2(k).YTickLabel = [];
    end
    if ~ismember(k, bottomRowIdx)
        ax2(k).XTickLabel = [];
    end
end

% Dummy colorbar for alignment
cb2 = colorbar;
cb2.Layout.Tile = 'east';
cb2.Visible = 'off';

% Match horizontal alignment
t2.InnerPosition(1) = t1.InnerPosition(1);
t2.InnerPosition(3) = t1.InnerPosition(3);

% Adjust outer position for shared X label
outerPos = t2.OuterPosition;
outerPos(2) = outerPos(2) - 0.05;
outerPos(4) = outerPos(4) + 0.05;
t2.OuterPosition = outerPos;

% Tick direction out
set(ax2, 'TickDir', 'out');

% Shared labels
t2.XLabel.String = '$t$';
t2.XLabel.Interpreter = 'latex';
t2.YLabel.String = '$y(t)$';
t2.YLabel.Interpreter = 'latex';
t2.Title.FontSize = 24;
t2.XLabel.FontSize = 24;
t2.YLabel.FontSize = 24;

% Export
set(fig2, 'PaperPositionMode', 'auto');
print(fig2, 'Figures_withDynamics/heatmap_sweepM_mleAbsolute_bottom.png', ...
      '-dpng', '-r300');

%% Plot how close the true values are to being contained in the interval V1
fig = figure(4);
fig.Position = [left, bottom, width, height];
t = tiledlayout(4,5);
t.TileIndexing = "columnMajor";
colormap(batlow)
clims = [min(intervalContained_cont1,[],'all'),max(intervalContained_cont1,[],'all')];
clims = [1e-4, 5]
levels = linspace(clims(1),clims(2),21);
for j = 1:length(indices)
    for i = 1:num_params
        % Reshape the data to grid
        hm = reshape(intervalContained_cont1(:,i,j),numGridpoints,numGridpoints);
        nexttile
        if i ~= 4
            imagesc(epsilon,sigma,hm,clims);
            yticks(0:0.25:0.5)
        else % plot the dynamics
            plot(time,model(allTrueVals(1,1:3,j),time),'LineWidth',2);
            xticks(0:10:60)
            yticks(0:7:29)
            ylim([0,29])
            grid on
        end
        % contourf(epsilon,sigma,hm,levels);
        set(gca,'YDir','normal');
        if j == 1
            if i ~= 4
                ylabel(strcat('$',params{i},'$'));
            else
                ylabel('PM GLUT4')
            end
        end
        if i == 1
            title(strcat("M = ",num2str(M_range(j))));
        end
    end
end
t.Title.String = 'True Value in Interval V1';
t.XLabel.String = 'Potential Delay (\Delta) - Bottom Row Time';
t.YLabel.String = '\sigma';
t.Title.FontSize = 24;
t.XLabel.FontSize = 24;
t.YLabel.FontSize = 24;
cb = colorbar;
cb.Layout.Tile = 'east';

%% Plot how close the true values are to being contained in the interval V2
%% === SETTINGS for square, matching tiles ===
tileSize_cm = 4;   % side length of each tile in cm
nCols = 5;
nRows_top = 3;     % top figure rows
nRows_bottom = 1;  % bottom figure rows

figWidth_cm  = tileSize_cm * nCols;
figHeight_top_cm    = tileSize_cm * nRows_top;
figHeight_bottom_cm = tileSize_cm * nRows_bottom;

%% === TOP FIGURE ===
fig1 = figure(5);
fig1.Units = 'centimeters';
fig1.Position = [2, 2, figWidth_cm, figHeight_top_cm];

t1 = tiledlayout(fig1, nRows_top, nCols);
t1.TileIndexing = "columnMajor";
t1.Padding = 'loose';
t1.TileSpacing = 'loose';
colormap([1, 1, 1; batlow])

% Colour limits
clims = [0, max(intervalContained_cont2,[],'all')];
levels = linspace(clims(1), clims(2), 21);

% Plot only the top 3 rows (i = 1:3)
for j = 1:length(indices)
    for i = 1:(num_params-1)
        hm = reshape(intervalContained_cont2(:,i,j), numGridpoints, numGridpoints);
        nexttile
        if i ~= 4
            imagesc(epsilon, sigma, hm, clims);
            yticks(0:0.25:0.5)
        end
        set(gca,'YDir','normal');
        if j == 1
            if i ~= 4
                ylabel(strcat('$',params{i},'$'));
            else
                ylabel('PM GLUT4')
            end
        end
        if i == 1
            title(strcat("$M = ", num2str(M_range(j)),'$'));
        end
    end
end

% Tick label suppression
ax = flipud(findall(t1, 'Type', 'axes'));
leftColIdx   = 1:nRows_top;
bottomRowIdx = nRows_top:nRows_top:(nRows_top*nCols);
for k = 1:numel(ax)
    if ~ismember(k, leftColIdx)
        ax(k).YTickLabel = [];
    end
    if ~ismember(k, bottomRowIdx)
        ax(k).XTickLabel = [];
    end
end
set(ax, 'TickDir', 'out');

% Shared labels
%t1.Title.String = 'Distance of True Value to CI';
t1.XLabel.String = '\Delta';
t1.YLabel.String = '\sigma_Y';
t1.Title.FontSize = 24;
t1.XLabel.FontSize = 24;
t1.YLabel.FontSize = 24;

% Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';

% Axes font/line
set(findobj(fig1,'type','axes'), 'FontSize', 16, 'LineWidth', 1);

% Export
set(fig1, 'PaperPositionMode', 'auto');
print(fig1, 'Figures_withDynamics/heatmap_sweepM_intervalDist_top.png', ...
      '-dpng', '-r300');

% === BOTTOM FIGURE ===
fig2 = figure(10);
fig2.Units = 'centimeters';
fig2.Position = [2, 2, fig1.Position(3), figHeight_bottom_cm];

t2 = tiledlayout(fig2, nRows_bottom, nCols, 'TileIndexing', 'columnMajor');
t2.Padding = t1.Padding;
t2.TileSpacing = t1.TileSpacing;

% Plot only the bottom row (i = 4)
for j = 1:length(indices)
    nexttile
    plot(time, model(allTrueVals(1,1:3,j), time), 'LineWidth', 2);
    xticks(0:20:60)
    yticks(0:7:29)
    ylim([0,29])
    grid on
end

% Tick label suppression
ax2 = flipud(findall(t2, 'Type', 'axes'));
leftColIdx   = 1:nRows_bottom;
bottomRowIdx = nRows_bottom:nRows_bottom:(nRows_bottom*nCols);
for k = 1:numel(ax2)
    if ~ismember(k, leftColIdx)
        ax2(k).YTickLabel = [];
    end
    if ~ismember(k, bottomRowIdx)
        ax2(k).XTickLabel = [];
    end
end

% Dummy colorbar for alignment
cb2 = colorbar;
cb2.Layout.Tile = 'east';
cb2.Visible = 'off';

% Match horizontal alignment
t2.InnerPosition(1) = t1.InnerPosition(1);
t2.InnerPosition(3) = t1.InnerPosition(3);

% Adjust outer position for shared X label
outerPos = t2.OuterPosition;
outerPos(2) = outerPos(2) - 0.05;
outerPos(4) = outerPos(4) + 0.05;
t2.OuterPosition = outerPos;

% Tick direction out
set(ax2, 'TickDir', 'out');

% Shared labels
t2.XLabel.String = '$t$';
t2.XLabel.Interpreter = 'latex';
t2.YLabel.String = '$P(t)$';
t2.YLabel.Interpreter = 'latex';
t2.Title.FontSize = 24;
t2.XLabel.FontSize = 24;
t2.YLabel.FontSize = 24;

% Export
set(fig2, 'PaperPositionMode', 'auto');
print(fig2, 'Figures_withDynamics/heatmap_sweepM_intervalDist_bottom.png', ...
      '-dpng', '-r300');

% %% Plot the interval width with contours for each rate on same plot
% fig = figure(1);
% fig.Position = [left, bottom, width, height/2];
% t = tiledlayout(1,4);
% t.TileIndexing = "columnMajor";
% colormap(batlow)
% clims = [min(intervalWidth,[],'all'),max(intervalWidth,[],'all')];
% levels = [0.02 0.02];
% for i = 1:num_params
%     nexttile
%     hold on
%     for j = 1:length(indices)
%         % Reshape the data to grid
%         hm = reshape(intervalWidth(:,i,j),numGridpoints,numGridpoints);
% 
%         %imagesc(epsilon,sigma,hm,clims);
%         contour(epsilon,sigma,hm,levels);
%         set(gca,'YDir','normal');
% 
%     end
%     title(params{i});
%     hold off
% end
% t.Title.String = 'Interval Width';
% t.XLabel.String = '\epsilon';
% t.YLabel.String = '\sigma';
% t.Title.FontSize = 24;
% t.XLabel.FontSize = 24;
% t.YLabel.FontSize = 24;
% cb = colorbar;
% cb.Layout.Tile = 'east';
% 
% %% Plot the MLE estimate with contours for each rate on same plot
% fig = figure(1);
% fig.Position = [left, bottom, width, height/2];
% t = tiledlayout(1,4);
% t.TileIndexing = "columnMajor";
% colormap(batlow)
% colours = crameri('batlow',length(indices));
% mlLims_bottom = min(abs(mlEstimate - allTrueVals),[],'all');
% mlLims_top = abs(mlEstimate - allTrueVals);
% mlLims_top(mlLims_top == Inf) = 0;
% clims = [mlLims_bottom,...
%     max(mlLims_top,[],'all','omit')];
% levels = [0.015 0.015];
% for i = 1:num_params
%     nexttile
%     hold on
%     for j = 1:length(indices)
%         % Reshape the data to grid
%         hm = reshape(mlEstimate(:,i,j),numGridpoints,numGridpoints);
%         % get relative error
%         if i <= 3 % param not sigma
%             hm = (hm - allTrueVals(1,i,j)) ./ allTrueVals(1,i,j);
%         else % param is sigma
%             hm = (hm - sigma') ./ sigma';
%         end
% 
%         %imagesc(epsilon,sigma,hm,clims);
%         contour(epsilon,sigma,hm,levels,'color',colours(j,:));
%         set(gca,'YDir','normal');
% 
%     end
%     title(params{i});
%     hold off
% end
% t.Title.String = 'Interval Width';
% t.XLabel.String = '\epsilon';
% t.YLabel.String = '\sigma';
% t.Title.FontSize = 24;
% t.XLabel.FontSize = 24;
% t.YLabel.FontSize = 24;
% cb = colorbar;
% cb.Layout.Tile = 'east';