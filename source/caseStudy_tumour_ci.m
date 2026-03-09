% Investigate an estimate of carrying capacity for Gompertzian tumour
% growth. Data is incomplete and the steady state is not observed
% All data has units of tumour cells x10^6 (or tumour volume in mm^3)
rng(100)
model = @(theta,t) theta(:,3) .* (theta(:,1)./theta(:,3)).^(exp(-theta(:,2).*t));

% Load tumour data
load tumourData.mat

tumourData = origdata(origdata.experiment == 1,:);
tumourData = tumourData(tumourData.mouse == 6,:);

time = tumourData.time;
tumourSize = tumourData.tumour;

lw = 6;
fs = 36;

%% Fit to the single mouse data to get parameters

ft = fittype(@(N0,r,L,t) model([N0,r,L],t), 'independent','t');
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',zeros(1,3), ...
               'StartPoint', [100, 1, 1000]);

[curve,gof,output] = fit(time,tumourSize,ft,fo);

residuals=output.residuals;

figure(1)
plot(curve,time,tumourSize);

% Plot residuals
figure(2)
plot(curve,time,tumourSize,"residuals")

% standard deviation of the residuals
sigma = std(residuals);

fittedParam = coeffvalues(curve);
scaledParam = fittedParam .* [1e6, 1, 1e6];
t = linspace(0,max(time),1001);
figure(3) % Scaled back to 10^6
hold on
scatter(time,tumourSize.*10^6,'rx');
plot(t,model(scaledParam,t))
hold off

%% Create synthetic data for 6 tumours
% Carrying capacity is given (tumour grown in constrained environment)
% Sample initial condition from uniform (like Kim et al measurements begin
% when a tumour is between given sizes
% Sample growth rate from normal distribution
% Want to infer the mean growth rate
% Do inference with perfect data
% Do inference with correlated time errors, i.e., each mouse is measured
% sequentially and it takes some mean time to take measurements for each
% mouse. Consider two cases, each day that measurements are made the mice
% are always measured in the same order, and each day the mice are measured
% in random order

numSamples = 6;

% sample param values for each mouse
theta = nan(numSamples, 3);
% Sample N0 initial value
theta(:,1) = unifrnd(100,120,[numSamples,1]);
% sample r growth rate
mur = 0.08; % mean growth rate
sigmar = 0.02; % growth rate std
theta(:,2) = normrnd(mur,sigmar,[numSamples,1]);
% sample carrying capacity L
% Will set all to be equal
L_true = 4000;
theta(:,3) = L_true .* ones(numSamples,1);

% Generate synthetic data
perfectData_example = model(theta,time')';
sampleTime = repmat(time,6,1);
perfectData_example = perfectData_example(:);

% plot the sample data
figure(4) 
scatter(sampleTime,perfectData_example)

%% Fit to the perfect Data

ft = fittype(@(r,L,t) model([mean(perfectData_example(1,:)),r,L],t), 'independent','t');
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',zeros(1,2), ...
               'StartPoint', [1, 1000]);
[curve_perfect,gof_perfect,output_perfect] = fit(sampleTime,perfectData_example,ft,fo);

width = 900;
height = 700;
left = 50;
bottom = 80;

fig = figure(5);
fig.Position = [left, bottom, width, height];
plot(curve_perfect,'k--',sampleTime,perfectData_example,'k^');
legend('off')
xlabel('Time (days)')
ylabel('Tumour Cells')

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth',2);
set(findobj(gcf,'type','line'),'LineWidth',3,'MarkerSize',10);



%% Generate correlated time error data 

% Generate time error data
% muX = 1/48 (average about 30 min to get measurement)
% sigmaX = 0.05 (seems to give a decent looking distribution
muX = 1/6;
sigmaX = 0.1;

% Visualise the error distribution
figure(3)
x = linspace(0,2,1001);
plot(x,lognpdf(x,log(muX.^2 ./ sqrt(muX.^2 + sigmaX.^2)),sqrt(log(1 + sigmaX.^2 ./ (muX.^2)))));

[orderedTime, randomOrderTime] = generate_correlatedTimes(time,muX,sigmaX,numSamples,'lognormal');

orderedSynthData = nan(length(time),numSamples);
randomOrderSynthData = orderedSynthData;
% run model
for i = 1:numSamples
    orderedSynthData(:,i) = model(theta(i,:),orderedTime(:,i));
    randomOrderSynthData(:,i) = model(theta(i,:),randomOrderTime(:,i));
end
    
%% Fit time error data

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',zeros(1,3), ...
               'StartPoint', [110, 1, 1000]);

ft = fittype(@(a,r,L,t) model([a,r,L],t), 'independent','t');
[curve_ordered,gof_ordered,output_ordered] = fit(sampleTime,orderedSynthData(:),ft,fo);

ft = fittype(@(a,r,L,t) model([a,r,L],t), 'independent','t');
[curve_randomOrder,gof_randomOrder,output_randomOrder] = fit(sampleTime,randomOrderSynthData(:),ft,fo);

ft = fittype(@(a,r,L,t) model([a,r,L],t), 'independent','t');
[curve_perfect,gof_perfect,output_perfect] = fit(sampleTime,perfectData_example,ft,fo);

curve_ordered; 
curve_randomOrder; 
curve_perfect; 

gof_ordered; 
gof_randomOrder; 
gof_perfect; 

fig = figure(6);
fig.Position = [left, bottom, width, height];
plot(curve_ordered,'b',sampleTime,orderedSynthData,'bo');
hold on
plot(curve_randomOrder,'r--',sampleTime,randomOrderSynthData,'rx');
plot(curve_perfect,'k-.',sampleTime,perfectData_example,'k^');
hold off

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth',2);
set(findobj(gcf,'type','line'),'LineWidth',3,'MarkerSize',10);

%% Fit time error data - set start point as mean of data

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',zeros(1,2), ...
               'StartPoint', [1, 1000]);

ft = fittype(@(r,L,t) model([mean(orderedSynthData(1,:)),r,L],t), 'independent','t');
[curve_ordered,gof_ordered,output_ordered] = fit(sampleTime,orderedSynthData(:),ft,fo);

ft = fittype(@(r,L,t) model([mean(randomOrderSynthData(1,:)),r,L],t), 'independent','t');
[curve_randomOrder,gof_randomOrder,output_randomOrder] = fit(sampleTime,randomOrderSynthData(:),ft,fo);

ft = fittype(@(r,L,t) model([mean(perfectData_example(1,:)),r,L],t), 'independent','t');
[curve_perfect,gof_perfect,output_perfect] = fit(sampleTime,perfectData_example,ft,fo);

curve_ordered; 
curve_randomOrder; 
curve_perfect; 

gof_ordered; 
gof_randomOrder; 
gof_perfect; 

%% -------- Figure 8  --------
% CI appearance
ciLevel   = 0.95;   % 95% confidence interval for mean response
fillAlpha = 0.20;   % transparent shading

fig = figure(8);
fig.Position = [left, bottom, width, height];

% Ordered measurements (blue solid)
hOrdAll = plot(curve_ordered,'b', sampleTime, orderedSynthData,'bo');
hold on
hOrd = hOrdAll(1);                      % first handle is the fitted curve
set(hOrd,'LineWidth',3);
xg_ord = get(hOrd,'XData'); xg_ord = xg_ord(:);
yci_ord = predint(curve_ordered, xg_ord, ciLevel, 'functional');
ciPatch(gca, xg_ord, yci_ord(:,1), yci_ord(:,2), 'b', fillAlpha);

% Random order measurements (red dashed)
hRandAll = plot(curve_randomOrder,'r--', sampleTime, randomOrderSynthData,'rx');
hRand = hRandAll(1);
set(hRand,'LineWidth',3);
xg_rand = get(hRand,'XData'); xg_rand = xg_rand(:);
yci_rand = predint(curve_randomOrder, xg_rand, ciLevel, 'functional');
ciPatch(gca, xg_rand, yci_rand(:,1), yci_rand(:,2), 'r', fillAlpha);

% Perfect data fit (black dash-dot)
hPerfAll = plot(curve_perfect,'k-.', sampleTime, perfectData_example,'k^');
hPerf = hPerfAll(1);
set(hPerf,'LineWidth',3);
xg_perf = get(hPerf,'XData'); xg_perf = xg_perf(:);
yci_perf = predint(curve_perfect, xg_perf, ciLevel, 'functional');
ciPatch(gca, xg_perf, yci_perf(:,1), yci_perf(:,2), 'k', fillAlpha);

% Legend with only the fitted lines
lgd = legend([hOrd, hRand, hPerf], ...
             {'Ordered Measurements','Random Order Measurements','Perfect Data'}, ...
             'Location','southeast');

hold off
xlabel('Time (days)')
ylabel('Tumour Cells')

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth',2);
set(findobj(gcf,'type','line'),'LineWidth',3,'MarkerSize',10);


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
    % keep shading behind markers/lines
    try
        uistack(h,'bottom');
    catch
        % ignore if uistack unavailable
    end
end
