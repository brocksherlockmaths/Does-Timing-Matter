%% Set up the parasitemia model

rng(5)
fs = 20;

% CI appearance 
ciLevel   = 0.95;   % 95% confidence interval for mean response
fillAlpha = 0.20;   

model = @(theta,t) theta(2).*theta(1).*exp(theta(4).*t) + theta(3)*exp(-theta(5).*t);

% Param values for synth data
V0 = 0.002;
f = 1;
B0 = 0.02;
g = 0.531;
k = 0.323;

theta = [V0,f,B0,g,k];

%% Run some sample data and plot

t = 0:2:24;

exampleData = model(theta,t);

scatter(t,exampleData)

%% Noisy Data
numSamples = 1;
sigma = 1/4;
exampleData_noise = normalGenerator_noTimeError(theta,t,sigma,model,numSamples);

scatter(t,exampleData_noise);
hold on
plot(t,exampleData);
hold off

%% Time error data

numSamples = 1;
sigma = 0;
epsilon = 1/4;
exampleData_timeError = normalGenerator_normalTimeError(theta,t,sigma,model,epsilon,numSamples);

scatter(t,exampleData_timeError);
hold on
plot(t,exampleData);
hold off

%% Fit model to data no time error
% Set the initial background value. Know everything about background and
% are only inferring parasitemia

ft = fittype(@(V0,B0,g,k,t) model([V0,f,B0,g,k],t), 'independent','t');
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',zeros(1,4));

fitData = exampleData_noise;

[curve1,gof1] = fit(t',fitData',ft,fo)

figure(1)
% Plot curve1 with its data and add CI band
hAll1 = plot(curve1,t,fitData); 
hold on
hLine1 = hAll1(1);
xg1 = get(hLine1,'XData'); xg1 = xg1(:);
yci1 = predint(curve1, xg1, ciLevel, 'functional'); 
ciPatch(gca, xg1, yci1(:,1), yci1(:,2), get(hLine1,'Color'), fillAlpha);

%% Fit model to data time error
fitData = exampleData_timeError;

[curve2,gof2] = fit(t',fitData',ft,fo)

% Plot curve2 on top (keep your style; default red below)
hLine2 = plot(curve2,'r'); 
xg2 = get(hLine2(1),'XData'); xg2 = xg2(:);
yci2 = predint(curve2, xg2, ciLevel, 'functional');
ciPatch(gca, xg2, yci2(:,1), yci2(:,2), get(hLine2(1),'Color'), fillAlpha);

hold off

%% Layout & export settings for later figures
width = 900;
height = 700;
left = 50;
bottom = 80;

%% Figure 2: true model vs curve2 fit (with CI) + inset (with CI)
fig = figure(2);
fig.Position = [left, bottom, width, height];
time = linspace(0,25,1001);

% Main axes
pTrue = plot(time,model(theta,time),'k'); 
hold on
hMain2 = plot(curve2,'r--'); % fitted line on main axes
plot(t,fitData,'rx')
xg_main2 = get(hMain2(1),'XData'); xg_main2 = xg_main2(:);
yci_main2 = predint(curve2, xg_main2, ciLevel, 'functional');
ciPatch(gca, xg_main2, yci_main2(:,1), yci_main2(:,2), get(hMain2(1),'Color'), fillAlpha);
hold off

legend('True Model', 'Fit Model','Data','Location','northwest')
xlabel('Time (days)')
ylabel('Parasite Numbers')

ylim([-20,1200])

% --- Inset (with CI) ---
insetWidth  = 0.35;
insetHeight = 0.35;
insetLeft   = 0.20;
insetBottom = 0.20;

insetAxes = axes('Position', [insetLeft insetBottom insetWidth insetHeight]);

zoomX = [0 5];   % time range
zoomY = [0 200]; % parasite number range

zoomIdx = time >= zoomX(1) & time <= zoomX(2) & ...
          model(theta,time) >= zoomY(1) & model(theta,time) <= zoomY(2);

plot(time(zoomIdx), model(theta,time(zoomIdx)), 'k', 'LineWidth', 2); hold on
% Inset fit line
plot(time(zoomIdx), feval(curve2,time(zoomIdx))', 'r--', 'LineWidth', 2);
% Inset CI (compute on zoom grid)
yci_inset2 = predint(curve2, time(zoomIdx)', ciLevel, 'functional');
ciPatch(insetAxes, time(zoomIdx)', yci_inset2(:,1), yci_inset2(:,2), [1 0 0], fillAlpha);
% Inset data points
maskPoints = t >= zoomX(1) & t <= zoomX(2) & ...
             fitData >= zoomY(1) & fitData <= zoomY(2);
plot(t(maskPoints), fitData(maskPoints), 'rx', 'LineWidth', 1.5);
ylim([0,0.06])
xlim([0,5])
hold off

% Style inset
set(insetAxes, 'FontSize', 8, 'Box', 'on'); % small font for tick labels

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth',2);
set(findobj(gcf,'type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(gcf,'type','axes'),'FontSize',fs,'LabelFontSizeMultiplier',1.8);


%% Fit only initial parasite load
ft = fittype(@(V0,g,t) model([V0,f,B0,g,k],t), 'independent','t');
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',zeros(1,2));

[curve3,gof3] = fit(t',fitData',ft,fo)

%% Figure 3: true model vs curve3 fit (with CI) + inset (with CI)
fig = figure(3);
fig.Position = [left, bottom, width, height];
time = linspace(0,25,1001);

% Main axes
plot(time,model(theta,time),'k');
hold on
hMain3 = plot(curve3,'r--');
plot(t,fitData,'rx')
xg_main3 = get(hMain3(1),'XData'); xg_main3 = xg_main3(:);
yci_main3 = predint(curve3, xg_main3, ciLevel, 'functional');
ciPatch(gca, xg_main3, yci_main3(:,1), yci_main3(:,2), get(hMain3(1),'Color'), fillAlpha);
hold off

legend('True Model', 'Fit Model','Data','Location','northwest')
xlabel('Time (days)')
ylabel('Parasite Numbers')

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth',2);
set(findobj(gcf,'type','line'),'LineWidth',3,'MarkerSize',10);

% --- Inset (with CI) ---
insetWidth  = 0.35;
insetHeight = 0.35;
insetLeft   = 0.20;
insetBottom = 0.20;

insetAxes = axes('Position', [insetLeft insetBottom insetWidth insetHeight]);

zoomX = [0 5];   % time range
zoomY = [0 200]; % parasite number range

zoomIdx = time >= zoomX(1) & time <= zoomX(2) & ...
          model(theta,time) >= zoomY(1) & model(theta,time) <= zoomY(2);

% Inset
plot(time(zoomIdx), model(theta,time(zoomIdx)), 'k', 'LineWidth', 2); hold on
plot(time(zoomIdx), feval(curve3,time(zoomIdx)), 'r--', 'LineWidth', 2);
yci_inset3 = predint(curve3, time(zoomIdx)', ciLevel, 'functional');
ciPatch(insetAxes, time(zoomIdx)', yci_inset3(:,1), yci_inset3(:,2), [1 0 0], fillAlpha);
maskPoints = t >= zoomX(1) & t <= zoomX(2) & ...
             fitData >= zoomY(1) & fitData <= zoomY(2);
plot(t(maskPoints), fitData(maskPoints), 'rx', 'LineWidth', 1.5);
hold off

% Style inset
set(insetAxes, 'FontSize', 8, 'Box', 'on');

set(findobj(gcf,'type','axes'),'FontSize',18,'LineWidth',2);
set(findobj(gcf,'type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(gcf,'type','axes'),'FontSize',fs,'LabelFontSizeMultiplier',1.8);

%% ---------- Local helper (robust CI patch) ----------
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
