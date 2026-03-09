%% Parasitemia V0 and g recovery sweeps (time-error only)
clear; clc; close all

%% Model definition
model = @(theta,t) theta(2).*theta(1).*exp(theta(4).*t) + theta(3).*exp(-theta(5).*t);

%% Fixed known parameters
f  = 1;
B0 = 0.02;
k  = 0.323;

%% Simulation settings
t_meas         = 0:2:24;               % nominal measurement times
numNoiseLevels = 20;
sigmaVals      = linspace(0,0.5,numNoiseLevels); % measurement noise SD
epsilonVals    = linspace(0,0.5,numNoiseLevels); % time error SD
numMC          = 100;                  % Monte Carlo replicates per grid cell

% Random draw functions
getV0 = @(sz) 0.002 * ones(sz);
getg  = @(sz) 0.531 * ones(sz);

% Storage arrays — V0
biasV0     = nan(numNoiseLevels,numNoiseLevels);
rmseV0     = nan(numNoiseLevels,numNoiseLevels);
varV0est   = nan(numNoiseLevels,numNoiseLevels);
relErrV0   = nan(numNoiseLevels,numNoiseLevels);
ciWidthV0  = nan(numNoiseLevels,numNoiseLevels);

% Storage arrays — g
biasg     = nan(numNoiseLevels,numNoiseLevels);
rmseg     = nan(numNoiseLevels,numNoiseLevels);
vargest   = nan(numNoiseLevels,numNoiseLevels);
relErrg   = nan(numNoiseLevels,numNoiseLevels);
ciWidthg  = nan(numNoiseLevels,numNoiseLevels);

%% Grid loop
for i = 1:numNoiseLevels
    for j = 1:numNoiseLevels
        sigma   = sigmaVals(i);
        epsilon = epsilonVals(j);

        estV0_all  = nan(numMC,1);
        estg_all   = nan(numMC,1);
        trueV0_all = nan(numMC,1);
        trueg_all  = nan(numMC,1);
        ciWidthsV0_all = nan(numMC,1);
        ciWidthsg_all  = nan(numMC,1);

        parfor m = 1:numMC
            % Draw parameters
            V0 = getV0(1);
            g  = getg(1);
            theta = [V0,f,B0,g,k];
            trueV0_all(m) = V0;
            trueg_all(m)  = g;

            % Generate synthetic data
            y_data = normalGeneratorMultiplicative_normalTimeError(theta,t_meas,sigma,model,epsilon,1);

            % Fit V0 and g jointly
            ft = fittype(@(V0_fit,g_fit,t) model([V0_fit,f,B0,g_fit,k],t), ...
                'independent','t','coefficients',{'V0_fit','g_fit'});
            fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0]);

            try
                curve = fit(t_meas', y_data', ft, fo);
                estV0_all(m) = curve.V0_fit;
                estg_all(m)  = curve.g_fit;

                % Get 95% CIs
                ci = confint(curve, 0.95); % 2xN array
                ciWidthsV0_all(m) = ci(2,1) - ci(1,1);
                ciWidthsg_all(m)  = ci(2,2) - ci(1,2);
            catch
                estV0_all(m)     = NaN;
                estg_all(m)      = NaN;
                ciWidthsV0_all(m)= NaN;
                ciWidthsg_all(m) = NaN;
            end
        end

        % --- Compute stats: V0 ---
        validV0 = ~isnan(estV0_all);
        if any(validV0)
            diffsV0         = estV0_all(validV0) - trueV0_all(validV0);
            biasV0(i,j)     = mean(diffsV0);
            rmseV0(i,j)     = sqrt(mean(diffsV0.^2));
            varV0est(i,j)   = var(estV0_all(validV0),1);
            relErrV0(i,j)   = mean(diffsV0 ./ trueV0_all(validV0));
            ciWidthV0(i,j)  = mean(ciWidthsV0_all(validV0));
        end

        % --- Compute stats: g ---
        validg = ~isnan(estg_all);
        if any(validg)
            diffsg         = estg_all(validg) - trueg_all(validg);
            biasg(i,j)     = mean(diffsg);
            rmseg(i,j)     = sqrt(mean(diffsg.^2));
            vargest(i,j)   = var(estg_all(validg),1);
            relErrg(i,j)   = mean(diffsg ./ trueg_all(validg));
            ciWidthg(i,j)  = mean(ciWidthsg_all(validg));
        end
    end
end

save("parasitemia_sweepsData_V0_g.mat");

%% --- Plotting ---
load("parasitemia_sweepsData_V0_g.mat");
width = 900; height = 700; left = 50; bottom = 80;
fs = 26; lw = 2;

% === V0 plots ===
maxAbsBias = max(abs(biasV0(:)), [], 'omitnan');
maxAbsRelerr = max(abs(relErrV0(:)), [], 'omitnan');

plots = { ...
    {biasV0,    'Bias in V_0 estimate',          'Figures/V0_bias_heatmap.png', true, maxAbsBias}, ...
    {rmseV0,    'RMSE of V_0 estimate',          'Figures/V0_rmse_heatmap.png', false, []}, ...
    {varV0est,  'Variance of V_0 estimates',     'Figures/V0_variance_heatmap.png', false, []}, ...
    {relErrV0,  'Mean Relative Error in V_0 estimates', 'Figures/V0_relError_heatmap.png', true, maxAbsRelerr}, ...
    {ciWidthV0, 'Mean 95% CI width for V_0',     'Figures/V0_CIwidth_heatmap.png', false, []} ...
};

for p = 1:numel(plots)
    fig = figure; fig.Position = [left,bottom,width,height];
    imagesc(epsilonVals,sigmaVals,plots{p}{1});
    set(gca,'YDir','normal','FontSize',fs,'LineWidth',lw);
    if plots{p}{4}
        caxis([-plots{p}{5}, plots{p}{5}]);
        colormap(crameri('-vik',256));
    else
        colormap(crameri('batlow',256));
    end
    colorbar;
    xlabel('\sigma_U','FontSize',fs);
    ylabel('\sigma_\epsilon','FontSize',fs);
    %title(plots{p}{2},'FontSize',fs);
    set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);
    exportgraphics(fig,plots{p}{3});
end

% === g plots ===
maxAbsBiasg = max(abs(biasg(:)), [], 'omitnan');
maxAbsRelerrg = max(abs(relErrg(:)), [], 'omitnan');

plots_g = { ...
    {biasg,    'Bias in g estimate',          'Figures/g_bias_heatmap.png', true, maxAbsBiasg}, ...
    {rmseg,    'RMSE of g estimate',          'Figures/g_rmse_heatmap.png', false, []}, ...
    {vargest,  'Variance of g estimates',     'Figures/g_variance_heatmap.png', false, []}, ...
    {relErrg,  'Mean Relative Error in g estimates', 'Figures/g_relError_heatmap.png', true, maxAbsRelerrg}, ...
    {ciWidthg, 'Mean 95% CI width for g',     'Figures/g_CIwidth_heatmap.png', false, []} ...
};

for p = 1:numel(plots_g)
    fig = figure; fig.Position = [left,bottom,width,height];
    imagesc(epsilonVals,sigmaVals,plots_g{p}{1});
    set(gca,'YDir','normal','FontSize',fs,'LineWidth',lw);
    if plots_g{p}{4}
        caxis([-plots_g{p}{5}, plots_g{p}{5}]);
        colormap(crameri('-vik',256));
    else
        colormap(crameri('batlow',256));
    end
    colorbar;
    xlabel('\sigma_U','FontSize',fs);
    ylabel('\sigma_\epsilon','FontSize',fs);
    %title(plots_g{p}{2},'FontSize',fs);
    set(findobj(gcf,'type','axes'),'FontSize',fs,'LineWidth',2,'LabelFontSizeMultiplier',1.4);
    exportgraphics(fig,plots_g{p}{3});
end

disp('Sweep complete for V0')