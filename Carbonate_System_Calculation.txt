%% Carbonate_System_Calculation.m
% -------------------------------------------------------------------------
% Calculation of seawater carbonate system parameters from measured pH
% (total scale), temperature, salinity, and total alkalinity (TA) using
% the CO2SYS package for MATLAB (van Heuven et al., 2011; Lewis &
% Wallace, 1998).
%
% Input:
%   input_carbonate_data.csv — comma-separated file with columns:
%     Group           — experimental treatment label (character)
%     Enclosure       — replicate identifier (numeric)
%     pH_total        — pH on the total scale (dimensionless)
%     TA_umol_kg      — total alkalinity (umol kg-1)
%     Temperature_C   — seawater temperature (degrees Celsius)
%     Salinity_PSU    — practical salinity (PSU)
%
% Output:
%   carbonate_system_all_enclosures.csv — per-enclosure results
%   carbonate_system_summary.csv        — group-level mean and SD
%   carbonate_system_fig.png / .pdf     — bar plots of key parameters
%
% Dependency:
%   CO2SYS.m v3 (https://github.com/jamesorr/CO2SYS-MATLAB)
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. Read measured parameters from file
inputFile = 'input_carbonate_data.csv';
T = readtable(inputFile);
n = height(T);

%% 2. CO2SYS configuration
PAR1TYPE  = 1;   % total alkalinity
PAR2TYPE  = 3;   % pH
pHscaleIN = 1;   % total scale
K1K2const = 4;   % Mehrbach refit by Dickson & Millero (1987)
KSO4const = 1;   % Dickson (1990)
KFconst   = 1;   % Dickson & Riley (1979)
TBconst   = 1;   % Uppstrom (1974)

%% 3. Run CO2SYS
[DATA, HEADERS] = CO2SYS( ...
    T.TA_umol_kg, T.pH_total, PAR1TYPE, PAR2TYPE, ...
    T.Salinity_PSU, T.Temperature_C, T.Temperature_C, 0, 0, ...
    zeros(n,1), zeros(n,1), ...
    pHscaleIN, K1K2const, KSO4const, KFconst, TBconst);

col = @(name) find(strcmp(HEADERS, name));

%% 4. Extract carbonate system parameters
T.pCO2_uatm       = DATA(:, col('pCO2in'));
T.DIC_umol_kg      = DATA(:, col('TCO2in'));
T.HCO3_umol_kg     = DATA(:, col('HCO3in'));
T.CO3_umol_kg      = DATA(:, col('CO3in'));
T.CO2_umol_kg      = DATA(:, col('CO2in'));
T.Omega_calcite    = DATA(:, col('OmegaCAin'));
T.Omega_aragonite  = DATA(:, col('OmegaARin'));
T.Revelle_factor   = DATA(:, col('RFin'));

%% 5. Group-level summary statistics
paramCols = {'pH_total','TA_umol_kg','pCO2_uatm','DIC_umol_kg', ...
             'HCO3_umol_kg','CO3_umol_kg','CO2_umol_kg', ...
             'Omega_calcite','Omega_aragonite','Revelle_factor'};

groups = unique(T.Group);
nG     = numel(groups);
nP     = numel(paramCols);

summaryRows = cell(nG, 1 + 2*nP);
for g = 1:nG
    idx = strcmp(T.Group, groups{g});
    summaryRows{g,1} = groups{g};
    for p = 1:nP
        vals = T.(paramCols{p})(idx);
        summaryRows{g, 1+p}    = mean(vals);
        summaryRows{g, 1+nP+p} = std(vals);
    end
end

summaryHeaders = ['Group', strcat(paramCols,'_mean'), strcat(paramCols,'_std')];
summaryTable   = cell2table(summaryRows, 'VariableNames', summaryHeaders);

%% 6. Export results
writetable(T,            'carbonate_system_all_enclosures.csv');
writetable(summaryTable, 'carbonate_system_summary.csv');

%% 7. Visualisation
plotParams = {'pCO2_uatm','DIC_umol_kg','HCO3_umol_kg', ...
              'CO3_umol_kg','Omega_aragonite','Omega_calcite'};
plotLabels = {'{\itp}CO_2 (\muatm)','DIC (\mumol kg^{-1})', ...
              'HCO_3^- (\mumol kg^{-1})','CO_3^{2-} (\mumol kg^{-1})', ...
              '\Omega_{Aragonite}','\Omega_{Calcite}'};

cmap = lines(nG);

figure('Units','centimeters','Position',[3 3 28 16]);
for sp = 1:numel(plotParams)
    subplot(2, 3, sp); hold on;
    m = zeros(nG,1);
    s = zeros(nG,1);
    for g = 1:nG
        idx = strcmp(T.Group, groups{g});
        m(g) = mean(T.(plotParams{sp})(idx));
        s(g) = std(T.(plotParams{sp})(idx));
    end
    b = bar(1:nG, m, 0.5, 'FaceColor','flat','EdgeColor','k','LineWidth',0.8);
    b.CData = cmap(1:nG,:);
    errorbar(1:nG, m, s, 'k.', 'LineWidth',1.2, 'CapSize',6);
    set(gca, 'XTick',1:nG, 'XTickLabel',groups, ...
        'FontSize',9, 'Box','off', 'TickDir','out');
    ylabel(plotLabels{sp}, 'FontSize',10);
    hold off;
end
sgtitle('Carbonate System Parameters', 'FontSize',12, 'FontWeight','bold');

print(gcf, 'carbonate_system_fig', '-dpng', '-r300');
print(gcf, 'carbonate_system_fig', '-dpdf');
