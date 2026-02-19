%% ------------
%% MAIN SCRIPT 
%% ------------

clc; clear; close all;

% Set paths
dataPath = 'D:\Documents\MATLAB\dFBA';
outDir = fullfile(dataPath, 'RP_Simulation');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Load model
modelFile = fullfile(dataPath, 'yeast-GEM.mat');
if exist(modelFile, 'file')
    tmp = load(modelFile);
    flds = fieldnames(tmp);
    model = tmp.(flds{1});
else
    error('Model file not found');
end

% Set solver
changeCobraSolver('glpk', 'LP', 0);

% Biomass reaction
biomassRxnID = 'r_2111';

% Conditions
conditions = {
    'Anaerobic_S_P', ...
    'Anaerobic_R_P', ...
    'Aerobic_S_P', ...
    'Aerobic_R_P'};

% Reaction mapping
fluxMap.qs = 'r_1714';     % Glucose
fluxMap.qco2 = 'r_1672';   % CO2
fluxMap.qo2 = 'r_1992';    % O2
fluxMap.qace = 'r_1634';   % Acetate
fluxMap.qcit = 'r_1687';   % Citrate
fluxMap.qlac = 'r_1546';   % Lactate
fluxMap.qpyr = 'r_2033';   % Pyruvate
fluxMap.qeth = 'r_1761';   % Ethanol
fluxMap.qgly = 'r_1808';   % Glycerol

% Pyruvate metabolism reactions
fluxMap.pyk = 'r_0962';    % Pyruvate kinase (PEP → pyruvate) - CDC19
fluxMap.pdc = 'r_0959';    % Pyruvate decarboxylase (main PDC)
fluxMap.pdh = 'r_0961';    % Pyruvate dehydrogenase (mitochondrial)

% Run dFBA for each condition
results = struct();
for c = 1:numel(conditions)
    condName = conditions{c};
    fprintf('\n=== %s ===\n', condName);
    
    xlsFile = fullfile(dataPath, [condName '.xlsx']);
    if ~exist(xlsFile, 'file')
        warning('File %s not found', xlsFile);
        continue;
    end
    
    % Read data
    T = readtable(xlsFile, 'Sheet', 'mmol');
    
    % Run dFBA
    res = run_dFBA(model, T, fluxMap, biomassRxnID, condName);
    results.(condName) = res;
end

%% FIGURE 1: MAIN COMPARISON SP VS RP WITH ALL BYPRODUCTS (NO LEGENDS)
fprintf('\n=== Generating Figure 1: Main Comparison SP vs RP ===\n');
createMainComparisonPlot(results, fluxMap, biomassRxnID, outDir);

%% FIGURE 2: PYRUVATE ROUTING ANALYSIS WITH CDC19 CONSTRAINTS (NO LEGENDS)
fprintf('\n=== Generating Figure 2: Pyruvate Routing Analysis ===\n');
analyzePyruvateRouting(results, fluxMap, outDir);

%% FIGURE 3: LEGEND FOR MAIN COMPARISON PLOT
fprintf('\n=== Generating Legend for Main Comparison ===\n');
createMainLegend(outDir);

%% FIGURE 4: LEGEND FOR PYRUVATE ROUTING ANALYSIS
fprintf('\n=== Generating Legend for Pyruvate Routing ===\n');
createPyruvateRoutingLegend(outDir);

fprintf('\n=== DONE ===\n');

%% --------------------------------
%% dFBA FUNCTION - WITH CONSTRAINTS
%% --------------------------------
function res = run_dFBA(model0, T, fluxMap, biomassRxnID, condName)

    time_s = T.Time;
    nT = height(T);
    nRxn = numel(model0.rxns);
    
    fluxMatrix = zeros(nT, nRxn);
    mu_predicted = zeros(nT, 1);
    status = cell(nT, 1);
    
    % Detect anaerobic/aerobic from condition name
    isAnaerobic = contains(lower(condName), 'anaerobic');
    
    % Feed phase detection (high glucose uptake > 5 mmol/gCDW/h)
    feed_phase = false(nT, 1);
    if ismember('qs', T.Properties.VariableNames)
        feed_phase = abs(T.qs) > 5;
    end
    
    % Loop through time points
    for k = 1:nT
        model = model0;
        
        %% 1. MAXIMISE BIOMASS
        model.c(:) = 0;
        biomassIdx = findRxnIDs(model, biomassRxnID);
        model.c(biomassIdx) = 1;
        
        %% 2. ATP MAINTENANCE
        atpm_idx = findRxnIDs(model, 'r_4046');
        if atpm_idx > 0
            model.lb(atpm_idx) = 1.5;
            model.ub(atpm_idx) = 2.5;
        end
        
        %% 3. GROWTH RATE LIMITS
        if isAnaerobic
            model.ub(biomassIdx) = 0.25;
            model.lb(biomassIdx) = 0;
            
            % Allow anaerobic supplements
            erg_idx = findRxnIDs(model, 'r_1897');
            if erg_idx > 0
                model.lb(erg_idx) = -0.005;
                model.ub(erg_idx) = 0;
            end
            
            oleate_idx = findRxnIDs(model, 'r_2134');
            if oleate_idx > 0
                model.lb(oleate_idx) = -0.002;
                model.ub(oleate_idx) = 0;
            end
        else
            model.ub(biomassIdx) = 0.45;
            model.lb(biomassIdx) = 0;
        end
        
        if feed_phase(k)
            if isAnaerobic
                model.ub(biomassIdx) = 0.25;
            else
                model.ub(biomassIdx) = 0.45;
            end
        else
            model.ub(biomassIdx) = 0.05;
        end
        
        %% 4. CONSTRAIN GLUCOSE UPTAKE
        if ismember('qs', T.Properties.VariableNames)
            qs = T.qs(k);
            idx = findRxnIDs(model, fluxMap.qs);
            if idx > 0
                model.lb(idx) = qs;
                model.ub(idx) = 0;
            end
        end
        
        %% 5. BLOCK O2 FOR ANAEROBIC
        if isAnaerobic && isfield(fluxMap, 'qo2')
            o2_idx = findRxnIDs(model, fluxMap.qo2);
            if o2_idx > 0
                model.lb(o2_idx) = 0;
                model.ub(o2_idx) = 0;
            end
        end
        
        %% 6. PYRUVATE KINASE CONSTRAINT
        if isfield(fluxMap, 'pyk')
            pyk_idx = findRxnIDs(model, fluxMap.pyk);
            if pyk_idx > 0 && ismember('qs', T.Properties.VariableNames)
                q_glc_abs = abs(T.qs(k));
                
                if feed_phase(k)
                    if isAnaerobic
                        % Anaerobic feed: 50% constraint
                        normal_pyk = q_glc_abs * 2.0;
                        model.ub(pyk_idx) = normal_pyk * 0.5;
                    else
                        % Aerobic feed: 80% constraint
                        normal_pyk = q_glc_abs * 2.0;
                        model.ub(pyk_idx) = normal_pyk * 0.8;
                    end
                else
                    % Post-feed: unconstrained
                    model.ub(pyk_idx) = 1000;
                end
                model.lb(pyk_idx) = 0;
            end
        end
        
        %% 7. Solve
        sol = optimizeCbModel(model, 'max');
        
        if sol.stat == 1
            mu_predicted(k) = sol.v(biomassIdx);
            fluxMatrix(k,:) = sol.v(:)';
            status{k} = 'optimal';
        else
            status{k} = 'infeasible';
        end
    end
    
    %% Package results
    res = struct();
    res.condition = condName;
    res.time_s = time_s;
    res.mu_predicted = mu_predicted;
    res.fluxMatrix = fluxMatrix;
    res.rxnIDs = model0.rxns;
    res.rxnNames = model0.rxnNames;
    res.status = status;
    res.feed_phase = feed_phase;
    res.isAnaerobic = isAnaerobic;
    
    % Store experimental data
    res.measured = struct();
    res.measured.time = time_s;
    res.measured.qs = T.qs;
    
    if ismember('mu', T.Properties.VariableNames)
        res.measured.mu = T.mu;
    end
    if ismember('qco2', T.Properties.VariableNames)
        res.measured.qco2 = T.qco2;
    end
    if ismember('qo2', T.Properties.VariableNames)
        res.measured.qo2 = T.qo2;
    end
    
    products = {'qeth','qgly','qace','qcit','qlac','qpyr'};
    for i = 1:length(products)
        if ismember(products{i}, T.Properties.VariableNames)
            res.measured.(products{i}) = T.(products{i});
        end
    end
end

%% -------------------------------------------------------------
%% FIGURE 1: MAIN COMPARISON PLOT - SP VS RP WITH ALL BYPRODUCTS
%% -------------------------------------------------------------
function createMainComparisonPlot(results, fluxMap, biomassRxnID, outDir)
    % Creates one figure with all conditions together
    % 12 panels showing all parameters with SP (solid) and RP (broken)
    % NO LEGENDS - legends are in separate figure
    
    % Define colors
    colors = struct();
    colors.Aerobic_S_P = [0.6294, 0.8078, 0.9804]; % Light blue
    colors.Aerobic_R_P = [0, 0, 0.6451];            % Dark blue
    colors.Anaerobic_S_P = [1.0000, 0.6471, 0];     % Light orange
    colors.Anaerobic_R_P = [0.8000, 0.5000, 0];     % Dark orange
    
    % Order of conditions for plotting
    cond_order = {'Aerobic_S_P', 'Aerobic_R_P', 'Anaerobic_S_P', 'Anaerobic_R_P'};
    
    % Parameters to plot (name, rxnID, ylabel, invert)
    params = {
        'Glucose Uptake', fluxMap.qs, 'Glucose uptake [mmol/gDW/h]', true;
        'Growth Rate', biomassRxnID, 'Growth rate [h^{-1}]', false;
        'CO_2 Production', fluxMap.qco2, 'CO_2 [mmol/gDW/h]', false;
        'O_2 Uptake', fluxMap.qo2, 'O_2 uptake [mmol/gDW/h]', true;
        'Ethanol', fluxMap.qeth, 'Ethanol [mmol/gDW/h]', false;
        'Glycerol', fluxMap.qgly, 'Glycerol [mmol/gDW/h]', false;
        'Acetate', fluxMap.qace, 'Acetate [mmol/gDW/h]', false;
        'Citrate', fluxMap.qcit, 'Citrate [mmol/gDW/h]', false;
        'Lactate', fluxMap.qlac, 'Lactate [mmol/gDW/h]', false;
        'Pyruvate', fluxMap.qpyr, 'Pyruvate [mmol/gDW/h]', false;
        'Pyruvate Kinase (PYK)', fluxMap.pyk, 'PYK flux [mmol/gDW/h]', false;
        'PDC & PDH Routing', '', 'PDC/PDH flux [mmol/gDW/h]', false;
    };
    
    % Set font
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultTextFontName', 'Helvetica');
    
    % Create figure
    fig = figure('Position', [50, 50, 1800, 1200], 'Color', 'white');
    
    % Get time points and feed phase from first condition
    first_cond = cond_order{1};
    if ~isfield(results, first_cond)
        first_cond = fieldnames(results);
        first_cond = first_cond{1};
    end
    time_s = results.(first_cond).time_s;
    feed_phase = results.(first_cond).feed_phase;
    
    % Plot each parameter
    for p = 1:12
        subplot(3, 4, p);
        hold on;
        box on;
        
        param_name = params{p, 1};
        rxn_id = params{p, 2};
        ylabel_str = params{p, 3};
        invert_plot = params{p, 4};
        
        % Special handling for last panel (PDC & PDH)
        if p == 12
            % Plot PDC and PDH for all conditions
            for c = 1:length(cond_order)
                condName = cond_order{c};
                if ~isfield(results, condName)
                    continue;
                end
                res = results.(condName);
                
                % Determine line style based on condition
                if contains(condName, '_S_')
                    line_style = '-';  % SP - solid
                else
                    line_style = '--'; % RP - broken
                end
                
                % Plot PDC
                pdc_idx = findRxnIDs(res.rxnIDs, fluxMap.pdc);
                if pdc_idx > 0
                    pdc_flux = res.fluxMatrix(:, pdc_idx);
                    plot(time_s, pdc_flux, 'Color', colors.(condName), ...
                        'LineStyle', line_style, 'LineWidth', 1.5);
                end
                
                % Plot PDH (dotted line to distinguish from PDC)
                pdh_idx = findRxnIDs(res.rxnIDs, fluxMap.pdh);
                if pdh_idx > 0
                    pdh_flux = res.fluxMatrix(:, pdh_idx);
                    plot(time_s, pdh_flux, 'Color', colors.(condName), ...
                        'LineStyle', ':', 'LineWidth', 1.5);
                end
            end
            
        else
            % Regular parameter - plot all conditions
            for c = 1:length(cond_order)
                condName = cond_order{c};
                if ~isfield(results, condName)
                    continue;
                end
                res = results.(condName);
                
                idx = findRxnIDs(res.rxnIDs, rxn_id);
                if idx > 0
                    flux = res.fluxMatrix(:, idx);
                    if invert_plot
                        flux = -flux;
                    end
                    
                    % Determine line style
                    if contains(condName, '_S_')
                        line_style = '-';  % SP - solid
                    else
                        line_style = '--'; % RP - broken
                    end
                    
                    plot(time_s, flux, 'Color', colors.(condName), ...
                        'LineStyle', line_style, 'LineWidth', 2);
                end
            end
        end
        
        % Set x-axis limits
        xlim([0 150]);
        
        % Calculate y-limits with padding
        children = get(gca, 'Children');
        y_min = inf;
        y_max = -inf;
        for i = 1:length(children)
            if strcmp(get(children(i), 'Type'), 'line')
                ydata = get(children(i), 'YData');
                y_min = min(y_min, min(ydata));
                y_max = max(y_max, max(ydata));
            end
        end
        
        if isfinite(y_min) && isfinite(y_max)
            y_padding = 0.1 * (y_max - y_min);
            if y_padding < 0.01
                y_padding = 0.1;
            end
            ylim([y_min - y_padding, y_max + y_padding]);
        end
        
        % Highlight feed phase
        if any(feed_phase)
            yl = ylim;
            feed_end = max(time_s(feed_phase));
            fill([0 feed_end feed_end 0], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
            line([feed_end feed_end], yl, 'Color', [0.5 0.5 0.5], ...
                'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
        end
        
        xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel(ylabel_str, 'FontSize', 10, 'FontWeight', 'bold');
        title(param_name, 'FontSize', 15, 'FontWeight', 'bold');
        
        % NO LEGEND - legends are in separate figure
        
        grid on;
    end
    
    % Add overall title
    sgtitle('Model Predictions: SP vs RP Comparison with Pyruvate Routing Analysis', ...
        'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica');
    
    % Save figure
    pngFile = fullfile(outDir, 'Main_Comparison_SP_RP.png');
    exportgraphics(fig, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved Main Comparison Plot: %s\n', pngFile);
    
    %% Reset default font
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultTextFontName', 'Helvetica');
end

%% ----------------------------------------------------------
%% FIGURE 2: PYRUVATE ROUTING ANALYSIS WITH CDC19 CONSTRAINTS
%% ----------------------------------------------------------
function analyzePyruvateRouting(results, fluxMap, outDir)
    % Creates a 2x2 figure showing pyruvate routing for each condition
    % NO LEGENDS - legends are in separate figure
        
    % Define colors for each condition
    colors = struct();
    colors.Aerobic_S_P = [0.6294, 0.8078, 0.9804]; % Light blue
    colors.Aerobic_R_P = [0, 0, 0.6451];            % Dark blue
    colors.Anaerobic_S_P = [1.0000, 0.6471, 0];     % Light orange
    colors.Anaerobic_R_P = [0.8000, 0.5000, 0];     % Dark orange
    
    % Order for subplots
    subplot_order = {'Aerobic_S_P', 'Aerobic_R_P', 'Anaerobic_S_P', 'Anaerobic_R_P'};
    
    % Set font
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultTextFontName', 'Helvetica');
    
    fig = figure('Position', [50, 50, 1400, 1000], 'Color', 'white');
    
    for c = 1:4
        condName = subplot_order{c};
        
        if ~isfield(results, condName)
            subplot(2, 2, c);
            text(0.5, 0.5, sprintf('%s\nData not available', condName), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, ...
                'FontName', 'Helvetica', 'FontWeight', 'bold');
            axis off;
            continue;
        end
        
        res = results.(condName);
        
        % Get indices for pyruvate-related reactions
        pyk_idx = findRxnIDs(res.rxnIDs, fluxMap.pyk);
        pdc_idx = findRxnIDs(res.rxnIDs, fluxMap.pdc);
        pdh_idx = findRxnIDs(res.rxnIDs, fluxMap.pdh);
        
        time_s = res.time_s;
        fluxMatrix = res.fluxMatrix;
        feed_phase = res.feed_phase;
        
        % Get fluxes
        pyk_flux = fluxMatrix(:, pyk_idx);
        pdc_flux = zeros(length(time_s), 1);
        pdh_flux = zeros(length(time_s), 1);
        
        if pdc_idx > 0
            pdc_flux = fluxMatrix(:, pdc_idx);
        end
        
        if pdh_idx > 0
            pdh_flux = fluxMatrix(:, pdh_idx);
        end
        
        % Plot for this condition
        subplot(2, 2, c);
        hold on;
        box on;
        
        % Plot fluxes with condition-specific colors
        plot(time_s, pyk_flux, 'Color', colors.(condName), 'LineStyle', '-', ...
            'LineWidth', 2.5);
        plot(time_s, pdc_flux, 'Color', colors.(condName), 'LineStyle', '--', ...
            'LineWidth', 2.5);
        plot(time_s, pdh_flux, 'Color', colors.(condName), 'LineStyle', ':', ...
            'LineWidth', 2.5);
        
        % Set x-axis limits
        xlim([0 150]);
        
        % Calculate y-limits with padding
        all_flux_data = [pyk_flux(:); pdc_flux(:); pdh_flux(:)];
        y_min = min(all_flux_data);
        y_max = max(all_flux_data);
        
        if y_max > y_min
            y_padding = 0.15 * (y_max - y_min);
            y_min = y_min - y_padding;
            y_max = y_max + y_padding;
        else
            y_min = y_min - 1;
            y_max = y_max + 1;
        end
        
        % Highlight feed phase
        if any(feed_phase)
            feed_end = max(time_s(feed_phase));
            fill([0 feed_end feed_end 0], [y_min y_min y_max y_max], ...
                 [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
            line([feed_end feed_end], [y_min y_max], ...
                 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
                 'LineStyle', '--', 'HandleVisibility', 'off');
        end
        
        ylim([y_min y_max]);
        
        xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
        ylabel('Flux [mmol/gDW/h]', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
        
        % Format condition name for title
        title_cond = strrep(condName, '_', ' ');
        title(title_cond, 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Helvetica');
        
        % NO LEGEND - legends are in separate figure
        
        grid on;
        
        % Calculate and display statistics in BOTTOM RIGHT corner
        if any(feed_phase)
            mean_pyk = mean(pyk_flux(feed_phase));
            mean_pdc = mean(pdc_flux(feed_phase));
            mean_pdh = mean(pdh_flux(feed_phase));
            
            if mean_pyk > 0
                pdc_percent = mean_pdc/mean_pyk*100;
                pdh_percent = mean_pdh/mean_pyk*100;
            else
                pdc_percent = 0;
                pdh_percent = 0;
            end
            
            text_str = {
                sprintf('Feed Phase Averages:');
                sprintf('PYK: %.2f', mean_pyk);
                sprintf('PDC: %.2f (%.0f%%)', mean_pdc, pdc_percent);
                sprintf('PDH: %.2f (%.0f%%)', mean_pdh, pdh_percent)
            };
            
            % Position text in top right corner
            text(0.6, 0.6, text_str, 'Units', 'normalized', ...
                 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
                 'FontSize', 10, 'FontName', 'Helvetica', ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                 'FontWeight', 'bold');
        end
    end

    sgtitle('Pyruvate Routing Analysis with CDC19 Constraints', ...
            'FontSize', 18, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    
    % Save plot
    pngFile = fullfile(outDir, 'Pyruvate_Routing_Analysis.png');
    exportgraphics(fig, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved Pyruvate Routing Analysis: %s\n', pngFile);
    
    %% Reset default font
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultTextFontName', 'Helvetica');
end

%% ------------------------------------------------------------------------
%% FIGURE 3: LEGEND FOR MAIN COMPARISON PLOT
%% ------------------------------------------------------------------------
function createMainLegend(outDir)
    % Creates a separate figure with legend for the main comparison plot
    
    figure('Position', [100, 100, 400, 300], 'Color', 'white');
    axis off;
    
    % Create dummy lines for legend
    hold on;
    
    % Aerobic SP
    plot(nan, nan, 'Color', [0.6294, 0.8078, 0.9804], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic SP');
    
    % Aerobic RP
    plot(nan, nan, 'Color', [0, 0, 0.6451], 'LineStyle', '--', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic RP');
    
    % Anaerobic SP
    plot(nan, nan, 'Color', [1.0000, 0.6471, 0], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic SP');
    
    % Anaerobic RP
    plot(nan, nan, 'Color', [0.8000, 0.5000, 0], 'LineStyle', '--', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic RP');
    
    % Feed phase indicator
    fill([0 1 1 0], [0 0 1 1], [0.8 0.8 0.8], 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName', 'Feed Phase');
    
    legend('Location', 'best', 'FontSize', 10, 'FontName', 'Helvetica', ...
        'FontWeight', 'bold', 'Box', 'on');
    
    title('Legend: Main Comparison Plot', 'FontSize', 14, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica');
    
    % Save legend
    pngFile = fullfile(outDir, 'Legend_Main_Comparison.png');
    exportgraphics(gcf, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved Main Legend: %s\n', pngFile);
end

%% ------------------------------------------------------------------------
%% FIGURE 4: LEGEND FOR PYRUVATE ROUTING ANALYSIS
%% ------------------------------------------------------------------------
function createPyruvateRoutingLegend(outDir)
    % Creates a separate figure with legend for pyruvate routing analysis
    
    figure('Position', [100, 100, 400, 350], 'Color', 'white');
    axis off;
    
    % Create dummy lines for legend
    hold on;
    
    % Line styles
    plot(nan, nan, 'Color', [0, 0, 0], 'LineStyle', '-', ...
        'LineWidth', 2.5, 'DisplayName', 'PYK Flux');
    
    plot(nan, nan, 'Color', [0, 0, 0], 'LineStyle', '--', ...
        'LineWidth', 2.5, 'DisplayName', 'PDC Flux');
    
    plot(nan, nan, 'Color', [0, 0, 0], 'LineStyle', ':', ...
        'LineWidth', 2.5, 'DisplayName', 'PDH Flux');
    
    % Colors
    text(0.1, 0.6, 'Color Coding:', 'FontSize', 11, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica');
    
    % Color samples
    plot(nan, nan, 'Color', [0.6294, 0.8078, 0.9804], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic SP');
    
    plot(nan, nan, 'Color', [0, 0, 0.6451], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic RP');
    
    plot(nan, nan, 'Color', [1.0000, 0.6471, 0], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic SP');
    
    plot(nan, nan, 'Color', [0.8000, 0.5000, 0], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic RP');
    
    % Feed phase
    fill([0 1 1 0], [0 0 1 1], [0.8 0.8 0.8], 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName', 'Feed Phase');
    
    % Statistics box
    rectangle('Position', [0.1, 0.05, 0.8, 0.15], 'EdgeColor', 'black', ...
        'LineWidth', 1);
    text(0.5, 0.12, 'Feed Phase Averages Box', 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    
    legend('Location', 'best', 'FontSize', 9, 'FontName', 'Helvetica', ...
        'FontWeight', 'bold', 'Box', 'on');
    
    title('Legend: Pyruvate Routing Analysis', 'FontSize', 14, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica');
    
    % Save legend
    pngFile = fullfile(outDir, 'Legend_Pyruvate_Routing.png');
    exportgraphics(gcf, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved Pyruvate Routing Legend: %s\n', pngFile);
end

%% Helper function to find reaction index
function idx = findRxnIDs(rxnList, rxnID)
    if iscell(rxnList)
        idx = find(strcmp(rxnList, rxnID));
    else
        idx = find(strcmp(rxnList.rxns, rxnID));
    end
    if isempty(idx)
        idx = 0;
    else
        idx = idx(1);
    end
end