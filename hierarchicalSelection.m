function [levels, centrality_dynamic, coverage_history] = hierarchicalSelection()  
    %% Main Function: Hierarchical Monitoring Network Construction  
    % Input: Coordinates, VOI, Coverage Area, Initial Centrality from data.xlsx  
    % Output: Level assignment matrix, Dynamic centrality, Coverage history  
      
    %% Data Reading and Preprocessing  
    try  
        data = readtable('data.xlsx');  
        x = data.x;  
        y = data.y;  
        VOI = data.VOI;  
        centrality_init = data.centrality_init;  
        coverage_areas = data.coverage_areas;  
    catch ME  
        error('Data reading failed: %s\nPlease ensure:\n1. data.xlsx exists in current directory\n2. Column names include: x, y, VOI, centrality_init, coverage_areas', ME.message);  
    end  
   
    %% Core Parameter Settings (Adjustable as needed)  
    params.alpha = 1;        % VOI weight  
    params.beta = 0;         % Centrality weight  
    params.dist_threshold = 7000; % Spatial adjacency threshold (meters)  
    params.converge_thres = 0.003; % Incremental convergence threshold (0.3%)  
    params.max_levels = 10;     % Maximum allowed levels  
    params.min_level_points = 1;% Minimum points per level  
    params.gain_threshold = 0.01; % Gain threshold (5%)  
    params.redundancy_threshold = 0.1; % Information redundancy threshold  
    params.max_nodes_for_full_calc = 500; % Maximum nodes for full centrality calculation  
   
    %% Initialization  
    n = length(x);  
    [VOI_std, centrality_std] = standardizeFeatures(VOI, centrality_init);  
    centrality_dynamic = centrality_init;  
    covered = false(1, n);  
    total_coverage = 0;  
    coverage_history = [];  
    level_metrics = struct('level', {}, 'voi_sum', {}, 'cent_sum', {}, 'coverage', {}, 'point_count', {}, 'redundancy', {});  
   
    %% 1. Enhanced Starting Point Selection (VOI × Centrality Product)  
    start_score = VOI_std .* centrality_std;  
    [~, start_idx] = max(start_score);  
    levels = -ones(n, 1);  
    levels(start_idx) = 0;  
    covered(start_idx) = true;  
    total_coverage = coverage_areas(start_idx);  
    coverage_history(end+1) = total_coverage;  
      
    % Record starting point metrics  
    level_metrics(1).level = 0;  
    level_metrics(1).voi_sum = VOI(start_idx);  
    level_metrics(1).cent_sum = centrality_init(start_idx);  
    level_metrics(1).coverage = coverage_areas(start_idx);  
    level_metrics(1).point_count = 1;  
    level_metrics(1).redundancy = 0; % Level 0 has no redundancy  
      
    prev_voi = level_metrics(1).voi_sum;  
    prev_cent = level_metrics(1).cent_sum;  
      
    %% Main Loop Improvement
    current_level = 1;
    remaining = setdiff(1:n, start_idx);
    last_valid_cent = centrality_dynamic(start_idx);
    converge_counter = 0; % Add convergence counter

    while ~isempty(remaining) && current_level <= params.max_levels
        % Dynamic weight adjustment (enhance centrality influence)
        gamma = params.alpha + (1-params.alpha)*(1 - 1/(current_level+1));
        
        % Generate candidate points (with spatial distribution constraints)
        [candidates, gains] = generateCandidates(remaining, levels, x, y, VOI_std, ...
                                               centrality_dynamic, coverage_areas, gamma, params);
        
        % Call the unified selectPoints function
        [selected, remaining] = selectPoints(candidates, gains, remaining, params);
        
        if isempty(selected)
            if current_level > 2
                % Try to supplement nodes from historical levels
                prev_level = current_level-1;
                candidates = find(levels == prev_level & ~covered);
                if ~isempty(candidates)
                    selected = candidates(1:min(2,length(candidates)));
                    remaining = setdiff(remaining, selected);
                end
            end
            if isempty(selected)
                break;
            end
        end
 
        % Update system state
        levels(selected) = current_level;
        covered(selected) = true;
        total_coverage = total_coverage + sum(coverage_areas(selected));
        coverage_history(end+1) = total_coverage;
          
        % Real-time centrality update  
        [centrality_dynamic, ~] = updateCentrality(levels, x, y, params);  
          
        % === Fix Point: Add current level metrics recording ===
        level_metrics(current_level+1).level = current_level;
        level_metrics(current_level+1).voi_sum = sum(VOI(selected));
        level_metrics(current_level+1).cent_sum = sum(centrality_dynamic(selected));
        level_metrics(current_level+1).coverage = sum(coverage_areas(selected));
        level_metrics(current_level+1).point_count = numel(selected);
          
        % Calculate information redundancy with previous level  
        if current_level > 0  
            prev_points = find(levels == current_level-1);  
            curr_points = selected;  
            redundancy = computeLevelRedundancy(prev_points, curr_points, VOI, centrality_dynamic);  
            level_metrics(current_level+1).redundancy = redundancy;  
            fprintf('Information redundancy between Level %d and Level %d: %.4f\n', current_level-1, current_level, redundancy);  
              
            % Redundancy check  
            if redundancy > params.redundancy_threshold  
                fprintf('Stopped at Level %d due to information redundancy %.4f exceeding threshold %.4f\n', current_level, redundancy, params.redundancy_threshold);  
                break;  
            end  
        else  
            level_metrics(current_level+1).redundancy = 0;  
        end  
          
        % Convergence detection  
        current_voi = level_metrics(current_level+1).voi_sum;  
        current_cent = level_metrics(current_level+1).cent_sum;  
        voi_delta = abs((current_voi - prev_voi)/max(prev_voi, 1e-8));  
        cent_delta = abs((current_cent - prev_cent)/max(prev_cent, 1e-8));  
          
        fprintf('Level %d: Selected %d points, VOI increment=%.2f%%, Centrality increment=%.2f%%\n', ...  
                current_level, numel(selected), voi_delta*100, cent_delta*100);  
          
        if voi_delta < params.converge_thres && cent_delta < params.converge_thres  
            converge_counter = converge_counter + 1;  
            if converge_counter > 3  
                fprintf('Stopped at Level %d due to incremental convergence\n', current_level);  
                break;  
            end  
        else  
            converge_counter = 0;  
        end  
          
        prev_voi = current_voi;  
        prev_cent = current_cent;  
        current_level = current_level + 1;  
    end  
   
    %% Post-processing Optimization  
    [levels, centrality_dynamic] = optimizeRedundancy(levels, VOI, centrality_dynamic, params);  
      
    %% Calculate extraction efficiency and contribution for each point (using level weights)  
    [efficiency, contribution] = calculatePointMetrics(levels, VOI, centrality_dynamic);  
      
    %% Calculate information redundancy between all levels  
    level_redundancy = computeAllLevelRedundancy(levels, VOI, centrality_dynamic);  
    
    %% ================= New: Calculate Information Transfer Efficiency ================= %%
    [transfer_efficiency, flow_intensity] = computeTransferEfficiency(levels, level_metrics);
    
    %% Visualization
    visualizeResults(x, y, levels, coverage_areas, level_metrics, efficiency, contribution, ...
                     level_redundancy, transfer_efficiency, flow_intensity);
    
    %% Save Results
    saveResults(levels, coverage_areas, VOI, centrality_dynamic, efficiency, contribution, ...
                level_metrics, level_redundancy, transfer_efficiency, flow_intensity);
      
    % Return results  
    coverage_history = coverage_history';  
end  

%% ================= Information Transfer Efficiency Calculation Function ================= %%
function [transfer_efficiency, flow_intensity] = computeTransferEfficiency(levels, level_metrics)
    % Calculate information transfer efficiency between levels
    % Input: Level assignment, Level metrics
    % Output: Transfer efficiency matrix, Flow intensity matrix
    
    % Get valid levels
    valid_levels = [level_metrics.level];
    num_levels = length(valid_levels);
    
    % Initialize output matrices
    transfer_efficiency = zeros(num_levels, num_levels);
    flow_intensity = zeros(num_levels, num_levels);
    
    % Calculate average VOI and average centrality for each level
    avg_voi = zeros(1, num_levels);
    avg_cent = zeros(1, num_levels);
    
    for i = 1:num_levels
        if level_metrics(i).point_count > 0
            avg_voi(i) = level_metrics(i).voi_sum / level_metrics(i).point_count;
            avg_cent(i) = level_metrics(i).cent_sum / level_metrics(i).point_count;
        end
    end
    
    % Calculate information transfer efficiency between levels
    for from_level = 1:num_levels
        for to_level = (from_level+1):num_levels
            % VOI transfer rate = Target level average VOI / Source level average VOI
            if avg_voi(from_level) > 0
                voi_transfer = avg_voi(to_level) / avg_voi(from_level);
            else
                voi_transfer = 0;
            end
            
            % Flow intensity = Target level average centrality / Source level average centrality
            if avg_cent(from_level) > 0
                flow_strength = avg_cent(to_level) / avg_cent(from_level);
            else
                flow_strength = 0;
            end
            
            % Comprehensive transfer efficiency (geometric mean)
            transfer_efficiency(from_level, to_level) = sqrt(voi_transfer * flow_strength);
            flow_intensity(from_level, to_level) = flow_strength;
        end
    end
end

%% ================= Enhanced Centrality Calculation ================= %%  
function [new_cent, adj_matrix] = updateCentrality(levels, x, y, params)
    % Always calculate weighted average of three centrality measures
    selected = find(levels ~= -1);
    if isempty(selected)
        new_cent = ones(size(levels));
        adj_matrix = [];
        return;
    end
 
    % Build adjacency matrix
    adj_matrix = buildAdjacencyMatrix(x(selected), y(selected), params.dist_threshold);
    n = size(adj_matrix,1);
 
    % 1. Degree Centrality (Corrected version)
    deg_cent = sum(adj_matrix, 2);
    deg_cent = deg_cent / max([deg_cent; 1]); % Safe division
 
    % 2. Betweenness Centrality (Brandes Algorithm)
    bet_cent = computeBetweenness(adj_matrix);
    bet_cent = bet_cent / max([bet_cent; 1]);
 
    % 3. Closeness Centrality (Floyd Algorithm)
    clos_cent = computeCloseness(adj_matrix);
    clos_cent = clos_cent / max([clos_cent; 1]);
 
    % Combined Centrality (Equal weight weighted average)
    combined_cent = (deg_cent + bet_cent + clos_cent) / 3;
 
    % Update to entire network
    new_cent = zeros(size(levels));
    new_cent(selected) = combined_cent;
end
 
function bet_cent = computeBetweenness(adj)
    % Brandes Algorithm Implementation
    n = size(adj,1);
    bet_cent = zeros(n,1);
    
    for s = 1:n
        dist = -ones(1,n);
        sigma = zeros(1,n);
        pred = cell(1,n);
        delta = zeros(1,n);
        
        dist(s) = 0;
        sigma(s) = 1;
        Q = s;
        
        while ~isempty(Q)
            v = Q(1);
            Q = Q(2:end);
            neighbors = find(adj(v,:));
            
            for w = neighbors
                if dist(w) < 0
                    dist(w) = dist(v) + 1;
                    Q = [Q, w];
                end
                if dist(w) == dist(v) + 1
                    sigma(w) = sigma(w) + sigma(v);
                    pred{w} = [pred{w}, v];
                end
            end
        end
        
        W = find(dist > 0);
        [~, idx] = sort(dist(W), 'descend');
        sorted_w = W(idx);
        
        for w = sorted_w
            for v = pred{w}
                delta(v) = delta(v) + (sigma(v)/sigma(w)) * (1 + delta(w));
            end
            if w ~= s
                bet_cent(w) = bet_cent(w) + delta(w);
            end
        end
    end
    
    % Fix single node case
    if n == 1
        bet_cent = ones(n,1);
    else
        bet_cent = bet_cent / ((n-1)*(n-2)); % Normalization
    end
end
 
function clos_cent = computeCloseness(adj)
    % Floyd-Warshall Algorithm Implementation
    n = size(adj,1);
    dist = inf(n);
    
    for i = 1:n
        dist(i,i) = 0;
        neighbors = find(adj(i,:));
        dist(i,neighbors) = 1;
    end
    
    for k = 1:n
        for i = 1:n
            for j = 1:n
                if dist(i,k) + dist(k,j) < dist(i,j)
                    dist(i,j) = dist(i,k) + dist(k,j);
                end
            end
        end
    end
    
    clos_cent = zeros(n,1);
    for i = 1:n
        valid_dists = dist(i, dist(i,:) < inf & dist(i,:) > 0);
        if ~isempty(valid_dists)
            clos_cent(i) = (n-1) / sum(valid_dists);
        else
            clos_cent(i) = 0;
        end
    end
    
    % Fix single node case
    if n == 1
        clos_cent = ones(n,1);
    end
end
  
%% ================= Level-based Contribution Calculation ================= %%  
function [efficiency, contribution] = calculatePointMetrics(levels, VOI, centrality)  
    % Calculate extraction efficiency and contribution for each monitoring point (using level weights)  
    n = numel(levels);  
    efficiency = zeros(n, 1);  
    contribution = zeros(n, 1);  
      
    % Get all valid levels  
    valid_levels = levels(levels >= 0);  
    unique_levels = unique(valid_levels);  
      
    % Calculate level weights (Level VOI sum / All levels VOI sum)  
    total_voi = sum(VOI(levels >= 0));  
    level_weights = zeros(max(unique_levels)+1, 1);  
      
    for lvl = unique_levels'  
        idx = (levels == lvl);  
        level_voi = sum(VOI(idx));  
        level_weights(lvl+1) = level_voi / total_voi;  
    end  
      
    % Calculate metrics for each point  
    for i = 1:n  
        if levels(i) >= 0  
            lvl = levels(i);  
              
            % Extraction efficiency: Lower level means higher efficiency  
            efficiency(i) = 1 / (lvl + 1);  
              
            % Contribution = (Point VOI / Level VOI sum) * Level weight  
            lvl_points = find(levels == lvl);  
            lvl_voi_sum = sum(VOI(lvl_points));  
            voi_contribution = VOI(i) / (lvl_voi_sum + eps);  
              
            contribution(i) = voi_contribution * level_weights(lvl+1);  
        else  
            efficiency(i) = 0;  
            contribution(i) = 0;  
        end  
    end  
end  
  
%% ================= Other Helper Functions ================= %%  
function [VOI_std, centrality_std] = standardizeFeatures(VOI, centrality)  
    % Feature standardization (preserving original distribution characteristics)  
    VOI_std = (VOI - min(VOI)) / (max(VOI) - min(VOI) + eps);  
    centrality_std = (centrality - min(centrality)) / (max(centrality) - min(centrality) + eps);  
end  
  
function adj_matrix = buildAdjacencyMatrix(x, y, threshold)  
    % Build spatial adjacency matrix (efficient vectorized implementation)  
    n = length(x);  
    coords = [x(:), y(:)];  
    dist_mat = pdist2(coords, coords);  
    adj_matrix = dist_mat <= threshold;  
    adj_matrix(logical(eye(n))) = 0; % Remove self-connections  
end  
  
function [candidates, gains] = generateCandidates(remaining, levels, x, y, VOI, centrality, coverage, gamma, params)
    % Generate candidate points (including spatial decay factor)
    n = numel(remaining);
    candidates = zeros(n,1);
    gains = zeros(n,1);
    
    selected = find(levels ~= -1);
    if isempty(selected)
        gains = ones(n,1);
        candidates = remaining;
        return;
    end
    
    for i = 1:n
        idx = remaining(i);
        % Coverage gain (Gaussian spatial decay)
        distances = pdist2([x(idx),y(idx)], [x(selected),y(selected)]);
        decay = exp(-(distances/params.dist_threshold).^2);
        cov_gain = coverage(idx) * max(decay);
        
        % Structural gain (considering network connections)
        connections = sum(distances <= params.dist_threshold);
        cent_gain = connections / numel(selected);
        
        % Comprehensive gain
        gains(i) = gamma * cov_gain + (1-gamma) * cent_gain;
        candidates(i) = idx;
    end
    
    [gains, sort_idx] = sort(gains, 'descend');
    candidates = candidates(sort_idx);
end

%% ====== Fixed selectPoints Function (Single Version) ====== %%
function [selected, remaining] = selectPoints(candidates, gains, remaining, params)
    % Simplified selection strategy: Select points meeting gain threshold
    selected = [];
    if isempty(gains)
        return;
    end
    
    min_gain = max(gains) * params.gain_threshold;
    
    % Select all points above threshold
    for i = 1:numel(candidates)
        if gains(i) < min_gain
            break;
        end
        selected = [selected; candidates(i)];
    end
    
    % Force minimum points requirement
    if numel(selected) < params.min_level_points
        needed = params.min_level_points - numel(selected);
        if numel(candidates) >= needed
            selected = [selected; candidates(1:needed)];
        end
    end
    
    % Update remaining points set
    remaining = setdiff(remaining, selected);
end

function redundancy = computeLevelRedundancy(prev_points, curr_points, VOI, centrality)  
    % Calculate information redundancy between adjacent levels  
    if isempty(prev_points) || isempty(curr_points)  
        redundancy = 0;  
        return;  
    end  
      
    % Extract features  
    prev_voi = VOI(prev_points);  
    curr_voi = VOI(curr_points);  
    prev_cent = centrality(prev_points);  
    curr_cent = centrality(curr_points);  
      
    % Calculate VOI feature similarity (Cosine Similarity)  
    mean_prev_voi = mean(prev_voi);  
    mean_curr_voi = mean(curr_voi);  
    std_prev_voi = std(prev_voi);  
    std_curr_voi = std(curr_voi);  
      
    voi_sim = (mean_prev_voi * mean_curr_voi) / (std_prev_voi * std_curr_voi + eps);  
      
    % Calculate centrality feature similarity  
    mean_prev_cent = mean(prev_cent);  
    mean_curr_cent = mean(curr_cent);  
    std_prev_cent = std(prev_cent);  
    std_curr_cent = std(curr_cent);  
      
    cent_sim = (mean_prev_cent * mean_curr_cent) / (std_prev_cent * std_curr_cent + eps);  
      
    % Calculate distribution overlap (Histogram Intersection)  
    [voi_overlap, ~] = histogramIntersection(prev_voi, curr_voi, 10);  
    [cent_overlap, ~] = histogramIntersection(prev_cent, curr_cent, 10);  
      
    % Comprehensive redundancy (Weighted average)  
    redundancy = 1 - 0.4*voi_sim - 0.3*cent_sim - 0.3*(voi_overlap + cent_overlap)/2;  
    redundancy = max(0, min(1, redundancy)); % Limit to [0,1] range  
end  
  
function [overlap, bins] = histogramIntersection(data1, data2, num_bins)  
    % Calculate histogram intersection of two datasets  
    combined = [data1(:); data2(:)];  
    bins = linspace(min(combined), max(combined), num_bins+1);  
      
    hist1 = histcounts(data1, bins, 'Normalization', 'probability');  
    hist2 = histcounts(data2, bins, 'Normalization', 'probability');  
      
    overlap = sum(min(hist1, hist2));  
end  
  
function level_redundancy = computeAllLevelRedundancy(levels, VOI, centrality)  
    % Calculate information redundancy between all adjacent levels  
    unique_levels = unique(levels(levels >= 0));  
    level_redundancy = [];  
      
    if numel(unique_levels) < 2  
        return;  
    end  
      
    % Calculate redundancy for each pair of adjacent levels  
    for i = 1:(numel(unique_levels)-1)  
        prev_level = unique_levels(i);  
        curr_level = unique_levels(i+1);  
          
        prev_points = find(levels == prev_level);  
        curr_points = find(levels == curr_level);  
          
        redundancy = computeLevelRedundancy(prev_points, curr_points, VOI, centrality);  
        level_redundancy(end+1, :) = [prev_level, curr_level, redundancy];  
    end  
end  
  
function [opt_levels, opt_cent] = optimizeRedundancy(levels, VOI, cent, params)  
    % Redundancy optimization (Level merging)  
    opt_levels = levels;  
    opt_cent = cent;  
      
    unique_levels = unique(levels(levels >= 0));  
    if numel(unique_levels) < 2  
        return;  
    end  
      
    for i = 2:numel(unique_levels)  
        prev = unique_levels(i-1);  
        curr = unique_levels(i);  
        prev_idx = (levels == prev);  
        curr_idx = (levels == curr);  
          
        % Calculate feature differences between levels  
        voi_diff = abs(mean(VOI(prev_idx)) - mean(VOI(curr_idx)));  
        cent_diff = abs(mean(cent(prev_idx)) - mean(cent(curr_idx)));  
          
        if voi_diff < 0.1 && cent_diff < 0.1  
            opt_levels(curr_idx) = prev;  
            fprintf('Merged Level %d to Level %d (Small feature difference: VOI=%.4f, Centrality=%.4f)\n', ...  
                    curr, prev, voi_diff, cent_diff);  
        end  
    end  
end  

%% ================= Enhanced Visualization Function ================= %%
function visualizeResults(x, y, levels, coverage, metrics, efficiency, contribution, ...
                          level_redundancy, transfer_efficiency, flow_intensity)  
    % Enhanced visualization (Five-panel plot)  
    figure('Position', [50 100 2500 300], 'Color', 'w');  
      
    % 1. Spatial Distribution Plot  
    subplot(1,5,1);  
    valid_levels = levels(levels >= 0);  
    unique_levels = unique(valid_levels);  
    colors = lines(length(unique_levels));  
      
    hold on;  
    for i = 1:length(unique_levels)  
        idx = levels == unique_levels(i);  
        scatter(x(idx), y(idx), 30, colors(i,:), 'filled');  
    end  
    title('Hierarchical Spatial Distribution');  
    xlabel('X Coordinate');  
    ylabel('Y Coordinate');  
    grid on;  
      
    % Create legend  
    legend_labels = arrayfun(@(l) sprintf('Level %d (%d points)', l, sum(levels == l)), unique_levels, 'UniformOutput', false);  
    legend(legend_labels, 'Location', 'best');  
      
    % 2. Extraction Efficiency Heatmap  
    subplot(1,5,2);  
    scatter(x, y, 50, efficiency, 'filled');  
    title('Extraction Efficiency Heatmap');  
    xlabel('X Coordinate');  
    ylabel('Y Coordinate');  
    colorbar;  
    colormap(jet);  
    grid on;  
      
    % 3. Contribution Heatmap  
    subplot(1,5,3);  
    scatter(x, y, 50, contribution, 'filled');  
    title('Contribution Heatmap');  
    xlabel('X Coordinate');  
    ylabel('Y Coordinate');  
    colorbar;  
    colormap(jet);  
    grid on;  
      
    % 4. Redundancy Analysis  
    subplot(1,5,4);  
    if ~isempty(level_redundancy)  
        bar(level_redundancy(:,3), 'FaceColor', [0.7, 0.2, 0.2]);  
        title('Inter-level Information Redundancy');  
        xlabel('Level Pair');  
        ylabel('Redundancy');  
        xticks(1:size(level_redundancy,1));  
        xticklabels(arrayfun(@(i) sprintf('%d-%d', level_redundancy(i,1), level_redundancy(i,2)), 1:size(level_redundancy,1), 'UniformOutput', false));  
        grid on;  
        ylim([0 1]);  
    else  
        text(0.5, 0.5, 'No Level Redundancy Data', 'HorizontalAlignment', 'center');  
        title('Inter-level Information Redundancy');  
        axis off;  
    end  
      
    % 5. Information Transfer Efficiency Analysis  
    subplot(1,5,5);  
    if ~isempty(transfer_efficiency) && ~isempty(flow_intensity)
        % Create level labels
        levels = 0:(size(transfer_efficiency,1)-1);
        level_labels = arrayfun(@num2str, levels, 'UniformOutput', false);
        
        % Draw heatmap
        imagesc(transfer_efficiency);
        title('Inter-level Information Transfer Efficiency');
        xlabel('Target Level');
        ylabel('Source Level');
        colorbar;
        colormap(jet);
        
        % Add text annotations
        for i = 1:size(transfer_efficiency,1)
            for j = 1:size(transfer_efficiency,2)
                if transfer_efficiency(i,j) > 0
                    text(j, i, sprintf('%.2f\n(%.2f)', transfer_efficiency(i,j), flow_intensity(i,j)), ...
                         'HorizontalAlignment', 'center', 'FontSize', 8);
                end
            end
        end
        
        set(gca, 'XTick', 1:length(levels));
        set(gca, 'XTickLabel', level_labels);
        set(gca, 'YTick', 1:length(levels));
        set(gca, 'YTickLabel', level_labels);
    else
        text(0.5, 0.5, 'No Transfer Efficiency Data', 'HorizontalAlignment', 'center');
        title('Inter-level Information Transfer Efficiency');
        axis off;
    end
    
    % Add total coverage information  
    total_coverage = sum(coverage);  
    final_coverage = cumsum([metrics.coverage]);  
    final_coverage = final_coverage(end);  
    annotation('textbox', [0.4, 0.05, 0.2, 0.05], 'String', ...  
        sprintf('Total Coverage: %.2f%%', (final_coverage/total_coverage)*100), ...  
        'FitBoxToText', 'on', 'BackgroundColor', 'white');  
end  

%% ================= Enhanced Results Saving Function ================= %%
function saveResults(levels, coverage, VOI, centrality, efficiency, contribution, ...
                     metrics, level_redundancy, transfer_efficiency, flow_intensity)  
    % Results saving (Enhanced format)  
    result = table(levels, coverage, VOI, centrality, efficiency, contribution, ...  
        'VariableNames', {'Level', 'CoverageArea', 'VOI', 'Centrality', 'Efficiency', 'Contribution'});  
    writetable(result, 'Hierarchical_Results.xlsx');  
      
    % Save point-level metrics statistics  
    point_metrics = table(...  
        mean(VOI), std(VOI), min(VOI), max(VOI), ...  
        mean(centrality), std(centrality), min(centrality), max(centrality), ...  
        mean(efficiency), std(efficiency), min(efficiency), max(efficiency), ...  
        mean(contribution), std(contribution), min(contribution), max(contribution), ...  
        'VariableNames', {'VOI_Mean', 'VOI_Std', 'VOI_Min', 'VOI_Max', ...  
                          'Centrality_Mean', 'Centrality_Std', 'Centrality_Min', 'Centrality_Max', ...  
                          'Efficiency_Mean', 'Efficiency_Std', 'Efficiency_Min', 'Efficiency_Max', ...  
                          'Contribution_Mean', 'Contribution_Std', 'Contribution_Min', 'Contribution_Max'});  
    writetable(point_metrics, 'Point_Level_Metrics_Statistics.xlsx');  
      
    % Save level metrics  
    levels_axis = [metrics.level]';  
    voi_sum = [metrics.voi_sum]';  
    cent_sum = [metrics.cent_sum]';  
    coverage_sum = [metrics.coverage]';  
    point_count = [metrics.point_count]';  
    redundancy = [metrics.redundancy]';  
      
    level_table = table(levels_axis, voi_sum, cent_sum, coverage_sum, point_count, redundancy, ...  
        'VariableNames', {'Level', 'TotalVOI', 'TotalCentrality', 'TotalCoverage', 'PointCount', 'Redundancy'});  
    writetable(level_table, 'Level_Metrics.xlsx');  
      
    % Save level redundancy  
    if ~isempty(level_redundancy)  
        redundancy_table = array2table(level_redundancy, 'VariableNames', {'Level_i', 'Level_j', 'Redundancy'});  
        writetable(redundancy_table, 'Level_Redundancy.xlsx');  
    end  
    
    % Save information transfer efficiency
    if ~isempty(transfer_efficiency)
        % Create row and column labels
        levels = 0:(size(transfer_efficiency,1)-1);
        row_labels = arrayfun(@(x) sprintf('Level_%d', x), levels, 'UniformOutput', false);
        col_labels = arrayfun(@(x) sprintf('Level_%d', x), levels, 'UniformOutput', false);
        
        % Convert efficiency matrix to table
        eff_table = array2table(transfer_efficiency, 'RowNames', row_labels, 'VariableNames', col_labels);
        writetable(eff_table, 'Information_Transfer_Efficiency.xlsx', 'WriteRowNames', true);
        
        % Save flow intensity
        flow_table = array2table(flow_intensity, 'RowNames', row_labels, 'VariableNames', col_labels);
        writetable(flow_table, 'Flow_Intensity.xlsx', 'WriteRowNames', true);
    end
end