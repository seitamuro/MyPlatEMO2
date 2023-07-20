classdef NSGAIIMY < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            % [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            %fitness = calcFitnessByConvhull(Population.objs);
            fitness = calcFitnessByDelaunay(Population.objs);

            mean_fitness = zeros(1, 25000);
            min_fitness = zeros(1, 25000);
            max_fitness = zeros(1, 25000);

            cnt = 1;

            mean_fitness(cnt) = mean(fitness);
            min_fitness(cnt) = min(fitness);
            max_fitness(cnt) = max(fitness);

            figure;
            h1 = plot(1, mean_fitness(cnt), "r-");
            hold on;
            h2 = plot(1,  min_fitness(cnt), "g-");
            h3 = plot(1, max_fitness(cnt), "b-");
            xlabel("Generation");
            ylabel("Average Fitness");
            legend("Mean fitness", "Min fitness", "Max fitness");

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,-fitness);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                allPopulation = [Population Offspring];
                %[~, sortIndex] = sort(calcFitnessByConvhull(allPopulation.objs),'descend');
                [~, sortIndex] = sort(calcFitnessByDelaunay(allPopulation.objs),'descend');
                Population = allPopulation(sortIndex(1:Problem.N));
                %fitness = calcFitnessByConvhull(Population.objs);
                fitness = calcFitnessByDelaunay(Population.objs);

                cnt = cnt + 1;

                mean_fitness(cnt) = mean(fitness);
                min_fitness(cnt) = min(fitness);
                max_fitness(cnt) = max(fitness);
 
                set(h1, 'XData', 1:cnt, 'YData', mean_fitness(1:cnt));
                set(h2, 'XData', 1:cnt, 'YData', min_fitness(1:cnt));
                set(h3, 'XData', 1:cnt, 'YData', max_fitness(1:cnt));
                drawnow;
            end
        end
    end
end

function changeSliderValue(sld, field, DTArc, PopulationArc, plot_fitness_convhull, plot_fitness_delaunay)
sld.Value = field.Value;    
updateSld(sld, DTArc, PopulationArc, plot_fitness_convhull);
updateSld(sld, DTArc, PopulationArc, plot_fitness_delaunay);
end

function updateSldProp(sld, cnt)
    if(cnt == 1)
        updateSldProp(sld, 2);
    else
        sld.Limits = [1 cnt];
        sld.Value = cnt;
    end
end

function updateSld(sld, DTArc, PopulationArc, plot_fitness_convhull, plot_fitness_delaunay)
    tiledlayout(2, 1);
    nexttile
    all = PopulationArc(int32(sld.Value));
    triplot(DTArc(int32(sld.Value)), all(:, 1), all(:, 2));
    nexttile
    plot(plot_fitness_convhull(int32(sld.Value)));
    plot(plot_fitness_delaunay(int32(sld.Value)));
end

function total = total_sum_edges(all)
    % ドロネー三角形分割する

    % Find unique rows while preserving the order
    [uniqueData, ~, ic] = unique(all, 'rows', 'stable');
    
    % Count the number of occurrences of each unique row
    counts = histc(ic, 1:size(uniqueData, 1));
    
    % Find the indices of rows that occur more than once
    duplicateIndices = find(counts > 1);
    
    % Retrieve the actual duplicate rows
    duplicateData = uniqueData(duplicateIndices, :);

    disp("DuplicateIndices");
    disp(duplicateIndices);
    disp("DuplicateData");
    disp(duplicateData);

    DT = delaunay(all); 
    warning('error', 'MATLAB:delaunay:DupPtsDelaunayWarnId');

    % 各三角形の辺の長さの合計を求める
    total = sum(sqrt(sum((all(DT(:, 1), :) - all(DT(:, end), :)).^2, 1)));
    for i=1:size(all, 2)
        total = total + sum(sqrt(sum((all(DT(:, i), :) - all(DT(:, i+1), :)).^2, 1)));
    end
end

function fitness = calcFitnessByConvhull(all)
    % allの凸包を計算し、その面積をtotal_areaに保存
    [~, total_area] = convhull(all);

    % fitnessを初期化
    fitness = zeros(length(all),1);
    
    for i = 1:length(all)
        % allから現在の要素を除去
        subset = all;
        subset(i,:) = [];
        
        % subsetの凸包の面積を計算
        [~, subset_area] = convhull(subset);

        % ドロネー三角形分割し、すべての辺の長さの総和を求める
        %subset_edges = total_sum_edges(subset);

        % fitnessにtotal_areaからsubset_areaを引いた値を代入
        fitness(i) = total_area - subset_area;

        % fitnessのtotal_edgesからsubset_edgesを引いた値を代入
        %fitness(i) = total_edges - subset_edges;
    end
end

function fitness = calcFitnessByDelaunay(all)
    % ドロネー三角形分割し、すべての辺の長さの総和を求める
    total_edges = total_sum_edges(all);

    % fitnessを初期化
    fitness = zeros(length(all),1);
    
    for i = 1:length(all)
        % allから現在の要素を除去
        subset = unique(all, 'rows','stable');

        % all(i)と一致する行をsubsetから削除
        idx = ismember(subset, all(i,:), 'rows');
        subset(idx, :) = [];

        disp("subset");
        disp(subset);
        disp("duplicate check");
        disp(length(subset));
        disp(length(unique(subset, "rows", "stable")));
        disp("idx");
        disp(idx);

        % subsetの凸包の面積を計算
        % [~, subset_area] = convhull(subset);

        % ドロネー三角形分割し、すべての辺の長さの総和を求める
        subset_edges = total_sum_edges(subset);

        % fitnessにtotal_areaからsubset_areaを引いた値を代入
        % fitness(i) = total_area - subset_area;

        % fitnessのtotal_edgesからsubset_edgesを引いた値を代入
        fitness(i) = total_edges - subset_edges;
    end
end
