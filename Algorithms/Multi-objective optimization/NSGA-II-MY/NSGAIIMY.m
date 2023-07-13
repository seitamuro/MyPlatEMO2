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
            fitness = calcFitnessByConvhull(Population.objs);

            [~, maxIdx] = max(fitness);
            plot_fitness_convhull = [mean(calcFitnessByConvhull(Population.objs))];
            plot_fitness_delaunay = [mean(calcFitnessByDelaunay(Population.objs))];
            DT = delaunay(Population.objs);

            cnt = 1;

            DTArc = containers.Map(cnt, DT);
            PopulationArc = containers.Map(cnt, Population.objs);
            %plot_fitness_convhull = containers.Map(1, [0 total_sum_edges(Population.objs)]);
            %plot_fitness_delaunay = containers.Map

            fig = uifigure();
            sld = uislider(fig, "Position", [150 150 120 3], "ValueChangedFcn", @(sld, event) updateSld(sld, DTArc, PopulationArc, plot_fitness_convhull, plot_fitness_delaunay));
            field = uieditfield(fig, "numeric", "Position", [150, 180, 100, 22], "ValueChangedFcn", @(field, ~) changeSliderValue(sld, field, DTArc, PopulationArc, plot_fitness_convhull, plot_fitness_delaunay));
            updateSldProp(sld, cnt);
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,-fitness);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                allPopulation = [Population Offspring];
                [~, sortIndex] = sort(calcFitnessByConvhull(allPopulation.objs),'descend');
                Population = allPopulation(sortIndex(1:Problem.N));
                fitness = calcFitnessByConvhull(Population.objs);

                cnt = cnt + 1;
                disp("----plot_fitness_convhull");
                disp(plot_fitness_convhull(cnt-1));
                disp("---mean");
                disp(mean(calcFitnessByConvhull(Population.objs)));
                disp(length(horzcat(plot_fitness_convhull(cnt-1), mean(calcFitnessByConvhull(Population.objs)))));
                disp(length(plot_fitness_convhull));
                plot_fitness_convhull(cnt) = horzcat(plot_fitness_convhull(cnt-1), mean(calcFitnessByConvhull(Population.objs)));
                plot_fitness_delaunay(cnt) = horzcat(plot_fitness_delaunay(cnt-1), mean(calcFitnessByDelaunay(Population.objs)));
                updateSldProp(sld, cnt);
                all = Population.objs;
                DTArc(cnt) = delaunay(all);
                PopulationArc(cnt) = all;
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
    DT = delaunay(all);  

    % 各三角形の辺の長さの合計を求める
    total = sum(sqrt(sum((all(DT(:, 1), :) - all(DT(:, end), :)).^2, 1)));
    for i=1:size(all, 2)
        total = total + sum(sqrt(sum((all(DT(:, i), :) - all(DT(:, i+1), :)).^2, 1)));
    end
end

function fitness = calcFitnessByConvhull(all)
    % allの凸包を計算し、その面積をtotal_areaに保存
    [~, total_area] = convhull(all);
    
    % ドロネー三角形分割し、すべての辺の長さの総和を求める
    %all = [1 0; 0 0; 0 1; 1 1; 2 2];
    %total_edges = total_sum_edges(all);

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
    % allの凸包を計算し、その面積をtotal_areaに保存
    % [~, total_area] = convhull(all);
    
    % ドロネー三角形分割し、すべての辺の長さの総和を求める
    all = [1 0; 0 0; 0 1; 1 1; 2 2];
    total_edges = total_sum_edges(all);

    % fitnessを初期化
    fitness = zeros(length(all),1);
    
    for i = 1:length(all)
        % allから現在の要素を除去
        subset = all;
        subset(i,:) = [];
        
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
