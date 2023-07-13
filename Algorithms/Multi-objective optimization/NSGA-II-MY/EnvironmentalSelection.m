function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = MyNDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function [FrontNo, MaxFNo] = MyNDSort(objs, cons, N)
    % FrontNoをInfで初期化
    FrontNo = inf(size(objs, 1), 1);
    P = [1 1; 1 0; 0 0; 0 1];

    for rank=1:N
        % iがInfで1:size(objs, 1)
        i = find(isinf(FrontNo));
        for i=1:length(i)
            % k = i以外の1:size(objs, 1)かつFrontNo(k)はinfである
            k = setdiff(1:size(objs, 1), i);
            k = k(isinf(FrontNo(k)));
            % もしFrontNoがinfなら
            if isinf(FrontNo(i))
                for j=1:size(objs, 2)
                    flag1 = all(objs(i, j) >= objs(k, j));
                    flag3 = all(objs(i, j) <= objs(k, j));

                    % もしobjs[i, j] >= objs[k, j] かつ いずれかのiが objs[i, j] > objs[k, j]
                    % もしくはobjs[i, j] <= objs[k, j]かつ いずれかのiが objs[i, j] < objs[k, j]
                    % ならFrontNo[i]はrank
                    if any(flag1 | flag3)
                        FrontNo(i) = rank;
                    end
                end
            end
        end
    end
    % FrontNoを転置する
    FrontNo = FrontNo';
    % MaxFNoをFrontNoの最大値にする
    MaxFNo = max(FrontNo);
end

function fitness = calcFitness(all)
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

        % fitnessにtotal_areaからsubset_areaを引いた値を代入
        disp(total_area);
        disp(subset_area);
        fitness(i) = total_area - subset_area;
    end
end
