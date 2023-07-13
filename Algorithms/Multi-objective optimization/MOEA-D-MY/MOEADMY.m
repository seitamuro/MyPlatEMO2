classdef MOEADMY < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Multiobjective evolutionary algorithm based on decomposition
% type --- 1 --- The type of aggregation function

%------------------------------- Reference --------------------------------
% Q. Zhang and H. Li, MOEA/D: A multiobjective evolutionary algorithm based
% on decomposition, IEEE Transactions on Evolutionary Computation, 2007,
% 11(6): 712-731.
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
            %% Parameter setting
            type = Algorithm.ParameterSet(1);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M); % Algorithm/Utility functionsにある。ランダムに均一に散らばった点を生成する
            % ターミナルからはaddPathでUnility functionsフォルダを追加することで呼び出すことが出来る
            % Problem.N Populationの数。個体の数
            % Problem.M 目的関数の数
            T = ceil(Problem.N/10); % 渡された数値を丸める

            %% Detect the neighbours of each solution
            B = pdist2(W,W); % 観測値の2つの集合管のペアワイズ距離
            % ペアワイズ距離とは2点間の距離のこと。指定したアルゴリズムで距離を計算してくれる。デフォルトはユークリッド距離
            % 格子状に各点との距離が計算され、行列にその結果が格納される
            [~,B] = sort(B,2);
            B = B(:,1:T);% ソートして上位10%を選択?

            %% Generate random population
            Population = Problem.Initialization(); % 個体を初期化
            Z = min(Population.objs,[],1); % 最も目的関数を最小化できているものを取得
            %Z = [1 1];

            %% Optimization
            while Algorithm.NotTerminated(Population)
                %% For each solution
                for i = 1 : Problem.N
                    % Choose the parents
                    P = B(i,randperm(size(B,2)));

                    % Generate an offspring
                    Offspring = OperatorGAhalf(Problem,Population(P(1:2)));

                    % Update the ideal point
                    Z = min(Z,Offspring.obj);

                    % Update the neighbours
                    % Tchebycheff approach
                    g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
                    Population(P(g_old>=g_new)) = Offspring;
                end
            end
        end
    end
end