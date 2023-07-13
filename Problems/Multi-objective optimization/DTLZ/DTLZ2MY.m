classdef DTLZ2MY < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Benchmark MOP proposed by Deb, Thiele, Laumanns, and Zitzler

%------------------------------- Reference --------------------------------
% K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, Scalable test problems
% for evolutionary multiobjective optimization, Evolutionary multiobjective
% Optimization. Theoretical Advances and Applications, 2005, 105-145.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    methods
        %% Default settings of the problem
        function Setting(obj)
            div = 4;
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 4; end
            obj.N = obj.D^div;
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        function Population = Initialization(obj,~)
        %Initialization - Generate multiple initial solutions.
        %
        %   P = obj.Initialization() randomly generates the decision
        %   variables of obj.N solutions and returns the SOLUTION objects.
        %
        %   P = obj.Initialization(N) generates N solutions.
        %
        %   This function is usually called at the beginning of algorithms.
        %
        %   Example:
        %       Population = Problem.Initialization()
        
            div = 4;
        	N = obj.N;
            PopDec = zeros(N,obj.D);
            Type   = arrayfun(@(i)find(obj.encoding==i),1:5,'UniformOutput',false);
            if ~isempty(Type{1})        % Real variables
                for i=1:N
                    k = i;
                    for j=1:obj.D
                        PopDec(i, j) = rem(k, div+1);
                        k = ceil(k / (div+1));
                    end
                end
                PopDec = PopDec / div;
                %disp(PopDec);
                %disp(size(PopDec));
            end
            if ~isempty(Type{2})        % Integer variables
                PopDec(:,Type{2}) = round(unifrnd(repmat(obj.lower(Type{2}),N,1),repmat(obj.upper(Type{2}),N,1)));
            end
            if ~isempty(Type{3})        % Label variables
                PopDec(:,Type{3}) = round(unifrnd(repmat(obj.lower(Type{3}),N,1),repmat(obj.upper(Type{3}),N,1)));
            end
            if ~isempty(Type{4})        % Binary variables
                PopDec(:,Type{4}) = logical(randi([0,1],N,length(Type{4})));
            end
            if ~isempty(Type{5})        % Permutation variables
                [~,PopDec(:,Type{5})] = sort(rand(N,length(Type{5})),2);
            end
            Population = obj.Evaluation(PopDec);
        end
        
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            g1      = abs(2 - PopDec(:,2));
            g2      = abs(sin(PopDec(:, 1) + pi/2));
            PopObj = [g1 g2];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end