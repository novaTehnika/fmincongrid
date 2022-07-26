function x = fmincongrid(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fmincongrid.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% March 9, 2022
%
% PURPOSE:
% This function mimics the built-in MATLAB function fmincon(), except
% performs a grid search as the method of finding the minimum of the
% funtion under the given constaints.
%
% USAGE:
% Given: x0,A,b,Aeq,beq,lb,ub,nonlcon as used by fmincon()
%  NOTE: x0 is not used
%
% Specify the number of grid points for each input
% options.numSteps = [number of steps for input 1, for input 2,... ];
%
% Specify type of grid: either linearly spaced or logscale with base 10
% options.gridType = ['linear','log10','linear'];
%
% x = fmincongrid(fun,[],A,b,Aeq,beq,lb,ub,nonlcon,options)
%
% FUNCTION M-FILES
%
% UPDATES
% 3/9/2022 - Creation.
%
% Copyright (C) 2022  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m = length(lb); % number of inputs to function

    % Error checks: input dimensions
    assert(m == length(lb) ...
        && m == length(ub) ...
        && m == length(options.numSteps),...
        'Argument dimension are not compatible') 
    assert(isstring(options.gridType),'Grid types must given as strings.')

    for j = 1:m
        if options.gridType(j) == 'log10'
            assert(lb(j) > 0 && ub(j) > 0, ...
                ['For a base 10 log grid, intput ',num2str(j),...
                ', the range must be positive definite.'])
        end
    end

    try
        options.nRefine;
    catch
        options.nRefine = 0;
    end

    % determine the number of grid points for each function input and the 
    % scaling parameter, b based on the grid type
    numSteps = zeros(m,1);
    bgrid = zeros(m,1);
    for j = 1:m
        numSteps(j) = options.numSteps(j);
        switch options.gridType(j)
                case 'linear'
                    bgrid(j) = (ub(j) - lb(j))/(numSteps(j) - 1);
                case 'log10'
                    bgrid(j) = (log10(ub(j)) - log10(lb(j)))/(numSteps(j) - 1);
            otherwise
                error(['The grid type for input ',num2str(j),...
                    ' is not specified properly.'])
        end
    end

    % determine the number of combinations and sub-combinations
    n=1; % number of combinations
    nSub = ones(m,1); % number of combinations in remaining dimensions
    for j1 = 1:m
        n = n*numSteps(j1);
        for j2 = j1+1:m
            nSub(j1) = nSub(j1)*numSteps(j2);
        end
    end
    
    % fill array of all index combinations
    index = ones(n,m);
    for i = 1:n
        for j = 1:m
            index(i,j) = mod( ceil(i/nSub(j))-1, numSteps(j)) + 1;
        end
    end

    % calculate function result for all combinations
    x = zeros(m,1);
    funval = zeros(n,1);
    for i = 1:n
        % fill input array for current combination of parameters
         for j = 1:m
            switch options.gridType(j)
                case 'linear'
                    x(j) = bgrid(j)*(index(i,j)-1) + lb(j);
                case 'log10'
                    x(j) = 10^(bgrid(j)*(index(i,j)-1) + log10(lb(j)));
            end
         end

        % Evaluate the function for the current combination of parameters
        funval(i) = fun(x);

        % Evaluate whether constraints are met 
        % default to true incase constraint is no used
        con_lin = 1;
        if ~isempty(A); con_lin = (A*x <= b); end

        con_lineq = 1;
        if ~isempty(Aeq); con_lineq = (Aeq*x == beq); end

        con_nonlin = 1;
        con_nonlineq = 1;
        if ~isempty(nonlcon)
            [c,ceq] = nonlcon(x);
            if ~isempty(c); 
                con = (c <= 0); 
                for j = 1:length(con)
                    con_nonlin = con_nonlin*con(j);
                end
            end
            if ~isempty(ceq); 
                con = (ceq == 0); end
                for j = 1:length(con)
                    con_nonlineq = con_nonlineq*con(j);
                end
        end

        % set function value to nan is an constraint is not met
        if ~(con_lin && con_lineq && con_nonlin && con_nonlineq)
            funval(i) = NaN;
        end

    end


    % Find optimal result
    [~, iOpt] = min(funval);

    % Calculate inputs corresponding to the optimal result and return
    for j = 1:m
        switch options.gridType(j)
            case 'linear'
                x(j) = bgrid(j)*(index(iOpt,j)-1) + lb(j);
            case 'log10'
                x(j) = 10^(bgrid(j)*(index(iOpt,j)-1) + log10(lb(j)));
        end
    end

    % Refine result through recursion
    if options.nRefine > 0
        optionsref = options;
        optionsref.nRefine = options.nRefine - 1;
        lbref = zeros(size(lb));
        ubref = zeros(size(ub));
        for j = 1:m
            switch options.gridType(j)
                case 'linear'
                    lbref(j) = bgrid(j)*(index(iOpt,j)-2) + lb(j);
                    lbref(j) = max(lbref(j),lb(j));
                    ubref(j) = bgrid(j)*(index(iOpt,j)) + lb(j);
                    ubref(j) = min(ubref(j),ub(j));
                case 'log10'
                    lbref(j) = 10^(bgrid(j)*(index(iOpt,j)-2) + log10(lb(j)));
                    lbref(j) = max(lbref(j),lb(j));
                    ubref(j) = 10^(bgrid(j)*(index(iOpt,j)) + log10(lb(j)));
                    ubref(j) = min(ubref(j),ub(j));
            end
        end
        x = fmincongrid(fun,x0,A,b,Aeq,beq,lbref,ubref,nonlcon,optionsref);
    end



end