function result = Opt_Optimization
%OPTIMIZATION provides the interface between the optimization solver and
%the objective function
%
% last revised: 2021/03/23

global p
global d
global stopOptimization

stopOptimization = 0;

%% Optimization solver options
    %% Fminsearch and Fmincon
    p.FminOptions = optimset('MaxIter',     p.MaxIter,...
                             'MaxFunEvals', p.MaxEval,...
                             'Display',     'iter');
                       
    %% fSSm
    if ( p.LogPara )
        p.fSSmOptions.log_var  = []; % 1:p.In.nOpt+p.Ex.nOpt;
    else
        p.fSSmOptions.log_var  = 1:p.In.nOpt+p.Ex.nOpt;
    end
    p.fSSmOptions.maxtime      = p.MaxTime;            % maximal time limit in seconds (e.g. 1 day = 60*60*24)
    p.fSSmOptions.maxeval      = p.MaxEval;            % maximal number of function evaluations (e.g. 1000000)
    p.fSSmOptions.local.solver = 0;
    
    if ( p.SmallRefset )
        nvar = p.In.nOpt + p.Ex.nOpt;
        p.fSSmOptions.ndiverse     = 4*nvar; % 2 = arbitrary factor, normally 10, changed to 4 here
        tmp = roots([1, -1, -3*nvar]);       % 3 = arbitrary factor, normally 10
        p.fSSmOptions.dim_refset   = ceil(tmp(tmp>0));
        if ( mod(p.fSSmOptions.dim_refset,2) )
            p.fSSmOptions.dim_refset = p.fSSmOptions.dim_refset + 1; % dim_refset should be an even number
        end
    end
    
    %% SRES
    p.SRESOptions.lambda = 100;                        % population size (number of offspring) (100 to 200)
    p.SRESOptions.G      = 50;                         % maximum number of generations
    p.SRESOptions.mu     = 1/7*p.SRESOptions.lambda;   % parent number (mu/lambda usually 1/7)
    p.SRESOptions.pf     = 0.45;                       % pressure on fitness in [0 0.5] try around 0.45
    p.SRESOptions.varphi = 1;                          % expected rate of convergence (usually 1)
    
result = [];    
    
%% Prepare optimization of the intracellular model
% Find position of initial conditions subject to optimization in p.In.Ic
if(isempty(p.In.IcOpt))
    p.In.PosIcOpt = [];
else
    p.In.PosIcOpt = zeros(length(p.In.IcOpt),1);
    for CountIc=1:length(p.In.IcOpt)
        if(sum(strcmp(p.In.IcOpt(CountIc), p.In.StateNames))==0)
            disp(char(strcat('Warning(Optimization): "', p.In.IcOpt(CountIc), '" is no initial condition in the intracellular model.')));
            return
        else
            p.In.PosIcOpt(CountIc) = find(strcmp(p.In.IcOpt(CountIc), p.In.StateNames));
        end
    end
end

% Set optimization bounds and initial guess for the intracellular model
p.In.x_L = zeros(p.In.nOpt,1);
p.In.x_U = zeros(p.In.nOpt,1);
p.In.x0  = zeros(p.In.nOpt,1);
for CountOptPa=1:length(p.In.PaOpt)
    if(sum(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames))==0)
        disp(strcat('Warning(Optimization): "', p.In.PaOpt(CountOptPa),...
            '" is no parameter in the intracellular model.'));
        return
    end
    p.In.x_L(CountOptPa) = p.In.ParaValues(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames))/p.In.PaDev;
    p.In.x_U(CountOptPa) = p.In.ParaValues(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames))*p.In.PaDev;
    p.In.x0(CountOptPa)  = p.In.ParaValues(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames));
end
if(~isempty(p.In.IcOpt))
    for CountOptIc=1:length(p.In.IcOpt)
        p.In.x_L(length(p.In.PaOpt)+CountOptIc) = p.In.Ic.STV(strcmp(p.In.IcOpt(CountOptIc), p.In.StateNames))/p.In.IcDev;
        p.In.x_U(length(p.In.PaOpt)+CountOptIc) = p.In.Ic.STV(strcmp(p.In.IcOpt(CountOptIc), p.In.StateNames))*p.In.IcDev;
        p.In.x0(length(p.In.PaOpt)+CountOptIc)  = p.In.Ic.STV(strcmp(p.In.IcOpt(CountOptIc), p.In.StateNames));
    end
end


% Find position of states subject to optimization in the intracellular model
for CountFit=1:length(p.In.States2Fit)
    if(sum(strcmp(p.In.States2Fit(CountFit), [p.In.StateNames; p.In.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.In.States2Fit(CountFit), '" is no state or variable in the intracellular model.'));
        return
    elseif(sum(strcmp(p.In.States2Fit(CountFit), d.In.ExpStateNames))==0)
        disp(strcat('Warning(Optimization): "', p.In.States2Fit(CountFit), '" has not been measured.'));
        return
    elseif ( isempty(p.In.States2Fit) )
        p.In.PosState2Fit = [];
    else
        p.In.PosState2Fit(CountFit, 1) = find(strcmp(p.In.States2Fit(CountFit), [p.In.StateNames; p.In.VariableNames]));
        p.In.PosState2Fit(CountFit, 2) = find(strcmp(p.In.States2Fit(CountFit), d.In.ExpStateNames));
    end
end

% Find index for normalization to inital state value
if(strcmpi(p.In.ObjFunCal, 'normtoinitialstate') || strcmpi(p.In.ObjFunCal, 'normtostate'))
    if(sum(strcmp(p.In.NormState, d.In.ExpStateNames))==0 ||...
       sum(strcmp(p.In.NormState, [p.In.StateNames; p.In.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.In.NormState, '" is no state in the intracellular model and/or data set.'));
        return
    else
        p.In.IndexNormSim = find(strcmp(p.In.NormState, [p.In.StateNames; p.In.VariableNames]));
        p.In.IndexNormExp = find(strcmp(p.In.NormState, d.In.ExpStateNames));
    end
end

%% Prepare optimization of the extracellular model
% Find position of initial conditions subject to optimization in p.Ex.Ic
if(isempty(p.Ex.IcOpt))   % was In -> changed
    p.Ex.PosIcOpt = [];
else
    p.Ex.PosIcOpt = zeros(length(p.Ex.IcOpt),1);
    for CountIc=1:length(p.Ex.IcOpt)
        if(sum(strcmp(p.Ex.IcOpt(CountIc), p.Ex.StateNames))==0)
            disp(char(strcat('Warning(Optimization): "', p.Ex.IcOpt(CountIc), '" is no initial condition in the extracellular model.')));
            return
        else
            p.Ex.PosIcOpt(CountIc) = find(strcmp(p.Ex.IcOpt(CountIc), p.Ex.StateNames));
        end
    end
end

% Set optimization bounds and initial guess for the extracellular model
p.Ex.x_L = zeros(p.Ex.nOpt,1);
p.Ex.x_U = zeros(p.Ex.nOpt,1);
p.Ex.x0  = zeros(p.Ex.nOpt,1);
for CountOptPa=1:length(p.Ex.PaOpt)
    if(sum(strcmp(p.Ex.PaOpt(CountOptPa), p.Ex.ParaNames))==0)
        disp(strcat('Warning(Optimization): "', p.Ex.PaOpt(CountOptPa),...
            '" is no parameter in the extracellular model.'));
        return
    end
    p.Ex.x_L(CountOptPa) = p.Ex.(char(p.Ex.PaOpt(CountOptPa)))/p.Ex.PaDev;
    p.Ex.x_U(CountOptPa) = p.Ex.(char(p.Ex.PaOpt(CountOptPa)))*p.Ex.PaDev;
    p.Ex.x0(CountOptPa)  = p.Ex.(char(p.Ex.PaOpt(CountOptPa)));
end
if(~isempty(p.Ex.IcOpt))
    for CountOptIc=1:length(p.Ex.IcOpt)
        p.Ex.x_L(length(p.Ex.PaOpt)+CountOptIc) = p.Ex.y0(strcmp(p.Ex.IcOpt(CountOptIc), p.Ex.StateNames))/p.Ex.IcDev;
        p.Ex.x_U(length(p.Ex.PaOpt)+CountOptIc) = p.Ex.y0(strcmp(p.Ex.IcOpt(CountOptIc), p.Ex.StateNames))*p.Ex.IcDev;
        p.Ex.x0(length(p.Ex.PaOpt)+CountOptIc)  = p.Ex.y0(strcmp(p.Ex.IcOpt(CountOptIc), p.Ex.StateNames));
    end
end


% Find position of states subject to optimization in the extracellular model
for CountFit=1:length(p.Ex.States2Fit)
    if(sum(strcmp(p.Ex.States2Fit(CountFit), [p.Ex.StateNames p.Ex.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.Ex.States2Fit(CountFit), '" is no state or variable in the extracellular model.'));
        return
    elseif(sum(strcmp(p.Ex.States2Fit(CountFit), d.Ex.ExpStateNames))==0)
        disp(strcat('Warning(Optimization): "', p.Ex.States2Fit(CountFit), '" has not been measured.'));
        return
    else
        p.Ex.PosState2Fit(CountFit, 1) = find(strcmp(p.Ex.States2Fit(CountFit), [p.Ex.StateNames p.Ex.VariableNames]));
        p.Ex.PosState2Fit(CountFit, 2) = find(strcmp(p.Ex.States2Fit(CountFit), d.Ex.ExpStateNames));
    end
end

% Find index for normalization to inital state value
if(strcmpi(p.Ex.ObjFunCal, 'normtoinitialstate') || strcmpi(p.Ex.ObjFunCal, 'normtostate'))
    if(sum(strcmp(p.Ex.NormState, d.Ex.ExpStateNames))==0 ||...
       sum(strcmp(p.Ex.NormState, [p.Ex.StateNames p.Ex.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.Ex.NormState, '" is no state in the extracellular model and/or data set.'));
        return
    else
        p.Ex.IndexNormSim = find(strcmp(p.Ex.NormState, [p.Ex.StateNames p.Ex.VariableNames]));
        p.Ex.IndexNormExp = find(strcmp(p.Ex.NormState, d.Ex.ExpStateNames));
    end
end

%% Prepare vector to realize offset in RNA concentrations for intracellular model
% as free viral RNAs in the seed virus might adhere to cells and thus cause
% a constant offset in measurements we implement such an offset in the
% objective function (*_ObjFun.m)
% RnaNames = {'Rm5', 'RcSeg5', 'RvSeg5'};
% p.In.OffsetIndex = zeros(length(RnaNames), 2);
% if ( ~isempty(d.In.ExpStateNames) && ~isempty(p.InToPop.states) )
%     for CountRnas = 1:length(RnaNames)
%     %     p.In.OffsetIndex(CountRnas, 1) = find(strcmp(RnaNames(CountRnas), [p.In.StateNames; p.In.VariableNames]));
%         p.In.OffsetIndex(CountRnas, 1) = find(strcmp(RnaNames(CountRnas), p.InToPop.states));
%         p.In.OffsetIndex(CountRnas, 2) = find(strcmp(RnaNames(CountRnas), d.In.ExpStateNames));
%     end
% end

vRNANames = {'RvSeg1', 'RvSeg9', 'RvSeg5', 'Rm1', 'Rm9', 'Rm5'};
% vRNANames = p.In.States2Fit;
p.In.OffsetIndex = zeros(length(vRNANames), 2);
if ( ~isempty(d.In.ExpStateNames) && ~isempty(p.In.States2Fit) )
    for CountRnas = 1:length(vRNANames)
    %     p.In.OffsetIndex(CountRnas, 1) = find(strcmp(RnaNames(CountRnas), [p.In.StateNames; p.In.VariableNames]));
        p.In.OffsetIndex(CountRnas, 1) = find(strcmp(vRNANames(CountRnas), p.In.States2Fit));
        p.In.OffsetIndex(CountRnas, 2) = find(strcmp(vRNANames(CountRnas), d.In.ExpStateNames));
    end
end 

vRNANames = {'VOnlyRel', 'VTot', 'DTot'};
% check if eligible states are getting fitted
tmp = strcmp(repmat(vRNANames',1,length(p.Ex.States2Fit)), repmat(p.Ex.States2Fit,length(vRNANames),1));
vRNANames(sum(tmp,2) == 0 ) = [];

p.Ex.OffsetIndex = zeros(length(vRNANames), 2); 
for CountRnas = 1:length(vRNANames)
    p.Ex.OffsetIndex(CountRnas, 1) = find(strcmp(vRNANames(CountRnas), p.Ex.States2Fit)); % or SimStatesIn
    p.Ex.OffsetIndex(CountRnas, 2) = find(strcmp(vRNANames(CountRnas), d.Ex.ExpStateNames));
end

%% Check whether parameter occuring in both models are subject to optimization
p.In.IdxKfus   = false(1,p.In.nPaOpt);
p.In.IdxKfus   = strcmp('kFus', p.In.PaOpt);
p.In.IdxFm     = false(1,p.In.nPaOpt);
p.In.IdxFm     = strcmp('Fm', p.In.PaOpt);
p.In.IdxKsynM  = false(1,p.In.nPaOpt);
p.In.IdxKsynM  = strcmp('kSynM', p.In.PaOpt);

p.Ex.IdxF      = false(1,p.Ex.nPaOpt);
p.Ex.IdxF      = strcmp('F', p.Ex.PaOpt);
p.Ex.IdxSIE    = false(1,p.Ex.nPaOpt);
p.Ex.IdxSIE    = strcmp('SIE_Time', p.Ex.PaOpt);

if(sum(strcmp('kFus', p.Ex.PaOpt)))
    disp('Warning(Optimization): If kFus should be optimized include it only in p.In.PaOpt instead of p.Ex.PaOpt.');
    return
elseif ( sum(strcmp('Fm', p.Ex.PaOpt)) )
    disp('Warning(Optimization): If Fm should be optimized include it only in p.In.PaOpt instead of p.Ex.PaOpt.');
    return
elseif ( sum(strcmp('kSynM', p.Ex.PaOpt)) )
    disp('Warning(Optimization): If kSynM should be optimized include it only in p.In.PaOpt instead of p.Ex.PaOpt.');
    return
    
end

ParasInBothModels = {'NHi', 'NLo', 'kAtt', 'kEqHi', 'kEqLo', 'kEn'}; %, 'kApoT', 'kApoI', 'ASS', 'AST'};
for CountOptPa=1:length(p.In.PaOpt)
    if(sum(strncmp(p.In.PaOpt(CountOptPa), ParasInBothModels, 3)))
        disp('Warning(Optimization): Optimization of entry parameters which occur in both the intracellular and extracellular model is not implemented (except for kFus).');
        return
    end
end
for CountOptPa=1:length(p.Ex.PaOpt)
    if(sum(strncmp(p.Ex.PaOpt(CountOptPa), ParasInBothModels, 3)))
        disp('Warning(Optimization): Optimization of entry parameters which occur in both the intracellular and extracellular model is not implemented (except for kFus).');
        return
    end
end

%% Start parameter estimation

% Create safety measures in case network connection is lost or strange error occurs
p.safety.Bestf = 1e3;
p.safety.xbest = ones(size([p.In.x0;  p.Ex.x0],1),1);

% Set upper and lower bounds for special parameters
        p.In.x_U(strcmp('FPar',[p.In.PaOpt p.In.IcOpt])) = 1;
        
        p.In.x_L(strcmp('Fadv',p.In.PaOpt)) = 1e-3;
        
        p.Ex.x_L(strcmp('AST',p.Ex.PaOpt)) = 5;
        p.Ex.x_U(strcmp('AST',p.Ex.PaOpt)) = 20;

        p.Ex.x_L(strcmp('kApoI',p.Ex.PaOpt)) = 5e-2; % prevent regular fSSm failure      

        p.Ex.x_U(strcmp('mu_MODIP',p.Ex.PaOpt)) = 1;
        
        if ( p.LogPara )
            p.In.x_L = log10(p.In.x_L);
            p.In.x_U = log10(p.In.x_U);
            p.In.x0  = log10(p.In.x0);
            p.Ex.x_L = log10(p.Ex.x_L);
            p.Ex.x_U = log10(p.Ex.x_U);
            p.Ex.x0  = log10(p.Ex.x0);
        end                 

        p.problem.x_0     = [p.In.x0;  p.Ex.x0];
        p.problem.x_L     = [p.In.x_L; p.Ex.x_L];
        p.problem.x_U     = [p.In.x_U; p.Ex.x_U];      
        
        % initialize parameter set recording matrix
        p.record = zeros(0,size(p.problem.x_0,1));            
        
switch lower(p.Solver)
    case {'fminsearch'}
        %% Fminsearch
        [xbest, result.ObjFun] = ...
            fminsearch(@Opt_ObjFun_parallel, [p.In.x0; p.Ex.x0], p.FminOptions);
        result.In.xbest = xbest(1:p.In.nOpt);
        result.Ex.xbest = xbest(p.In.nOpt+1:end);

    case {'fmincon'}
        [xbest, result.ObjFun] = ...
            fmincon(@Opt_ObjFun_parallel, [p.In.x0; p.Ex.x0],[],[],[],[],...
            [p.In.x_L; p.Ex.x_L],[p.In.x_U; p.Ex.x_U],[],p.FminOptions);
        result.In.xbest = xbest(1:p.In.nOpt);
        result.Ex.xbest = xbest(p.In.nOpt+1:end); 

    case {'fssm'}
        %% fSSm
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));    

        p.problem.int_var = 0;
        p.problem.bin_var = 0;
        p.problem.c_L     = [];
        p.problem.c_U     = [];
        p.problem.N       = numel(d.In.Time) + numel(d.Ex.Time);

        p.problem.f       = 'Opt_ObjFun_parallel';
        [result.fSSm] = fssm_kernel(p.problem, p.fSSmOptions);           

        xbest  = result.fSSm.xbest;

        result.In.xbest  = xbest(1:p.In.nOpt);
        result.Ex.xbest  = xbest(p.In.nOpt+1:end);
        result.ObjFunVal = result.fSSm.fbest;
        
    case{'cmaes'}     

        %%% toolbox options
        N = numel(p.problem.x_0);
        
        % set sigma, confident search radius? ([0 1])
        max_insigma = ones(N,1)*p.SearchRadiusConfidence;
%         max_insigma = ones(N,1)*0.05;
%         max_insigma = [];  
        
        p.problem.insigma = max_insigma; % high influence of sigma in the optmization time

        % set optimizer options
        p.CMAESoptions              = cmaes;
        p.CMAESoptions.PopSize      = (4 + floor(3*log(N)));
        p.CMAESoptions.IncPopSize   = 2;
        p.CMAESoptions.LBounds      = [p.In.x_L; p.Ex.x_L];
        p.CMAESoptions.UBounds      = [p.In.x_U; p.Ex.x_U];
        p.CMAESoptions.MaxIter      = 1e8; 1e3*(N+5)^2/sqrt(N)*100;
        p.CMAESoptions.EvalParallel = 'off';
        p.CMAESoptions.TolX         = 1e-25; %1e-30*max(insigma);  % optional
%         p.CMAESoptions.TolFun       = 1e-18;  % optional
%         p.CMAESoptions.Restarts     = 10;     % optional, restart optimization when stuck and increase p.CMAESoptions.PopSize=p.CMAESoptions.PopSize*p.CMAESoptions.IncPopSize
        p.CMAESoptions.DiagonalOnly = 0; 1; % optional, 1 for large N
        p.CMAESoptions.CMA.active   = 2; % optional, 2 faster
        p.CMAESoptions.LogModulo    = 0; % optional
        p.CMAESoptions.DispModulo   = 5; % result display frequency, number of iteration not FunEvals
        p.CMAESoptions.Noise.on     = 0; % optional

        p.CMAESoptions.StopOnWarnings            = 'off';  % ''no''==''off''==0, ''on''==''yes''==1 ';
        p.CMAESoptions.StopOnEqualFunctionValues = 0; % '2 + N/3  % number of iterations';    
        p.CMAESoptions.StopOnStagnation          = 'off';
        
        inputsforobjfun             = {}; %In1 In2 ... Inn}; %optional
                        
%       [xbest, result.ObjFunVal, counteval, stopflag, out, result.bestever]
        [xbest, result.ObjFunVal, ~, ~, ~, result.bestever] = ...
             cmaes('Opt_ObjFun_parallel', p.problem.x_0, p.problem.insigma, p.CMAESoptions, inputsforobjfun);

        result.In.xbest  = xbest(1:p.In.nOpt)';
        result.Ex.xbest  = xbest(p.In.nOpt+1:end)';
       
    case {'sres'}
        %% SRES
        [xbest, result.SRES.Statistics,...
            result.SRES.Gm, result.SRES.conv_curve] =...
            SRES(@Opt_ObjFun,'min', [p.In.x_L; p.Ex.x_L; p.In.x_U p.Ex.x_U],...
            p.SRESOptions.lambda, p.SRESOptions.G,...
            p.SRESOptions.mu, p.SRESOptions.pf, p.SRESOptions.varphi);
        result.In.xbest = xbest(1:p.In.nOpt);
        result.Ex.xbest = xbest(p.In.nOpt+1:end);

    otherwise
        warning(char(['Optimization - No optimization algorithm "',p.solver,'" found!']));
        result = [];
end

