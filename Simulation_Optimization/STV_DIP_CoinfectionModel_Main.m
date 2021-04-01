%%% -------------------------------------------------------------- %%%
%%% Evaluation, Optimization and Simulation of a Multi-scale model %%%
%%% of STV and DIP co-infection                                    %%%
%%% -------------------------------------------------------------- %%%

% last revised: 2021/03/23

%% Housekeeping commands I
clear all; 
close all; 
clear mex; 
clc;
global p d

% define compiler (example, if not already integrated)
% setenv('MW_MINGW64_LOC','C:\...\mingw64');

% check for IQM toolbox (formerly known as SB Toolbox 2)
if ( ~exist('IQMsimulate','file') )
    error(sprintf(['Please install the IQM toolbox as it is required for model simulation!\n', ...
                   'Check: https://iqmtools.intiquan.com/']));
end 

tmp = dir('*.mexa64'); % delete MEX files to prevent occasional crashes
for i = 1:length(tmp);  delete(tmp(i).name);  end

ModelPath = pwd; Index = find(filesep==ModelPath);
ModelPath = strcat(ModelPath(1:Index(end)), 'ModelAndExperiments', filesep);
SearchPath = path; path(ModelPath, path);
clear Index

TSTART = tic;

%% Start of complete simulation run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Info.StartTime = clock;

p.Ex = ExLev_ParameterDeclaration(p); % load extracellular parameter names and base values
p.Ex.MOIMODIPConditions = [0, repmat(1e-3, 1, 4), repmat(3, 1, 4), repmat(30, 1, 4); ...
                           0, repmat([0, 1e-3, 3, 30], 1, 3)]';                      

%% Main options
p.Task           = 'simulate';  % task: 'simulate', 'optimize', 'plot'

%%% experiment 1 is a mock experiment without STV or COI, so all simualtions
%%% start with condition 2 and the last one is 13
p.SimCond      = 2:13;    % if 'simulate' is set: only simulate + plot these init. conditions
p.OptCond      = 2:13;  % if 'optimize' is set: only estimate to data for these init. conditions, but plot all at the end
p.LoadResultPara = 'FinalParameterSet'; 

p.ParallelComp = 1; % 0 - disable parallel computing (still creates pool, but uses only 1 worker)
                    % 1 - enable parallel computing
p.NumWorkers   = 1; % number of workers

%% DIP MSM options
p.Ex.T0      = 2.2e6;    %cells/ml  initial number of cells in the extracellular model
p.Ex.IV_0    = 0;        %cells/ml  initial number of STV-inf. cells in the extracellular model
p.Ex.ICO_0   = 0;        %cells/ml  initial number of co-inf. cells in the extracellular model

% Define other parameters
p.In.ReleaseAccuracy = 0.95; % define how much of the intracellular release should be considered, smaller value = faster simulation
p.Apo.Increase       = 'sigmoidal'; % choose apoptosis rate increase: 'none', 'linear', 'sigmoidal', 'hill*' (*=1-4) 
p.DIPdeNovo          = 0;    % de novo generation of DIPs, not implemented
p.mRNAPrimTrans      = 1;    % use primary transcription in DIP-infected cells for DI and S5 mRNA
p.In.kSynV_Threshold = 1e-3; % MODIP/MOI ratio below which kSynV is not reduced

p.SimTime    = 61;          %h  maximum simulation time of extracellular model
p.In.SimTime = p.SimTime;   %h  maximum simulation time of intracellular model (initially set to extracellular model)
p.In.h       = 0.1;         %h  step-size used to solve the extracellular model with Euler's-Method   
p.Ex.h       = 0.1;         %h  step-size used to solve the extracellular model with Euler's-Method   
p.Ex.hs      = 0.02;        %h  small step-size, used when regular step-size induces oscillating behavior around low values

p.Ex.tspan = 0 : p.Ex.h : p.SimTime;    %h  time points used to solve the extracellular model by Euler's-Method 
p.In.tspan = 0 : p.In.h : p.In.SimTime; %h  time points used to solve the intracellular model

%% Optimization options
p.MaxTime = 48*60*60; %maximum run time of optimization solver in seconds

% intracellular model
p.In.PaOpt = {'kFus', 'Fm', 'Kr', 'kSynV', 'kRelRed', 'kRel', 'KvRel', 'Fadv'}; %parameters subject to optimization in the intracellular model          
p.In.IcOpt = {'FPar'};    %initial conditions subject to optimization in the intracellular model, set Vex as the first variable
p.In.PaDev = 1000;        %upper and lower bounds of parameters deviate by */ PaDev from initial parameter guess
p.In.IcDev = 10;          %upper and lower bounds of IC deviate by */ IcDev from IC
p.In.States2Fit = {'RvSeg1', 'RvSeg9', 'RvSeg5', 'Rm1', 'Rm9', 'Rm5'}; %states subject to fitting for the intracellular model
p.In.RelStateWeight = ones(size(p.In.States2Fit));     %relative weight of the deviation of each fitted state in the objective function

% extracellular model
p.Ex.PaOpt = {'mu_MODIP', 'kLys', 'kApoT', 'kApoI', 'AST', 'V1', 'V2'}; %parameters subject to optimization in the extracellular model          
p.Ex.IcOpt = {};          %initial conditions subject to optimization in the extracellular model
p.Ex.PaDev = 1000;        %upper and lower bounds of parameters deviate by */ PaDev from initial parameter guess
p.Ex.IcDev = 2;           %upper and lower bounds of IC deviate by */ IcDev from IC
p.Ex.States2Fit = {'T', 'I', 'Ia', 'Ta', 'VOnlyRel', 'VTot', 'DTot'}; %states subject to fitting for the extracellular model
p.Ex.RelStateWeight = ones(size(p.Ex.States2Fit));     %relative weight of the deviation of each fitted state in the objective function

p.In.nPaOpt = length(p.In.PaOpt);                          %number of parameters subject to optimization
p.In.nOpt   = length(p.In.PaOpt)+length(p.In.IcOpt);       %number of parameters and initial conditions subject to optimization 
p.Ex.nPaOpt = length(p.Ex.PaOpt);                          %number of parameters subject to optimization
p.Ex.nOpt   = length(p.Ex.PaOpt)+length(p.Ex.IcOpt);       %number of parameters and initial conditions subject to optimization 

% Optimization minor options
p.Solver  = 'CMAES';             %optimization solver: 'CMAES', 'fSSm', 'fminsearch', 'fmincon' (if available)
p.MaxIter = 1e10;                %maximum iterations
p.MaxEval = 1e9;                 %maximum function evaluations
p.SearchRadiusConfidence = 0.2;  % parameter for CMA-ES algorithm, defines search range (confidence radius)

p.LogPara        = 1;            % use logarithmic parameters for estimation (converted internally)
p.SmallRefset    = 1;            % reduce Refset size for faster (but less robust) estimation   

p.In.ObjFunCal = 'log';                %'normal': model and experiment are compared quantitatively 
p.Ex.ObjFunCal = 'log/normal';         %'Log': model states and experimental values are logarithmized prior to comparison
                                       %'log-extra_0_rule': simulation number that go to 0 are handled differently, calculate percentage of ExpVal to max(ExpVal)
                                       %'NormToMax': model and experimental values are normalized to maximum prior to comparison
                                       %'NormToState': model and experimental values are normalized to a specific state given by p.NormState prior to comparison
                                       %'NormToInitial': model and experimental values are normalized to inital values of each state prior to comparison
                                       %'NormToInitialState: model and experimental values are normalized to initial value of a specific state given by p.NormState prior to comparison
p.In.NormState = '';            %state variables and experiment values are normalized to the inital value of this state prior to their comparison when p.ObjFunCal = 'NormToInitialState'
p.Ex.NormState = ''; 

%% Model and experimental files
p.In.ExpName   = ['MOI_MODIP_data', filesep, 'InLev_MOIMODIPscatter_*']; % 'old_ISX_data/Data_ISX_intracellular_36h';           %experiment name for the intracellular level, 'old_ISX_data/Data_ISX_intracellular_36h'
p.Ex.ExpName   = ['MOI_MODIP_data', filesep, 'ExLev_MOIMODIPscatter_*'];           %experiment name for the extracellular level, 'old_ISX_data/Data_ISX_extracellular'    

p.In.ModelName = 'InLev_STV_DIP_Coinfection.txt';      %name of intracellular model
p.Ex.ModelName = @ExLev_EquationSystem;      %function handle to m-file containin the extracellular model    

%% Initial parameter guesses

% Intracellular model

p.In.ChangedParaNames  = {'Kr', 'kSynV', 'kRelRed', 'kRel', 'KvRel', 'Fadv'}; %names  of parameters that ought to be changed in the intracellular model compared to the parameter declaration file
p.In.ChangedParaValues = [7.8e3, 20.1,     4.1e-4,  6.15e3,   1.8,    0.32]; 

p.In.ChangedIcNames    = {'FPar'}; %names  of initial conditions that ought to be changed in the intracellular model compared to the parameter declaration file
p.In.ChangedIcValues   = [3.6e-3];

% Extracellular model
p.Ex.ChangedParaNames  = {'mu_MODIP', 'kLys', 'kApoT', 'kApoI', 'AST', 'V1', 'V2'}; %names  of parameters that ought to be changed in the extracellular model compared to the parameter declaration file
p.Ex.ChangedParaValues = [  0.63,      0.16,  1.18e-2,  0.27,   6.65,   5.2,  0.1]; 

p.Ex.ChangedIcNames    = {}; %names  of initial conditions that ought to be changed in the extracellular model compared to the parameter declaration file
p.Ex.ChangedIcValues   = [];     

% Parameters included on both scales
p.ChangedParaNames     = {'kFus', 'Fm', 'kSynM'}; %names  of parameters that ought to be changed in both the intracellular and extracellular model (kFus, Fm, kSynM)    
p.ChangedParaValues    = [58.3,   0.12,  1.73e5];  

%% Load model and apply base + changed parameters
% load intracellular models
InModel = IQMmodel(fullfile(ModelPath, p.In.ModelName)); 

p.In.FPar0 = []; % parameter for fraction of infectious virions released (FIVR)

if ( isempty(p.LoadResultPara) )  
    error('Please define a base parameter set!');
else
    % Get base parameters from previous estimation results (excluding initial conditions)
    try        
        tmp = load(['results/', p.LoadResultPara]);    
        if ( strcmp(p.Task,'plot') )
            load(['results/', p.LoadResultPara])
            p.OptCond
            VisualizeResults(result, d, p);
            return;
        end
    catch
        error('Failed to load previously estimated parameters as initial values. Check the result file name! (if some plots appear: error during plotting)');
    end    
    
    % intracellular parameters and parameters used on both levels
    p.BaseParaNames = {'kFus', 'Fm', 'kSynM'};
    p.BaseParaValues = [];   
    p.In.BaseParaNames  = {'Fadv', 'kEn', 'kFus', 'kSynV', 'kSynC', 'kSynM', 'kBindM1', 'kDegM', 'kRel', 'KvRel', 'Kr', 'kRelRed', 'Fm'}; % parameters that were under estimation at some point
    p.In.BaseParaValues = [];
    for i = 1 : length(p.In.BaseParaNames)
        
        if ( strcmp('kFus', p.In.BaseParaNames{i}) )
            if ( sum(strcmp('kFus', tmp.p.In.PaOpt)) )  
                p.BaseParaValues(1) = tmp.Opt.In.xbest(strcmp('kFus', tmp.p.In.PaOpt));
            else
                p.BaseParaValues(1) = tmp.p.Ex.kFus;
            end
        elseif ( strcmp('Fm', p.In.BaseParaNames{i}) )
            if ( sum(strcmp('Fm', tmp.p.In.PaOpt)) )  
                p.BaseParaValues(2) = tmp.Opt.In.xbest(strcmp('Fm', tmp.p.In.PaOpt));
            else
                p.BaseParaValues(2) = tmp.p.Ex.Fm;
            end
        elseif ( strcmp('kSynM', p.In.BaseParaNames{i}) )
            if ( sum(strcmp('kSynM', tmp.p.In.PaOpt)) )  
                p.BaseParaValues(3) = tmp.Opt.In.xbest(strcmp('kSynM', tmp.p.In.PaOpt));
            else
                p.BaseParaValues(3) = tmp.p.Ex.kSynM;
            end
        end

        if ( sum(strcmp(p.In.BaseParaNames{i}, tmp.p.In.PaOpt)) )
            p.In.BaseParaValues(strcmp(p.In.BaseParaNames{i},p.In.BaseParaNames)) = tmp.Opt.In.xbest(strcmp(p.In.BaseParaNames{i}, tmp.p.In.PaOpt));
        else
            p.In.BaseParaValues(strcmp(p.In.BaseParaNames{i},p.In.BaseParaNames)) = tmp.p.In.ParaValues(strcmp(p.In.BaseParaNames{i},tmp.p.In.ParaNames));
        end

    end

    % extracellular parameters
    p.Ex.BaseParaNames  = {'kApoT', 'kApoI', 'kDegV', 'kLys', 'ASS', 'AST', 'mu_MODIP', 'SIE_Time', 'V1', 'V2'}; % parameters that were under estimation at some point
    if ( strcmp('kSynV_base',tmp.p.Ex.PaOpt) )
         p.Ex.BaseParaNames = [p.Ex.BaseParaNames, {'kSynV_base'}];
    end
    
    p.Ex.BaseParaValues = [];    
    for i = 1 : length(p.Ex.BaseParaNames)
        if ( sum(strcmp(p.Ex.BaseParaNames{i}, tmp.p.Ex.PaOpt)) )
            p.Ex.BaseParaValues(strcmp(p.Ex.BaseParaNames{i}, p.Ex.BaseParaNames)) = tmp.Opt.Ex.xbest(strcmp(p.Ex.BaseParaNames{i}, tmp.p.Ex.PaOpt));
        else
            p.Ex.BaseParaValues(strcmp(p.Ex.BaseParaNames{i}, p.Ex.BaseParaNames)) = tmp.p.Ex.(p.Ex.BaseParaNames{i});
        end
    end 

    % initial condition
    if ( sum(strcmp('FPar', tmp.p.In.IcOpt)) )
        p.In.FPar0 = tmp.Opt.In.xbest(length(tmp.p.In.PaOpt) + strcmp('FPar',tmp.p.In.IcOpt));
    else
        p.In.FPar0 = tmp.p.In.Ic.COI(strcmp('FPar', tmp.p.In.StateNames));
    end    

end

% apply base parameters to the intracellular model
InModel = IQMparameters(InModel, horzcat(p.BaseParaNames, p.In.BaseParaNames), [p.BaseParaValues p.In.BaseParaValues]);

% apply base parameters to the extracellular model
h.Ex.BaseParaNames  = horzcat(p.BaseParaNames, p.Ex.BaseParaNames);
h.Ex.BaseParaValues = [p.BaseParaValues, p.Ex.BaseParaValues];
for CountPara = 1:length(h.Ex.BaseParaNames)
    if(isfield(p.Ex, h.Ex.BaseParaNames(CountPara)))
        p.Ex.(char(h.Ex.BaseParaNames(CountPara))) = h.Ex.BaseParaValues(CountPara);
    else
        disp(strcat('Warning(Main): Parameter "', h.Ex.BaseParaNames(CountPara),...
            '" is no parameter in the model.'));
%         return
    end
end

% apply changed parameters from guesses to the intracellular model
InModel = IQMparameters(InModel, horzcat(p.ChangedParaNames, p.In.ChangedParaNames), ...
                                 [p.ChangedParaValues p.In.ChangedParaValues]);

% apply changed parameters from guesses to the extracellular model
h.Ex.ChangedParaNames  = horzcat(p.ChangedParaNames, p.Ex.ChangedParaNames);
h.Ex.ChangedParaValues = [p.ChangedParaValues, p.Ex.ChangedParaValues];
for CountPara = 1:length(h.Ex.ChangedParaNames)
    if(isfield(p.Ex, h.Ex.ChangedParaNames(CountPara)))
        p.Ex.(char(h.Ex.ChangedParaNames(CountPara))) = h.Ex.ChangedParaValues(CountPara);
    elseif ( strcmp('kSupp', h.Ex.ChangedParaNames(CountPara)) )
        p.Ex.kSupp = h.Ex.ChangedParaValues(CountPara);
    elseif ( strcmp('SIE_Time', h.Ex.ChangedParaNames(CountPara)) )
        p.Ex.SIE_Time = h.Ex.ChangedParaValues(CountPara);        
    else
        disp(strcat('Warining(Main): Parameter "', h.Ex.ChangedParaNames(CountPara),...
            '" is no parameter in the model.'));
        return
    end
end

if ( sum(strcmp(p.In.ChangedIcNames, 'FPar')) )
    p.In.FPar0 = p.In.ChangedIcValues(strcmp(p.In.ChangedIcNames, 'FPar'));
end

p.Ex.kApoI_D = p.Ex.kApoT; % specific apoptosis rate for DIP-only infected cells, use kApoT for now
p.Ex.kDegD   = p.Ex.kDegV; % specific degradation rate for DIPs, use kDegV for now

%% load experimental data and merge with model
%intracellular model
if ( sum(p.In.ExpName == '*') )
    d.In.Time           = [];
    d.In.ExpStateNames  = [];
    d.In.ExpStateValues = zeros(0,0,0);   
    
    tmp = find(p.In.ExpName == filesep, 1, 'last');
    SubFolder = p.In.ExpName(1:tmp);
    D = dir(fullfile(ModelPath, SubFolder)); 

    Id_FileName = p.In.ExpName(tmp+1:end-1); % remove '*'
    Id_Length   = length(Id_FileName);
    for i = 1 : size(D,1)
        % skip invalid files
        if ( length(D(i).name) < Id_Length );  continue;  end
        if ( ~strcmp(D(i).name(1:Id_Length), Id_FileName) );  continue;  end
        
        % get MOI MODIP conditions from filename
        tmp = regexp(D(i).name(1:end-4), '_', 'split');
        if ( length(tmp) < 3 );  error('Error during loading of data.');  end
        
        if ( strcmpi(tmp{3}, 'mock') ) % mock infection
            tmp_IC = [0 0];
        elseif ( length(tmp) == 3 )
            tmp_IC = [str2double(regexprep(tmp{3}, 'MOI', '')), 0];
        else
            tmp_IC = [str2double(regexprep(tmp{3}, 'MOI', '')), str2double(regexprep(tmp{4}, 'MODIP', ''))];
        end
        
        for j = 1 : 2
            if ( tmp_IC(j) < 0 );  tmp_IC(j) = 10^tmp_IC(j);  end % transform -3 to 10^-3
        end
        
        % find corresponding position in IC order
        IC_pos = find(sum(repmat(tmp_IC,size(p.Ex.MOIMODIPConditions,1),1) == p.Ex.MOIMODIPConditions,2) == 2);

        % get data from file
        InMeasurements = IQMmeasurement(fullfile(ModelPath, SubFolder, D(i).name));
        [Tmp_Time, Tmp_ExpStateNames, Tmp_ExpStateValues] = IQMmeasurementdata(InMeasurements);     

        d.In.ExpStateValues(:,:,IC_pos) = Tmp_ExpStateValues;
        
        % time + state info should be the same, if not -> error
        if ( isempty(d.In.Time) && isempty(d.In.ExpStateNames) )
            d.In.Time = Tmp_Time;
            d.In.ExpStateNames = Tmp_ExpStateNames;
        else
            if ( sum(d.In.Time == Tmp_Time) < numel(d.In.Time) || ~isempty(setdiff(d.In.ExpStateNames, Tmp_ExpStateNames)) )
                error('Error during data loading.')
            end
        end
    end
else
    InMeasurements = IQMmeasurement(fullfile(ModelPath, strcat(p.In.ExpName, '.csv')));
    InExperiment   =  IQMexperiment(fullfile(ModelPath, strcat(p.In.ExpName, '.exp')));
    [d.In.Time, d.In.ExpStateNames, d.In.ExpStateValues] = IQMmeasurementdata(InMeasurements);
end

%extracellular model
if ( sum(p.Ex.ExpName == '*') )
    d.Ex.Time           = [];
    d.Ex.ExpStateNames  = [];
    d.Ex.ExpStateValues = zeros(0,0,0);   
    
    tmp = find(p.Ex.ExpName == filesep, 1, 'last');
    SubFolder = p.Ex.ExpName(1:tmp);
    D = dir(fullfile(ModelPath, SubFolder));

    Id_FileName = p.Ex.ExpName(tmp+1:end-1); % remove '*'
    Id_Length   = length(Id_FileName);
    for i = 1 : size(D,1)
        % skip invalid files
        if ( length(D(i).name) < Id_Length );  continue;  end
        if ( ~strcmp(D(i).name(1:Id_Length), Id_FileName) );  continue;  end
        
        % get MOI MODIP conditions from filename
        tmp = regexp(D(i).name(1:end-4), '_', 'split');
        if ( length(tmp) < 3 );  error('Error during data loading.');  end
        
        if ( strcmpi(tmp{3}, 'mock') ) % mock infection
            tmp_IC = [0 0];
        elseif ( length(tmp) == 3 )
            tmp_IC = [str2double(regexprep(tmp{3}, 'MOI', '')), 0];
        else
            tmp_IC = [str2double(regexprep(tmp{3}, 'MOI', '')), str2double(regexprep(tmp{4}, 'MODIP', ''))];
        end
        
        for j = 1 : 2
            if ( tmp_IC(j) < 0 );  tmp_IC(j) = 10^tmp_IC(j);  end % transform -3 to 10^-3
        end
        
        % find corresponding position in IC order
        IC_pos = find(sum(repmat(tmp_IC,size(p.Ex.MOIMODIPConditions,1),1) == p.Ex.MOIMODIPConditions,2) == 2);

        
        % get data from file
        ExMeasurements = IQMmeasurement(fullfile(ModelPath, SubFolder, D(i).name));
        [Tmp_Time, Tmp_ExpStateNames, Tmp_ExpStateValues] = IQMmeasurementdata(ExMeasurements);     
        
        d.Ex.ExpStateValues(:,:,IC_pos) = Tmp_ExpStateValues;
        
        % time + state info should be the same, if not -> error
        if ( isempty(d.Ex.Time) && isempty(d.Ex.ExpStateNames) )
            d.Ex.Time = Tmp_Time;
            d.Ex.ExpStateNames = Tmp_ExpStateNames;
        else
            if ( sum(d.Ex.Time == Tmp_Time) < numel(d.Ex.Time) || ~isempty(setdiff(d.Ex.ExpStateNames, Tmp_ExpStateNames)) )
                error('Error during data loading')
            end
        end
    end
else
    ExMeasurements = IQMmeasurement(fullfile(ModelPath, strcat(p.Ex.ExpName, '.csv')));
    ExExperiment   =  IQMexperiment(fullfile(ModelPath, strcat(p.Ex.ExpName, '.exp')));
    [d.Ex.Time, d.Ex.ExpStateNames, d.Ex.ExpStateValues] = IQMmeasurementdata(ExMeasurements);
end

% MOI 10^-3 + MODIP 0 or 10^-3: total FL titer measured shortly after
% infection, next time points gave no results, unlikely viruses were
% already produced in this capacity
d.Ex.ExpStateValues(3,8,2) = NaN;
d.Ex.ExpStateValues(2,8,3) = NaN;
fprintf('Info: Measurement data at early time point for the total number\n      of FL removed for low MOI and low MODIP conditions!\n\n');

for CountTime=1:length(d.Ex.Time)
    if(isempty(find(p.Ex.tspan==d.Ex.Time(CountTime), 1)))
        disp(strcat('Warining(Main): Choose p.Ex.h such that a vector 0:p.Ex.h:max([p.Ex.Time; d.In.Time]) contains all measurement time points.'));
    else
        p.Ex.TimeIndexData(CountTime) = find(p.Ex.tspan==d.Ex.Time(CountTime), 1);
    end
end

%% Prepare simulation of the intracellular model
p.In.VariableNames                = IQMvariables(InModel);
[p.In.StateNames, ~, p.In.Ic.STV] = IQMstates(InModel);
[p.In.ParaNames, p.In.ParaValues] = IQMparameters(InModel);

p.In.Ic.COI = p.In.Ic.STV;
try
    IQMmakeMEXmodel(InModel,     'InModel_MexFile');    

    p.CompileFlag = 1;  % compilation successful
catch ModelCompileError
    p.InModel = InModel;
    
    p.CompileFlag = 0;  % compilation failed
    fprintf(['Warning: Model could not be compiled into a MEX-file.\n', ...
             '         Simulation will take longer than using a compiler and\n', ...
             '         may generate additional warnings during simulation.\n', ...
             '         The usage of a compiler is advised.\n\n']);
end

% changed IC for FPar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( ~isempty(p.In.FPar0) )
    p.In.Ic.STV(strcmp('FPar', p.In.StateNames)) = p.In.FPar0;
    p.In.Ic.COI(strcmp('FPar', p.In.StateNames)) = p.In.FPar0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              

%% Prepare simulation of the extracellular model
result.Ex.states    = {'T',   'Ta',                                 ...
                              'Ia', 'V', 'VAttHi', 'VAttLo', 'VEn', ...
                       'I_D',       'D', 'DAttHi', 'DAttLo', 'DEn', ...
                                    'VTot', 'DTot', 'Par',          ...
                       'Cviab', 'Dead', ...
                       'VRelTot', 'VOnlyRel', 'DOnlyRel',...
                       'Rm5_ID', 'Rm9_ID'};
result.Ex.variables = {'I_V', 'I_V_SIE', 'I_CO', 'I', 'Log10Tcid50', 'Log10Ha'};

p.Ex.StateNames = result.Ex.states; p.Ex.VariableNames = result.Ex.variables;

for i = 1 : length(result.Ex.states)
    p.Idx.(result.Ex.states{i}) = strcmp(result.Ex.states{i},result.Ex.states);
end
for i = 1 : length(result.Ex.variables)
    p.Idx.(result.Ex.variables{i}) = strcmp(result.Ex.variables{i},result.Ex.variables);
end    

%% Calculate the integral of kApoI(tau), exp(-kApoI(tau)) 
% modify apoptosis rates and integrals
[p.Ex.ApoI, p.Ex.IntKapoi] = getApoI(p);

%% Simulate the intracellular model once to obtain reaction names
if(p.CompileFlag)        
    result.In.COI     = IQMPsimulate('InModel_MexFile', [0 1]); 
else
    result.In.COI     = IQMsimulate(p.InModel,'ode23s',[0 1], p.In.Ic.COI);
end

% divide variables and reactions for CoinfectionModel (as the reactions are also declared as variables)
for j = 1 : length(result.In.COI.variables)
    if ( strcmp('r', result.In.COI.variables{j}(1)) )
        break; % get j
    end
end
result.In.COI.reactions      = result.In.COI.variables(j:end);
result.In.COI.reactionvalues = result.In.COI.variablevalues(:,j:end);        

%% prepare saving of intracellular dynamics converted to population level                  
p.InToPop.StateVarNames = [p.In.StateNames; p.In.VariableNames];
p.InToPop.RecordedStates =  {'RcSeg1', 'RcSeg5', 'RcSeg9', ...
                             'RvSeg1', 'RvSeg5', 'RvSeg9', ...
                             'Rm1', 'Rm5', 'Rm9', ...
                             'VpNuc1', 'VpNucM11', 'VpCytM11',...
                             'VpNuc5', 'VpNucM15', 'VpCytM15',...
                             'VpNuc9', 'VpNucM19', 'VpCytM19',...
                             'vRNPnuc1', 'vRNPnuc5', 'vRNPnuc9',...
                             'VpNucPercentage5', 'vRNPnuc_DAPI',...
                             'Vcyt', 'Vrel', 'Prel', 'P_Rdrp', 'P_M1', 'P_Nep'};                 
                                                 
if ( ~strcmp(p.Task,'simulate') )      
    p.InToPop.states      = p.In.States2Fit;
    
    p.Ex.FullDistribution = cell(length(d.In.Time),1);
    
    p.InToPop.ExpTime = d.In.Time;
    for i = 1 : length(p.InToPop.states)
        p.InToPop.STV.statevalues{i} = zeros(length(p.In.tspan),length(d.In.Time));
        p.InToPop.COI.statevalues{i} = zeros(length(p.In.tspan),length(d.In.Time)); 
    end            
else
    p.InToPop.states        = p.InToPop.RecordedStates;      
    p.Ex.PopTick = p.SimTime/p.Ex.h;
    p.Ex.FullDistribution = cell(p.Ex.PopTick+1,1);

    p.InToPop.ExpTime = 0:p.SimTime/p.Ex.PopTick:p.SimTime;
    for i = 1 : length(p.InToPop.states)
        p.InToPop.STV.statevalues{i} = zeros(length(p.In.tspan),p.Ex.PopTick+1);
        p.InToPop.COI.statevalues{i} = zeros(length(p.In.tspan),p.Ex.PopTick+1); 
    end   
end  

for i = 1 : length(p.InToPop.states)
    p.InToPop.Idx.(p.InToPop.states{i}) = strcmp(p.InToPop.states{i},p.InToPop.states);
end

%% Solve the system using Euler's Method  
p.In.ReactionNames = result.In.COI.reactions;

switch (lower(p.Task))
    case {'simulate'}                     
%%% PARSIM
        NumSim = length(p.SimCond);
        for i = 1 : NumSim
            g(i) = p; % replacement structure for parallel simulation
        end        
        
        % open pool of workers
        CreatePool = 0;
        if ( ~p.ParallelComp )
            p.ParallelComp = 1;
            p.NumWorkers   = 1;
        end
        if ( p.ParallelComp && ~isempty(gcp('nocreate')) )
            tmp = gcp;
            if ( tmp.NumWorkers ~= p.NumWorkers )
                delete(gcp);
                CreatePool = 1;
            end
        elseif ( p.ParallelComp && isempty(gcp('nocreate')) )
            CreatePool = 1;
        end
        
        if ( CreatePool )
            if ( p.NumWorkers > 6 )
                error('Number of workers is too high! (but could be used when settings are changed)')
            else
                PoolObj = parpool(p.NumWorkers); 
            end
        end        
        
        % allocate results matrices
        result.Ex.statevalues    = zeros(length(p.Ex.tspan), length(result.Ex.states), size(p.Ex.MOIMODIPConditions,1));
        result.Ex.variablevalues = zeros(length(p.Ex.tspan), length(result.Ex.variables), size(p.Ex.MOIMODIPConditions,1));
        
        p.InToPop.populationvalues = zeros(p.Ex.PopTick+1,length(p.InToPop.states), size(p.Ex.MOIMODIPConditions,1));    

        for i = 1 : NumSim
            
            SimRun = p.SimCond(i);
            
            % define necessary intracellular simulation time to capture all measurement time points
            DataTimeLim_In = find(~isnan(d.In.ExpStateValues(:,1,SimRun)),1,'last'); % last time points measured in individual experiment  
            g(i).In.SimTime = d.Ex.Time(DataTimeLim_In) + g(i).In.h; % increase length by one tick to capture full dynamics of initially infected cells
            g(i).In.tspan = 0 : g(i).In.h : g(i).In.SimTime;               

            % effective MOI-dependent release rates for the intracellular model    
            g(i).DynMOI.STV.RrelV    = cell(length(g(i).Ex.tspan),1);
            g(i).DynMOI.STV.RrelD    = cell(length(g(i).Ex.tspan),1);
            g(i).DynMOI.STV.RrelDTot = cell(length(g(i).Ex.tspan),1);
            g(i).DynMOI.STV.RrelP    = cell(length(g(i).Ex.tspan),1);
            g(i).DynMOI.COI.RrelV    = cell(length(g(i).Ex.tspan),1);
            g(i).DynMOI.COI.RrelD    = cell(length(g(i).Ex.tspan),1);
            g(i).DynMOI.COI.RrelDTot = cell(length(g(i).Ex.tspan),1);
            g(i).DynMOI.COI.RrelP    = cell(length(g(i).Ex.tspan),1);               
            
            % initial conditions extracellular level
            g(i).Ex.Moi   = g(i).Ex.MOIMODIPConditions(SimRun,1);
            g(i).Ex.Modip = g(i).Ex.MOIMODIPConditions(SimRun,2);  
            
            g(i).Ex.y0 = zeros(1,length(result.Ex.states));

            g(i).Ex.y0(g(i).Idx.T)    = d.Ex.ExpStateValues(1, strcmp('T',d.Ex.ExpStateNames), SimRun);
            g(i).Ex.y0(g(i).Idx.Ta)   = d.Ex.ExpStateValues(1, strcmp('Ta',d.Ex.ExpStateNames), SimRun);
            g(i).Ex.y0(g(i).Idx.Ia)   = d.Ex.ExpStateValues(1, strcmp('Ia',d.Ex.ExpStateNames), SimRun);

            g(i).Ex.y0(g(i).Idx.V)    = g(i).Ex.Moi   * d.Ex.ExpStateValues(1, strcmp('Cviab',d.Ex.ExpStateNames), SimRun);
            g(i).Ex.y0(g(i).Idx.D)    = g(i).Ex.Modip * d.Ex.ExpStateValues(1, strcmp('Cviab',d.Ex.ExpStateNames), SimRun);
            g(i).Ex.y0(g(i).Idx.Par)  = 0; 
            
            % differentiate type of infected cell based on Poisson distribution of the initial MOI/MODIP conditions
            if ( SimRun < 4 )
                % set infected cells for very low MOI+MODIP conditions to 0, otherwise intracellular dynamics start too early
                g(i).Ex.IV_0          = 0; 
                g(i).Ex.ICO_0         = 0;
                g(i).Ex.y0(g(i).Idx.I_D) = 0;
            else
                tmp = poisspdf(0,[g(i).Ex.Moi, g(i).Ex.Modip]); % get Poisson values for MOI/MODIP
                
                InfPercentages   = [(1-tmp(1)) * tmp(2), ...     % I_V
                                    tmp(1)     * (1-tmp(2)), ... % I_D
                                    (1-tmp(1)) * (1-tmp(2))];    % I_CO
                
                g(i).Ex.IV_0  = InfPercentages(1)/sum(InfPercentages) * d.Ex.ExpStateValues(1, strcmp('I',d.Ex.ExpStateNames), SimRun);
                g(i).Ex.ICO_0 = InfPercentages(3)/sum(InfPercentages) * d.Ex.ExpStateValues(1, strcmp('I',d.Ex.ExpStateNames), SimRun);
                g(i).Ex.y0(g(i).Idx.I_D) = InfPercentages(2)/sum(InfPercentages) * d.Ex.ExpStateValues(1, strcmp('I',d.Ex.ExpStateNames), SimRun);
                
                % exception for MOI 0 + MODIP 0
                if ( sum(InfPercentages) == 0 )
                    g(i).Ex.IV_0             = 0;
                    g(i).Ex.ICO_0            = 0;
                    g(i).Ex.y0(g(i).Idx.I_D) = 0;
                end
            end

            % define MODIP/MOI-dependent kSynV
            if ( g(i).Ex.Modip > g(i).In.kSynV_Threshold ) 
                g(i).In.kSynV = g(i).In.ParaValues(strcmp(g(i).In.ParaNames, 'kSynV')) / (g(i).Ex.V1*(g(i).Ex.Modip/g(i).Ex.Moi).^g(i).Ex.V2);
            else
                g(i).In.kSynV = g(i).In.ParaValues(strcmp(g(i).In.ParaNames, 'kSynV'));
            end      
            
            % calculate low MOI infection dynamics to save computation time
            % output: g(i).In.LowMOI.STV/COI.RrelV/D/P, g(i).In.LowMOI.AllValuesSTV/COI
            g(i).xIn = [];            
            g(i) = LowMOI_PreSim(g(i));         
            
        end

        %% parallel calculation of infection dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h         = cell(1,2);
        resEx_tmp = [];
        
        fprintf('Simulating ... \n\n');
        
        parfor i = 1 : NumSim  
%             try
                % Simulate the extracellular model
                [h{i}, resEx_tmp(i).time, resEx_tmp(i).statevalues, resEx_tmp(i).variablevalues(:, g(i).Idx.I_V),...
                  resEx_tmp(i).variablevalues(:, g(i).Idx.I_V_SIE), resEx_tmp(i).variablevalues(:, g(i).Idx.I_CO),...
                  resEx_tmp(i).variablevalues(:, g(i).Idx.I), ~, ~] = ExLev_Euler(g(i)); 
%             catch 
%                 fprintf('Simulation run for condition %i failed.\n\n',g(i).SimCond(i));
%             end
        end         
      
        %% get all results in the correct place
        PopVal_tmp = p.InToPop.populationvalues;
        for k = 1 : NumSim
          
            p = h{k};  
            SimRun = p.SimCond(k);  

            result.Ex.time = resEx_tmp(k).time;
            result.Ex.statevalues(:,:,SimRun) = resEx_tmp(k).statevalues;
            result.Ex.variablevalues(:,1 : find(p.Idx.I==1),SimRun) = resEx_tmp(k).variablevalues;
                      
            result.Ex.variablevalues(:, strcmp('Log10Tcid50', result.Ex.variables), SimRun) =...
                log10(max(1, result.Ex.statevalues(:, p.Idx.V, SimRun)));   
            result.Ex.variablevalues(:, strcmp('Log10Ha', result.Ex.variables), SimRun) =...
                log10(max(1, result.Ex.statevalues(:, p.Idx.Par, SimRun)));              
            
            % calculate values at intracellular measurement time points 
            SumCells = result.Ex.statevalues(:,strcmp('Cviab',result.Ex.states),SimRun);

            for i = 1 : p.Ex.PopTick + 1                 
                for j = 1 : length(p.InToPop.states)
                    
                    min_tmp1 = min(size(p.Ex.FullDistribution{i},1),size(p.InToPop.STV.statevalues{j},1));

                    tmp = p.Ex.FullDistribution{i}(1:min_tmp1, 1) .* ...
                          p.InToPop.STV.statevalues{j}(1:min_tmp1, i);

                    min_tmp2 = min(size(p.Ex.FullDistribution{i},1),size(p.InToPop.COI.statevalues{j},1));
                    
                    if ( min_tmp1 ~= min_tmp2 )
                        if ( min_tmp1 < min_tmp2 )
                            tmp = [tmp; zeros(min_tmp2-min_tmp1,1)];
                        else
                            error('Error during simulation.');
                        end
                    end
                    tmp = tmp + p.Ex.FullDistribution{i}(1:min_tmp2, 2) .* ...
                          p.InToPop.COI.statevalues{j}(1:min_tmp2, i);

                    PopVal_tmp(i, j, SimRun) = sum(tmp) ./ SumCells(i);
                end           
            end  

            % add DIP-infected cell mRNA to others      
            if ( p.mRNAPrimTrans )      
                PopVal_tmp(:, strcmp('Rm5', p.InToPop.states), SimRun) = ...
                  PopVal_tmp(:, strcmp('Rm5', p.InToPop.states), SimRun) + ...
                  result.Ex.statevalues(:, strcmp('Rm5_ID', result.Ex.states), SimRun) .* ... 
                  result.Ex.statevalues(:,strcmp('I_D',result.Ex.states), SimRun) ./ SumCells;              

                PopVal_tmp(:, strcmp('Rm9', p.InToPop.states), SimRun) = ...
                  PopVal_tmp(:, strcmp('Rm9', p.InToPop.states), SimRun) + ...
                  result.Ex.statevalues(:, strcmp('Rm9_ID', result.Ex.states), SimRun) .* ... 
                  result.Ex.statevalues(:,strcmp('I_D',result.Ex.states), SimRun) ./ SumCells;      
            end
            
        end
        p.InToPop.populationvalues = PopVal_tmp;

        VisualizeResults(result, d, p); 

    case {'optimize'}

        % open pool of workers
        CreatePool = 0;
        if ( ~p.ParallelComp )
            p.ParallelComp = 1;
            p.NumWorkers   = 1;
        end
        if ( p.ParallelComp && ~isempty(gcp('nocreate')) )
            tmp = gcp;
            if ( tmp.NumWorkers ~= p.NumWorkers )
                delete(gcp);
                CreatePool = 1;
            end
        elseif ( p.ParallelComp && isempty(gcp('nocreate')) )
            CreatePool = 1;
        end
        
        if ( CreatePool )
            if ( p.NumWorkers > 6 )
                error('Number of workers is too high! (but could be used when settings are changed)')
            else
                PoolObj = parpool(p.NumWorkers); 
            end
        end

        % save time vectors + values before optimization
        p.SimTime_save    = p.SimTime;
        p.In.SimTime_save = p.In.SimTime; 
        p.Ex.tspan_save   = p.Ex.tspan;
        p.In.tspan_save   = p.In.tspan;

        % initial conditions extracellular level
        p.Ex.Moi   = p.Ex.MOIMODIPConditions(p.OptCond(1),1);
        p.Ex.Modip = p.Ex.MOIMODIPConditions(p.OptCond(1),2);  

        p.Ex.y0 = zeros(1,length(result.Ex.states));
        p.Ex.y0(p.Idx.T)    = p.Ex.T0;
%         p.Ex.y0(p.Idx.Ia) = 0;
        p.Ex.y0(p.Idx.V)    = p.Ex.Moi*p.Ex.y0(p.Idx.T);
        p.Ex.y0(p.Idx.D)    = p.Ex.Modip*p.Ex.y0(p.Idx.T);       
        p.Ex.y0(p.Idx.Par)  = 0; %20.6165*p.Ex.y0(p.Idx.V); % 638121526; % *= 20.6165*V0 *= 8.7414e6 * MOI

        disp('Optimizing...')
        fprintf('Start:         %s\n',datestr(p.Info.StartTime));
        if ( p.MaxTime >= 24*60*60 );  tmp = datevec(p.MaxTime./(60*60*24)) - [0 1 0 0 0 0]; 
        else                           tmp = datevec(p.MaxTime./(60*60*24));
        end
        fprintf('Projected end: %s\n',datestr(p.Info.StartTime + tmp));

        %% Optimize the multiscale model 
        p.TotalOptTime     = 0;
        p.TotalTimeCounter = 1;
        
        p.FailedRuns       = zeros(1,3);
        p.Error_099        = zeros(1,2);
        
        p.ParaOV           = zeros(0, numel(p.In.PaOpt)+numel(p.In.IcOpt)+numel(p.Ex.PaOpt)+numel(p.Ex.IcOpt));
        p.SuccessVector    = zeros(0,1);    
        
        Opt = Opt_Optimization;
        if(isempty(Opt));  disp('Error!');  return
        end
        
        % save optimization results to prevent data loss due to simualtion errors
        save('results/tmp_file.mat', 'd', 'Opt', 'p');
        
        % close pool of workers
        delete(gcp);

        %% re-set time span + scale parameters
        p.SimTime    = p.SimTime_save;
        p.In.SimTime = p.In.SimTime_save; 
        p.Ex.tspan   = p.Ex.tspan_save;
        p.In.tspan   = p.In.tspan_save;        
        
        if ( p.LogPara ) % re-transform logarithmic parameters and estimation bounds
            Opt.In.xbest  = 10.^Opt.In.xbest;
            Opt.Ex.xbest  = 10.^Opt.Ex.xbest;     

            p.problem.x_L = 10.^(p.problem.x_L);
            p.problem.x_U = 10.^(p.problem.x_U);
            p.problem.x_0 = 10.^(p.problem.x_0);                  
        end

        % Display estimated parameter values
        AdjWhiteSpace = cell(length(p.problem.x_0),2);          % adjustable white space, get values in right position
        tmp = getPotency(p.problem.x_0,2);  
        for i = 1 : length(p.problem.x_0)
            IdvWhiteSpace = max(tmp(:,2)) - tmp(i,2);
            if ( IdvWhiteSpace > max(tmp(:,2)) );  IdvWhiteSpace = max(tmp(:,2));  end
            for j = 1 : IdvWhiteSpace;  AdjWhiteSpace{i,1} = [AdjWhiteSpace{i,1} ' '];  end
        end
        tmp = getPotency([Opt.In.xbest, Opt.Ex.xbest],2);
        for i = 1 : length([Opt.In.xbest, Opt.Ex.xbest])
            IdvWhiteSpace = max(tmp(:,2)) - tmp(i,2);
            if ( IdvWhiteSpace > max(tmp(:,2)) );  IdvWhiteSpace = max(tmp(:,2));  end
            for j = 1 : IdvWhiteSpace;  AdjWhiteSpace{i,2} = [AdjWhiteSpace{i,2} ' '];  end
        end  
        OptNames = [p.In.PaOpt p.In.IcOpt p.Ex.PaOpt p.Ex.IcOpt];   % shorten names of parameters if they are too long
        for i = 1 : length(OptNames)
            if ( length(OptNames{i}) > 7 );  OptNames{i} = OptNames{i}(1:7);  end
        end

        fprintf('\n+++++++++++ Parameter estimation +++++++++++\n');
        fprintf('---initial---        ---estimated---\n');
        aws = 1;
        for i = 1 : size(p.In.PaOpt,2);  fprintf('  %.2e    %s\t%.2e\n',p.problem.x_0(aws),OptNames{aws},Opt.In.xbest(i));  aws = aws + 1; end
        for i = 1 : size(p.In.IcOpt,2);  fprintf('  %.2e    %s\t%.2e\n',p.problem.x_0(aws),OptNames{aws},Opt.In.xbest(size(p.In.PaOpt,2)+i));  aws = aws + 1; end
        fprintf('++++++++++++++++++++++++++++++++++++++++++++\n');
        for i = 1 : size(p.Ex.PaOpt,2);  fprintf('  %.2e    %s\t%.2e\n',p.problem.x_0(aws),OptNames{aws},Opt.Ex.xbest(i));  aws = aws + 1;  end        
        for i = 1 : size(p.Ex.IcOpt,2);  fprintf('  %.2e    %s\t%.2e\n',p.problem.x_0(aws),OptNames{aws},Opt.Ex.xbest(size(p.Ex.PaOpt,2)+i));  aws = aws + 1;  end  
        fprintf('++++++++++++++++++++++++++++++++++++++++++++\n\n');  

        % Simulate the optimized intracellular model
        p.In.Ic.STV(p.In.PosIcOpt) = Opt.In.xbest(p.In.nPaOpt+1:end);
        p.In.Ic.COI(p.In.PosIcOpt) = Opt.In.xbest(p.In.nPaOpt+1:end);

        p.xIn = Opt.In.xbest;

        if ( p.Ex.Moi <= 1 )
            p.In.Ic.STV(strcmp('Vex', p.In.StateNames)) = 1;
            p.In.Ic.COI(strcmp('Vex', p.In.StateNames)) = 1;
        else
            p.In.Ic.STV(strcmp('Vex', p.In.StateNames)) =  p.Ex.Moi;
            p.In.Ic.COI(strcmp('Vex', p.In.StateNames)) =  p.Ex.Moi;
        end
        if ( p.Ex.Modip <= 1 )
            p.In.Ic.COI(strcmp('Dex', p.In.StateNames)) = 1; %virions/cell  initial number of DIPs in the intracellular model  
        else
            p.In.Ic.COI(strcmp('Dex', p.In.StateNames)) = p.Ex.Modip;
        end      
        p.In.Ic.STV(strcmp('Dex', p.In.StateNames)) = 0;            

        clear mex
        h.In.PaBest = Opt.In.xbest(1:p.In.nPaOpt);
        if(p.CompileFlag)
            result.In.STV     = IQMPsimulate('InModel_MexFile', [0 1], p.In.Ic.STV, [p.In.PaOpt, 'DF'], [h.In.PaBest, 0]); 
            result.In.COI     = IQMPsimulate('InModel_MexFile', [0 1], p.In.Ic.COI, [p.In.PaOpt, 'DF'], [h.In.PaBest, 1]);             
        else          
            tmp = p.InModel;
            tmp = IQMparameters(tmp, [p.In.PaOpt, 'DF'], [h.In.PaBest, 0]);            
            result.In.STV     = IQMsimulate(tmp,'ode23s',[0 1], p.In.Ic.STV); 
            
            tmp = p.InModel;
            tmp = IQMparameters(tmp, [p.In.PaOpt, 'DF'], [h.In.PaBest, 1]);                
            result.In.COI     = IQMsimulate(tmp,'ode23s',[0 1], p.In.Ic.COI);    
        end

        % divide variables and reactions for CoinfectionModel (as the reactions are also declared as variables)
        for i = 1:length(result.In.STV.variables)
            if ( strcmp('r', result.In.STV.variables{i}(1)) )
                p.In.VarReaBreak = i;
                break;
            end
        end
        result.In.STV.reactions      = result.In.STV.variables(i:end);
        result.In.STV.reactionvalues = result.In.STV.variablevalues(:,i:end);
        result.In.STV.variables      = result.In.STV.variables(1:i-1);
        result.In.STV.variablevalues = result.In.STV.variablevalues(:,1:i-1);

        result.In.COI.reactions      = result.In.COI.variables(i:end);
        result.In.COI.reactionvalues = result.In.COI.variablevalues(:,i:end);
        result.In.COI.variables      = result.In.COI.variables(1:i-1);
        result.In.COI.variablevalues = result.In.COI.variablevalues(:,1:i-1);       

        %% Simulate the optimized extracellular model
        p.SimCond = 2 : size(p.Ex.MOIMODIPConditions);
        
        result.Ex.statevalues      = zeros(length(p.Ex.tspan), length(result.Ex.states), size(p.Ex.MOIMODIPConditions,1));
        result.Ex.variablevalues   = zeros(length(p.Ex.tspan), length(result.Ex.variables), size(p.Ex.MOIMODIPConditions,1));
        
        p.Ex.PopTick               = p.SimTime/p.Ex.h;        
        p.InToPop.states           = p.InToPop.RecordedStates;
        p.InToPop.populationvalues = zeros(p.Ex.PopTick+1,length(p.InToPop.states), size(p.Ex.MOIMODIPConditions,1));  
        
        for i = 1 : length(p.InToPop.states)
            p.InToPop.Idx.(p.InToPop.states{i}) = strcmp(p.InToPop.states{i},p.InToPop.states);
        end        

        fprintf('\nSimulate infections with best estimation results ...\n\n');
        for SimRun = p.SimCond

            % update intracellular parameters (special case for kSynV, rest is applied automatically in _EffectiveMOI.m)
            p.In.ParaValues(strcmp('kSynV', p.In.ParaNames)) = Opt.In.xbest(strcmp('kSynV', p.In.PaOpt));

            % update extracellular IC and parameters (needed for LowMOI_PreSim)
            p.Ex.y0(p.Ex.PosIcOpt) = Opt.Ex.xbest(p.Ex.nPaOpt+1:end);
            for CountPara = 1:p.Ex.nPaOpt
                p.Ex.(char(p.Ex.PaOpt(CountPara))) = Opt.Ex.xbest(CountPara);
            end
            if(sum(strcmp('kFus', p.In.PaOpt)));   p.Ex.kFus  = Opt.In.xbest(p.In.IdxKfus);   end      
            if(sum(strcmp('Fm', p.In.PaOpt)));     p.Ex.Fm    = Opt.In.xbest(p.In.IdxFm);   end
            if(sum(strcmp('kSynM', p.In.PaOpt)));  p.Ex.kSynM = Opt.In.xbest(p.In.IdxKsynM);   end

            % define necessary intracellular simulation time to capture all measurement time points
            DataTimeLim_In = find(~isnan(d.In.ExpStateValues(:,1,SimRun)),1,'last'); % last time points measured in individual experiment  
            p.In.SimTime = d.Ex.Time(DataTimeLim_In) + p.In.h; % increase length by one tick to capture full dynamics of initially infected cells
            p.In.tspan = 0 : p.In.h : p.In.SimTime;                            
            
            p.Ex.Moi   = p.Ex.MOIMODIPConditions(SimRun,1);
            p.Ex.Modip = p.Ex.MOIMODIPConditions(SimRun,2);                              
            
            % define MODIP/MOI-dependent kSynV
            if ( p.Ex.Modip > p.In.kSynV_Threshold )
                p.In.kSynV = p.In.ParaValues(strcmp(p.In.ParaNames, 'kSynV')) / (p.Ex.V1*(p.Ex.Modip/p.Ex.Moi).^p.Ex.V2);
            else
                p.In.kSynV = p.In.ParaValues(strcmp(p.In.ParaNames, 'kSynV'));
            end

            % calculate low MOI infection dynamics to save computation time
            % output: p.In.LowMOI.STV/COI.RrelV/D/P, p.In.LowMOI.AllValuesSTV/COI
            p = LowMOI_PreSim(p);

            % initial conditions extracellular level
            p.Ex.y0 = zeros(1,length(result.Ex.states));
            
            p.Ex.y0(p.Idx.T)    = d.Ex.ExpStateValues(1, strcmp('T',d.Ex.ExpStateNames), SimRun);
            p.Ex.y0(p.Idx.Ta)   = d.Ex.ExpStateValues(1, strcmp('Ta',d.Ex.ExpStateNames), SimRun);
            p.Ex.y0(p.Idx.Ia)   = d.Ex.ExpStateValues(1, strcmp('Ia',d.Ex.ExpStateNames), SimRun);

            p.Ex.y0(p.Idx.V)    = p.Ex.Moi   * d.Ex.ExpStateValues(1, strcmp('Cviab',d.Ex.ExpStateNames), SimRun);
            p.Ex.y0(p.Idx.D)    = p.Ex.Modip * d.Ex.ExpStateValues(1, strcmp('Cviab',d.Ex.ExpStateNames), SimRun);
            p.Ex.y0(p.Idx.Par)  = 0;           
            
            % differentiate type of infected cell based on Poisson distribution of the initial MOI/MODIP conditions
            if ( SimRun < 4 )
                % set infected cells for very low MOI+MODIP conditions to 0, otherwise intracellular dynamics start too early
                p.Ex.IV_0          = 0; 
                p.Ex.ICO_0         = 0;
                p.Ex.y0(p.Idx.I_D) = 0;
            else
                tmp = poisspdf(0,[p.Ex.Moi, p.Ex.Modip]); % get Poisson values for MOI/MODIP
                
                InfPercentages   = [(1-tmp(1)) * tmp(2), ...     % I_V
                                    tmp(1)     * (1-tmp(2)), ... % I_D
                                    (1-tmp(1)) * (1-tmp(2))];    % I_CO
                
                p.Ex.IV_0  = InfPercentages(1)/sum(InfPercentages) * d.Ex.ExpStateValues(1, strcmp('I',d.Ex.ExpStateNames), SimRun);
                p.Ex.ICO_0 = InfPercentages(3)/sum(InfPercentages) * d.Ex.ExpStateValues(1, strcmp('I',d.Ex.ExpStateNames), SimRun);
                p.Ex.y0(p.Idx.I_D) = InfPercentages(2)/sum(InfPercentages) * d.Ex.ExpStateValues(1, strcmp('I',d.Ex.ExpStateNames), SimRun);  
                
                % exception for MOI 0 + MODIP 0
                if ( sum(InfPercentages) == 0 )
                    p.Ex.IV_0          = 0;
                    p.Ex.ICO_0         = 0;
                    p.Ex.y0(p.Idx.I_D) = 0;
                end                
            end            

            %%% Calculate the integral of kApoI(tau), exp(-kApoI(tau)) 
            % modify apoptosis rates and integrals
            [p.Ex.ApoI, p.Ex.IntKapoi] = getApoI(p);            
            
            % effective MOI-dependent release rates for the intracellular model
            p.DynMOI.STV.RrelV    = cell(length(p.Ex.tspan),1);
            p.DynMOI.STV.RrelD    = cell(length(p.Ex.tspan),1);
            p.DynMOI.STV.RrelDTot = cell(length(p.Ex.tspan),1);
            p.DynMOI.STV.RrelP    = cell(length(p.Ex.tspan),1);
            p.DynMOI.COI.RrelV    = cell(length(p.Ex.tspan),1);
            p.DynMOI.COI.RrelD    = cell(length(p.Ex.tspan),1);
            p.DynMOI.COI.RrelDTot = cell(length(p.Ex.tspan),1);
            p.DynMOI.COI.RrelP    = cell(length(p.Ex.tspan),1); 

            % prepare saving of intracellular dynamics converted to population level      
            p.InToPop.states        = p.InToPop.RecordedStates; 
            p.InToPop.StateVarNames = [p.In.StateNames; p.In.VariableNames];
            p.Ex.PopTick            = p.SimTime/p.Ex.h;
            p.InToPop.ExpTime       = 0:p.SimTime/p.Ex.PopTick:p.SimTime;
            for i = 1 : length(p.InToPop.states)
                p.InToPop.STV.statevalues{i} = zeros(length(p.Ex.tspan),p.Ex.PopTick+1);
                p.InToPop.COI.statevalues{i} = zeros(length(p.Ex.tspan),p.Ex.PopTick+1); 
            end                     
            
            try
                [p, result.Ex.time, result.Ex.statevalues(:, :, SimRun), result.Ex.variablevalues(:, p.Idx.I_V, SimRun),...
                  result.Ex.variablevalues(:, p.Idx.I_V_SIE, SimRun), result.Ex.variablevalues(:, p.Idx.I_CO, SimRun),...
                  result.Ex.variablevalues(:, p.Idx.I, SimRun), FullIC, p.eFPar] = ExLev_Euler(p);    
            catch
                fprintf('Error during simulation of initial condition #%i\n',SimRun);
                continue;
            end

            result.Ex.variablevalues(:, strcmp('Log10Tcid50', result.Ex.variables), SimRun) =...
                log10(max(1, result.Ex.statevalues(:, p.Idx.V, SimRun)));   
            result.Ex.variablevalues(:, strcmp('Log10Ha', result.Ex.variables), SimRun) =...
                log10(max(1, result.Ex.statevalues(:, p.Idx.Par, SimRun)));                            
           
            % calculate values at intracellular measurement time points    
            SumCells = result.Ex.statevalues(:,strcmp('Cviab',result.Ex.states),SimRun); 
                     
            for i = 1 : p.Ex.PopTick + 1                 
                for j = 1 : length(p.InToPop.states)
                    
                    min_tmp1 = min(size(p.Ex.FullDistribution{i},1),size(p.InToPop.STV.statevalues{j},1));
                    tmp = p.Ex.FullDistribution{i}(1:min_tmp1, 1) .* ...
                          p.InToPop.STV.statevalues{j}(1:min_tmp1,i);

                    min_tmp2 = min(size(p.Ex.FullDistribution{i},1),size(p.InToPop.COI.statevalues{j},1));                   
                    if ( min_tmp1 ~= min_tmp2 )
                        if ( min_tmp1 < min_tmp2 )
                            tmp = [tmp; zeros(min_tmp2-min_tmp1,1)];
                        else
                            error('Error during simulation');
                        end
                    end
                    tmp = tmp + p.Ex.FullDistribution{i}(1:min_tmp2, 2) .* ...
                          p.InToPop.COI.statevalues{j}(1:min_tmp2,i);
                      
                    p.InToPop.populationvalues(i, j, SimRun) = sum(tmp) ./ SumCells(i);
                end   
            end
            p.InToPop.populationvalues(1, :, SimRun) = zeros(1,size(p.InToPop.populationvalues,2));            
             
            % add DIP-infected cell mRNA to others    
            if ( p.mRNAPrimTrans )
                p.InToPop.populationvalues(:, strcmp('Rm5', p.InToPop.states), SimRun) = ...
                  p.InToPop.populationvalues(:, strcmp('Rm5', p.InToPop.states), SimRun) + ...
                  result.Ex.statevalues(:, strcmp('Rm5_ID', result.Ex.states), SimRun) .* ... 
                  result.Ex.statevalues(:,strcmp('I_D',result.Ex.states), SimRun) ./ SumCells;              

                p.InToPop.populationvalues(:, strcmp('Rm9', p.InToPop.states), SimRun) = ...
                  p.InToPop.populationvalues(:, strcmp('Rm9', p.InToPop.states), SimRun) + ...
                  result.Ex.statevalues(:, strcmp('Rm9_ID', result.Ex.states), SimRun) .* ... 
                  result.Ex.statevalues(:,strcmp('I_D',result.Ex.states), SimRun) ./ SumCells;   
            end
          
        end
        result.Opt = Opt;

        VisualizeResults(result, d, p);

        %% Check if parameters got close to lower or upper bounds
        tmp = [Opt.In.xbest Opt.Ex.xbest];
        BoundLimit = zeros(length(tmp),2);
        for i = 1 : length(tmp)
            if ( tmp(i) < 1.1 * p.problem.x_L(i) );      BoundLimit(i,1) = 1;
            elseif ( tmp(i) > 0.9 * p.problem.x_U(i) );  BoundLimit(i,2) = 1;
            end    
        end

        fprintf('\n++++++++++ Estimation bound limitation ++++++++++');
        for i = 1 : length(tmp)
            if ( sum(sum(BoundLimit)) == 0)
                fprintf('\n\nnone');
                break;
            end

            if ( BoundLimit(i,1) )
                fprintf('\n\n%s\t is at lower bound: %.2e <- %.2e',OptNames{i}, p.problem.x_L(i), tmp(i));
            elseif ( BoundLimit(i,2) )
                fprintf('\n\n%s\t is at upper bound: %.2e -> %.2e',OptNames{i}, tmp(i), p.problem.x_U(i));
            end
        end
        fprintf('\n\n+++++++++++++++++++++++++++++++++++++++++++++++++\n\n');

        %% Display estimation information
        p.Info.EndTime = clock;
        fprintf('\n++++++++++ Additional information ++++++++++\n');
        fprintf('Start:      %s\n',datestr(p.Info.StartTime));
        fprintf('End:        %s\n',datestr(p.Info.EndTime));
        if ( datenum(p.Info.EndTime - p.Info.StartTime) < 0.0417 )
            fprintf('Duration:   %i min\n',ceil(1440*datenum(p.Info.EndTime - p.Info.StartTime)));
        else
            tmp = 24*datenum(p.Info.EndTime - p.Info.StartTime);
            fprintf('Duration:   %i h + %i min\n',floor(tmp),ceil(60*(tmp-floor(tmp))));
        end
        try
            fprintf('FunEvals:   %i\n',result.Opt.fSSm.numeval);
        catch ME
        end
        fprintf('ObjFunVal:  %i\n',result.Opt.ObjFunVal/1e6);
        fprintf('++++++++++++++++++++++++++++++++++++++++++++\n\n');        
        
otherwise
    disp('Warnin(Main): No such task.')                

end

%% Save final estimation results
if ( strcmp(p.Task,'optimize') && p.MaxTime > 300)
    SaveDir = dir([pwd '/results']);    % target directory
    tmp = 0;
    SaveNum = 1;                        % result number

    % find unused filename
    while 1

        for k = 1 : size(SaveDir,1)
            if ( strcmp(sprintf('PE_%s_result_%i.mat',datestr(date(),29,'local'),SaveNum),...
                        SaveDir(k).name) )
                SaveNum = SaveNum + 1;
                break;
            elseif ( k == size(SaveDir,1) )
                tmp = 1;
            end
        end

        if ( tmp );  break;  end
    end

    save(sprintf('results/PE_%s_result_%i.mat',datestr(date(),29,'local'),SaveNum), 'd', 'Opt', 'p', 'result');
end    
    
%% end of complete simulation run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% simulation time
TEND = toc(TSTART);
fprintf('\nSimulation took %.1f minutes.\n\n',TEND/60);

%% Housekeeping commands II
clear Count* h
path(SearchPath)
return
































