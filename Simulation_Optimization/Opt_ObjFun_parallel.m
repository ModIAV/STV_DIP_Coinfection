function ObjFunVal = Opt_ObjFun_parallel(x,varargin)
%OBJFUN determines the deviation of measurements and simulation by using W 
% as the weighting matrix and returns a scalar for optimization routines

% last revised: 2021/03/23

global p
global d
global stopOptimization % stop variable

clear mex

TSTART = tic;

NumOpt = length(p.OptCond);

% convert x into a row vector if CMEAS optimizer is used
if ( size(x,1) > 1 && size(x,2) > 1 )
    x = x(:,1)';
elseif ( size(x,1) > 1 )
    x = x';
else
%     error('ObjFun input has wrong format!')
end

p.ParaOV(end+1,:) = x;
p.SuccessVector(end+1) = 0;

% transform parameters if necessary
if ( p.LogPara )
    p.xIn = 10.^(x(1:p.In.nOpt));
    p.xEx = 10.^(x(p.In.nOpt+1:end));
else
    p.xIn = x(1:p.In.nOpt);
    p.xEx = x(p.In.nOpt+1:end);
end

%% Intracellular model new
%adopt new initial conditions
p.In.Ic.STV(p.In.PosIcOpt) = p.xIn(p.In.nPaOpt+1:p.In.nOpt);
p.In.Ic.COI(p.In.PosIcOpt) = p.xIn(p.In.nPaOpt+1:p.In.nOpt);

%% variables predefined for parfor loop
DataTimeLim_In = zeros(length(p.OptCond),1);
DataTimeLim_Ex = zeros(length(p.OptCond),1);
for i = 1 : NumOpt
    g(i) = p; % replacement structure for parallel simulation
end

% Prepare simulation
for i = 1 : NumOpt

    SimRun = p.OptCond(i);

    % define necessary simulation time to capture all measurement time points
    DataTimeLim_In(i) = find(~isnan(d.In.ExpStateValues(:,1,SimRun)),1,'last'); % last time points measured in individual experiment
    DataTimeLim_Ex(i) = find(~isnan(d.Ex.ExpStateValues(:,1,SimRun)),1,'last'); % last time points measured in individual experiment

    g(i).SimTime    = d.Ex.Time(DataTimeLim_Ex(i));    
    g(i).In.SimTime = d.Ex.Time(DataTimeLim_In(i)) + g(i).In.h; % increase length by one tick to capture full dynamics of initially infected cells
    
    g(i).Ex.tspan = 0 : g(i).Ex.h : g(i).SimTime;    
    g(i).In.tspan = 0 : g(i).In.h : g(i).In.SimTime;      
    
    % initial conditions extracellular level
    g(i).Ex.Moi   = g(i).Ex.MOIMODIPConditions(SimRun,1);
    g(i).Ex.Modip = g(i).Ex.MOIMODIPConditions(SimRun,2);  

    g(i).Ex.y0 = zeros(1,length(g(i).Ex.StateNames));
    
    g(i).Ex.y0(g(i).Idx.T)    = d.Ex.ExpStateValues(1, strcmp('T',d.Ex.ExpStateNames), SimRun);
    g(i).Ex.y0(g(i).Idx.Ta)   = d.Ex.ExpStateValues(1, strcmp('Ta',d.Ex.ExpStateNames), SimRun);
    g(i).Ex.y0(g(i).Idx.Ia)   = d.Ex.ExpStateValues(1, strcmp('Ia',d.Ex.ExpStateNames), SimRun);

    g(i).Ex.y0(g(i).Idx.V)    = g(i).Ex.Moi   * d.Ex.ExpStateValues(1, strcmp('Cviab',d.Ex.ExpStateNames), SimRun);
    g(i).Ex.y0(g(i).Idx.D)    = g(i).Ex.Modip * d.Ex.ExpStateValues(1, strcmp('Cviab',d.Ex.ExpStateNames), SimRun);

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
    
    %% Extracellular model
    %adopt new initial conditions
    g(i).Ex.y0(g(i).Ex.PosIcOpt) = g(i).xEx(g(i).Ex.nPaOpt+1:g(i).Ex.nOpt);

    %adopt intracellular parameter changes (special case for kSynV, rest is applied automatically in _EffectiveMOI.m)
    g(i).In.ParaValues(strcmp('kSynV', g(i).In.ParaNames)) = g(i).xIn(strcmp('kSynV', g(i).In.PaOpt));
    
    %adopt extracellular parameter changes
    for CountPara = 1:g(i).Ex.nPaOpt
        g(i).Ex.(char(g(i).Ex.PaOpt(CountPara))) = g(i).xEx(CountPara);
    end
    
    %adopt parameter changes on both levels
    if(sum(strcmp('kFus',  g(i).In.PaOpt)));   g(i).Ex.kFus   = g(i).xIn(g(i).In.IdxKfus);  end  
    if(sum(strcmp('Fm',    g(i).In.PaOpt)));   g(i).Ex.Fm     = g(i).xIn(g(i).In.IdxFm);  end
    if(sum(strcmp('kSynM', g(i).In.PaOpt)));   g(i).Ex.kSynM  = g(i).xIn(g(i).In.IdxKsynM);  end
    
    if ( g(i).Ex.kApoT > g(i).Ex.kApoI ) % baseline apoptosis rate is smaller than infected cell rate     
        p.TotalOptTime  = p.TotalOptTime + toc(TSTART);
        ObjFunVal = 1e8;
        ObjFunVal = NaN;
        return
    end    

    %%% Calculate the integral of kApoI(tau), exp(-kApoI(tau)) 
    % modify apoptosis rates and integrals
    [g(i).Ex.ApoI, g(i).Ex.IntKapoi] = getApoI(g(i));  

    % effective MOI-dependent release rates for the intracellular model
    g(i).DynMOI.STV.RrelV    = cell(length(g(i).Ex.tspan),1);
    g(i).DynMOI.STV.RrelD    = cell(length(g(i).Ex.tspan),1);
    g(i).DynMOI.COI.RrelDTot = cell(length(g(i).Ex.tspan),1); 
    g(i).DynMOI.STV.RrelP    = cell(length(g(i).Ex.tspan),1);
    g(i).DynMOI.COI.RrelV    = cell(length(g(i).Ex.tspan),1);
    g(i).DynMOI.COI.RrelD    = cell(length(g(i).Ex.tspan),1);
    g(i).DynMOI.COI.RrelDTot = cell(length(g(i).Ex.tspan),1); 
    g(i).DynMOI.COI.RrelP    = cell(length(g(i).Ex.tspan),1);  
    
    % define MODIP/MOI-dependent kSynV
    if ( g(i).Ex.Modip > g(i).In.kSynV_Threshold ) 
        g(i).In.kSynV = g(i).In.ParaValues(strcmp(g(i).In.ParaNames, 'kSynV')) / (g(i).Ex.V1*(g(i).Ex.Modip/g(i).Ex.Moi).^g(i).Ex.V2);
    else
        g(i).In.kSynV = g(i).In.ParaValues(strcmp(g(i).In.ParaNames, 'kSynV'));
    end
    
    % calculate low MOI infection dynamics to save computation time
    % output: p.In.LowMOI.STV/COI.RrelV/D/P, p.In.LowMOI.AllValuesSTV/COI  
    try
        g(i) = LowMOI_PreSim(g(i));    
    catch   
        p.FailedRuns(2) = p.FailedRuns(2) + 1;
        p.TotalOptTime  = p.TotalOptTime + toc(TSTART);
        ObjFunVal = NaN;
        return        
    end
  
    % prepare saving of intracellular dynamics converted to population level     
    g(i).InToPop.states        = g(i).In.States2Fit;
    g(i).InToPop.StateVarNames = [g(i).In.StateNames; g(i).In.VariableNames];
    g(i).InToPop.ExpTime       = d.In.Time;
    for j = 1 : length(g(i).InToPop.states)
        g(i).InToPop.STV.statevalues{j} = zeros(length(g(i).In.tspan),length(d.In.Time));
        g(i).InToPop.COI.statevalues{j} = zeros(length(g(i).In.tspan),length(d.In.Time)); 
    end  

end

%% parallel calculation of infection dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h         = cell(1,2);
resEx_tmp = [];

ErrorFlag = zeros(NumOpt,1);
parfor i = 1 : NumOpt  
    try       
        
        % Simulate the extracellular model
        [h{i}, resEx_tmp(i).time, resEx_tmp(i).statevalues, resEx_tmp(i).variablevalues(:, g(i).Idx.I_V),...
          resEx_tmp(i).variablevalues(:, g(i).Idx.I_V_SIE), resEx_tmp(i).variablevalues(:, g(i).Idx.I_CO),...
          resEx_tmp(i).variablevalues(:, g(i).Idx.I), ~, ~] = ExLev_Euler(g(i)); 
    catch 
        ErrorFlag(i,1) = 1;
    end
end

% catch simulation errors
if ( sum(ErrorFlag) ) 
    p.FailedRuns(3) = p.FailedRuns(3) + 1;
    p.TotalOptTime  = p.TotalOptTime + toc(TSTART); 
    ObjFunVal = NaN; 
    return
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% evaluate simulation results
TmpFunVal = zeros(NumOpt,1);
for i = 1 : NumOpt

    p = h{i};
    resEx = resEx_tmp(i);
    
    SimRun = p.OptCond(i);
    
    %% calculate values at intracellular measurement time points    
    SumCells = resEx.statevalues(:,strcmp('Cviab',p.Ex.StateNames));
             
    p.InToPop.populationvalues = zeros(length(d.In.Time),length(p.InToPop.states)); 
    for j = 1 : length(d.In.Time)
        if ( d.In.Time(j) > p.Ex.tspan(end) );  break;  end % cancel when maximum sim time for specific condition was reached

        for k = 1 : length(p.InToPop.states) 
            MaxTime = min(round(d.In.Time(j)/p.In.h)+1, size(p.InToPop.STV.statevalues{k},1));   
            STVstate = zeros(length(p.Ex.tspan_save),1);     
            STVstate(1:MaxTime) = p.InToPop.STV.statevalues{k}(1:MaxTime,j);            
            tmp = p.Ex.FullDistribution{round(d.In.Time(j)/p.In.h)+1}(:,1) .* ...
                  STVstate(1:round(d.In.Time(j)/p.In.h)+1);

            MaxTime = min(round(d.In.Time(j)/p.In.h)+1, size(p.InToPop.COI.statevalues{k},1));   
            COIstate = zeros(length(p.Ex.tspan),1);
            COIstate(1:MaxTime) = p.InToPop.COI.statevalues{k}(1:MaxTime,j);            
            tmp = tmp + p.Ex.FullDistribution{round(d.In.Time(j)/p.In.h)+1}(:,2) .* ...
                  COIstate(1:round(d.In.Time(j)/p.In.h)+1);            

            p.InToPop.populationvalues(j,k) = sum(tmp) ./ SumCells(round(d.In.Time(j)/p.In.h)+1);
        end           
    end    
    
    % calculate percentage of vRNP in nucleus
    if ( sum(strcmp('vRNPnuc_DAPI', p.In.States2Fit)) )
        Idx_tmp = [ find(p.InToPop.Idx.vRNPnuc1 + p.InToPop.Idx.vRNPnuc5 + p.InToPop.Idx.vRNPnuc9);
                    find(p.InToPop.Idx.RvSeg1 + p.InToPop.Idx.RvSeg5 + p.InToPop.Idx.RvSeg9) ];

        p.InToPop.populationvalues(:, p.InToPop.Idx.vRNPnuc_DAPI) = ...
          100 * ( sum(p.InToPop.populationvalues(:, Idx_tmp(1,:)), 2) ) ...
             ./ ( sum(p.InToPop.populationvalues(:, Idx_tmp(2,:)), 2) );       

        p.InToPop.populationvalues(1, p.InToPop.Idx.vRNPnuc_DAPI) = d.In.ExpStateValues(1, p.InToPop.Idx.vRNPnuc_DAPI, SimRun); % disregard difference in first time point -> pseudo offset      
    end
    
    % only take fitted states for further analysis
    SimValuesIn = p.InToPop.populationvalues(:,1:length(p.In.States2Fit));
    
    %%% add DIP-infected cell mRNA to others
    if ( p.mRNAPrimTrans )
        [~, tmp] = find(d.In.Time==resEx.time);

        SimValuesIn(1:length(tmp), strcmp('Rm5', p.In.States2Fit)) = ...
          SimValuesIn(1:length(tmp), strcmp('Rm5', p.In.States2Fit)) + ...
          resEx.statevalues(tmp, strcmp('Rm5_ID', p.Ex.StateNames)) .* ...
          resEx.statevalues(tmp, strcmp('I_D', p.Ex.StateNames)) ./ SumCells(tmp); 

        SimValuesIn(1:length(tmp), strcmp('Rm9', p.In.States2Fit)) = ...
          SimValuesIn(1:length(tmp), strcmp('Rm9', p.In.States2Fit)) + ...
          resEx.statevalues(tmp, strcmp('Rm9_ID', p.Ex.StateNames)) .* ...
          resEx.statevalues(tmp, strcmp('I_D', p.Ex.StateNames)) ./ SumCells(tmp);  
    end
  
    %%% OFFSET - intracellular level
    % To account for free viral RNAs in the seed virus which might adhere to
    % cells and cause a constant offset in measurements, we add the 0 hpi
    % simulation value of each RNA species to the subsequent time points    
    if ( ~isempty(d.In.ExpStateNames) && ~isempty(p.In.States2Fit) )        
        
        AddOffset = d.In.ExpStateValues(1, p.In.OffsetIndex(:,2), p.OptCond(i));
        if ( sum(isnan(AddOffset)) > 0 )
            if ( sum( sum(isnan(d.In.ExpStateValues(:, p.In.OffsetIndex(:,2), p.OptCond(i)))) == size(d.In.ExpStateValues,1) ) == 0 )
                AddOffset(isnan(AddOffset)) = 0;
            end
        end   

        SimValuesIn(:,p.In.OffsetIndex(:,1)) = SimValuesIn(:,p.In.OffsetIndex(:,1)) +...
            repmat(AddOffset, length(d.In.Time), 1);         
    end    

    %% calculate objective function values

    if ( ~isempty(d.In.ExpStateValues) )

        MeasuredValuesIn = d.In.ExpStateValues(:,p.In.PosState2Fit(:,2),p.OptCond(i));

        % To yield an HA titer cells need to produce more than 2e7 virions/ml,
        % which means that with a cell concentration of roughly 1.2e6/ml each cell
        % must at least produce 16.7 virions to yield a measureable HA titer
        a = SimValuesIn(:,strcmp('Vrel', p.In.States2Fit));
        a(a<16.7) = 0;
        SimValuesIn(:,strcmp('Vrel', p.In.States2Fit)) = a;

        % Define weighting matrix for extracellular model + logarithmization
        p.In.W = ones(length(d.In.Time),length(p.In.States2Fit));
        switch lower(p.In.ObjFunCal)
            case {'normal'}
                p.In.W = 1./repmat(max(d.In.ExpStateValues(:, p.In.PosState2Fit(:,2),SimRun)), length(d.In.Time), 1); %normalize to maximum of measurement

            case {'log'}        
                p.In.W = 1./log10(repmat(max(d.In.ExpStateValues(:, p.In.PosState2Fit(:,2),SimRun)), length(d.In.Time), 1));
                
                SimValuesIn      = log10(SimValuesIn);
                MeasuredValuesIn = log10(MeasuredValuesIn);
                SimValuesIn(SimValuesIn==-Inf)           = 0;
                MeasuredValuesIn(MeasuredValuesIn==-Inf) = 0;               
            otherwise
                error('Please define a supported option for data normalization!');                
        end
        p.In.W(isnan(d.In.ExpStateValues(:, p.In.PosState2Fit(:,2),SimRun))) = 0;
      
        DevStatesIn = SimValuesIn-MeasuredValuesIn;
        DevStatesIn(isnan(DevStatesIn)) = 0;            
        
        ObjFunValIn = sum(sum((p.In.W .* DevStatesIn).^2).*p.In.RelStateWeight);
    else
        ObjFunValIn = 0;
    end

    %% Extracellular model
    %calculate the deviation between simulation and measurements
    SimulationEx     = [resEx.statevalues resEx.variablevalues];

    SimValuesEx      = SimulationEx(p.Ex.TimeIndexData(1:DataTimeLim_Ex(i)), p.Ex.PosState2Fit(:,1));
    MeasuredValuesEx =           d.Ex.ExpStateValues(1:DataTimeLim_Ex(i), p.Ex.PosState2Fit(:,2),SimRun);
    
    %%% OFFSET - extracellular level
    % To account for free viral RNAs in the seed virus which might adhere to
    % cells and cause a constant offset in measurements, we add the 0 hpi
    % simulation value of each RNA species to the subsequent time points      
    AddOffset = d.Ex.ExpStateValues(2, p.Ex.OffsetIndex(:,2), SimRun);
    AddOffset(isnan(AddOffset)) = 0;
    AddOffset = AddOffset .* ~isnan(d.Ex.ExpStateValues(3, p.Ex.OffsetIndex(:,2), SimRun));    
    
    SimValuesEx(:,p.Ex.OffsetIndex(:,1)) = SimValuesEx(:,p.Ex.OffsetIndex(:,1)) +...
        repmat(AddOffset, size(SimValuesEx,1), 1);      
    
    % Define weighting matrix for extracellular model + logarithmization
    p.Ex.W = ones(length(d.Ex.Time),length(p.Ex.States2Fit));
    switch lower(p.Ex.ObjFunCal)
        case {'normal'}
            p.Ex.W = 1./repmat(max(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(:,2),SimRun)), length(d.Ex.Time), 1); %normalize to maximum of measurement

        case {'log'}        
            p.Ex.W = 1./log10(repmat(max(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(:,2),SimRun)), length(d.Ex.Time), 1));
             
            SimValuesEx      = log10(SimValuesEx);
            MeasuredValuesEx = log10(MeasuredValuesEx);
            MeasuredValuesEx(MeasuredValuesEx==-Inf) = 0;
            SimValuesEx(SimValuesEx==-Inf)           = 0; 
            
        case {'log/normal'}  

            CellPop = {'T', 'Ta', 'I', 'Ia'};
            CellPopPos = [];
            for kk = 1 : length(p.Ex.States2Fit)
                if ( sum(strcmp(p.Ex.States2Fit(kk),CellPop)) )
                    CellPopPos(end+1,1) = kk;
                end
            end
            NoCellPopPos = setdiff(1:length(p.Ex.States2Fit), CellPopPos);
            
            p.Ex.W_normal = 1./repmat(max(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(CellPopPos,2), SimRun)), ...
                                      length(d.Ex.Time), 1);

            p.Ex.W_log    = 1./log10(repmat(max(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(NoCellPopPos,2), SimRun)), ...
                                            length(d.Ex.Time), 1));
                                     
            p.Ex.W_normal(isnan(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(CellPopPos,2), SimRun))) = 0;
            p.Ex.W_log(isnan(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(NoCellPopPos,2), SimRun))) = 0;                      
             
            SimValuesEx_normal      = SimValuesEx(:,CellPopPos);
            SimValuesEx_log         = log10(SimValuesEx(:,NoCellPopPos));
            
            MeasuredValuesEx_normal = MeasuredValuesEx(:,CellPopPos);
            MeasuredValuesEx_log    = log10(MeasuredValuesEx(:,NoCellPopPos));
            
            MeasuredValuesEx_log(MeasuredValuesEx_log==-Inf) = 0;
            SimValuesEx_log(SimValuesEx_log==-Inf)           = 0 ; 

        otherwise
            error('Please define a supported option for data normalization!');                
    end
    p.Ex.W(isnan(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(:,2),SimRun))) = 0;     

    % calculate difference
    DevStatesEx = SimValuesEx-MeasuredValuesEx;
    DevStatesEx(isnan(MeasuredValuesEx)) = 0;
      
  
    if ( strcmp(p.Ex.ObjFunCal,'log/normal') ) 
        % calculate difference - normal
        DevStatesEx = SimValuesEx_normal - MeasuredValuesEx_normal;
        DevStatesEx(isnan(MeasuredValuesEx_normal)) = 0;

        ObjFunValEx_normal = sum(sum((p.Ex.W_normal(1:DataTimeLim_Ex(i),:) .* DevStatesEx).^2).*p.Ex.RelStateWeight(CellPopPos));

        % calculate difference - log
        DevStatesEx = SimValuesEx_log - MeasuredValuesEx_log;
        DevStatesEx(isnan(MeasuredValuesEx_log)) = 0;

        ObjFunValEx = ObjFunValEx_normal + ...
                      sum(sum((p.Ex.W_log(1:DataTimeLim_Ex(i),:) .* DevStatesEx).^2).*p.Ex.RelStateWeight(NoCellPopPos));              

    else
        ObjFunValEx = sum(sum((p.Ex.W(1:DataTimeLim_Ex(i),:) .* DevStatesEx).^2).*p.Ex.RelStateWeight);
    end
   
    %% Combine both objective function values
    %c calculate total numel
    if ( ~isempty(d.In.ExpStateValues) )
        TmpFunVal(i,1) = TmpFunVal(i,1) + ObjFunValIn/(numel(MeasuredValuesIn)-sum(sum(isnan(MeasuredValuesIn))));
    end
    TmpFunVal(i,1) = TmpFunVal(i,1) + ObjFunValEx/(numel(MeasuredValuesEx)-sum(sum(isnan(MeasuredValuesEx)))); 

end

ObjFunVal = sum(TmpFunVal(~isnan(TmpFunVal)));

% save best result in case of crash (does not work during parallel computation)
if ( ObjFunVal < p.safety.Bestf )
    p.safety.Bestf = ObjFunVal;
    p.safety.xbest = x;
    
    % record parameter sets for video purposes 
    p.record(end+1,:) = x;
end

ObjFunVal = ObjFunVal*1e6;
% fprintf('\nNew ObjFun: %.2f, reached in %.1f s\n',ObjFunVal,toc(TSTART)); 

%% optimization time calculation + display
p.TotalOptTime = p.TotalOptTime + toc(TSTART);
if ( p.TotalOptTime > p.MaxTime )
   stopOptimization = 1; 
elseif ( p.TotalOptTime > p.TotalTimeCounter*3600 )
    fprintf('### Current simulation time: %i h ###\n', p.TotalTimeCounter);
    p.TotalTimeCounter = p.TotalTimeCounter + 1;
end

p.SuccessVector(end) = 1;


