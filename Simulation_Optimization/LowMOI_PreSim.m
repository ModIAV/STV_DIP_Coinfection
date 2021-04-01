function [p] = LowMOI_PreSim(p)
% LowMOI_PreSim calculates low MOI and MODIP infection dynamics that can be
% assumed as the result of all infection with low viral input
%
% last revised: 2021/03/23

if ( strcmp('simulate',p.Task) )
    p.xIn = zeros(1,p.In.nPaOpt);
    for i = 1 : p.In.nPaOpt
        p.xIn(i) = p.In.ParaValues(strcmp(p.In.PaOpt{i},p.In.ParaNames));
    end
end

%% define indices for easier access
Idx.Vex    = strcmp('Vex',p.In.StateNames);
Idx.VattHi = strcmp('VattHi',p.In.StateNames);
Idx.VattLo = strcmp('VattLo',p.In.StateNames);
Idx.Ven    = strcmp('Ven',p.In.StateNames);
Idx.Vfus   = strcmp('Vfus',p.In.StateNames);
Idx.Vcyt   = strcmp('Vcyt',p.In.StateNames);
Idx.Dex    = strcmp('Dex',p.In.StateNames);
Idx.DattHi = strcmp('DattHi',p.In.StateNames);
Idx.DattLo = strcmp('DattLo',p.In.StateNames);
Idx.Den    = strcmp('Den',p.In.StateNames);
Idx.Dfus   = strcmp('Dfus',p.In.StateNames);
Idx.Dcyt   = strcmp('Dcyt',p.In.StateNames);
Idx.FPar   = strcmp('FPar',p.In.StateNames);

%% pre-simulation            
p.In.LowMOI.STV.RrelV = zeros(1, p.In.SimTime/p.In.h+1);
p.In.LowMOI.STV.RrelD = zeros(1, p.In.SimTime/p.In.h+1);
p.In.LowMOI.STV.RrelDTot = zeros(1, p.In.SimTime/p.In.h+1);
p.In.LowMOI.STV.RrelP = zeros(1, p.In.SimTime/p.In.h+1);
p.In.LowMOI.COI.RrelV = zeros(1, p.In.SimTime/p.In.h+1);
p.In.LowMOI.COI.RrelD = zeros(1, p.In.SimTime/p.In.h+1);
p.In.LowMOI.COI.RrelDTot = zeros(1, p.In.SimTime/p.In.h+1);
p.In.LowMOI.COI.RrelP = zeros(1, p.In.SimTime/p.In.h+1);

p.In.LowMOI.AllValuesSTV = {}; % cell(length(p.In.LowMOI.ICvalues),1);
p.In.LowMOI.AllValuesCOI = {}; % cell(length(p.In.LowMOI.ICvalues),length(p.In.LowMOI.ICvalues));

SimTimeVec = 0:p.Ex.h:p.In.SimTime; 

%%% STV-infected cells
% set initial conditions
IC = p.In.Ic.STV;

IC(strcmp('Vex',p.In.StateNames))  = 0;
IC(strcmp('Vcyt',p.In.StateNames)) = 1;
IC(strcmp('Dex',p.In.StateNames))  = 0;
IC(strcmp('Dcyt',p.In.StateNames)) = 0;

% simulate inntracellular model
if(p.CompileFlag)
    result.In.COI     = IQMPsimulate('InModel_MexFile', SimTimeVec, IC, ...
                                     [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 0]);
                                 
else
    tmp = p.InModel;
    tmp = IQMparameters(tmp, [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 0]);

    result.In.COI     = IQMsimulate(tmp, 'ode23s', SimTimeVec, IC);  
end            

% set all DIP related values to 0, values can reach 1e-12 due to numerical instabilities during SBPD simulate (?)
tmp = [21:9:84, 104, 106, 108]; % vector of intracellular states related to DIPs
result.In.COI.statevalues(:,tmp) = zeros(size(result.In.COI.statevalues,1), length(tmp));
tmp = [20, 29, 39:9:93, 107:9:170, 172, 200, 203, 205] + 19; % vector of intracellular states related to DIPs (+19 for offset with variables + reactions)
result.In.COI.reactionvalues(:,tmp) = zeros(size(result.In.COI.reactionvalues,1), length(tmp));                  

p.In.LowMOI.STV.RrelV(1,:) = result.In.COI.variablevalues(:,strcmp('rRel',result.In.COI.variables));
p.In.LowMOI.STV.RrelD(1,:) = 0 * result.In.COI.variablevalues(:,strcmp('rRelD',result.In.COI.variables));
p.In.LowMOI.STV.RrelDTot(1,:) = 0 * result.In.COI.variablevalues(:,strcmp('rRelDIPtot',result.In.COI.variables));
p.In.LowMOI.STV.RrelP(1,:) = result.In.COI.variablevalues(:,strcmp('rRelP',result.In.COI.variables));

p.In.LowMOI.AllValuesSTV = [result.In.COI.statevalues, result.In.COI.variablevalues]; % save state values for p.InToPop calculation
    
%%% COI-infected cells
% set initial conditions
IC = p.In.Ic.COI;

IC(strcmp('Vex',p.In.StateNames))  = 0;
IC(strcmp('Vcyt',p.In.StateNames)) = 1;
IC(strcmp('Dex',p.In.StateNames))  = 0;
IC(strcmp('Dcyt',p.In.StateNames)) = 1;

% simulate intracellular model
if(p.CompileFlag)  
    try
        result.In.COI     = IQMPsimulate('InModel_MexFile', SimTimeVec, IC, ...
                                         [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 1]);
    catch
        IC(Idx.Vex) = IC(Idx.Vex)*0.999;
        IC(Idx.Dex) = IC(Idx.Dex)*0.999;        
        result.In.COI     = IQMPsimulate('InModel_MexFile', SimTimeVec, IC, ...
                                         [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 1]);                                     
        p.Error_099(2) = p.Error_099(2) + 1;                                     
    end
else
    tmp = p.InModel;
    tmp = IQMparameters(tmp, [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 1]);
    
    result.In.COI     = IQMsimulate(tmp, 'ode23s', SimTimeVec, IC);                                
end       
    
p.In.LowMOI.COI.RrelV(1,:) = result.In.COI.variablevalues(:,strcmp('rRel',result.In.COI.variables));
p.In.LowMOI.COI.RrelD(1,:) = result.In.COI.variablevalues(:,strcmp('rRelD',result.In.COI.variables));
p.In.LowMOI.COI.RrelDTot(1,:) = result.In.COI.variablevalues(:,strcmp('rRelDIPtot',result.In.COI.variables));
p.In.LowMOI.COI.RrelP(1,:) = result.In.COI.variablevalues(:,strcmp('rRelP',result.In.COI.variables));

p.In.LowMOI.AllValuesCOI = [result.In.COI.statevalues, result.In.COI.variablevalues]; % save state values for p.InToPop calculation

%% allow virus particle release only after at least 1 full particle was formed
%  -> Vrel(Vrel < 1) = 0
AllReleases =  [p.In.LowMOI.STV.RrelV;...
                p.In.LowMOI.STV.RrelD;...
                p.In.LowMOI.STV.RrelDTot;...
                p.In.LowMOI.STV.RrelP;...
                p.In.LowMOI.COI.RrelV;...
                p.In.LowMOI.COI.RrelD;...
                p.In.LowMOI.COI.RrelDTot;...
                p.In.LowMOI.COI.RrelP];

% find time points at which cumualtive V/D/PRel > 1
CSum = cumsum(AllReleases,2);
FirstRelPos = zeros(8,1);

for i = 1 : size(CSum,1)
    if ( sum(CSum(i,:) > 1) )       
        FirstRelPos(i) = find(CSum(i,:)>1, 1, 'first');
    end
end

% handle special case when no STV infected cells exist
if ( isempty(p.In.LowMOI.STV.RrelV) )
    CSum = [zeros(size(CSum));CSum];    
    FirstRelPos = [zeros(4,1); FirstRelPos];
end

% adjust release vectors
if ( FirstRelPos(1) > 0 )
    p.In.LowMOI.STV.RrelV(1:FirstRelPos(1)-1)    = 0;
    p.In.LowMOI.STV.RrelV(max(FirstRelPos(1),1)) = CSum(1,max(FirstRelPos(1),1));
else
    p.In.LowMOI.STV.RrelV(:) = 0;
end

if ( FirstRelPos(2) > 0 )        
p.In.LowMOI.STV.RrelD(1:FirstRelPos(2)-1)    = 0;
p.In.LowMOI.STV.RrelD(max(FirstRelPos(2),1)) = CSum(2,max(FirstRelPos(2),1));
else
    p.In.LowMOI.STV.RrelD(:) = 0;
end   

if ( FirstRelPos(3) > 0 )        
p.In.LowMOI.STV.RrelDTot(1:FirstRelPos(3)-1)    = 0;
p.In.LowMOI.STV.RrelDTot(max(FirstRelPos(3),1)) = CSum(3,max(FirstRelPos(3),1));
else
    p.In.LowMOI.STV.RrelDTot(:) = 0;
end   

if ( FirstRelPos(4) > 0 )
p.In.LowMOI.STV.RrelP(1:FirstRelPos(4)-1)    = 0;
p.In.LowMOI.STV.RrelP(max(FirstRelPos(4),1)) = CSum(4,max(FirstRelPos(4),1));
else
    p.In.LowMOI.STV.RrelP(:) = 0;
end    

if ( FirstRelPos(5) > 0 )
p.In.LowMOI.COI.RrelV(1:FirstRelPos(5)-1)    = 0;
p.In.LowMOI.COI.RrelV(max(FirstRelPos(5),1)) = CSum(5,max(FirstRelPos(5),1));
else
    p.In.LowMOI.COI.RrelV(:) = 0;
end    

if ( FirstRelPos(6) > 0 )
p.In.LowMOI.COI.RrelD(1:FirstRelPos(6)-1)    = 0;
p.In.LowMOI.COI.RrelD(max(FirstRelPos(6),1)) = CSum(6,max(FirstRelPos(6),1));
else
    p.In.LowMOI.COI.RrelD(:) = 0;
end   

if ( FirstRelPos(7) > 0 )        
p.In.LowMOI.COI.RrelDTot(1:FirstRelPos(7)-1)    = 0;
p.In.LowMOI.COI.RrelDTot(max(FirstRelPos(7),1)) = CSum(7,max(FirstRelPos(7),1));
else
    p.In.LowMOI.STV.RrelDTot(:) = 0;
end   

if ( FirstRelPos(8) > 0 )
p.In.LowMOI.COI.RrelP(1:FirstRelPos(8)-1)    = 0;
p.In.LowMOI.COI.RrelP(max(FirstRelPos(8),1)) = CSum(8,max(FirstRelPos(8),1));
else
    p.In.LowMOI.COI.RrelP(:) = 0;
end    

return

