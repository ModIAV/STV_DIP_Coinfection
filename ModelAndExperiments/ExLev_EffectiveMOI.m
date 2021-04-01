function [p] = ExLev_EffectiveMOI(p, y, I_V, Iter)
%Intergrals calculates the integrals of virus production, infected
%   cells and the apoptosis rate of infected cells over all cell ages tau.
%
%
%   last revised: 2021/03/09

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

%% create empty release vectors
p.DynMOI.STV.RrelV{Iter+1,1} = zeros(1, length(p.In.tspan));
p.DynMOI.STV.RrelD{Iter+1,1} = zeros(1, length(p.In.tspan));
p.DynMOI.STV.RrelDTot{Iter+1,1} = zeros(1, length(p.In.tspan));
p.DynMOI.STV.RrelP{Iter+1,1} = zeros(1, length(p.In.tspan));
p.DynMOI.COI.RrelV{Iter+1,1} = zeros(1, length(p.In.tspan));
p.DynMOI.COI.RrelD{Iter+1,1} = zeros(1, length(p.In.tspan));
p.DynMOI.COI.RrelDTot{Iter+1,1} = zeros(1, length(p.In.tspan));
p.DynMOI.COI.RrelP{Iter+1,1} = zeros(1, length(p.In.tspan));

%% define low MOI pre-simulation cut-offs
CutoffSTV = 0.6;
CutoffDIP = 0.2;

%% define new ICs
p.DynMOI.STV.Ic = p.In.Ic.STV;
p.DynMOI.COI.Ic = p.In.Ic.COI;

if ( Iter == 1 && ( p.Ex.IntI_FirstWave(end,1) + p.Ex.IntI_FirstWave(end,3) > 0 ) ...
               && ( sum(diff(p.Ex.IntI_FirstWave(:,1))>0) > 0 ...
                    || sum(diff(p.Ex.IntI_FirstWave(:,3))>0) > 0 ) ...
   )
   
    TargetCells = sum(p.Ex.y_FirstWave(:,p.Idx.T|p.Idx.Ta|p.Idx.I_D),2) + p.Ex.IntI_FirstWave(:,1); % cells V and D can bind to, excludes I_CO, superinfection exclusion 
    
    %%% STV-infected cell
    if ( p.Ex.IntI_FirstWave(end,1) > 0 && sum(diff(p.Ex.IntI_FirstWave(:,1))>0) > 0 )    
       % calculate contribution factor
        Cfac = diff(p.Ex.IntI_FirstWave(:,1));
        Cfac(Cfac<0) = 0;
        Cfac = [Cfac/sum(Cfac); 0];
        
        % Ex, Fus + Cyt
        tmp = p.Ex.y_FirstWave(:,p.Idx.V) ./ TargetCells;
        tmp(isinf(tmp)) = 0;  % compute average number for infection vased on # of infected cells
        tmp(isinf(tmp)) = 0;  % prevent 0/0 result
        
        p.DynMOI.STV.Ic(Idx.Vex)   = tmp' * Cfac;   
        p.DynMOI.STV.Ic(Idx.Vfus)  = 1;
        p.DynMOI.STV.Ic(Idx.Vcyt)  = 1;    
        
        % AttHi, AttLo + En
        tmp = [p.Ex.y_FirstWave(:,p.Idx.VAttHi)   ./ TargetCells, ...
                 p.Ex.y_FirstWave(:,p.Idx.VAttLo) ./ TargetCells, ...
                   p.Ex.y_FirstWave(:,p.Idx.VEn)  ./ TargetCells];
        tmp(isinf(tmp)) = 0;  % compute average number for infection vased on # of infected cells
        tmp(isinf(tmp)) = 0;

        p.DynMOI.STV.Ic(Idx.VattHi)   = tmp(:,1)' * Cfac;
        p.DynMOI.STV.Ic(Idx.VattLo)   = tmp(:,2)' * Cfac;
        p.DynMOI.STV.Ic(Idx.Ven)      = tmp(:,3)' * Cfac;
        
        p.DynMOI.STV.Ic(Idx.Dex)      = 0;
    else
        p.DynMOI.STV.Ic(Idx.Dex)      = 0;          
    end
    
    %%% co-infected cell
    if ( p.Ex.IntI_FirstWave(end,3) > 0 && sum(diff(p.Ex.IntI_FirstWave(:,3))>0) > 0)
        % calculate contribution factor
        Cfac = diff(p.Ex.IntI_FirstWave(:,3));
        Cfac(Cfac<0) = 0;
        Cfac = [Cfac/sum(Cfac); 0];
        
        % Ex, Fus + Cyt
        tmp = p.Ex.y_FirstWave(:,p.Idx.V) ./ TargetCells;
        tmp(isinf(tmp)) = 0;  % compute average number for infection vased on # of infected cells
        tmp(isinf(tmp)) = 0;  % prevent 0/0 result
        
        p.DynMOI.COI.Ic(Idx.Vex)   = tmp' * Cfac;   
        p.DynMOI.COI.Ic(Idx.Vfus)  = 1;
        p.DynMOI.COI.Ic(Idx.Vcyt)  = 1; 

        tmp = p.Ex.y_FirstWave(:,p.Idx.D) ./ TargetCells;
        tmp(isinf(tmp)) = 0;  % compute average number for infection vased on # of infected cells    
        tmp(isinf(tmp)) = 0;
        
        p.DynMOI.COI.Ic(Idx.Dex)   = tmp' * Cfac;
        p.DynMOI.COI.Ic(Idx.Dfus)  = 1;
        p.DynMOI.COI.Ic(Idx.Dcyt)  = 1;        
        
        % AttHi, AttLo + En
        tmp = [p.Ex.y_FirstWave(:,p.Idx.VAttHi)   ./ TargetCells, ...
                 p.Ex.y_FirstWave(:,p.Idx.VAttLo) ./ TargetCells, ...
                   p.Ex.y_FirstWave(:,p.Idx.VEn)  ./ TargetCells];
        tmp(isinf(tmp)) = 0;  % compute average number for infection vased on # of infected cells
        tmp(isinf(tmp)) = 0;

        p.DynMOI.COI.Ic(Idx.VattHi)   = tmp(:,1)' * Cfac;
        p.DynMOI.COI.Ic(Idx.VattLo)   = tmp(:,2)' * Cfac;
        p.DynMOI.COI.Ic(Idx.Ven)      = tmp(:,3)' * Cfac;

        tmp = [p.Ex.y_FirstWave(:,p.Idx.DAttHi)   ./ TargetCells, ...
                 p.Ex.y_FirstWave(:,p.Idx.DAttLo) ./ TargetCells, ...
                   p.Ex.y_FirstWave(:,p.Idx.DEn)  ./ TargetCells];
        tmp(isinf(tmp)) = 0;  % compute average number for infection vased on # of infected cells
        tmp(isinf(tmp)) = 0;

        p.DynMOI.COI.Ic(Idx.DattHi)   = tmp(:,1)' * Cfac;
        p.DynMOI.COI.Ic(Idx.DattLo)   = tmp(:,2)' * Cfac;
        p.DynMOI.COI.Ic(Idx.Den)      = tmp(:,3)' * Cfac;      
    end     
        
else
    TargetCells = sum(y(p.Idx.T|p.Idx.Ta|p.Idx.I_D)) + I_V;
    % cells V and D can bind to, excludes I_CO, superinfection exclusion    
    
    %%% STV-infected cell
    if ( y(p.Idx.T) > 0 )

        p.DynMOI.STV.Ic(Idx.Vex)  = y(p.Idx.V)/TargetCells;
        p.DynMOI.STV.Ic(Idx.Vfus) = 1;
        p.DynMOI.STV.Ic(Idx.Vcyt) = 1;     

        p.DynMOI.STV.Ic(Idx.Dex)  = 0;

        if ( Iter > 0 )
            p.DynMOI.STV.Ic(Idx.VattHi) = y(p.Idx.VAttHi)/TargetCells;
            p.DynMOI.STV.Ic(Idx.VattLo) = y(p.Idx.VAttLo)/TargetCells;
            p.DynMOI.STV.Ic(Idx.Ven)    = y(p.Idx.VEn)   /TargetCells;
        end

    else
        p.DynMOI.STV.Ic = zeros(size(p.DynMOI.STV.Ic));

    end

    %%% co-infected cell
    if ( y(p.Idx.I_D) + I_V > 0 )

        p.DynMOI.COI.Ic(Idx.Vex)   = y(p.Idx.V)/TargetCells;
        p.DynMOI.COI.Ic(Idx.Vfus)  = 1;
        p.DynMOI.COI.Ic(Idx.Vcyt)  = 1; 

        p.DynMOI.COI.Ic(Idx.Dex)   = y(p.Idx.D)/TargetCells;
        p.DynMOI.COI.Ic(Idx.Dfus)  = 1;
        p.DynMOI.COI.Ic(Idx.Dcyt)  = 1;        

        if ( Iter > 0 )
            p.DynMOI.COI.Ic(Idx.VattHi) = y(p.Idx.VAttHi)/TargetCells;
            p.DynMOI.COI.Ic(Idx.VattLo) = y(p.Idx.VAttLo)/TargetCells;
            p.DynMOI.COI.Ic(Idx.Ven)    = y(p.Idx.VEn)   /TargetCells;

            p.DynMOI.COI.Ic(Idx.DattHi) = y(p.Idx.DAttHi)/TargetCells;
            p.DynMOI.COI.Ic(Idx.DattLo) = y(p.Idx.DAttLo)/TargetCells;
            p.DynMOI.COI.Ic(Idx.Den)    = y(p.Idx.DEn)   /TargetCells;
        end

    else
        p.DynMOI.COI.Ic = zeros(size(p.DynMOI.COI.Ic));        

    end
end

%% calculate dynamics considering current MOI / virus distribution
%%% STV-infected cell
if ( y(p.Idx.T) > 0 || Iter < 3 )
    
    % only simulate IC model when MOI > Cutoff, otherwise use pre-simualtion for low MOI infection
    if ( sum(p.DynMOI.STV.Ic(Idx.Vex|Idx.VattHi|Idx.VattLo|Idx.Ven)) > CutoffSTV ...
         || Iter < 3 )
        
        % apply upper limit for initial conditions
        if ( sum(p.DynMOI.STV.Ic > 1e5) )

            MaxNum = max(p.DynMOI.STV.Ic(Idx.Vex|Idx.VattHi|Idx.VattLo|Idx.Ven));

            p.DynMOI.STV.Ic(Idx.Vex)    = 1e5 * p.DynMOI.STV.Ic(Idx.Vex)    / MaxNum;
            p.DynMOI.STV.Ic(Idx.VattHi) = 1e5 * p.DynMOI.STV.Ic(Idx.VattHi) / MaxNum;
            p.DynMOI.STV.Ic(Idx.VattLo) = 1e5 * p.DynMOI.STV.Ic(Idx.VattLo) / MaxNum;
            p.DynMOI.STV.Ic(Idx.Ven)    = 1e5 * p.DynMOI.STV.Ic(Idx.Ven)    / MaxNum;     
        end       

        if ( p.CompileFlag )       
            try      
                result.DynMOI.STV = IQMPsimulate('InModel_MexFile', p.In.tspan, p.DynMOI.STV.Ic, ...
                                                 [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 0]);                                          
            catch
                error('Simulation with IQMPsimulate produced an error.')
            end            
        else
            tmp = p.InModel;
            tmp = IQMparameters(tmp, [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 0]);

            result.DynMOI.STV     = IQMsimulate(tmp, 'ode23s', p.In.tspan, p.DynMOI.STV.Ic);     
        end

        % divide variables and reactions for CoinfectionModel (as the reactions are also declared as variables)
        for i = 1 : length(result.DynMOI.STV.variables)
            if ( strcmp('r', result.DynMOI.STV.variables{i}(1)) )
                tmp = i;
                break;
            end
        end
        result.DynMOI.STV.reactions      = result.DynMOI.STV.variables(i:end);
        result.DynMOI.STV.reactionvalues = result.DynMOI.STV.variablevalues(:,i:end);
        result.DynMOI.STV.variables      = result.DynMOI.STV.variables(1:i-1);
        result.DynMOI.STV.variablevalues = result.DynMOI.STV.variablevalues(:,1:i-1);    

        % set all DIP related values to 0, values can reach 1e-12 due to numerical instabilities during SBsPD simulate (?)
        tmp = [21:9:84, 104, 106, 108]; % vector of intracellular states related to DIPs
        result.DynMOI.STV.statevalues(:,tmp) = zeros(size(result.DynMOI.STV.statevalues,1), length(tmp));
        tmp = [20, 29, 39:9:93, 107:9:170, 172, 200, 203, 205]; % vector of intracellular states related to DIPs
        result.DynMOI.STV.reactionvalues(:,tmp) = zeros(size(result.DynMOI.STV.reactionvalues,1), length(tmp));    

        % define release rates for the specific MOI conditions
        p.DynMOI.STV.RrelV{Iter+1,1}    = result.DynMOI.STV.reactionvalues(:,strcmp('rRel', p.In.ReactionNames))';
        p.DynMOI.STV.RrelD{Iter+1,1}    = 0 * result.DynMOI.STV.reactionvalues(:,strcmp('rRelD', p.In.ReactionNames))'; % p.DIPdeNovo * result.DynMOI.STV.reactionvalues(:,strcmp('rRelD', p.In.ReactionNames))';
        p.DynMOI.STV.RrelDTot{Iter+1,1} = 0 * result.DynMOI.STV.reactionvalues(:,strcmp('rRelDIPtot', p.In.ReactionNames))'; % p.DIPdeNovo * result.DynMOI.STV.reactionvalues(:,strcmp('rRelDIPtot', p.In.ReactionNames))';
        p.DynMOI.STV.RrelP{Iter+1,1}    = result.DynMOI.STV.reactionvalues(:,strcmp('rRelP', p.In.ReactionNames))';

        % failsave to prevent negative virus particle release values due to
        % numerical inaccuricies during SBPD simulate
        if  ( sum(p.DynMOI.STV.RrelP{Iter+1,1}<-1e-1) > 0 )
            if ( sum(p.DynMOI.STV.RrelP{Iter+1,1}(16:end) < -1e-1) > 0 )
                error('SBPDsimulate produced negative values for virus particle release!')
            end
        end 
        if ( sum(sum(result.DynMOI.STV.statevalues<-1)) > 0 )
            error('SBPDsimulate produced negative values during simulation!')
        end        

        % set negative values to zero (some parameter combinations induce negative release (?))
        p.DynMOI.STV.RrelV{Iter+1,1}(p.DynMOI.STV.RrelV{Iter+1,1}<0) = 0;
        p.DynMOI.STV.RrelD{Iter+1,1}(p.DynMOI.STV.RrelD{Iter+1,1}<0) = 0;
        p.DynMOI.STV.RrelDTot{Iter+1,1}(p.DynMOI.STV.RrelDTot{Iter+1,1}<0) = 0;
        p.DynMOI.STV.RrelP{Iter+1,1}(p.DynMOI.STV.RrelP{Iter+1,1}<0) = 0;

        AllValuesSTV = [result.DynMOI.STV.statevalues, result.DynMOI.STV.variablevalues];

    else     
        % use pre-simulated results for low MOI infections
        p.DynMOI.STV.RrelV{Iter+1,1} = p.In.LowMOI.STV.RrelV;
        p.DynMOI.STV.RrelD{Iter+1,1} = p.In.LowMOI.STV.RrelD; % p.DIPdeNovo * result.DynMOI.STV.reactionvalues(:,strcmp('rRelD', p.In.ReactionNames))';
        p.DynMOI.STV.RrelDTot{Iter+1,1} = p.In.LowMOI.STV.RrelDTot; % p.DIPdeNovo * result.DynMOI.STV.reactionvalues(:,strcmp('rRelDIPtot', p.In.ReactionNames))';
        p.DynMOI.STV.RrelP{Iter+1,1} = p.In.LowMOI.STV.RrelP;        
        
        AllValuesSTV = p.In.LowMOI.AllValuesSTV; 
    end        
        
    % save RNA dynamics + particle release for parameter optimization   
    ShiftedTime = round(p.InToPop.ExpTime/p.In.h)+1-Iter;
    if ( Iter > 0 );  tmp = find(ShiftedTime<=0,1,'last');
    else              tmp = 0;
    end
    ShiftedTime = ShiftedTime(tmp+1:end);
    ShiftedTime(ShiftedTime>p.In.SimTime/p.In.h) = [];

    for i = 1 : length(p.InToPop.states)  
        
        
        p.InToPop.STV.statevalues{i}(Iter+1,tmp+1:tmp+length(ShiftedTime)) = ...
          AllValuesSTV(ShiftedTime,strcmp(p.InToPop.states{i},p.InToPop.StateVarNames));                     
    end         

end

%%% co-infected cell
if ( I_V > 0 || y(p.Idx.I_D) > 0 || Iter < 3 )
    
    % only simulate IC model when MOI/MODIP > Cutoff, otherwise use pre-simualtion for low MOI infection
    if ( sum(p.DynMOI.COI.Ic(Idx.Vex|Idx.VattHi|Idx.VattLo|Idx.Ven)) > CutoffSTV ...
         || sum(p.DynMOI.COI.Ic(Idx.Dex|Idx.DattHi|Idx.DattLo|Idx.Den)) > CutoffDIP ...
         || Iter < 3 )

        % apply upper limit for initial conditions, but keep ratio for Vex/Dex
        if ( sum(p.DynMOI.COI.Ic > 1e5) )      

            MaxNum = max(p.DynMOI.COI.Ic(Idx.Vex|Idx.VattHi|Idx.VattLo|Idx.Ven| ...
                                         Idx.Dex|Idx.DattHi|Idx.DattLo|Idx.Den));

            p.DynMOI.COI.Ic(Idx.Vex)    = 1e5 * p.DynMOI.COI.Ic(Idx.Vex)    / MaxNum;
            p.DynMOI.COI.Ic(Idx.VattHi) = 1e5 * p.DynMOI.COI.Ic(Idx.VattHi) / MaxNum;
            p.DynMOI.COI.Ic(Idx.VattLo) = 1e5 * p.DynMOI.COI.Ic(Idx.VattLo) / MaxNum;
            p.DynMOI.COI.Ic(Idx.Ven)    = 1e5 * p.DynMOI.COI.Ic(Idx.Ven)    / MaxNum;   

            p.DynMOI.COI.Ic(Idx.Dex)    = 1e5 * p.DynMOI.COI.Ic(Idx.Dex)    / MaxNum;
            p.DynMOI.COI.Ic(Idx.DattHi) = 1e5 * p.DynMOI.COI.Ic(Idx.DattHi) / MaxNum;
            p.DynMOI.COI.Ic(Idx.DattLo) = 1e5 * p.DynMOI.COI.Ic(Idx.DattLo) / MaxNum;
            p.DynMOI.COI.Ic(Idx.Den)    = 1e5 * p.DynMOI.COI.Ic(Idx.Den)    / MaxNum;        

        end   
        
        if ( p.CompileFlag )              
            try
                result.DynMOI.COI = IQMPsimulate('InModel_MexFile', p.In.tspan, p.DynMOI.COI.Ic, ...
                                                 [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 1]);  
            catch   
                try % second try to avoid rare error causing a crash with very specific initial conditions
                    p.DynMOI.COI.Ic(Idx.Vex) = p.DynMOI.COI.Ic(Idx.Vex)*0.999;
                    p.DynMOI.COI.Ic(Idx.Dex) = p.DynMOI.COI.Ic(Idx.Dex)*0.999;
                    result.DynMOI.COI = IQMPsimulate('InModel_MexFile', p.In.tspan, p.DynMOI.COI.Ic, ...
                                                     [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 1]);                            
                    p.Error_099(1) = p.Error_099(1) + 1;
                catch
                    error('Simulation with IQMPsimulate produced an error.')
                end
            end
        else
            tmp = p.InModel;
            tmp = IQMparameters(tmp, [p.In.PaOpt, 'kSynV', 'DF'], [p.xIn(1:p.In.nPaOpt), p.In.kSynV, 1]);

            result.DynMOI.COI     = IQMsimulate(tmp, 'ode23s', p.In.tspan, p.DynMOI.COI.Ic);    
        end

        % divide variables and reactions for CoinfectionModel (as the reactions are also declared as variables)
        for i = 1 : length(result.DynMOI.COI.variables)
            if ( strcmp('r', result.DynMOI.COI.variables{i}(1)) )
                tmp = i;
                break;
            end
        end
        result.DynMOI.COI.reactions      = result.DynMOI.COI.variables(i:end);
        result.DynMOI.COI.reactionvalues = result.DynMOI.COI.variablevalues(:,i:end);
        result.DynMOI.COI.variables      = result.DynMOI.COI.variables(1:i-1);
        result.DynMOI.COI.variablevalues = result.DynMOI.COI.variablevalues(:,1:i-1);

        % define release rates for the specific MOI/MODIP conditions
        p.DynMOI.COI.RrelV{Iter+1,1}    = result.DynMOI.COI.reactionvalues(:,strcmp('rRel', p.In.ReactionNames))';
        p.DynMOI.COI.RrelD{Iter+1,1}    = result.DynMOI.COI.reactionvalues(:,strcmp('rRelD', p.In.ReactionNames))';
        p.DynMOI.COI.RrelDTot{Iter+1,1} = result.DynMOI.COI.reactionvalues(:,strcmp('rRelDIPtot', p.In.ReactionNames))';
        p.DynMOI.COI.RrelP{Iter+1,1}    = result.DynMOI.COI.reactionvalues(:,strcmp('rRelP', p.In.ReactionNames))';    
        
        % set negative values to zero (some parameter combinations induce negative release (?))
        p.DynMOI.COI.RrelV{Iter+1,1}(p.DynMOI.COI.RrelV{Iter+1,1}<0) = 0;
        p.DynMOI.COI.RrelD{Iter+1,1}(p.DynMOI.COI.RrelD{Iter+1,1}<0) = 0;
        p.DynMOI.COI.RrelDTot{Iter+1,1}(p.DynMOI.COI.RrelDTot{Iter+1,1}<0) = 0;
        p.DynMOI.COI.RrelP{Iter+1,1}(p.DynMOI.COI.RrelP{Iter+1,1}<0) = 0;    

        
        AllValuesCOI = [result.DynMOI.COI.statevalues, result.DynMOI.COI.variablevalues];
    else              
        % use pre-simulated results for low MOI infections
        p.DynMOI.COI.RrelV{Iter+1,1} = p.In.LowMOI.COI.RrelV;
        p.DynMOI.COI.RrelD{Iter+1,1} = p.In.LowMOI.COI.RrelD;
        p.DynMOI.COI.RrelDTot{Iter+1,1} = p.In.LowMOI.COI.RrelDTot;
        p.DynMOI.COI.RrelP{Iter+1,1} = p.In.LowMOI.COI.RrelP;        
        
        AllValuesCOI = p.In.LowMOI.AllValuesCOI;       
    end
    
    % save RNA dynamics + particle release for parameter optimization
    ShiftedTime = round(p.InToPop.ExpTime/p.In.h)+1-Iter;
    if ( Iter > 0 );  tmp = find(ShiftedTime<=0,1,'last');
    else              tmp = 0;
    end
    ShiftedTime = ShiftedTime(tmp+1:end);
    ShiftedTime(ShiftedTime>p.In.SimTime/p.In.h) = [];

    for i = 1 : length(p.InToPop.states)
        p.InToPop.COI.statevalues{i}(Iter+1,tmp+1:tmp+length(ShiftedTime)) = ...
          AllValuesCOI(ShiftedTime,strcmp(p.InToPop.states{i},p.InToPop.StateVarNames));                     
    end   
        
end

%% allow virus particle release only after at least 1 full particle was formed
%  -> Vrel(Vrel < 1) = 0
AllReleases =  [p.DynMOI.STV.RrelV{Iter+1,1};...
                p.DynMOI.STV.RrelD{Iter+1,1};...
                p.DynMOI.STV.RrelDTot{Iter+1,1};...
                p.DynMOI.STV.RrelP{Iter+1,1};...
                p.DynMOI.COI.RrelV{Iter+1,1};...
                p.DynMOI.COI.RrelD{Iter+1,1};...
                p.DynMOI.COI.RrelDTot{Iter+1,1};...
                p.DynMOI.COI.RrelP{Iter+1,1}];

% find time points at which cumualtive V/D/PRel > 1
CSum = cumsum(AllReleases,2);
FirstRelPos = zeros(8,1);

for i = 1 : size(CSum,1)
    if ( sum(CSum(i,:) > 1) )       
        FirstRelPos(i) = find(CSum(i,:)>1, 1, 'first');
    end
end

% handle special case when no STV infected cells exist
if ( isempty(p.DynMOI.STV.RrelV{Iter+1,1}) )
    CSum = [zeros(size(CSum));CSum];    
    FirstRelPos = [zeros(4,1); FirstRelPos];
end

% adjust release vectors
if ( y(p.Idx.T) > 0 || Iter < 3 )
    if ( FirstRelPos(1) > 0 )
        p.DynMOI.STV.RrelV{Iter+1,1}(1:FirstRelPos(1)-1)    = 0;
        p.DynMOI.STV.RrelV{Iter+1,1}(max(FirstRelPos(1),1)) = CSum(1,max(FirstRelPos(1),1));
    else
        p.DynMOI.STV.RrelV{Iter+1,1}(:) = 0;
    end

    if ( FirstRelPos(2) > 0 )        
    p.DynMOI.STV.RrelD{Iter+1,1}(1:FirstRelPos(2)-1)    = 0;
    p.DynMOI.STV.RrelD{Iter+1,1}(max(FirstRelPos(2),1)) = CSum(2,max(FirstRelPos(2),1));
    else
        p.DynMOI.STV.RrelD{Iter+1,1}(:) = 0;
    end   
    
    if ( FirstRelPos(3) > 0 )        
    p.DynMOI.STV.RrelDTot{Iter+1,1}(1:FirstRelPos(3)-1)    = 0;
    p.DynMOI.STV.RrelDTot{Iter+1,1}(max(FirstRelPos(3),1)) = CSum(3,max(FirstRelPos(3),1));
    else
        p.DynMOI.STV.RrelDTot{Iter+1,1}(:) = 0;
    end     

    if ( FirstRelPos(4) > 0 )
    p.DynMOI.STV.RrelP{Iter+1,1}(1:FirstRelPos(4)-1)    = 0;
    p.DynMOI.STV.RrelP{Iter+1,1}(max(FirstRelPos(4),1)) = CSum(4,max(FirstRelPos(4),1));
    else
        p.DynMOI.STV.RrelP{Iter+1,1}(:) = 0;
    end    
end

if ( I_V > 0 || y(p.Idx.I_D) > 0 || Iter < 3 )
    if ( FirstRelPos(5) > 0 )
    p.DynMOI.COI.RrelV{Iter+1,1}(1:FirstRelPos(5)-1)    = 0;
    p.DynMOI.COI.RrelV{Iter+1,1}(max(FirstRelPos(5),1)) = CSum(5,max(FirstRelPos(5),1));
    else
        p.DynMOI.COI.RrelV{Iter+1,1}(:) = 0;
    end    

    if ( FirstRelPos(6) > 0 )
    p.DynMOI.COI.RrelD{Iter+1,1}(1:FirstRelPos(6)-1)    = 0;
    p.DynMOI.COI.RrelD{Iter+1,1}(max(FirstRelPos(6),1)) = CSum(6,max(FirstRelPos(6),1));
    else
        p.DynMOI.COI.RrelD{Iter+1,1}(:) = 0;
    end   
    
    if ( FirstRelPos(7) > 0 )
    p.DynMOI.COI.RrelDTot{Iter+1,1}(1:FirstRelPos(7)-1)    = 0;
    p.DynMOI.COI.RrelDTot{Iter+1,1}(max(FirstRelPos(7),1)) = CSum(7,max(FirstRelPos(7),1));
    else
        p.DynMOI.COI.RrelDTot{Iter+1,1}(:) = 0;
    end      

    if ( FirstRelPos(8) > 0 )
    p.DynMOI.COI.RrelP{Iter+1,1}(1:FirstRelPos(8)-1)    = 0;
    p.DynMOI.COI.RrelP{Iter+1,1}(max(FirstRelPos(8),1)) = CSum(8,max(FirstRelPos(8),1));
    else
        p.DynMOI.COI.RrelP{Iter+1,1}(:) = 0;
    end    
end









