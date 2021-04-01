function [IntRrelI, IntI, IntApoI, RowMatI, z] = ExLev_Integrals(n, y, MatI, p, varargin)
%Intergrals calculates the integrals of virus production, infected
%   cells and the apoptosis rate of infected cells over all cell ages tau.
%
%   last revised: 2021/03/09

% initialize values
RowMatI     = zeros(n,2);         % vector containing concentration of I_V and I_CO cells at current time step "t" for different ages "tau"
                                  % (LineMatI(1,:): tau = t; LineMatI(end,:): tau = 0)                             
                                  
if ( isempty(varargin) ) 
    SIE_step = p.Ex.SIE_Time/p.Ex.h + 1; % time step until infected cell prevent superinfection
    InSim_step = p.In.SimTime/p.Ex.h + 1; % time step until intracellular cells are simulated
    ExSim_step = p.SimTime/p.Ex.h + 1; % time step for full simulation
else
    RedFact = varargin{1};                      % reduction factor for step size   
    SIE_step = p.Ex.SIE_Time/(p.Ex.h*RedFact) + 1; % use original p.Ex.h in case it was reduced for SmallerSteps calculation
    InSim_step = p.In.SimTime/(p.Ex.h*RedFact) + 1; % time step until intracellular cells are simulated
    ExSim_step = p.SimTime/(p.Ex.h*RedFact) + 1; % time step for full simulation
    
    % set the last apoptosis rate values that are used to the current value when smaller ste-size is used
    Iter0   = varargin{2};    
    p.Ex.ApoI = [p.Ex.ApoI(1:Iter0), repmat(p.Ex.ApoI(Iter0), 1, RedFact)]; 
end

if ( SIE_step < n )
    I_V     = sum(MatI(end-SIE_step+1:end,1));
else
    I_V     = sum(MatI(:,1));
end
I_CO        = sum(MatI(:,2)); % number of I_CO cells in previous time step

%% calculate infection forces
TargetCellsV = y(p.Idx.T) + y(p.Idx.Ta) + y(p.Idx.I_D);
TargetCellsD = y(p.Idx.T) + y(p.Idx.Ta) + I_V; 

if ( TargetCellsV > 0 )        
    InfForceV_T  = p.Ex.Finf*p.Ex.kFus*y(p.Idx.VEn)*y(p.Idx.T) ...   % rInfV*T/(T+Ta+I_D)
                   /TargetCellsV;            

    InfForceV_ID = p.Ex.Finf*p.Ex.kFus*y(p.Idx.VEn)*y(p.Idx.I_D)...  % rInfV*I_D/(T+Ta+I_D)
                   /TargetCellsV;
         
else
    InfForceV_T  = 0;
    InfForceV_ID = 0;  
end

if ( TargetCellsD > 0 )  
    InfForceD    = p.Ex.Finf*p.Ex.kFus*y(p.Idx.DEn)...               % rInfD/(T+Ta+I_V+Ia_V)   
                   /TargetCellsD;   
               
    InfForceD_IV = p.Ex.Finf*p.Ex.kFus*y(p.Idx.DEn)*I_V ...          % rInfD*I_V/(T+Ta+I_V)
                   /TargetCellsD; 
           
else
    InfForceD    = 0;
    InfForceD_IV = 0;
end

% check if infection force is too high and would lead to overconsumption of available cells
% -> scale forces down to prevent overshooting
if ( y(p.Idx.T) > 0 )
    InfForceV = max(0, p.Ex.Finf*p.Ex.kFus*y(p.Idx.VEn)/TargetCellsV);
    if ( InfForceV + InfForceD > 1 / p.Ex.h )
        InfForceV_T = InfForceV/(InfForceV+InfForceD)/p.Ex.h * y(p.Idx.T);
    end
end

if ( y(p.Idx.I_D) > 0 )
    InfForceV = max(0, p.Ex.Finf*p.Ex.kFus*y(p.Idx.VEn)/TargetCellsV);
    if ( InfForceV > 1 / p.Ex.h )
        InfForceV_ID = y(p.Idx.I_D) / p.Ex.h;
    end
end

if ( I_V > 0 )
    if ( InfForceD > 1 / p.Ex.h )
        InfForceD    = 1 / p.Ex.h; % if D would infect 100+% of I_V cells -> they all disappear + become I_CO
        InfForceD_IV = I_V / p.Ex.h; % maximum of I_V -> I_CO can be taken from 100% of I_V cells
    end
end

%% calculate integrals with Euler's method

% create new sub-population, newly infected STV- and COI-cells ( further down in RowMatI = younger cell population)
RowMatI(end,1) = InfForceV_T * p.Ex.h;                   % I_V cells
RowMatI(end,2) = (InfForceV_ID + InfForceD_IV) * p.Ex.h; % I_CO cells

% STV cells
RemovalTerm = zeros(n-1,1);
if ( n > SIE_step )    
    RemovalTerm(1:n-SIE_step,1)     = max(-1, - (p.Ex.ApoI(n:-1:SIE_step+1))           * p.Ex.h);
    RemovalTerm(n-SIE_step+1:end,1) = max(-1, - (p.Ex.ApoI(SIE_step:-1:2) + InfForceD) * p.Ex.h);
else
    RemovalTerm(:,1)                = max(-1, - (p.Ex.ApoI(n:-1:2) + InfForceD)        * p.Ex.h);
end
RowMatI(1:end-1,1) = MatI(:,1) + RemovalTerm .* MatI(:,1);

% COI cells
RemovalTerm = [];
RemovalTerm(:,1) = max(-1, - (p.Ex.ApoI(n:-1:2)) * p.Ex.h);
RowMatI(1:end-1,2) = MatI(:,2) + RemovalTerm .* MatI(:,2); 

% virus particle release
ReleaseCells = {p.DynMOI.STV.RrelV; p.DynMOI.STV.RrelD; ...
                p.DynMOI.STV.RrelDTot; p.DynMOI.STV.RrelP;...
                p.DynMOI.COI.RrelV; p.DynMOI.COI.RrelD; ...
                p.DynMOI.COI.RrelDTot; p.DynMOI.COI.RrelP};
IntRrel = zeros(1,8);            
for i = 1 : 8    
    if ( sum( cellfun(@sum, ReleaseCells{i}) ) == 0 )    
        IntRrel(i) = 0;
        continue;
    end
    
    if ( n > InSim_step - 1 )
        tmpMat = ReleaseCells{i}(n-InSim_step+1 : min(n, ExSim_step));
    else                                                   
        tmpMat = ReleaseCells{i}(1:n);
    end 

    tmpMat = cell2mat(tmpMat);
    
    % add empty rows as zeros     
    tmpMat = [tmpMat; zeros(min(n, InSim_step)-size(tmpMat,1), InSim_step)]; 

    j = size(tmpMat,1);
    tmpDiag = tmpMat(j*(j-1)+1 : -j+1 : j);    

    if ( n > j )   
        if ( i < 5 )
            tmpSum = RowMatI(n-j+1:end,1) .* tmpDiag';
        else
            tmpSum = RowMatI(n-j+1:end,2) .* tmpDiag';
        end
    else
        if ( i < 5 )
            tmpSum = RowMatI(:,1) .* tmpDiag';
        else
            tmpSum = RowMatI(:,2) .* tmpDiag';
        end
    end

    IntRrel(i) = sum(tmpSum);

end

IntRrelIV_V = IntRrel(1);
IntRrelIV_D = IntRrel(2);
IntRrelIV_DTot = IntRrel(3);
IntRrelIV_P = IntRrel(4);

IntRrelICO_V = IntRrel(5);
IntRrelICO_D = IntRrel(6);
IntRrelICO_DTot = IntRrel(7);
IntRrelICO_P = IntRrel(8);

% cell apoptosis
IntApoI_V  = sum(RowMatI(:,1) .* [p.Ex.ApoI(size(RowMatI,1):-1:2)'; 0]);
IntApoI_CO = sum(RowMatI(:,2) .* [p.Ex.ApoI(size(RowMatI,1):-1:2)'; 0]);

% infected cell infection time categories
if ( sum(RowMatI(:,1)) > 1e-5 )   
    if ( SIE_step < n )
        IntI_V     = sum(RowMatI(end-SIE_step+1:end,1));
        IntI_V_SIE = sum(RowMatI(1:end-SIE_step,1)); % Iv cells that prevent superinfection (are already infected longer than p.Ex.SIE_Time)
    else
        IntI_V     = sum(RowMatI(:,1));
        IntI_V_SIE = 0;
    end
else
    IntI_V      = 0;
    IntI_V_SIE  = 0;
end

if ( sum(RowMatI(:,2)) > 1e-20 )
    IntI_CO = sum(RowMatI(:,2));
else
    IntI_CO = 0;
end

if ( isempty(IntI_V) || isempty(IntI_CO) )
    error('Integral calculation step failed!')
end

% collect results
IntI     = [IntI_V, IntI_V_SIE, IntI_CO];
IntRrelI = [IntRrelIV_V,  IntRrelIV_D,  IntRrelIV_DTot,  IntRrelIV_P, ...
            IntRrelICO_V, IntRrelICO_D, IntRrelICO_DTot, IntRrelICO_P];
IntApoI  = [IntApoI_V, IntApoI_CO];




























