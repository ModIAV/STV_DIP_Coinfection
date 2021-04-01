function dydt = ExLev_EquationSystem(t, y, IntRrelI, IntI, IntApoI, p)
%EquationSystem contains the model equations of the age-segregated 
%   model for influenza virus infection at the extracellular level
%
%   last revised: 2021/03/09

%% retrieve state variables                               
T      = y(p.Idx.T);
Ta     = y(p.Idx.Ta);

Ia     = y(p.Idx.Ia);
V      = y(p.Idx.V);
VAttHi = y(p.Idx.VAttHi);
VAttLo = y(p.Idx.VAttLo);
VEn    = y(p.Idx.VEn);

I_D    = y(p.Idx.I_D);
D      = y(p.Idx.D);
DAttHi = y(p.Idx.DAttHi);
DAttLo = y(p.Idx.DAttLo);
DEn    = y(p.Idx.DEn);

% VTot   = y(p.Idx.VTot);
% DTot   = y(p.Idx.DTot);
% Par    = y(p.Idx.Par);
% 
% Cviab  = y(p.Idx.Cviab);
% Dead   = y(p.Idx.Dead);
% 
% VRelTot= y(p.Idx.VRelTot);

VOnlyRel = y(p.Idx.VOnlyRel);
DOnlyRel = y(p.Idx.DOnlyRel);

Rm5_ID = y(p.Idx.Rm5_ID);
Rm9_ID = y(p.Idx.Rm9_ID);

%% values from "Integrals" file
IntI_V      = IntI(1);
IntI_V_SIE  = IntI(2);
IntI_CO     = IntI(3);

IntRrelIV_V  = IntRrelI(1);
IntRrelIV_D  = IntRrelI(2);
IntRrelIV_DTot  = IntRrelI(3);
IntRrelIV_P  = IntRrelI(4);
IntRrelICO_V = IntRrelI(5);
IntRrelICO_D = IntRrelI(6);
IntRrelICO_DTot = IntRrelI(7);
IntRrelICO_P = IntRrelI(8);

IntApoI_V  = IntApoI(1);
IntApoI_CO = IntApoI(2);        
        
%% set parameter values
p.Ex.kDegD   = p.Ex.kDegV;
p.Ex.kApoI_D = p.Ex.kApoT;

%% algebraic equations
if ( p.Ex.Moi + p.Ex.Modip > 6 )
    mu_factor = p.Ex.mu_MODIP;
else
    mu_factor = 1;
end
mu = mu_factor * max(0, p.Ex.MuMax/p.Ex.Tmax * (p.Ex.Tmax - (T + IntI_V + IntI_V_SIE + I_D + IntI_CO)));

% number of free binding sites, specific for each STVs and DIPs (assuming superinfection protection)
BHiV = p.Ex.BtotHi*(T + Ta + I_D) ...                 % 1) total # of sites that V could bind to
       - (T + Ta)*( max(0, DAttHi/(T + Ta + IntI_V) ...
                    + VAttHi/(T + Ta + I_D)) ) ...    % 2) V and D particles on T and Ta that block sites
       - (I_D)   *( max(0, VAttHi/(T + Ta + I_D)) );  % 3) V particles on I_D that block sites   
BHiV = max(0,BHiV);

BLoV = p.Ex.BtotLo*(T + Ta + I_D) ...
       - (T + Ta)*( max(0, DAttLo/(T + Ta + IntI_V) ...
                        + VAttLo/(T + Ta + I_D)) ) ...
       - (I_D)   *( max(0, VAttLo/(T + Ta + I_D)) ); 
BLoV = max(0,BLoV);

BHiD = p.Ex.BtotHi*(T + Ta + IntI_V) ...                 % 1) total # of sites that D could bind to
       - (T + Ta)*( max(0, DAttHi/(T + Ta + IntI_V) ...
                        + VAttHi/(T + Ta + I_D)) ) ...   % 2) V and D particles on T and Ta that block sites
       - (IntI_V)*( max(0, DAttHi/(T + Ta + IntI_V)) );   % 3) D particles on I_V that block sites 
BHiD = max(0,BHiD);   

BLoD = p.Ex.BtotLo*(T + Ta + IntI_V) ...
       - (T + Ta)*( max(0, DAttLo/(T + Ta + IntI_V) ...
                        + VAttLo/(T + Ta + I_D)) ) ...
       - (IntI_V)*( max(0, DAttLo/(T + Ta + IntI_V)) ); 
BLoD = max(0,BLoD);  

if ( T + Ta > 0 )
    rLysTa = p.Ex.kLys*Ta / (T + Ta + IntI_V + I_D);
else
    rLysTa = 0;
end


if ( T + Ta + I_D > 0 )    
    rInfV = p.Ex.Finf*p.Ex.kFus*VEn/(T + Ta + I_D);
else
    rInfV = 0;
end

if ( T + Ta + IntI_V > 0 )
    rInfD = p.Ex.Finf*p.Ex.kFus*DEn/(T + Ta + IntI_V);  
else
    rInfD = 0;
end

% check if infection force is too high and would lead to overconsumption of available cells
% -> scale forces down to prevent overshooting (and undershooting)
ScaleFlag = 0;
if ( (T + Ta + I_D) > 0 && (mu - rInfV - rInfD - p.Ex.kApoT) * p.Ex.h < -1 )
    rInfV_save = rInfV; % keep original value to apply for attached + enveloped virus particles
    ScaleFlag = 1;
    
    p_InfV = rInfV / (rInfV + rInfD);
    rInfV = p_InfV     * (1/p.Ex.h + mu - p.Ex.kApoT);
    rInfD = (1-p_InfV) * (1/p.Ex.h + mu - p.Ex.kApoT);
end


%% differential equations
dydt = zeros(size(y));

%%% cell populations
dydt(p.Idx.T)      = (mu              -  rInfV      - rInfD      - p.Ex.kApoT)  *T;
dydt(p.Idx.Ta)     = p.Ex.kApoT*T     - (rInfV      + rInfD      + p.Ex.kLys)   *Ta;
%dydt(p.Idx.I_V)  -- calculated in "Integrals"

if ( ~ScaleFlag )
    dydt(p.Idx.I_D)    = rInfD*T + mu*I_D - (rInfV                     + p.Ex.kApoI_D)*I_D;
else % use original rInfV as it does not conflict with rInfD when infecting I_D
    dydt(p.Idx.I_D)    = rInfD*T + mu*I_D - (min(rInfV_save, 1/p.Ex.h) + p.Ex.kApoI_D)*I_D;
end

%dydt(p.Idx.I_CO) -- calculated in "Integrals"
dydt(p.Idx.Ia)     = IntApoI_V + p.Ex.kApoI_D*I_D + IntApoI_CO + (rInfV+rInfD)*Ta - p.Ex.kLys*Ia;
%dydt(p.Idx.Cviab) -- calculated in "Euler"
dydt(p.Idx.Dead)   = p.Ex.kLys*(Ta + Ia);

%%% attachment to the cell surface
if ( ScaleFlag && rInfV > 0 ) % if rInfV was scaled, readjust virus particle absorption for infection    
    InfRatio = rInfV*T / (rInfV*T + min(rInfV_save, 1/p.Ex.h)*I_D);
    if ( isnan(InfRatio) )
        rInfV = 0;
    else
        rInfV    = InfRatio*rInfV + (1-InfRatio)*min(rInfV_save, 1/p.Ex.h);
    end
end

dydt(p.Idx.VAttHi) = p.Ex.kAttcHi*BHiV*V - (p.Ex.kDisHi + p.Ex.kEn)*VAttHi - (rInfV + rLysTa)*VAttHi; % could be made more complicated by also considering removal of Vs when a I_V cell is infected by a D (if you assume still Vs on an I_V cell)
dydt(p.Idx.VAttLo) = p.Ex.kAttcLo*BLoV*V - (p.Ex.kDisLo + p.Ex.kEn)*VAttLo - (rInfV + rLysTa)*VAttLo;

dydt(p.Idx.DAttHi) = p.Ex.kAttcHi*BHiD*D - (p.Ex.kDisHi + p.Ex.kEn)*DAttHi - (rInfD + rLysTa)*DAttHi;
dydt(p.Idx.DAttLo) = p.Ex.kAttcLo*BLoD*D - (p.Ex.kDisLo + p.Ex.kEn)*DAttLo - (rInfD + rLysTa)*DAttLo;

%%% endocytosis and fusion
dydt(p.Idx.VEn)    = p.Ex.kEn*(VAttHi+VAttLo) - p.Ex.kFus*VEn             - (rInfV + rLysTa)*VEn;

dydt(p.Idx.DEn)    = p.Ex.kEn*(DAttHi+DAttLo) - p.Ex.kFus*DEn             - (rInfD + rLysTa)*DEn;

%%% virus concentrations                 
dydt(p.Idx.V)      = IntRrelIV_V + IntRrelICO_V - p.Ex.kDegV*V ... 
                     + p.Ex.kDisHi*VAttHi + p.Ex.kDisLo*VAttLo - (p.Ex.kAttcHi*BHiV + p.Ex.kAttcLo*BLoV)*V;
               
dydt(p.Idx.D)      = IntRrelIV_D + IntRrelICO_D - p.Ex.kDegD*D ...
                     + p.Ex.kDisHi*DAttHi + p.Ex.kDisLo*DAttLo - (p.Ex.kAttcHi*BHiD + p.Ex.kAttcLo*BLoD)*D;                  
                 
dydt(p.Idx.DTot)   = IntRrelIV_DTot + IntRrelICO_DTot;               

dydt(p.Idx.Par)    = IntRrelIV_P + IntRrelICO_P;

dydt(p.Idx.VTot)   = dydt(p.Idx.Par) - dydt(p.Idx.DTot); % calculate change in VTot based on Par and DTot

dydt(p.Idx.VRelTot)= IntRrelIV_V + IntRrelICO_V;

dydt(p.Idx.VOnlyRel) = IntRrelIV_V + IntRrelICO_V - p.Ex.kDegV*VOnlyRel;
dydt(p.Idx.DOnlyRel) = IntRrelIV_D + IntRrelICO_D - p.Ex.kDegD*DOnlyRel;

%%% mRNA in DIP-infected cells
if ( p.mRNAPrimTrans )
    kDegM         = p.In.ParaValues(strcmp('kDegM',p.In.ParaNames));
    Lm5           = p.In.ParaValues(strcmp('Lm5',p.In.ParaNames));
    Lm9           = p.In.ParaValues(strcmp('Lm9',p.In.ParaNames));

    dydt(p.Idx.Rm5_ID) =           p.Ex.kSynM/Lm5 * p.Ex.Modip - kDegM * Rm5_ID; 
    dydt(p.Idx.Rm9_ID) = p.Ex.Fm * p.Ex.kSynM/Lm9 * p.Ex.Modip - kDegM * Rm9_ID; 
end

if ( sum(isnan(dydt)) )
    error('The Euler time step "p.Ex.h" needs to be reduced for calculation!');
end
