function [p, t, y, I_V_tot, I_V_SIE, I_CO, I, InfIC, eFPar] = ExLev_Euler(p)
%Euler solves the equation system in the function file specified by 
%   the handle p.Ex.Model at the time points p.Ex.tspan using the step size
%   p.Ex.h and the intial condition p.Ex.y0
%
%   last revised: 2021/03/09

% allocate result collection vectors
NumSteps = length(p.Ex.tspan)-1;

I_V      = [p.Ex.IV_0;  zeros(NumSteps,1)];
I_V_SIE  = zeros(NumSteps+1,1);
I_V_tot  = [p.Ex.IV_0;  zeros(NumSteps,1)];
I_CO     = [p.Ex.ICO_0; zeros(NumSteps,1)];
I        = [p.Ex.IV_0 + p.Ex.ICO_0 + p.Ex.y0(p.Idx.I_D); zeros(NumSteps,1)];
y        = [p.Ex.y0;    zeros(NumSteps,length(p.Ex.y0))];
eFPar    = [0.26;  zeros(NumSteps,1)];

% calculate release rates for cell already infected at cultivation start (e.g. in high MOI conditions)
p = ExLev_EffectiveMOI(p, p.Ex.y0, I_V(1), 0); % Iter = 0

% initial values for the integrals
IntI     = [I_V(1), I_V_SIE(1), I_CO(1)];

IntRrelI = [IntI(1)*p.DynMOI.STV.RrelV{1,1}, IntI(1)*p.DynMOI.STV.RrelD{1,1}, ...
            IntI(1)*p.DynMOI.STV.RrelDTot{1,1}, IntI(1)*p.DynMOI.STV.RrelP{1,1}, ...
            IntI(2)*p.DynMOI.COI.RrelV{1,1}, IntI(2)*p.DynMOI.COI.RrelD{1,1}, ...
            IntI(2)*p.DynMOI.COI.RrelDTot{1,1}, IntI(2)*p.DynMOI.COI.RrelP{1,1}];

IntApoI  = [IntI(1)*p.Ex.ApoI(1), IntI(3)*p.Ex.ApoI(1)];

% cell that contains the concentration of I_V and I_CO cells for all infection ages
RowMatI = [IntI(1), IntI(3)];
p.Ex.FullDistribution{1} = RowMatI;

y(1,p.Idx.Cviab) = sum(y(1,p.Idx.T|p.Idx.Ta|p.Idx.Ia|p.Idx.I_D)) + sum(IntI);

% matrices that contain the initial conditions for cells infected at specific times
InfIC = cell(1,2);
InfIC{1} = [[p.Ex.y0(p.Idx.V)/p.Ex.y0(p.Idx.T) 0 0 0];...
            zeros(length(p.Ex.tspan)-1,4)];  % Vex, VAttHi, VAttLo, Ven during cell infection
InfIC{2} = [[p.Ex.y0(p.Idx.V)/p.Ex.y0(p.Idx.T) 0 0 0 p.Ex.y0(p.Idx.D)/p.Ex.y0(p.Idx.T) 0 0 0];...
            zeros(length(p.Ex.tspan)-1,8)];  % V|Dex, V|DAttHi, V|DAttLo, V|Den during cell infection

h_ReductionFactor = p.Ex.h / p.Ex.hs;        % factor describing difference btw. big and small step
h_save = p.Ex.h;
        
% if ( strcmp('simulate',p.Task) );  fprintf('Simulating ...   0 %%');  end 

for Iter = 1 : NumSteps 

%     % display simulation progress
%     if ( mod(Iter,50) == 0 && strcmp('simulate',p.Task) )
%         fprintf('\b\b\b\b\b%3.0f %%',Iter/NumSteps*100);
%     end
    
    %%% calculate dynamics using the smaller time step %%%%%%%%%%%%%%%%%%%%
    p.Ex.h = p.Ex.hs;

    Iter_ss = Iter;
    y_ss    = [y(1:Iter,:); zeros(h_ReductionFactor,size(y,2))]; 
    MatI_ss = RowMatI;

    for i = 1 : h_ReductionFactor  
        
        % calculate the dynamics of the extracelullar system
        dydt = feval(p.Ex.ModelName, (Iter-1)*h_save, y_ss(Iter_ss,:), IntRrelI, IntI, IntApoI, p); 

        y_ss(Iter_ss+1,:) = max(zeros(size(p.Ex.y0)), y_ss(Iter_ss,:) + p.Ex.h*dydt);
        y_ss(Iter_ss+1, y_ss(Iter_ss+1,:)<1e-17) = 0;

        if ( y_ss(Iter_ss+1, p.Idx.I_D) < 1e-5 )
            y_ss(Iter_ss+1, p.Idx.I_D) = 0;
        end            
        
        [IntRrelI, IntI, IntApoI, RowMatI] = ExLev_Integrals(Iter_ss+1, ...
                                    y_ss(Iter_ss,:), MatI_ss, p, h_ReductionFactor, Iter); % for newly created infected cell slice: use 0 release at initiation time point

        MatI_ss = RowMatI;    

        Iter_ss = Iter_ss + 1;   

        % back up state values for infection calculation in case everything was infected in first time step
        if ( Iter == 1 )
            if ( i == 1 )
                p.Ex.IntI_FirstWave = [I_V(1), I_V_SIE(1), I_CO(1); IntI];
            else
                p.Ex.IntI_FirstWave = [p.Ex.IntI_FirstWave; IntI];                
            end
            
            if ( i == h_ReductionFactor )
                p.Ex.y_FirstWave = y_ss;
            end
        end

    end 
    
    y(Iter+1,:) = y_ss(end,:);
    RowMatI = [RowMatI(1:Iter,:); sum(RowMatI(Iter+1:end,:),1)];    

    p.Ex.h = h_save;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % remove values in release vectors when cell reached intracellular release limit
    if ( Iter > p.In.SimTime/p.Ex.h )
        p.DynMOI.STV.RrelV{Iter - p.In.SimTime/p.Ex.h} = [];
        p.DynMOI.STV.RrelD{Iter - p.In.SimTime/p.Ex.h} = [];
        p.DynMOI.STV.RrelDTot{Iter - p.In.SimTime/p.Ex.h} = [];
        p.DynMOI.STV.RrelP{Iter - p.In.SimTime/p.Ex.h} = [];

        p.DynMOI.COI.RrelV{Iter - p.In.SimTime/p.Ex.h} = [];
        p.DynMOI.COI.RrelD{Iter - p.In.SimTime/p.Ex.h} = [];
        p.DynMOI.COI.RrelDTot{Iter - p.In.SimTime/p.Ex.h} = [];
        p.DynMOI.COI.RrelP{Iter - p.In.SimTime/p.Ex.h} = [];
    end
    
    I_V(Iter+1)     = IntI(1);
    I_V_SIE(Iter+1) = IntI(2);
    I_CO(Iter+1)    = IntI(3);  
    
    I_V_tot(Iter+1) = IntI(1) + IntI(2);
    I(Iter+1)       = IntI(1) + IntI(2) + IntI(3) + y(Iter+1, p.Idx.I_D);

    y(Iter+1,p.Idx.Cviab) = sum(y(Iter+1,p.Idx.T|p.Idx.Ta|p.Idx.Ia)) + I(Iter+1);

    p.Ex.FullDistribution{Iter+1} = RowMatI;

    IntRrelIV_V  = IntRrelI(1);  IntRrelICO_V = IntRrelI(4);
    IntRrelIV_P  = IntRrelI(3);  IntRrelICO_P = IntRrelI(6);
    
    if ( y(Iter+1,p.Idx.VRelTot) >= 1 )
        eFPar(Iter+1) = (IntRrelIV_V + IntRrelICO_V) / (IntRrelIV_P + IntRrelICO_P);
    else
        eFPar(Iter+1) = NaN;
    end   
    
    %% calculate virus particle release rates for new effective MOI
    %  (when in the previous step infectable material was still present)
    if ( (y(Iter,p.Idx.T) > 0   && y(Iter,p.Idx.VEn) > 0) || ...
         (y(Iter,p.Idx.I_D) > 0 && y(Iter,p.Idx.VEn) > 0) || ...
         (I_V(Iter) > 0         && y(Iter,p.Idx.DEn) > 0) || ...
         Iter < 3 )     

        % MOI conditions for next step are only calculated according to the large step-size
        p = ExLev_EffectiveMOI(p, y(Iter,:), I_V(Iter), Iter);
    end

    %% prevent overloading of high-affinity binding sites (e.g. due to numerical inaccuracies)
    ParRatio = y(Iter+1,p.Idx.V)/(y(Iter+1,p.Idx.V) + y(Iter+1,p.Idx.D));
    if ( p.Ex.BtotHi*sum(y(Iter+1,p.Idx.T|p.Idx.Ta|p.Idx.I_D)) - y(Iter+1,p.Idx.VAttHi) < 0 )
        %save attached virus concentration
        a = y(Iter+1, p.Idx.VAttHi);    
        %calculate attached virus concentration that saturates all binding sites  
        y(Iter+1, p.Idx.VAttHi) = ParRatio * p.Ex.BtotHi*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)) ...
                                  + p.Ex.BtotHi*sum(y(Iter+1,p.Idx.I_D));                           
        %transfer the excessive virions back to the non-attached state
        y(Iter+1, p.Idx.V)      = y(Iter+1, p.Idx.V) + a - y(Iter+1, p.Idx.VAttHi);            
    end    

    if ( p.Ex.BtotHi*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)+I_V(Iter+1)) - y(Iter+1,p.Idx.DAttHi) < 0 )
        %save attached virus concentration
        a = y(Iter+1, p.Idx.DAttHi);
        %calculate attached virus concentration that saturates all binding sites        
        y(Iter+1, p.Idx.DAttHi) = (1-ParRatio) * p.Ex.BtotHi*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)) ...
                                  + p.Ex.BtotHi*I_V(Iter+1);
        %transfer the excessive virions back to the non-attached state
        y(Iter+1, p.Idx.D)      = y(Iter+1, p.Idx.D) + a - y(Iter+1, p.Idx.DAttHi); 
    end

    %prevent overloading of low-affinity binding sites (e.g. due to numerical inaccuracies)
    if ( p.Ex.BtotLo*sum(y(Iter+1,p.Idx.T|p.Idx.Ta|p.Idx.I_D)) - y(Iter+1,p.Idx.VAttLo) < 0 )
        %save attached virus concentration
        a = y(Iter+1, p.Idx.VAttLo);
        %calculate attached virus concentration that saturates all binding sites  
        y(Iter+1, p.Idx.VAttLo) = ParRatio * p.Ex.BtotLo*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)) ...
                                  + p.Ex.BtotLo*sum(y(Iter+1,p.Idx.I_D));                           
        %transfer the excessive virions back to the non-attached state
        y(Iter+1, p.Idx.V)      = y(Iter+1, p.Idx.V) + a - y(Iter+1, p.Idx.VAttLo);            
    end       

    if ( p.Ex.BtotLo*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)+I_V(Iter+1)) - y(Iter+1,p.Idx.DAttLo) < 0 )
        %save attached virus concentration
        a = y(Iter+1, p.Idx.DAttLo);
        %calculate attached virus concentration that saturates all binding sites  
        y(Iter+1, p.Idx.DAttLo) = (1-ParRatio) * p.Ex.BtotLo*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)) ...
                                  + p.Ex.BtotLo* I_V(Iter+1);                        
        %transfer the excessive virions back to the non-attached state
        y(Iter+1, p.Idx.D)      = y(Iter+1, p.Idx.D) + a - y(Iter+1, p.Idx.DAttLo);            
    end  

    if ( p.Ex.BtotHi*sum(y(Iter+1,p.Idx.T|p.Idx.Ta|p.Idx.I_D)+I_V(Iter+1)) ...
         - y(Iter+1,p.Idx.VAttHi) - y(Iter+1,p.Idx.DAttHi) < 0 )  
     
        if ( y(Iter+1,p.Idx.VAttHi) + y(Iter+1,p.Idx.DAttHi) < 1e-3 )
            % if error occurs when concentration are very low -> set them to 0
            y(Iter+1,p.Idx.VAttHi) = 0;
            y(Iter+1,p.Idx.DAttHi) = 0;
        else
            fprintf('%.2f\n',(Iter+1)*p.Ex.h);   
            error('Binding site overload!');            
        end
        
    end    

    %% check results for negative values
    if ( sum(y(Iter+1,:) < 0) )
        fprintf('%.2f\n',(Iter+1)*p.Ex.h);
        error('Negative results!');
    end    

end

t = p.Ex.tspan;


























