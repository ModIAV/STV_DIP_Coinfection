 function VisualizeResults(result, d, p)

% last revised: 2021/03/31 

% close all;
%% Retrieve states
% SimTimeIn   = result.In.COI.time;
SimStatesIn = p.InToPop.states;
SimValuesIn = p.InToPop.populationvalues;
SimValuesIn(isnan(SimValuesIn)) = 0;

SimTimeEx   = result.Ex.time;
SimStatesEx = [result.Ex.states, result.Ex.variables];
SimValuesEx = [result.Ex.statevalues result.Ex.variablevalues];

for i = 1 : length(SimStatesIn)
    SimIndex.(SimStatesIn{i}) = strcmp(SimStatesIn{i}, SimStatesIn);
end
for i = 1 : length(SimStatesEx)
    SimIndex.(SimStatesEx{i}) = strcmp(SimStatesEx{i}, SimStatesEx);
end
for i = 1 : length(d.In.ExpStateNames)
    ExpIndex.(d.In.ExpStateNames{i}) = strcmp(d.In.ExpStateNames{i}, d.In.ExpStateNames);
end
for i = 1 : length(d.Ex.ExpStateNames)
    ExpIndex.(d.Ex.ExpStateNames{i}) = strcmp(d.Ex.ExpStateNames{i}, d.Ex.ExpStateNames);
end

PlotLetters = 1;

%% plotting

%%% OFFSET - extracellular level
% virus was washed from the culture after infection, however, particles
% may have gotten stuck to the cell membrane and remain in the culture
% and cause a constant offset in measurements -> implementation of such an offset for plotting

vRNANames = {'VOnlyRel', 'VTot', 'DTot'};
p.Ex.OffsetIndex = zeros(length(vRNANames), 2); 
if ( ~isempty(d.Ex.ExpStateNames) )
    for CountRnas = 1:length(vRNANames)
        p.Ex.OffsetIndex(CountRnas, 1) = find(strcmp(vRNANames(CountRnas), p.Ex.StateNames)); % or SimStatesIn
        p.Ex.OffsetIndex(CountRnas, 2) = find(strcmp(vRNANames(CountRnas), d.Ex.ExpStateNames));
    end    

    for SimRun = p.SimCond % apply to simulation results              

        AddOffset = d.Ex.ExpStateValues(2, p.Ex.OffsetIndex(:,2), SimRun);
        AddOffset(isnan(AddOffset)) = 0;

        AddOffset = AddOffset .* ~isnan(d.Ex.ExpStateValues(3, p.Ex.OffsetIndex(:,2), SimRun));

        SimValuesEx(:,p.Ex.OffsetIndex(:,1),SimRun) = SimValuesEx(:,p.Ex.OffsetIndex(:,1),SimRun) +...
            repmat(AddOffset, size(SimValuesEx,1), 1);             
    end
end       

%% extracellular level (titers, cell populations)

%%% total number of infectious + virus particles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [1 6 35 18],...
               'Name', 'MOI/MODIP over multiple passages');

h.FigH = 0.22;
h.FigW = 0.17;
Row_counter = 1;
x_pos = 1;
           
for i = 2 : 13
    axes;
%     subplot(3,4,i-1)

    h.p1 = plot(SimTimeEx, SimValuesEx(:, SimIndex.VOnlyRel, i), 'LineWidth', 1.5, 'Color', 'b');
    hold on;
    h.p2 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.VOnlyRel, i), 'bx', 'LineWidth', 1, 'MarkerSize', 6);

    h.p3 = plot(SimTimeEx, SimValuesEx(:, SimIndex.VTot, i), 'LineWidth', 1.5, 'Color', 'k');
    h.p4 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.VTot, i), 'd', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 6);          

    h.p5 = plot(SimTimeEx, SimValuesEx(:, SimIndex.DTot, i), 'LineWidth', 1.5, 'Color', [0 0.5 0]);
    h.p6 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.DTot, i), 's', 'Color', [0 0.5 0], 'LineWidth', 1, 'MarkerSize', 6);        

    set(gca, 'FontName', 'Arial', 'FontSize', 10, ...
             'XLim', [0 33], ...
             'YLim', [1e1 1e12], 'YTick', 10.^(2:2:20), 'YScale', 'log',...
             'Position', [0.11 + (x_pos-1)*0.22, 0.15 + (3-Row_counter)*0.28, h.FigW ,h.FigH]);         
         
    x_pos = x_pos + 1;
    if ( ~mod(i-1,4) )
        Row_counter = Row_counter + 1;
        x_pos = 1;
    end

    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-15, 5e6, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('number of virus\nparticles (virions/mL)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end    
    
    if ( i == 13 )
        h.p01 = plot(-1, -1, 'bx-', 'LineWidth', 1, 'MarkerSize', 6);
        h.p02 = plot(-1, -1, 'kd-', 'LineWidth', 1, 'MarkerSize', 6);
        h.p03 = plot(-1, -1, 's-', 'Color', [0 0.5 0], 'LineWidth', 1, 'MarkerSize', 6);
        h.L = legend([h.p01, h.p02, h.p03], 'infectious STV', 'total STV', 'total DIP');
        set(h.L, 'Position', [0.47 0.0001 0.1 0.1], 'Orientation', 'horizontal',...
                 'FontName', 'Arial', 'FontSize', 10);
    end    
end

%%% cell populations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-dead cells    
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [1 1 30 15],...
               'Name', 'MOI/MODIP over multiple passages');

for i = 2 : 13
    subplot(3,4,i-1); hold on;

    if ( i == 4 || i == 5 )        
        plot([33 33], [eps 1e20], 'k--', 'LineWidth', 1.5)
    end

    h.p1 = plot(SimTimeEx, SimValuesEx(:, SimIndex.T, i), 'LineWidth', 1.5, 'Color', 'b'); 
    h.p2 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.T, i), 'bx', 'LineWidth', 1, 'MarkerSize', 6);   

    h.p3 = plot(SimTimeEx, SimValuesEx(:, SimIndex.I, i), 'LineWidth', 1.5, 'Color', [0 0.5 0]);
    h.p4 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.I, i), 'o', 'Color', [0 0.5 0], 'LineWidth', 1, 'MarkerSize', 6);   

    h.p5 = plot(SimTimeEx, SimValuesEx(:, SimIndex.Ia, i), 'LineWidth', 1.5, 'Color', [0 0 0]);
    h.p6 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.Ia, i), 'd', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 6);        

    h.p7 = plot(SimTimeEx, SimValuesEx(:, SimIndex.Ta, i), 'LineWidth', 1.5, 'Color', [0.9 0.42 0]);
    h.p8 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.Ta, i), 's', 'Color', [0.9 0.42 0], 'LineWidth', 1, 'MarkerSize', 6);               

    set(gca, 'XLim', [0 33], ...
             'YLim', [0 4e6], 'YTick', 0:1e6:4e6);

    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i == 6 )
        legend([h.p1, h.p3, h.p5, h.p7], 'T', 'I', 'Ia', 'Ta', 'Location', 'NorthEast')
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-18, 1.5e6, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('cell concentration\n(cells/mL)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end    

    box on;
end  

%%% cell population fraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-dead cells       
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [10 1 30 15],...
               'Name', 'MOI/MODIP over multiple passages');

h.FigH = 0.22;
h.FigW = 0.17;
Row_counter = 1;
x_pos = 1;           
           
for i = 2 : 13

    Frac_T  = 100 * SimValuesEx(:, SimIndex.T, i)  ./ SimValuesEx(:, SimIndex.Cviab, i);
    Frac_I  = 100 * SimValuesEx(:, SimIndex.I, i)  ./ SimValuesEx(:, SimIndex.Cviab, i);
    Frac_Ta = 100 * SimValuesEx(:, SimIndex.Ta, i) ./ SimValuesEx(:, SimIndex.Cviab, i);
    Frac_Ia = 100 * SimValuesEx(:, SimIndex.Ia, i) ./ SimValuesEx(:, SimIndex.Cviab, i);           

    subplot(3,4,i-1); hold on;

    if ( i == 4 || i == 5 )        
        plot([33 33], [eps 1e20], 'k--', 'LineWidth', 1.5)
    end        

    h.p1 = plot(SimTimeEx, Frac_T, 'LineWidth', 1.5, 'Color', 'b');        
    h.p2 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.pT, i), 'bx', 'LineWidth', 1, 'MarkerSize', 6);   

    h.p3 = plot(SimTimeEx, Frac_I, 'LineWidth', 1.5, 'Color', [0 0.5 0]);
    h.p4 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.pI, i), 'd', 'Color', [0 0.5 0], 'LineWidth', 1, 'MarkerSize', 6);   

    h.p5 = plot(SimTimeEx, Frac_Ia, 'LineWidth', 1.5, 'Color', [0 0 0]);
    h.p6 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.pIa, i), 's', 'Color', [0 0 0], 'LineWidth', 1, 'MarkerSize', 6);        

    h.p7 = plot(SimTimeEx, Frac_Ta, 'LineWidth', 1.5, 'Color', [0.9 0.42 0]);
    h.p8 = plot(d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.pTa, i), 'o', 'Color', [0.9 0.42 0], 'LineWidth', 1, 'MarkerSize', 6); 

    set(gca, 'FontName', 'Arial', 'FontSize', 10, ...
             'XLim', [0 33], ...
             'YLim', [0 100], 'YTick', 0:20:100,...
             'Position', [0.11 + (x_pos-1)*0.22, 0.15 + (3-Row_counter)*0.28, h.FigW ,h.FigH]);         
         
    x_pos = x_pos + 1;
    if ( ~mod(i-1,4) )
        Row_counter = Row_counter + 1;
        x_pos = 1;
    end

    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
%             xlabel('time post infection* (h)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
    end

    if ( i < 6 )
%         title(['MOI ', num2str(p.Ex.MOIMODIPConditions(i,1))]); % , 
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-15, 50, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('fraction of cells (%%)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end   
    
    
    if ( i == 13 )
        h.p01 = plot(-1, -1, 'bx-', 'LineWidth', 1, 'MarkerSize', 6);
        h.p02 = plot(-1, -1, 'd-', 'Color', [0 0.5 0], 'LineWidth', 1, 'MarkerSize', 6);
        h.p03 = plot(-1, -1, 'ks-', 'LineWidth', 1, 'MarkerSize', 6);
        h.p04 = plot(-1, -1, 'o-', 'Color', [0.9 0.42 0], 'LineWidth', 1, 'MarkerSize', 6);
        h.L = legend([h.p01, h.p02, h.p03, h.p04], 'uninfected', 'infected', 'infected, apoptotic', 'uninfected, apoptotic');
        set(h.L, 'Position', [0.47 0.0001 0.1 0.1], 'Orientation', 'horizontal',...
                 'FontName', 'Arial', 'FontSize', 10);
    end      

    box on;
end      

%% intracellular RNAs (S5, FL, DI)

%%% OFFSET - intracellular level
% as free viral RNAs in the seed virus might adhere to cells and thus cause
% a constant offset in measurements we implement such an offset for plotting
try
    vRNANames = {'RvSeg1', 'RvSeg9', 'RvSeg5', 'Rm1', 'Rm9', 'Rm5', 'vRNPnuc1', 'vRNPnuc9', 'vRNPnuc5'};
    p.In.OffsetIndex = zeros(length(vRNANames), 2); 
    if ( ~isempty(d.In.ExpStateNames) && ~isempty(p.InToPop.states) )
        for CountRnas = 1:length(vRNANames)
            p.In.OffsetIndex(CountRnas, 1) = find(strcmp(vRNANames(CountRnas), p.InToPop.states)); % or SimStatesIn
            p.In.OffsetIndex(CountRnas, 2) = find(strcmp(vRNANames(CountRnas), d.In.ExpStateNames));
        end

        for SimRun = p.SimCond % apply to simulation results  

            AddOffset = d.In.ExpStateValues(1, p.In.OffsetIndex(:,2), SimRun);
            AddOffset(isnan(AddOffset)) = 0;
            
            SimValuesIn(:,p.In.OffsetIndex(:,1),SimRun) = SimValuesIn(:,p.In.OffsetIndex(:,1),SimRun) +...
                repmat(AddOffset, size(SimValuesIn,1), 1);             
        end
    end 
catch
    vRNANames = {'RvSeg1', 'RvSeg9', 'RvSeg5'};
    p.In.OffsetIndex = zeros(length(vRNANames), 2); 
    if ( ~isempty(d.In.ExpStateNames) && ~isempty(p.InToPop.states) )
        for CountRnas = 1:length(vRNANames)
            p.In.OffsetIndex(CountRnas, 1) = find(strcmp(vRNANames(CountRnas), p.InToPop.states)); % or SimStatesIn
            p.In.OffsetIndex(CountRnas, 2) = find(strcmp(vRNANames(CountRnas), d.In.ExpStateNames));
        end

        for SimRun = p.SimCond % apply to simulation results  

            AddOffset = d.In.ExpStateValues(1, p.In.OffsetIndex(:,2), SimRun);
            AddOffset(isnan(AddOffset)) = 0;
            
            SimValuesIn(:,p.In.OffsetIndex(:,1),SimRun) = SimValuesIn(:,p.In.OffsetIndex(:,1),SimRun) +...
                repmat(AddOffset, size(SimValuesIn,1), 1);             
        end
    end        
end

if ( sum(sum(sum(isnan(SimValuesIn)))) )
    error('Error during application of intracellular offset!');
end

%%% vRNA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [20 10 30 15],...
               'Name', 'Seg5 vRNA levels');

for i = 2 : 13
    subplot(3,4,i-1)

    % find experimental data time limit
    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.RvSeg5, i)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime);

    % vRNA
    h.p1 = plot(SimPlotTime, SimValuesIn(1:length(SimPlotTime),SimIndex.RvSeg5,i), ...
                'LineWidth', 1.5, 'Color', 'b');
    hold on;
    h.p2 = plot(d.In.Time, d.In.ExpStateValues(:, ExpIndex.RvSeg5, i), ...
                'bo', 'LineWidth', 1, 'MarkerSize', 6);    

    plot([-1 1e3],[1e0 1e0],'k-','LineWidth',1)

    set(gca, 'XLim', [0 33], ...
             'YLim', [1e-2 1e6], 'YTick', 10.^(-10:2:10), 'YScale', 'log');

    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-18, 1e2, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('vRNA levels\n(vRNA/cell)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end      

    box on;
end  

%%% FL
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [20 8 30 15],...
               'Name', 'FL vRNA levels');

for i = 2 : 13
    subplot(3,4,i-1)

    % find experimental data time limit
    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.RvSeg1, i)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime); 

    h.p1 = plot(SimPlotTime, SimValuesIn(1:length(SimPlotTime),SimIndex.RvSeg1,i), ...
                'LineWidth', 1.5, 'Color', [0 0.5 0]);
    hold on;
    h.p2 = plot(d.In.Time, d.In.ExpStateValues(:, ExpIndex.RvSeg1, i), ...
                'bo', 'LineWidth', 1, 'MarkerSize', 6, 'Color', [0 0.5 0]);

    plot([-1 1e3],[1e0 1e0],'k-','LineWidth',1)

    set(gca, 'XLim', [0 33], ...
             'YLim', [1e-2 1e6], 'YTick', 10.^(-10:2:10), 'YScale', 'log');
                   
    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-18, 1e2, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('vRNA levels\n(vRNA/cell)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end      

    box on;
end   

%%% DI
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [20 6 30 15],...
               'Name', 'DI vRNA levels');

for i = 2 : 13
    subplot(3,4,i-1)

    % find experimental data time limit
    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.RvSeg9, i)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime);  

    h.p1 = plot(SimPlotTime, SimValuesIn(1:length(SimPlotTime),SimIndex.RvSeg9,i), ...
                'LineWidth', 1.5, 'Color', [0.9 0.42 0]);
    hold on;
    h.p2 = plot(d.In.Time, d.In.ExpStateValues(:, ExpIndex.RvSeg9, i), ...
                'bo', 'LineWidth', 1, 'MarkerSize', 6, 'Color', [0.9 0.42 0]);

    plot([-1 1e3],[1e0 1e0],'k-','LineWidth',1)

    set(gca, 'XLim', [0 33], ...
             'YLim', [1e-2 1e6], 'YTick', 10.^(-10:2:10), 'YScale', 'log');
            
    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-18, 1e2, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('vRNA levels\n(vRNA/cell)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end      

    box on;
end   

%%% mRNA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [5 10 30 15],...
               'Name', 'S5 mRNA levels');

for i = 2 : 13
    subplot(3,4,i-1)

    % find experimental data time limit
    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.Rm5, i)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime);

    % vRNA
    h.p1 = plot(SimPlotTime, SimValuesIn(1:length(SimPlotTime),SimIndex.Rm5,i), ...
                'LineWidth', 1.5, 'Color', 'b');
    hold on;
    h.p2 = plot(d.In.Time, d.In.ExpStateValues(:, ExpIndex.Rm5, i), ...
                'bo', 'LineWidth', 1, 'MarkerSize', 6);   
            

    plot([-1 1e3],[1e0 1e0],'k-','LineWidth',1)

    set(gca, 'XLim', [0 33], ...
             'YLim', [1e-2 1e6], 'YTick', 10.^(-10:2:10), 'YScale', 'log');
 
    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-18, 1e2, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('mRNA levels\n(mRNA/cell)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end      

    box on;
end  

%%% FL
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [5 8 30 15],...
               'Name', 'FL mRNA levels');

for i = 2 : 13
    subplot(3,4,i-1)

    % find experimental data time limit
    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.Rm1, i)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime); 

    h.p1 = plot(SimPlotTime, SimValuesIn(1:length(SimPlotTime),SimIndex.Rm1,i), ...
                'LineWidth', 1.5, 'Color', [0 0.5 0]);
    hold on;
    h.p2 = plot(d.In.Time, d.In.ExpStateValues(:, ExpIndex.Rm1, i), ...
                'bo', 'LineWidth', 1, 'MarkerSize', 6, 'Color', [0 0.5 0]);

    plot([-1 1e3],[1e0 1e0],'k-','LineWidth',1)

    set(gca, 'XLim', [0 33], ...
             'YLim', [1e-2 1e6], 'YTick', 10.^(-10:2:10), 'YScale', 'log');
                  
    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-18, 1e2, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('mRNA levels\n(mRNA/cell)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end

    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end  
    
    box on;
end   

%%% DI
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [5 6 30 15],...
               'Name', 'DI mRNA levels');

for i = 2 : 13
    subplot(3,4,i-1)

    % find experimental data time limit
    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.Rm9, i)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime);  

    h.p1 = plot(SimPlotTime, SimValuesIn(1:length(SimPlotTime),SimIndex.Rm9,i), ...
                'LineWidth', 1.5, 'Color', [0.9 0.42 0]);
    hold on;
    h.p2 = plot(d.In.Time, d.In.ExpStateValues(:, ExpIndex.Rm9, i), ...
                'bo', 'LineWidth', 1, 'MarkerSize', 6, 'Color', [0.9 0.42 0]);

    plot([-1 1e3],[1e0 1e0],'k-','LineWidth',1)

    set(gca, 'XLim', [0 33], ...
             'YLim', [1e-2 1e6], 'YTick', 10.^(-10:2:10), 'YScale', 'log');
             
    if ( i == 4 || i == 5 )
        set(gca, 'XLim', [0 80])
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 )
        text(-18, 1e2, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        ylabel(sprintf('mRNA levels\n(mRNA/cell)'), 'FontName', 'Arial', 'FontSize', 10);  
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end
    
    if ( PlotLetters )
        X = get(gca, 'XLim'); Y = get(gca, 'YLim');
        text(-0.2*X(2), 1.049*Y(2), char(63+i), 'Color', 'k',...
             'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold'); 
    end      

    box on;
end   

return

%% vRNP export %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

%%% vRNP % in nucleus
h.Fig1 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [20 12 30 15],...
               'Name', 'vRNP % in nuc');

for i = 2 : 13
    subplot(3,4,i-1)

    % find experimental data time limit
    MaxTime = find(isnan(d.In.ExpStateValues(:, ExpIndex.vRNPnuc_DAPI, i)) == 0, 1, 'last');
    SimPlotTime = 0:p.Ex.h:d.In.Time(MaxTime);  

    h.p1 = plot(SimPlotTime, SimValuesIn(1:length(SimPlotTime),SimIndex.vRNPnuc_DAPI,i), ...
                'LineWidth', 1.5, 'Color', 'b');
    hold on;
    h.p2 = plot(d.In.Time, d.In.ExpStateValues(:, ExpIndex.vRNPnuc_DAPI, i), ...
                'o', 'LineWidth', 1, 'MarkerSize', 6, 'Color', 'b');        

    plot([-1 1e3],[1e0 1e0],'k-','LineWidth',1)

    set(gca, 'XLim', [0 7], 'XTick', 0:3:9, ...
             'YLim', [0 100], 'YTick', 0:20:100);       
                   
    if ( i < 6 )
        plot([7 7], [-10 1e6], 'k--', 'LineWidth', 1.5)
        set(gca, 'XLim', [0 21], 'XTick', 0:6:24)
        xlbl = xticklabels;
        for j=1:length(xlbl)                 % generate the new ones for those locations...
            xlbl(j)={num2str(str2double(xlbl(j)),'\\bf%d')};
        end
        xticklabels(xlbl);
    end
    if ( i < 6 )
        title(sprintf('MODIP %s', num2str(p.Ex.MOIMODIPConditions(i,2)))); 
    end
    if ( mod(i-1,4) == 1 ) % -3
        ylabel(sprintf('vRNP in nuc (%%)'), 'FontName', 'Arial', 'FontSize', 10);  
        
        if ( i < 6 )
            text(-10, 50, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
                 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        else
            text(-3.35, 50, sprintf('MOI\n%s', num2str(p.Ex.MOIMODIPConditions(i,1))), ...
                 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
        end
        
    end
    if ( i > 9 )
        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
    end

    box on;
end 

%%% comparison like in excel file
h.Fig2 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [20 1 25 8],...
               'Name', 'Misc. stats');

vRNA_Names = {'FL vRNA', 'DI vRNA', 'S5 vRNA'};    
colors = [0 0 0; 0 0.5 0; 0.9 0.42 0; 0 0 1];
symbols = {'s', 'o', 'd', 'p'};
linestyles = {'-','-.','--',':'};

for moi = 1 : 3

    subplot(1,3,4-moi)

    for modip = 1 : 4    

        IC_num = (moi-1)*4 + modip + 1;

        % find experimental data time limit
        MaxTime = find(isnan(d.Ex.ExpStateValues(:, ExpIndex.VTot, IC_num)) == 0, 1, 'last');
        SimPlotTime = 0:p.Ex.h:d.Ex.Time(MaxTime); 

        hold on;

        h.p5 = plot(d.Ex.Time, d.In.ExpStateValues(:, ExpIndex.vRNPnuc_DAPI, IC_num), ...
                    sprintf('%s',symbols{modip}), 'Color', colors(modip,:), ...
                    'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', colors(modip,:));          

        h.(['m', num2str(modip)]) = plot(SimPlotTime,  SimValuesIn(1:length(SimPlotTime),SimIndex.vRNPnuc_DAPI, IC_num), ...
                    linestyles{modip},'LineWidth', 2, 'Color', colors(modip,:)); 

    end

    set(gca, 'FontName', 'Arial', 'FontSize', 12) 
    set(gca, 'YLim', [0 100], 'YTick', 0:20:100)

    if ( moi > 1 )
        set(gca, 'XLim', [0 24], 'XTick', 0:6:48);
    else
        plot([24 24], [eps 1e50], 'k--', 'LineWidth', 2);
        set(gca, 'XLim', [0 63], 'XTick', 0:12:120)
    end

    xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 12);
    ylabel(sprintf('vRNP in nuc (%%)'), 'FontName', 'Arial', 'FontSize', 12);  

    title(sprintf('%s', sprintf('MOI %s', num2str(p.Ex.MOIMODIPConditions(IC_num,1)))), 'FontName', 'Arial', 'FontSize', 12);             

    if ( moi == 3 )
        legend([h.m1, h.m2, h.m3, h.m4], 'MODIP 0', 'MODIP 0.001', 'MODIP 3', 'MODIP 30','Location', 'NE')
    end

    box on;


end

return
