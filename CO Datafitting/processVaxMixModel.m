% Gather output from the covidVaxMixModel

%% Set up
load('T_extend_oct9.mat')
    T=T_extend_oct9;
    T=groupedData(T);
    T.Properties.DimensionNames={'Row','Variables'};
    T.Properties.VariableNames={'days','hTotalCO'};
    T.Properties.IndependentVariableName='days';

endOfData = 1095;
step = 14; tend = 1095; %endOfData + 6*step;
t_indices =[24:14:822 830:4:868];%switch data frequency from daily to weekly occurs at entry 830
tef=T.days(t_indices)'+1;
tL = length(tef);

load('bMatrixOCT1.mat')
load('bMatrixOCT2.mat')
load('bMatrixOCT3.mat')
load('bMatrixOCT4.mat')
bMatrix = [bMatrixOCT1;bMatrixOCT2;bMatrixOCT3;bMatrixOCT4];%bMatrixOCT1_IG;%

while size(bMatrix, 2) < length(tef)
    bMatrix = [bMatrix, bMatrix(end)];
end

% grab parameter indices
pk = [(1:15)' reshape(17:61, [15 3])];
sInd = [pk(1, :); pk(9, :)]; % row 1 = S, row 2 = SV
eInd = [pk(2, :); pk(10, :)];
iInd = [pk(3, :); pk(11, :)]; 
aInd = [pk(4, :); pk(12, :)];
rInd = [pk(5, :); pk(13, :)];
hInd = [pk(7, :); pk(14, :)];
vInd = [pk(6, :); pk(15, :)]; % row 1 = V, row 2 = B
dInd = pk(8, :);
bInd = [63 64 65 68 69 70 71]; % [b bAuv bIuv bAsviv bIsviv bAsi bIsi]
NV = 67; % total vaccinated (SV + IV + ...)

%% Run simulations

epVec = fliplr(0:0.05:1); %%dont have starting eps right now
tVec = 0:tend;
outputMat = zeros(34, length(0:tend), length(epVec));
for i = 1:length(epVec)
    
    if i == 1
        m = covidVaxMixModel(tef, bMatrix(i, :), epVec(i));
        csObj = getconfigset(m, 'active');
        csObj.SolverOptions.MaxStep = 1;
        set(csObj, 'StopTime', tend, 'SolverType', 'Sundials');
        verify(m);
        sbioaccelerate(m, csObj);
        csObj = getconfigset(m, 'active');
    end
    
    tic
    [t1, x1, n] = sbiosimulate(m, csObj);
    disp('Simulation Time'), toc
    [tu1, itu] = unique(t1);
    xu1 = interp1(tu1, x1(itu, :), tVec);

    outputMat(1:8, :, i) = xu1(:, [sInd(1, :) sInd(2, :)])'; % S and SV by age
    outputMat(9:16, :, i) = xu1(:, [iInd(1, :) iInd(2, :)])'; % I and IV by age
    outputMat(17:24, :, i) = xu1(:, [aInd(1, :) aInd(2, :)])'; % A and AV by age
    outputMat(25:31, :, i) = xu1(:, bInd)'; % betas
    outputMat(32, :, i) = xu1(:, 16); % hosp total
    outputMat(33, :, i) = sum(xu1(:, hInd(1, :)), 2); % unvaccinated hosp total
    outputMat(34, :, i) = sum(xu1(:, hInd(2, :)), 2); % vaccinated hosp total
    
    % update epsilon value for next simulation
    if i < length(epVec)
        m.parameter(1).Value = epVec(i + 1);
        pVec = 3983:3:numel(m.parameter);
        for j = 1:tL
            m.parameter(pVec(j)).Value = bMatrix(i + 1, j);
        end
    end
    
end

%% Calculate exposures and make plots

expMat = zeros(4, length(tVec), length(epVec)); expSmooth = zeros(1, length(tVec));
for i = 1:length(epVec)
    
    sumS = [sum(outputMat(1:4, :, i)); sum(outputMat(5:8, :, i))];
    sumI = [sum(outputMat(9:12, :, i)); sum(outputMat(13:16, :, i))];
    sumA = [sum(outputMat(17:20, :, i)); sum(outputMat(21:24, :, i))];
    expMat(1, :, i) = outputMat(31, :, i).*sumS(1, :).*sumI(1, :) + outputMat(30, :, i).*sumS(1, :).*sumA(1, :); % S + I
    expMat(2, :, i) = outputMat(27, :, i).*sumS(1, :).*sumI(2, :) + outputMat(26, :, i).*sumS(1, :).*sumA(2, :); % S + IV
    expMat(3, :, i) = outputMat(27, :, i).*sumS(2, :).*sumI(1, :) + outputMat(26, :, i).*sumS(2, :).*sumA(1, :); % SV + I
    expMat(4, :, i) = outputMat(29, :, i).*sumS(2, :).*sumI(2, :) + outputMat(28, :, i).*sumS(2, :).*sumA(2, :); % SV + IV
    
    sw = 14;
    for j = 1:4
        for k = (sw + 1):(length(tVec) - sw)
            expSmooth(k) = mean(expMat(j, (k - sw):(k + sw), i));
        end
        expMat(j, (sw + 1):(length(tVec) - sw), i) = expSmooth((sw + 1):(length(tVec) - sw));
    end
    
    expMat(5, :, i) = sum(expMat(1:4, :, i)); % total exposures
    
end

%%
xticker = [31 60 91 121 152 182 213 244 274 305 335 366 397 425 456 486 517 547 578 609 639 670 700 731 762 790 821 851 882 912 943 973 1004 1034 1065] + 1;
xlabeler = {'2/1', '3/1', '4/1', '5/1', '6/1', '7/1', '8/1', '9/1', '10/1', '11/1', '12/1', '1/1', '2/1', '3/1', '4/1', '5/1', '6/1', '7/1', '8/1', '9/1', '10/1', '11/1', '12/1', '1/1', '2/1', '3/1', '4/1', '5/1', '6/1', '7/1','8/1','9/1','10/1','11/1','12/1'};
sDate = 336; eDate = 1070; 



%%
% % plot hospitalizations by epsilon
f1 = figure;%('Position', get(0, 'Screensize'));
set(f1, 'color', 'w');
% yyaxis left
plot(T.days, T.hTotalCO, 'ko', 'MarkerFaceColor', 'k')
hold on
for i = 2:3%length(epVec)

    plot(tVec, outputMat(32, :, i), 'Color',  [1 1 1]*0.9, 'LineWidth', 2)

end
for i = 1%:length(epVec)

    plot(tVec, outputMat(32, :, i), 'Color',  [0.8500 0.3250 0.0980]*0.9, 'LineWidth', 2)

end

xlabel('Date'), ylabel('Hospitalized Individuals')
xlim([31 eDate]), xticks(xticker), xticklabels(xlabeler)
ax = gca; ax.FontSize = 30;
f2 = getframe(f1);
rectangle('Position',[813 0 90 2500])
rectangle('Position',[990 0 52 2500])

%%
Usource=expMat(1,:,1)+expMat(3,:,1);
Vsource=expMat(2,:,1)+expMat(4,:,1);
source_count=[Usource;Vsource];
ll=length(tVec);
f1 = figure
set(f1, 'color', 'w');
hold on

bar(tVec, source_count,"stacked")

legend('Total','$\mathcal{I}_{U\rightarrow}$', '$\mathcal{I}_{V\rightarrow}$','Interpreter','latex', 'Location', 'northeast')
xlabel('Date'), ylabel('Incidence per day')
xlim([31 eDate]), %ylim([0 100]), 
xticks(xticker), xticklabels(xlabeler),
ax = gca; ax.FontSize = 30;

rectangle('Position',[813 0 90 14000])
rectangle('Position',[990 0 52 14000])
