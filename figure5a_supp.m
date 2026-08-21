%% Figure 5a: Proportion of incidence per day by source vaccination status

eps = [0.25 0.5 0.75 1]; ppV = 0.85; pImm = 0.3; tspan=[0 200];
results = cell(1, 2);
for i = 1:4
    [t, y, pVec] = seir_simulator(eps(i), ppV, pImm,tspan);
    totalInf = y(:, 3) + y(:, 8);
    uvInf = y(:, 3);
    vInf = y(:, 8);
    
    expSuIu = pVec(2)*(1 - pVec(10))*y(:, 1).*(y(:, 3)/pVec(9)) + pVec(2)*pVec(10)*y(:, 1).*(y(:, 3)/pVec(1));
    expSuIv = pVec(2)*pVec(10)*y(:, 1).*y(:, 8)/pVec(1);
    expSvIu = pVec(2)*pVec(10)*y(:, 6).*y(:, 3)/pVec(1);
    expSvIv = pVec(2)*(1 - pVec(10))*y(:, 6).*(y(:, 8)/pVec(8)) + pVec(2)*pVec(10)*y(:, 6).*(y(:, 8)/pVec(1));
    totalExp = expSuIu + expSuIv + expSvIu + expSvIv;

    cell{i} = [t totalInf uvInf vInf expSuIu expSuIv expSvIu expSvIv totalExp];
    
end

%%
f1 = figure;
set(f1, 'color', 'w');

hold on
rMat = cell{4}; % eps = 1 
%plot(rMat(:, 1), 100*sum(rMat(:, [5 7]), 2)./rMat(:, 9), 'Color', [0 0.4470 0.7410], 'LineWidth', 1)
plot(rMat(:, 1), 100*sum(rMat(:, [6 8]), 2)./rMat(:, 9),'Color',  [0.8500 0.3250 0.0980], 'LineWidth', 1)
rMat = cell{3}; % eps = 0.75
%plot(rMat(:, 1), 100*sum(rMat(:, [5 7]), 2)./rMat(:, 9), 'Color', [0 0.4470 0.7410], 'LineStyle', '-.', 'LineWidth', 1)
plot(rMat(:, 1), 100*sum(rMat(:, [6 8]), 2)./rMat(:, 9), 'Color',  [0.8500 0.3250 0.0980],'LineStyle', '-.', 'LineWidth', 1)
rMat = cell{2}; % eps = 0.5 
%plot(rMat(:, 1), 100*sum(rMat(:, [5 7]), 2)./rMat(:, 9), 'Color', [0 0.4470 0.7410],'LineStyle', ':', 'LineWidth', 1)
plot(rMat(:, 1), 100*sum(rMat(:, [6 8]), 2)./rMat(:, 9),'Color',  [0.8500 0.3250 0.0980],'LineStyle', ':', 'LineWidth', 1)
rMat = cell{1}; % eps = 0.25 
%plot(rMat(:, 1), 100*sum(rMat(:, [5 7]), 2)./rMat(:, 9), 'Color', [0 0.4470 0.7410],'LineStyle', '--', 'LineWidth', 1)
plot(rMat(:, 1), 100*sum(rMat(:, [6 8]), 2)./rMat(:, 9),'Color',  [0.8500 0.3250 0.0980],'LineStyle', '--', 'LineWidth', 1)
plot(rMat(:,1), 50*ones(length(rMat(:,1))), 'k-','LineWidth',2)

legend('VsI, \epsilon = 1', ...
    'VsI, \epsilon = 0.75', ...
    'VsI, \epsilon = 0.5', ...
    'VsI, \epsilon = 0.25', 'Location', 'East')
xlabel('Time (days)'), ylabel('% Incidence per day'), ylim([0, 100])
ax = gca; ax.FontSize = 48;

f2 = getframe(f1);