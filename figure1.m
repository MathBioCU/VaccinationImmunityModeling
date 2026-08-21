%% Figure 1: Simulation of incidence by vax status
eps = 0.7; ppV = 0.85; pImm = 0.3; tspan=[0 500];

    [t, y, pVec] = seir_simulator(eps, ppV, pImm,tspan);
    totalInf = y(:, 3) + y(:, 8);
    uvInf = y(:, 3);
    vInf = y(:, 8);
    expSuIu = pVec(2)*(1 - pVec(10))*y(:, 1).*(y(:, 3)/pVec(9)) + pVec(2)*pVec(10)*y(:, 1).*(y(:, 3)/pVec(1));
    expSuIv = pVec(2)*pVec(10)*pVec(11)*y(:, 1).*y(:, 8)/pVec(1);
    expSvIu = pVec(2)*pVec(10)*y(:, 6).*y(:, 3)/pVec(1);
    expSvIv = pVec(2)*pVec(11)*(1 - pVec(10))*y(:, 6).*(y(:, 8)/pVec(8)) + pVec(2)*pVec(11)*pVec(10)*y(:, 6).*(y(:, 8)/pVec(1));
    totalExp = expSuIu + expSuIv + expSvIu + expSvIv;

    results = [t totalInf uvInf vInf expSuIu expSuIv expSvIu expSvIv totalExp];
        g50ExpV=(expSuIv+expSvIv)./totalExp;
        ind_50=find(g50ExpV>0.5,1);
        if isempty(ind_50)==1
            day_50_2=NaN;
        else
            day_50_2=t(ind_50);
        end
        re=totalExp(2:end)./totalInf(1:end-1)/pVec(4);
        ind_re=min(find(re<1,1));
        if isempty(ind_re)==1
            day_re_2=NaN;
        else
            day_re_2=t(ind_re);
        end
        rect2=[min(day_re_2,day_50_2) 0 abs(day_50_2-day_re_2) max(totalExp)];
        rect3=[min(day_re_2,day_50_2) 0 abs(day_50_2-day_re_2) 4];

%% plot
f1 = figure;
set(f1, 'color', 'w');
Umat=rMat(:, 5)+rMat(:,7);
Vmat=rMat(:, 6)+rMat(:,8);
Imat=[Umat'; Vmat'];
rMat = results;
subplot(3,4,[1 2 3 4 5 6 7 8])
hold on
plot(rMat(:, 1), rMat(:, 9), 'k', 'LineWidth', 2)
bar(rMat(:, 1)', Imat,"stacked")%, 'FaceColor', [0 0.4470 0.7410])
% plot(rMat(:, 1), rMat(:, 5)+rMat(:,7), 'Color', [0 0.4470 0.7410], 'LineWidth', 2,'LineStyle','--')
% plot(rMat(:, 1), rMat(:, 6)+rMat(:,8), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2,'LineStyle',':')
rectangle('Position',rect2,'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.5)
legend('Total Incidence','Unvaccinated-source Incidence', 'Vaccinated-source Incidence', 'Location', 'Northeast')
ylabel('Incidence per day') 
xlim([50,150])
ax = gca; ax.FontSize = 20;

subplot(3,4,[9 10 11 12])
hold on
plot(rMat(1:end-1, 1), rMat(2:end,9)./rMat(1:end-1,2)/pVec(4),'k','LineWidth', 2)
plot(rMat(1:end-1, 1),ones(1,length(rMat(1:end-1, 1))),'Color',[1 1 1]*0.5,'LineStyle','--')
rectangle('Position',rect3,'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.5)
alpha(0.5)
ylabel('R_e'), ylim([0, 4])
xlabel('Time (days)'), xlim([50,150])
ax = gca; ax.FontSize = 20;