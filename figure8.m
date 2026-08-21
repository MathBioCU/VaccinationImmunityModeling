%% Figure 8: Simulation of incidence by vax status for Type III with various masking implementations

eps = 0.7; ppV = 0.85; pImm = 0.3; tspan=[0 500];
maskeff=0.5;

%% No Masking
mask_u=0; mask_v=0;

    [t, y, pVec] = seir_simulator_masking(eps, ppV, pImm,tspan,mask_u,mask_v);
    totalInf = y(:, 3) + y(:, 8);
    uvInf = y(:, 3);
    vInf = y(:, 8);
    expSuIu = (1-pVec(11))*(pVec(2)*(1 - pVec(10))*y(:, 1).*(y(:, 3)/pVec(9)) + pVec(2)*pVec(10)*y(:, 1).*(y(:, 3)/pVec(1)));
    expSuIv = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 1).*y(:, 8)/pVec(1);
    expSvIu = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 6).*y(:, 3)/pVec(1);
    expSvIv = (1-pVec(12))*(pVec(2)*(1 - pVec(10))*y(:, 6).*(y(:, 8)/pVec(8)) + pVec(2)*pVec(10)*y(:, 6).*(y(:, 8)/pVec(1)));
    totalExp = expSuIu + expSuIv + expSvIu + expSvIv;

    results = [t totalInf uvInf vInf expSuIu expSuIv expSvIu expSvIv totalExp];

%plot
figure(1);
hold on
%peak
ind_pk=find(totalExp==max(totalExp)); Ipk=totalExp(ind_pk); tpk=t(ind_pk);

rMat = results;
hold on
plot(rMat(:, 1), rMat(:, 9), 'k', 'LineWidth', 2)
plot(rMat(:, 1), rMat(:, 5)+rMat(:,7), 'Color', [0 0.4470 0.7410], 'LineWidth', 2,'LineStyle','--')
plot(rMat(:, 1), rMat(:, 6)+rMat(:,8), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2,'LineStyle',':')
plot(tpk,Ipk,'ko')
legend('Total','UsI', 'VsI', 'Location', 'Northeast')
ylabel('Incidence per day') 
xlabel('Time (days)')
xlim([0 200])
ylim([0 2000])
ax = gca; ax.FontSize = 48;

%% Unvax Only Mask
mask_u=maskeff; mask_v=0;

    [t, y, pVec] = seir_simulator_masking(eps, ppV, pImm,tspan,mask_u,mask_v);
    totalInf = y(:, 3) + y(:, 8);
    uvInf = y(:, 3);
    vInf = y(:, 8);
    expSuIu = (1-pVec(11))*(pVec(2)*(1 - pVec(10))*y(:, 1).*(y(:, 3)/pVec(9)) + pVec(2)*pVec(10)*y(:, 1).*(y(:, 3)/pVec(1)));
    expSuIv = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 1).*y(:, 8)/pVec(1);
    expSvIu = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 6).*y(:, 3)/pVec(1);
    expSvIv = (1-pVec(12))*(pVec(2)*(1 - pVec(10))*y(:, 6).*(y(:, 8)/pVec(8)) + pVec(2)*pVec(10)*y(:, 6).*(y(:, 8)/pVec(1)));
    totalExp = expSuIu + expSuIv + expSvIu + expSvIv;

    results = [t totalInf uvInf vInf expSuIu expSuIv expSvIu expSvIv totalExp];

%plot
figure(2);
hold on

rMat = results;
hold on
plot(rMat(:, 1), rMat(:, 9), 'k', 'LineWidth', 2)
plot(rMat(:, 1), rMat(:, 5)+rMat(:,7), 'Color', [0 0.4470 0.7410], 'LineWidth', 2,'LineStyle','--')
plot(rMat(:, 1), rMat(:, 6)+rMat(:,8), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2,'LineStyle',':')
plot(tpk,Ipk,'ko')
ylabel('Incidence per day') 
xlabel('Time (days)')
xlim([0 200])
ylim([0 2000])
ax = gca; ax.FontSize = 48;

%% Vax Only Mask
mask_u=0; mask_v=maskeff;

    [t, y, pVec] = seir_simulator_masking(eps, ppV, pImm,tspan,mask_u,mask_v);
    totalInf = y(:, 3) + y(:, 8);
    uvInf = y(:, 3);
    vInf = y(:, 8);
    expSuIu = (1-pVec(11))*(pVec(2)*(1 - pVec(10))*y(:, 1).*(y(:, 3)/pVec(9)) + pVec(2)*pVec(10)*y(:, 1).*(y(:, 3)/pVec(1)));
    expSuIv = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 1).*y(:, 8)/pVec(1);
    expSvIu = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 6).*y(:, 3)/pVec(1);
    expSvIv = (1-pVec(12))*(pVec(2)*(1 - pVec(10))*y(:, 6).*(y(:, 8)/pVec(8)) + pVec(2)*pVec(10)*y(:, 6).*(y(:, 8)/pVec(1)));
    totalExp = expSuIu + expSuIv + expSvIu + expSvIv;

    results = [t totalInf uvInf vInf expSuIu expSuIv expSvIu expSvIv totalExp];

%plot
figure(3);
hold on

rMat = results;
hold on
plot(rMat(:, 1), rMat(:, 9), 'k', 'LineWidth', 2)
plot(rMat(:, 1), rMat(:, 5)+rMat(:,7), 'Color', [0 0.4470 0.7410], 'LineWidth', 2,'LineStyle','--')
plot(rMat(:, 1), rMat(:, 6)+rMat(:,8), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2,'LineStyle',':')
plot(tpk,Ipk,'ko')
ylabel('Incidence per day') 
xlabel('Time (days)')
xlim([0 200])
ylim([0 2000])
ax = gca; ax.FontSize = 48;

%% All Masking
mask_u=maskeff; mask_v=maskeff;

    [t, y, pVec] = seir_simulator_masking(eps, ppV, pImm,tspan,mask_u,mask_v);
    totalInf = y(:, 3) + y(:, 8);
    uvInf = y(:, 3);
    vInf = y(:, 8);
    expSuIu = (1-pVec(11))*(pVec(2)*(1 - pVec(10))*y(:, 1).*(y(:, 3)/pVec(9)) + pVec(2)*pVec(10)*y(:, 1).*(y(:, 3)/pVec(1)));
    expSuIv = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 1).*y(:, 8)/pVec(1);
    expSvIu = ((1-0.5*pVec(11)-0.5*pVec(12)))*pVec(2)*pVec(10)*y(:, 6).*y(:, 3)/pVec(1);
    expSvIv = (1-pVec(12))*(pVec(2)*(1 - pVec(10))*y(:, 6).*(y(:, 8)/pVec(8)) + pVec(2)*pVec(10)*y(:, 6).*(y(:, 8)/pVec(1)));
    totalExp = expSuIu + expSuIv + expSvIu + expSvIv;

    results = [t totalInf uvInf vInf expSuIu expSuIv expSvIu expSvIv totalExp];

%plot
figure(4);
hold on

rMat = results;
hold on
plot(rMat(:, 1), rMat(:, 9), 'k', 'LineWidth', 2)
plot(rMat(:, 1), rMat(:, 5)+rMat(:,7), 'Color', [0 0.4470 0.7410], 'LineWidth', 2,'LineStyle','--')
plot(rMat(:, 1), rMat(:, 6)+rMat(:,8), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2,'LineStyle',':')
plot(tpk,Ipk,'ko')
ylabel('Incidence per day') 
xlabel('Time (days)')
xlim([100 300])
ylim([0 2000])
ax = gca; ax.FontSize = 48;
