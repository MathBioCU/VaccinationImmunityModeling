%% Figure 5b: Contour of assortivity vs time impacting outbreak
eps_range = 0:0.001:1; ppV_range = 0:0.001:1; pImm = 0.3; tspan=[0 500];

contour_val=zeros(length(ppV_range),length(eps_range));

for i=1:1:length(ppV_range)
    for j=1:1:length(eps_range)
        ppV=ppV_range(i); 
        eps=eps_range(j);
        [t, y, pVec] = seir_simulator(eps, ppV, pImm,tspan);
        totalInf = y(:, 3) + y(:, 8);
        expSuIu = pVec(2)*(1 - pVec(10))*y(:, 1).*(y(:, 3)/pVec(9)) + pVec(2)*pVec(10)*y(:, 1).*(y(:, 3)/pVec(1));
        expSuIv = pVec(2)*pVec(10)*pVec(11)*y(:, 1).*y(:, 8)/pVec(1);
        expSvIu = pVec(2)*pVec(10)*y(:, 6).*y(:, 3)/pVec(1);
        expSvIv = pVec(2)*pVec(11)*(1 - pVec(10))*y(:, 6).*(y(:, 8)/pVec(8)) + pVec(2)*pVec(11)*pVec(10)*y(:, 6).*(y(:, 8)/pVec(1));
        totalExp = expSuIu + expSuIv + expSvIu + expSvIv;
        g50ExpV=(expSuIv+expSvIv)./totalExp;
        ind_50=find(g50ExpV>0.5,1);
        if isempty(ind_50)==1
            day_50=NaN;
        else
            day_50=t(ind_50);
        end
        if isnan(day_50)
            contour_val(i,j)=500;
        else
            contour_val(i,j)=day_50;
        end
    end
end

%% plot
f1 = figure;
set(f1, 'color', 'w');
hold on

contourf(100*ppV_range,100*eps_range,contour_val','LineWidth',2)
colorbar()
colormap(turbo)
xlabel('% Vaccinated'), xlim([0,100])
ylabel('Assortivity'), ylim([0,100])
ax = gca; ax.FontSize = 48;
