%% Figure 7a: Outbreak Type by assortivity and vaccination coverage
eps = 0.7; ppV_range = 0:0.001:1; pImm = 0.3; tspan=[0 500];
tau_range=0:1:500;

contour_val=zeros(length(ppV_range),length(tau_range));

for i=1:1:length(ppV_range)
    for j=1:1:length(tau_range)
        ppV=ppV_range(i); 
        tau=tau_range(j);
        [t, y, pVec] = seir_simulator2(eps, ppV, pImm,tau,tspan);
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
        re=totalExp(2:end)./totalInf(1:end-1)/pVec(4);
        if re(1)>1
            ind_int=0;
            ind_re=min(find(re<1,1));
        else
            ind_int=min(max(find(re>1,1)));
            ind_re=min(find(re(ind_int+1:end)<1,1));
        end
        if isempty(ind_re)==1
            day_re=NaN;
        else
            day_re=t(ind_re+ind_int);
        end
        if isnan(day_50)
            contour_val(i,j)=0;
        else
            diff=day_re-day_50;
            if diff>0
                contour_val(i,j)=1;
            else
                contour_val(i,j)=-1;
            end
        end
    end
end

%% plot
f1 = figure;
set(f1, 'color', 'w');
hold on

contourf(100*ppV_range,tau_range,contour_val','LineWidth',2)
map = [
    0.7 0.45 0.88
    0 0.45 0.74 
    0.85 0.33 0.1
    ];
colormap(map)
xlabel('% Vaccinated'), xlim([0,100])
ylabel('\tau_b'), ylim([0,500])
ax = gca; ax.FontSize = 30;