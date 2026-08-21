%% Figure 4: Total cases for varying initial populatio immunity
eps = 0.7; ppV_range = 0:0.01:1; pImm_range = 0:0.01:1; tspan=[0 200];
% 
contour_val=zeros(length(ppV_range),length(pImm_range));

for i=1:1:length(ppV_range)
    for j=1:1:length(pImm_range)
        ppV=ppV_range(i); 
        pImm=pImm_range(j);
        [t, y, pVec] = seir_simulator(eps, ppV, pImm,tspan);
         cumm_cases=sum(y(:,3)+y(:,8));
        contour_val(i,j)=cumm_cases;
    end
end
%% plot
f1 = figure;
set(f1, 'color', 'w');

contourf(100*ppV_range,100*pImm_range,contour_val')
colorbar()
colormap(turbo)
xlabel('% Vaccinated'), xlim([0,100])
ylabel('% Infection Acquired Immunity'), ylim([0,100])
title('Total Cases')
ax = gca; ax.FontSize = 20;

f2 = getframe(f1);