function [t, y, pVec] = seir_simulator_masking_rand(eps, ppV, pImm,tspan,mask_u,mask_v,rand_p)

    N = 100000; % population size
    beta = 0.75; % transmission rate
    alpha = 4; % latency period
    gamma = 1/9; % 1/(infectious period)
    tau1 = 365; % duration of prior-infection acquired immunity (in days)
    tau2 = tau1; % duration of vaccine plus prior-infection acquired immunity (in days)
    tau3 = 213; % duration of vaccine acquired immunity (in days)
%     ppV = 0.7; % percent of population vaccinated
%     eps = 0.5; % degree of associativity (random mixing eps = 1)
%     pImm = 0.3; % percent of population with prior-infection immunity initially
    pVec = [N beta alpha gamma tau1 tau2 tau3 ppV*N (1 - ppV)*N eps mask_u, mask_v, rand_p];
    
    
    y0 = [((1 - pImm)*pVec(9) - 1); 0; 1; pImm*pVec(9); (1 - pImm)*pVec(8); 0; 0; 0; pImm*pVec(8)]; 
%     y0 = [(pVec(9) - 1); 0; 1; 0; pVec(8); 0; 0; 0; 0];
    options = odeset('MaxStep', 1);
%     [t, y] = ode45(@seir2, tspan, y0, [], pVec);
    [t, y] = ode45(@(t, y) seir2(t, y, pVec), tspan, y0, options);
    
%     figure
%     plot(t, y(:, 1), t, y(:, 2), t, y(:, 3), t, y(:, 4))
%     legend('S', 'E', 'I', 'R')
%     
%     figure
%     plot(t, y(:, 5), t, y(:, 6), t, y(:, 7), t, y(:, 8), t, y(:, 9))
%     legend('V', 'SV', 'EV', 'IV', 'RV')
%     
%     figure
    totalInf = y(:, 3) + y(:, 8);
%     plot(t, totalInf, t, y(:, 3), t, y(:, 8))
%     legend('Total Infections', 'Unvaccinated Infectious', 'Vaccinated Infectious')
%     
%     figure
%     plot(t, 100*y(:, 3)./totalInf, t, 100*y(:, 8)./totalInf)
%     legend('Percent Unvaccinated Infectious', 'Percent Vaccinated Infectious')
%     
%     F = [0 beta*y0(1)*((1 - eps)/pVec(9) + eps/N) 0 beta*y0(1)*(eps/N); 0 0 0 0; ...
%         0 beta*y0(6)*(eps/N) 0 beta*y0(6)*((1 - eps)/pVec(8) + eps/N); 0 0 0 0]; 
%     V = [1/alpha 0 0 0; -1/alpha gamma 0 0; 0 0 1/alpha 0; 0 0 -1/alpha gamma];
%     evals = eig(F/V); R0 = evals(1); disp(['R0 = ' num2str(R0)])
%     figure
% %     plot(t, R0*(y(:, 1) + y(:, 6))/N)
%     plot(t, R0*y(:, 1)/pVec(9) + R0*y(:, 6)/pVec(8))
%     hold on
% %     plot(t, R0*y(:, 1)/N, t, R0*y(:, 6)/N)
%     plot(t, R0*y(:, 1)/pVec(9), t, R0*y(:, 6)/pVec(8))
%     legend('Reff', 'Reff Unvaccinated', 'Reff Vaccinated')

    pUnv = sum(y(:, 3))/sum(totalInf); % percent of all cases among unvaccinated

end

% SEIR model with vaccination
% y0 = [(pVec(9) - 1); 0; 1; 0; pVec(8); 0; 0; 0; 0]; % SU, EU, IU, RU, V, SV, EV, IV, RV
function f = seir2(t, y, pVec)
    
    f = zeros(9, 1);
    f(1) = -pVec(2)*y(1)*(1-pVec(11))*pVec(13)*((1 - pVec(10))*(y(3)/pVec(9)) + pVec(10)*((y(8)/(1-pVec(11))*pVec(13)*((1-0.5*pVec(11)-0.5*pVec(12))) + y(3))/pVec(1))) + (1/pVec(5))*y(4);
    f(2) = pVec(2)*y(1)*(1-pVec(11))*pVec(13)*((1 - pVec(10))*(y(3)/pVec(9)) + pVec(10)*((y(8)/(1-pVec(11))*pVec(13)*((1-0.5*pVec(11)-0.5*pVec(12))) + y(3))/pVec(1))) - (1/pVec(3))*y(2);
    f(3) = (1/pVec(3))*y(2) - pVec(4)*y(3);
    f(4) = pVec(4)*y(3) - (1/pVec(5))*y(4);
    f(5) = -(1/pVec(7))*y(5);
    f(6) = -pVec(2)*y(6)*(1-pVec(12))*pVec(13)*((1 - pVec(10))*(y(8)/pVec(8)) + pVec(10)*((y(8) + y(3)/(1-pVec(12))*pVec(13)*((1-0.5*pVec(11)-0.5*pVec(12))))/pVec(1))) + (1/pVec(7))*y(5) + (1/pVec(6))*y(9);
    f(7) = pVec(2)*y(6)*(1-pVec(12))*pVec(13)*((1 - pVec(10))*(y(8)/pVec(8)) + pVec(10)*((y(8) + y(3)/(1-pVec(12))*pVec(13)*((1-0.5*pVec(11)-0.5*pVec(12))))/pVec(1))) - (1/pVec(3))*y(7);
    f(8) = (1/pVec(3))*y(7) - pVec(4)*y(8);
    f(9) = pVec(4)*y(8) - (1/pVec(6))*y(9);

end