function [t, y, pVec] = seir_simulator2(eps, ppV, pImm,tau3,tspan)

    N = 100000; % population size
    effT=1;
    beta = 0.75; % transmission rate
    alpha = 4; % latency period
    gamma = 1/9; % 1/(infectious period)
    tau1 = 365; % duration of prior-infection acquired immunity (in days)
    tau2 = tau1; % duration of vaccine plus prior-infection acquired immunity (in days)
%    tau3 = 213; % duration of vaccine acquired immunity (in days)
%     ppV = 0.7; % percent of population vaccinated
%     eps = 0.5; % degree of associativity (random mixing eps = 1)
%     pImm = 0.3; % percent of population with prior-infection immunity initially
    pVec = [N beta alpha gamma tau1 tau2 tau3 ppV*N (1 - ppV)*N eps effT];
    
    
    y0 = [((1 - pImm)*pVec(9) - 1); 0; 1; pImm*pVec(9); (1 - pImm)*pVec(8); 0; 0; 0; pImm*pVec(8)]; 
    options = odeset('MaxStep', 1);
    [t, y] = ode45(@(t, y) seir2(t, y, pVec), tspan, y0, options);
    totalInf = y(:, 3) + y(:, 8);

end

% SEIR model with vaccination
function f = seir2(t, y, pVec)
    
    f = zeros(4, 1);
    f(1) = -pVec(2)*y(1)*((1 - pVec(10))*(y(3)/pVec(9)) + pVec(10)*((pVec(11)*y(8) + y(3))/pVec(1))) + (1/pVec(5))*y(4);
    f(2) = pVec(2)*y(1)*((1 - pVec(10))*(y(3)/pVec(9)) + pVec(10)*((pVec(11)*y(8) + y(3))/pVec(1))) - (1/pVec(3))*y(2);
    f(3) = (1/pVec(3))*y(2) - pVec(4)*y(3);
    f(4) = pVec(4)*y(3) - (1/pVec(5))*y(4);
    f(5) = -(1/pVec(7))*y(5);
    f(6) = -pVec(2)*y(6)*((1 - pVec(10))*(pVec(11)*y(8)/pVec(8)) + pVec(10)*((pVec(11)*y(8) + y(3))/pVec(1))) + (1/pVec(7))*y(5) + (1/pVec(6))*y(9);
    f(7) = pVec(2)*y(6)*((1 - pVec(10))*(pVec(11)*y(8)/pVec(8)) + pVec(10)*((pVec(11)*y(8) + y(3))/pVec(1))) - (1/pVec(3))*y(7);
    f(8) = (1/pVec(3))*y(7) - pVec(4)*y(8);
    f(9) = pVec(4)*y(8) - (1/pVec(6))*y(9);

end