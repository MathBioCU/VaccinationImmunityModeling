function m = covidVaxMixModel(tVec, bVec, epsilon)
    % Sabina L. Altus and David M. Bortz (April 2022)
    
    %% Set up the model with Colorado population data
    tVec = [10, tVec]; % switch dates for beta
    tVec = tVec - 1; % shift since events executed after SEIR
    evtBd = 1; % force events to occur on correct day
    ageCat = {'0to19','20to39','40to64','65toinf'};
    Na = length(ageCat);
    rNum = 1;
    
    popByAge = [1411161     1697671     1818147      878522     5805501];
    load('currentLOS.mat'), hospLOS2 = zeros(1, 6, Na);
    load('currentVaccData.mat', 'lphaVacMult', 'lphaDailyB')
    vaccMult = zeros(size(lphaDailyB, 1), 1, Na); boostMat = zeros(size(lphaDailyB, 1), 1, Na);
    for i = 1:Na
        hospLOS2(1, :, i) = [hospLOSmetro(1, :, i) 0.6*hospLOSmetro(1, 5, i)];
        vaccMult(:, 1, i) = sum(lphaVacMult(:, :, i), 2);
        boostMat(:, 1, i) = sum(lphaDailyB(:, :, i), 2);
    end
    
    m = sbiomodel('covidVaxMixModel');   

    %% Parameters
    
    % Global constants
    p(1) = addparameter(m, 'epsilon', epsilon); % degree of between group mixing (random mixing: eps = 1)
    p(2) = addparameter(m, 'lambda', 1.395); % difference in infectiousness: symptomatic/asymptomatic
    p(3) = addparameter(m, 'alpha', 2.5); % incubation period
    p(4) = addparameter(m, 'gamma', 1/5); % recovery rate
    p(5) = addparameter(m, 'I0', 1); % initial number of infected individuals
    p(6) = addparameter(m, 'N', popByAge(5)); % total population
    p(7) = addparameter(m, 'N1', popByAge(1)); % age1 population
    p(8) = addparameter(m, 'N2', popByAge(2)); % age2 population
    p(9) = addparameter(m, 'N3', popByAge(3)); % age3 population
    p(10) = addparameter(m, 'N4', popByAge(4)); % age4 population
    
    % Parameter value change dates
    p(11) = addparameter(m, 'ioBegin', 24 - 1); % initial infection begins day 24 = 01/24/2020
    p(12) = addparameter(m, 'ioBd', 24 + evtBd - 1);
    p(22) = addparameter(m, 'tsh', 213 - 1); % switch date for hospitalization l.o.s. parameters (July 31)
    p(23) = addparameter(m, 'tscc1', 170 - 1); % first switch date for hosp + cc parameters (June 18)
    p(24) = addparameter(m, 'tscc2', 273 - 1); % second switch date for hosp + cc parameters (September 29)
    p(25) = addparameter(m, 'vDuration', 1/213);%CHANGE 8/22/2024 % 1/duration of immunity from vaccine (1/vd) pre-Delta
    set(p(25), 'ConstantValue', false);
    p(26) = addparameter(m, 'recoImm', 1/365);%CHANGE 8/22/2024 % 1/duration of immunity for recovery
    p(27) = addparameter(m, 'boostImm', 1/213);%CHANGE 8/22/2024 % 1/duration of immunity from booster

    % Time Dependent Parameters
    p(30) = addparameter(m, 'b', 0); % transmission rate
    set(p(30), 'ConstantValue', false);
    p(31) = addparameter(m, 'betaAuv', p(30).Value*p(1).Value/p(6).Value);
    set(p(31), 'ConstantValue', false);
    p(32) = addparameter(m, 'betaIuv', p(2).Value*p(31).Value);
    set(p(32), 'ConstantValue', false);
    p(33) = addparameter(m, ['bOptDay' num2str(tVec(1) + 1)], 0.01);
    
    % Age specific parameters
    pS = [0.110023 0.35705 0.561205 0.774879]; % fraction of symptomatic cases
    dnh = [0.000013 0.0000822 0.000622705 0.027971448]; % death rate for non-hospitalized patients
    chosp = [0.0386    0.0622    0.0863    0.1343]; % combined rate hospitalized (hosp + cc)
    chlos = [5.8303 6.366447985 10.54149451 10.53665902]; % combined hospitalization length of stay
    dch = [0.005504587 0.014716188 0.057513214 0.160570555]; % death rate for hospitalized patients

    % Hospitalization parameter values after switch
    losSwitches = [95 183 275 398 548] - 1; 
    chlosb = [5.5747 5.231230631 8.486187615 8.136581647]; % B: July 31, tsh
    chospb = [0.0386    0.0637    0.0850    0.1343]; % B: June 18, tscc1 (updated 11/08/2021)
    chospc = [0.0239    0.0398    0.0563    0.1343]; % C: September 29, tscc2 (updated 11/08/2021)
    
    p(34) = addparameter(m, 'hospMult', 1); % initial hospitalization rate multiplier (no variants)
    set(p(34), 'ConstantValue', false);
    p(35) = addparameter(m, 'tVariant1', 412); % B.1.1.7 (02/15/2021)
    p(36) = addparameter(m, 'hospMult1', 1.2243); % B.1.1.7 multiplier
    p(37) = addparameter(m, 'tVariant2', 501); % B.1.1.7 and CA (05/15/21)
    p(38) = addparameter(m, 'hospMult2', 1.5598); % B.1.1.7 multiplier
    p(39) = addparameter(m, 'tVariant3', 562); % Delta (07/15/2021)
    p(40) = addparameter(m, 'hospMult3', 1.6810); % Delta multiplier
    
    % Immunity rate multipliers
    p(41) = addparameter(m, 'iMult', 0.925); % symptomatic immunity rate multiplier
    p(42) = addparameter(m, 'aMult', 0.85); % symptomatic immunity rate multiplier
    p(43) = addparameter(m, 'aRecov', p(4).Value*p(42).Value); % asymptomatic recovery rate
    p(44) = addparameter(m, 'aSus', p(4).Value*(1 - p(42).Value)); % rate of asymptomatic to susceptible
    p(45) = addparameter(m, 'vDurationDelta', 1/213); %CHANGE 8/22/2024 % 1/duration of immunity from vaccine (1/vd)
    
    % Omicron
    p(46) = addparameter(m, 'tVariant4', 562); % Omicron (01/01/2022)
    p(47) = addparameter(m, 'hospMult4', 0.5); % Omicron multiplier
    
    % Beta terms for mixing
    p(48) = addparameter(m, 'NV', 1); % number vaccinated
    set(p(48), 'ConstantValue', false);
    p(49) = addparameter(m, 'betaAsviv', p(30).Value*(((1 - p(1).Value)/p(48).Value) + p(1).Value/p(6).Value));
    set(p(49), 'ConstantValue', false);
    p(50) = addparameter(m, 'betaIsviv', p(2).Value*p(49).Value);
    set(p(50), 'ConstantValue', false);
    p(51) = addparameter(m, 'betaAsi', p(30).Value*(((1 - p(1).Value)/(p(6).Value - p(48).Value)) + p(1).Value/p(6).Value));
    set(p(51), 'ConstantValue', false);
    p(52) = addparameter(m, 'betaIsi', p(2).Value*p(51).Value);
    set(p(52), 'ConstantValue', false);
    
    % Hospitalization length of stay
    j = 55;
    for i = 1:length(losSwitches)
        p(j) = addparameter(m, ['tLOS' num2str(i)], losSwitches(i));
        p(j + 1) = addparameter(m, ['hRecovP' num2str(i) 'A1'], (p(41).Value*(1 - dch(1)))/chlosb(1));
        p(j + 2) = addparameter(m, ['hRecovP' num2str(i) 'A2'], (p(41).Value*(1 - dch(2)))/hospLOS2(rNum, i, 2));
        p(j + 3) = addparameter(m, ['hRecovP' num2str(i) 'A3'], (p(41).Value*(1 - dch(3)))/hospLOS2(rNum, i, 3));
        p(j + 4) = addparameter(m, ['hRecovP' num2str(i) 'A4'], (p(41).Value*(1 - dch(4)))/hospLOS2(rNum, i, 4));
        p(j + 5) = addparameter(m, ['hDeathP' num2str(i) 'A1'], dch(1)/chlosb(1));
        p(j + 6) = addparameter(m, ['hDeathP' num2str(i) 'A2'], dch(2)/hospLOS2(rNum, i, 2));
        p(j + 7) = addparameter(m, ['hDeathP' num2str(i) 'A3'], dch(3)/hospLOS2(rNum, i, 3));
        p(j + 8) = addparameter(m, ['hDeathP' num2str(i) 'A4'], dch(4)/hospLOS2(rNum, i, 4));
        p(j + 9) = addparameter(m, ['hSusP' num2str(i) 'A1'], ((1 - p(41).Value)/p(41).Value)*p(j + 1).Value);
        p(j + 10) = addparameter(m, ['hSusP' num2str(i) 'A2'], ((1 - p(41).Value)/p(41).Value)*p(j + 2).Value);
        p(j + 11) = addparameter(m, ['hSusP' num2str(i) 'A3'], ((1 - p(41).Value)/p(41).Value)*p(j + 3).Value);
        p(j + 12) = addparameter(m, ['hSusP' num2str(i) 'A4'], ((1 - p(41).Value)/p(41).Value)*p(j + 4).Value);
        j = j + 13;
    end
    
    for i = 1:Na % p(3) = alpha, p(4) = gamma
        p(j) = addparameter(m, ['iRecovery' num2str(i)], p(4).Value*p(41).Value*(1 - p(34).Value*chosp(i) - dnh(i)));
        set(p(j), 'ConstantValue', false);
        p(j + 1) = addparameter(m, ['iRecoveryB' num2str(i)], p(4).Value*p(41).Value*(1 - p(34).Value*chospb(i) - dnh(i)));
        p(j + 2) = addparameter(m, ['iRecoveryC' num2str(i)], p(4).Value*p(41).Value*(1 - p(34).Value*chospc(i) - dnh(i)));
        
        p(j + 3) = addparameter(m, ['hRecov' num2str(i)], (p(41).Value*(1 - dch(i)))/chlos(i));
        set(p(j + 3), 'ConstantValue', false);
        
        p(j + 4) = addparameter(m, ['iDeath' num2str(i)], dnh(i)*p(4).Value);
        p(j + 5) = addparameter(m, ['hDeath' num2str(i)], dch(i)/chlos(i));
        set(p(j + 5), 'ConstantValue', false);
        
        p(j + 6) = addparameter(m, ['hHospitalized' num2str(i)], p(34).Value*chosp(i)*p(4).Value);
        set(p(j + 6), 'ConstantValue', false);
        p(j + 7) = addparameter(m, ['hHospitalizedB' num2str(i)], p(34).Value*chospb(i)*p(4).Value);
        p(j + 8) = addparameter(m, ['hHospitalizedC' num2str(i)], p(34).Value*chospc(i)*p(4).Value);
        
        p(j + 9) = addparameter(m, ['eInfectious' num2str(i)], pS(i)/p(3).Value);
        p(j + 10) = addparameter(m, ['eAsymptomatic' num2str(i)], (1 - pS(i))/p(3).Value);
        
        p(j + 11) = addparameter(m, ['chosp' num2str(i)], chospc(i));
        p(j + 12) = addparameter(m, ['dnh' num2str(i)], dnh(i));
        p(j + 13) = addparameter(m, ['iSus' num2str(i)], ((1 - p(41).Value)/p(41).Value)*p(j).Value);
        set(p(j + 13), 'ConstantValue', false);
        p(j + 14) = addparameter(m, ['hSus' num2str(i)], ((1 - p(41).Value)/p(41).Value)*p(j + 3).Value);
        set(p(j + 14), 'ConstantValue', false);
        
        p(j + 15) = addparameter(m, ['hvHospitalized' num2str(i)], chospc(i)*p(4).Value);
        set(p(j + 15), 'ConstantValue', false);
        
        j = j + 16;
    end
    
    % Vaccinations
    for i = 1:Na
        p(j) = addparameter(m, ['NbarA' num2str(i)], popByAge(i));
        set(p(j), 'ConstantValue', false);
        p(j + 1) = addparameter(m, ['vaccNumA' num2str(i)], 0);
        set(p(j + 1), 'ConstantValue', false);
        p(j + 2) = addparameter(m, ['vaccRateA' num2str(i)], p(j + 1).Value/p(j).Value);
        set(p(j + 2), 'ConstantValue', false);
        p(j + 3) = addparameter(m, ['NboostA' num2str(i)], 10);
        set(p(j + 3), 'ConstantValue', false);
        p(j + 4) = addparameter(m, ['boostNumA' num2str(i)], 0);
        set(p(j + 4), 'ConstantValue', false);
        p(j + 5) = addparameter(m, ['boostRateA' num2str(i)], p(j + 4).Value/p(j + 3).Value);
        set(p(j + 5), 'ConstantValue', false);
        j = j + 6;
    end
    
    for i = 1:size(vaccMult, 1)
        for k = 1:Na
            p(j) = addparameter(m, ['vaccNum' num2str(i) 'A' num2str(k)], vaccMult(i, rNum, k));
            p(j + 1) = addparameter(m, ['boostNum' num2str(i) 'A' num2str(k)], boostMat(i, rNum, k));
            j = j + 2;
        end
    end
    
    % Parameters for Optimization
    for i = 1:length(tVec)
        p(j) = addparameter(m, ['tOptDay' num2str(tVec(i) + 1)], tVec(i));
        p(j + 1) = addparameter(m, ['tOptDay' num2str(tVec(i) + 1) 'bd'], tVec(i) + evtBd);
        j = j + 2;
        if i >= 2
            p(j) = addparameter(m, ['bOptDay' num2str(tVec(i) + 1)], bVec(i - 1));
            j = j + 1;
        end
    end
    
    %% Species and Compartments
    
    for i = 1:Na
        CO{i} = addcompartment(m, ageCat{i});
        S{i} = addspecies(CO{i}, 'S');    % Susceptible
        E{i} = addspecies(CO{i}, 'E');    % Exposed
        I{i} = addspecies(CO{i}, 'I');    % (Symptomatic) Infectious  
        A{i} = addspecies(CO{i}, 'A');    % Asymptomatic Infectious
        R{i} = addspecies(CO{i}, 'R');    % Recovered
        V{i} = addspecies(CO{i}, 'V');    % Vaccinated
        H{i} = addspecies(CO{i}, 'H');    % Hospitalized (including in ICU)
        D{i} = addspecies(CO{i}, 'D');    % Deceased
        SV{i} = addspecies(CO{i}, 'SV');    % Susceptible & Vaccinated
        EV{i} = addspecies(CO{i}, 'EV');    % Exposed & Vaccinated
        IV{i} = addspecies(CO{i}, 'IV');    % (Symptomatic) Infectious & Vaccinated
        AV{i} = addspecies(CO{i}, 'AV');    % Asymptomatic Infectious & Vaccinated
        RV{i} = addspecies(CO{i}, 'RV');    % Recovered & Vaccinated
        HV{i} = addspecies(CO{i}, 'HV');    % Hospitalized & Vaccinated
        B{i} = addspecies(CO{i}, 'B');      % Boosted
    end
    hosptot = addspecies(CO{1}, 'hosptot');

    %% Reactions
    
    for i = 1:Na
            
        % H Recovery
        H2R{i} = addreaction(m, [ageCat{i} '.H -> ' ageCat{i} '.R']);
        H2R{i}.Name = [ageCat{i} ':H2R'];
        H2R{i}.addkineticlaw('MassAction');
        setparameter(H2R{i}.KineticLaw, 'Forward Rate Parameter', ['hRecov' num2str(i)]);

        % Hospitalization Death
        H2D{i} = addreaction(m, [ageCat{i} '.H -> ' ageCat{i} '.D']);
        H2D{i}.Name = [ageCat{i} ':H2D'];
        H2D{i}.addkineticlaw('MassAction');
        setparameter(H2D{i}.KineticLaw, 'Forward Rate Parameter', ['hDeath' num2str(i)]);
        
        % HV Recovery
        HV2RV{i} = addreaction(m, [ageCat{i} '.HV -> ' ageCat{i} '.RV']);
        HV2RV{i}.Name = [ageCat{i} ':HV2RV'];
        HV2RV{i}.addkineticlaw('MassAction');
        setparameter(HV2RV{i}.KineticLaw, 'Forward Rate Parameter', ['hRecov' num2str(i)]);

        % HV Death
        HV2D{i} = addreaction(m, [ageCat{i} '.HV -> ' ageCat{i} '.D']);
        HV2D{i}.Name = [ageCat{i} ':HV2D'];
        HV2D{i}.addkineticlaw('MassAction');
        setparameter(HV2D{i}.KineticLaw, 'Forward Rate Parameter', ['hDeath' num2str(i)]);

        % Susceptible to Vaccinated
        S2V{i} = addreaction(m, [ageCat{i} '.S -> ' ageCat{i} '.V']);
        S2V{i}.Name = [ageCat{i} ': S2V'];
        S2V{i}.addkineticlaw('MassAction');
        setparameter(S2V{i}.KineticLaw, 'Forward Rate Parameter', ['vaccRateA' num2str(i)]);
        
        % Vaccinated to Boosted
        V2B{i} = addreaction(m, [ageCat{i} '.V -> ' ageCat{i} '.B']);
        V2B{i}.Name = [ageCat{i} ': V2B'];
        V2B{i}.addkineticlaw('MassAction');
        setparameter(V2B{i}.KineticLaw, 'Forward Rate Parameter', ['boostRateA' num2str(i)]);
        
        % Susceptible Vaccinated to Boosted
        SV2B{i} = addreaction(m, [ageCat{i} '.SV -> ' ageCat{i} '.B']);
        SV2B{i}.Name = [ageCat{i} ': SV2B'];
        SV2B{i}.addkineticlaw('MassAction');
        setparameter(SV2B{i}.KineticLaw, 'Forward Rate Parameter', ['boostRateA' num2str(i)]);
        
        % Susceptible Recovered to Boosted
        RV2B{i} = addreaction(m, [ageCat{i} '.RV -> ' ageCat{i} '.B']);
        RV2B{i}.Name = [ageCat{i} ': RV2B'];
        RV2B{i}.addkineticlaw('MassAction');
        setparameter(RV2B{i}.KineticLaw, 'Forward Rate Parameter', ['boostRateA' num2str(i)]);

        % Recovered to Vaccinated
        R2V{i} = addreaction(m, [ageCat{i} '.R -> ' ageCat{i} '.V']);
        R2V{i}.Name = [ageCat{i} ': R2V'];
        R2V{i}.addkineticlaw('MassAction');
        setparameter(R2V{i}.KineticLaw, 'Forward Rate Parameter', ['vaccRateA' num2str(i)]);

        % Recovered to Susceptible
        R2S{i} = addreaction(m, [ageCat{i} '.R -> ' ageCat{i} '.S']);
        R2S{i}.Name = [ageCat{i} ': R2S'];
        R2S{i}.addkineticlaw('MassAction');
        setparameter(R2S{i}.KineticLaw, 'Forward Rate Parameter', 'recoImm');
        
        % RV to Susceptible
        RV2SV{i} = addreaction(m, [ageCat{i} '.RV -> ' ageCat{i} '.SV']);
        RV2SV{i}.Name = [ageCat{i} ': RV2SV'];
        RV2SV{i}.addkineticlaw('MassAction');
        setparameter(RV2SV{i}.KineticLaw, 'Forward Rate Parameter', 'recoImm');

        % Vaccinated to Susceptible
        V2SV{i} = addreaction(m, [ageCat{i} '.V -> ' ageCat{i} '.SV']);
        V2SV{i}.Name = [ageCat{i} ': V2SV'];
        V2SV{i}.addkineticlaw('MassAction');
        setparameter(V2SV{i}.KineticLaw, 'Forward Rate Parameter', 'vDuration');
        
        % Boosted to Susceptible
        B2SV{i} = addreaction(m, [ageCat{i} '.B -> ' ageCat{i} '.SV']);
        B2SV{i}.Name = [ageCat{i} ': B2SV'];
        B2SV{i}.addkineticlaw('MassAction');
        setparameter(B2SV{i}.KineticLaw, 'Forward Rate Parameter', 'boostImm');
      
        % Symptomatic to Susceptible
        I2S{i} = addreaction(m, [ageCat{i} '.I -> ' ageCat{i} '.S']);
        I2S{i}.Name = [ageCat{i} ': I2S'];
        I2S{i}.addkineticlaw('MassAction');
        setparameter(I2S{i}.KineticLaw, 'Forward Rate Parameter', ['iSus' num2str(i)]);
        
        % Asymptomatic to Susceptible
        A2S{i} = addreaction(m, [ageCat{i} '.A -> ' ageCat{i} '.S']);
        A2S{i}.Name = [ageCat{i} ': A2S'];
        A2S{i}.addkineticlaw('MassAction');
        setparameter(A2S{i}.KineticLaw, 'Forward Rate Parameter', 'aSus');
        
        % Hospitalized to Susceptible
        H2S{i} = addreaction(m, [ageCat{i} '.H -> ' ageCat{i} '.S']);
        H2S{i}.Name = [ageCat{i} ': H2S'];
        H2S{i}.addkineticlaw('MassAction');
        setparameter(H2S{i}.KineticLaw, 'Forward Rate Parameter', ['hSus' num2str(i)]);

        % Asymptomatic Recovery
        A2R{i} = addreaction(m, [ageCat{i} '.A -> ' ageCat{i} '.R']);
        A2R{i}.Name = [ageCat{i} ': A2R'];
        A2R{i}.addkineticlaw('MassAction');
        setparameter(A2R{i}.KineticLaw, 'Forward Rate Parameter', 'aRecov');
        
        % Asymptomatic Vacc Recovery
        AV2RV{i} = addreaction(m, [ageCat{i} '.AV -> ' ageCat{i} '.RV']);
        AV2RV{i}.Name = [ageCat{i} ': AV2RV'];
        AV2RV{i}.addkineticlaw('MassAction');
        setparameter(AV2RV{i}.KineticLaw, 'Forward Rate Parameter', 'aRecov');

        % (Symptomatic) Infectious Recovery
        I2R{i} = addreaction(m, [ageCat{i} '.I -> ' ageCat{i} '.R']);
        I2R{i}.Name = [ageCat{i} ': I2R'];
        I2R{i}.addkineticlaw('MassAction');
        setparameter(I2R{i}.KineticLaw, 'Forward Rate Parameter', ['iRecovery' num2str(i)]);
        
        % (Symptomatic) Infectious Vaccinated Recovery
        IV2RV{i} = addreaction(m, [ageCat{i} '.IV -> ' ageCat{i} '.RV']);
        IV2RV{i}.Name = [ageCat{i} ': IV2RV'];
        IV2RV{i}.addkineticlaw('MassAction');
        setparameter(IV2RV{i}.KineticLaw, 'Forward Rate Parameter', ['iRecovery' num2str(i)]);

        % Non-Hospitalized Death
        I2D{i} = addreaction(m, [ageCat{i} '.I -> ' ageCat{i} '.D']);
        I2D{i}.Name = [ageCat{i} ': I2D'];
        I2D{i}.addkineticlaw('MassAction');
        setparameter(I2D{i}.KineticLaw, 'Forward Rate Parameter', ['iDeath' num2str(i)]);

        % Exposed to (Symptomatic) Infectious
        E2I{i} = addreaction(m, [ageCat{i} '.E -> ' ageCat{i} '.I']);
        E2I{i}.Name = [ageCat{i} ': E2I'];
        E2I{i}.addkineticlaw('MassAction');
        setparameter(E2I{i}.KineticLaw, 'Forward Rate Parameter', ['eInfectious' num2str(i)]);

        % Exposed to Asymptomatic Infectious
        E2A{i} = addreaction(m, [ageCat{i} '.E -> ' ageCat{i} '.A']);
        E2A{i}.Name = [ageCat{i} ': E2A'];
        E2A{i}.addkineticlaw('MassAction');
        setparameter(E2A{i}.KineticLaw, 'Forward Rate Parameter', ['eAsymptomatic' num2str(i)]);
        
        % Exposed Vaccinated to (Symptomatic) Infectious Vaccinated
        EV2IV{i} = addreaction(m, [ageCat{i} '.EV -> ' ageCat{i} '.IV']);
        EV2IV{i}.Name = [ageCat{i} ': EV2IV'];
        EV2IV{i}.addkineticlaw('MassAction');
        setparameter(EV2IV{i}.KineticLaw, 'Forward Rate Parameter', ['eInfectious' num2str(i)]);

        % Exposed Vaccinated to Asymptomatic Infectious Vaccinated
        EV2AV{i} = addreaction(m, [ageCat{i} '.EV -> ' ageCat{i} '.AV']);
        EV2AV{i}.Name = [ageCat{i} ': EV2AV'];
        EV2AV{i}.addkineticlaw('MassAction');
        setparameter(EV2AV{i}.KineticLaw, 'Forward Rate Parameter', ['eAsymptomatic' num2str(i)]);

        % I2H
        I2H{i} = addreaction(m, [ageCat{i} '.I -> ' ageCat{i} '.H']);
        I2H{i}.Name = [ageCat{i} ':I2H'];
        I2H{i}.addkineticlaw('MassAction');
        setparameter(I2H{i}.KineticLaw, 'Forward Rate Parameter', ['hHospitalized' num2str(i)]);
        
        % IV2HV
        IV2HV{i} = addreaction(m, [ageCat{i} '.IV -> ' ageCat{i} '.HV']);
        IV2HV{i}.Name = [ageCat{i} ':IV2HV'];
        IV2HV{i}.addkineticlaw('MassAction');
        setparameter(IV2HV{i}.KineticLaw, 'Forward Rate Parameter', ['hvHospitalized' num2str(i)]);

        % Nonlinear terms: Susceptible to Exposed
        for j = 1:Na

            % Exposed via Symptomatic Infectious Individual
            EbyI{i} = addreaction(m, [ageCat{i} '.S + ' ageCat{j} '.I -> ' ageCat{j} '.I + ' ageCat{i} '.E']);
            EbyI{i}.Name = [ageCat{i} '.S - ' ageCat{j} '.I: Infection'];
            EbyI{i}.addkineticlaw('MassAction');
            setparameter(EbyI{i}.KineticLaw, 'Forward Rate Parameter', 'betaIsi');

            % Exposed via Asymptomatic Infectious Individual
            EbyA{i} = addreaction(m, [ageCat{i} '.S + ' ageCat{j} '.A -> ' ageCat{j} '.A + ' ageCat{i} '.E']);
            EbyA{i}.Name = [ageCat{i} '.S - ' ageCat{j} '.A: Infection'];
            EbyA{i}.addkineticlaw('MassAction');
            setparameter(EbyA{i}.KineticLaw, 'Forward Rate Parameter', 'betaAsi');
            
            % Exposed Vaccinated via Symptomatic Infectious Individual
            EVbyI{i} = addreaction(m, [ageCat{i} '.SV + ' ageCat{j} '.I -> ' ageCat{j} '.I + ' ageCat{i} '.EV']);
            EVbyI{i}.Name = [ageCat{i} '.SV - ' ageCat{j} '.I: Infection'];
            EVbyI{i}.addkineticlaw('MassAction');
            setparameter(EVbyI{i}.KineticLaw, 'Forward Rate Parameter', 'betaIuv');

            % Exposed Vaccinataed via Asymptomatic Infectious Individual
            EVbyA{i} = addreaction(m, [ageCat{i} '.SV + ' ageCat{j} '.A -> ' ageCat{j} '.A + ' ageCat{i} '.EV']);
            EVbyA{i}.Name = [ageCat{i} '.SV - ' ageCat{j} '.A: Infection'];
            EVbyA{i}.addkineticlaw('MassAction');
            setparameter(EVbyA{i}.KineticLaw, 'Forward Rate Parameter', 'betaAuv');
            
            % Exposed via Symptomatic Infectious Vaccinated Individual
            EbyIV{i} = addreaction(m, [ageCat{i} '.S + ' ageCat{j} '.IV -> ' ageCat{j} '.IV + ' ageCat{i} '.E']);
            EbyIV{i}.Name = [ageCat{i} '.S - ' ageCat{j} '.IV: Infection'];
            EbyIV{i}.addkineticlaw('MassAction');
            setparameter(EbyIV{i}.KineticLaw, 'Forward Rate Parameter', 'betaIuv');

            % Exposed via Asymptomatic Infectious Vaccinated Individual
            EbyAV{i} = addreaction(m, [ageCat{i} '.S + ' ageCat{j} '.AV -> ' ageCat{j} '.AV + ' ageCat{i} '.E']);
            EbyAV{i}.Name = [ageCat{i} '.S - ' ageCat{j} '.AV: Infection'];
            EbyAV{i}.addkineticlaw('MassAction');
            setparameter(EbyAV{i}.KineticLaw, 'Forward Rate Parameter', 'betaAuv');
            
            % Exposed Vaccinated via Symptomatic Vaccinated Infectious Individual
            EVbyIV{i} = addreaction(m, [ageCat{i} '.SV + ' ageCat{j} '.IV -> ' ageCat{j} '.IV + ' ageCat{i} '.EV']);
            EVbyIV{i}.Name = [ageCat{i} '.SV - ' ageCat{j} '.IV: Infection'];
            EVbyIV{i}.addkineticlaw('MassAction');
            setparameter(EVbyIV{i}.KineticLaw, 'Forward Rate Parameter', 'betaIsviv');

            % Exposed Vaccinataed via Asymptomatic Vaccinated Infectious Individual
            EVbyAV{i} = addreaction(m, [ageCat{i} '.SV + ' ageCat{j} '.AV -> ' ageCat{j} '.AV + ' ageCat{i} '.EV']);
            EVbyAV{i}.Name = [ageCat{i} '.SV - ' ageCat{j} '.AV: Infection'];
            EVbyAV{i}.addkineticlaw('MassAction');
            setparameter(EVbyAV{i}.KineticLaw, 'Forward Rate Parameter', 'betaAsviv');

        end
        
    end

    %% Initial Conditions 
    
    for i = 1:Na
        m.Compartment(i).Species(1).InitialAmount = popByAge(i);
    end

    %% Events
    
    % Initial infection
    evt1 = addevent(m, '(time >= ioBegin) && (time < ioBd)', {});
        set(evt1, 'Name', 'Initial infection (day 24 = 01/24/2020)');
        set(evt1, 'EventFcns', {['[' ageCat{1} '].S = ([' ageCat{1} '].S - 1)'], ...
            ['[' ageCat{1} '].I = I0']});
    
    % Changes in hospitalization parameters
    for i = 1:length(losSwitches)
        evt2(i) = addevent(m, ['time >= tLOS' num2str(i)], {}); % July 31
        set(evt2(i), 'Name', ['Change in hospitalization length of stay ' num2str(i)]);
        set(evt2(i), 'EventFcns', {['hRecov1 = hRecovP' num2str(i) 'A1'], ...
            ['hRecov2 = hRecovP' num2str(i) 'A2'], ...
            ['hRecov3 = hRecovP' num2str(i) 'A3'], ...
            ['hRecov4 = hRecovP' num2str(i) 'A4'], ...
            ['hSus1 = ((1 - iMult)/iMult)*hRecovP' num2str(i) 'A1'], ...
            ['hSus2 = ((1 - iMult)/iMult)*hRecovP' num2str(i) 'A2'], ...
            ['hSus3 = ((1 - iMult)/iMult)*hRecovP' num2str(i) 'A3'], ...
            ['hSus4 = ((1 - iMult)/iMult)*hRecovP' num2str(i) 'A4'], ...
            ['hDeath1 = hDeathP' num2str(i) 'A1'], ...
            ['hDeath2 = hDeathP' num2str(i) 'A2'], ...
            ['hDeath3 = hDeathP' num2str(i) 'A3'], ...
            ['hDeath4 = hDeathP' num2str(i) 'A4']});
    end
    
    for i = 1:Na
        evt3(i) = addevent(m, 'time >= tscc1', {}); % June 18
        set(evt3(i), 'Name', ['First change in fraction hospitalized ' ageCat{i}]);
        set(evt3(i), 'EventFcns', {['iRecovery' num2str(i) ' = iRecoveryB' num2str(i)], ...
            ['iSus' num2str(i) ' = ((1 - iMult)/iMult)*iRecoveryB' num2str(i)], ...
            ['hHospitalized' num2str(i) ' = hHospitalizedB' num2str(i)]});
    end

    for i = 1:Na
        evt4(i) = addevent(m, 'time >= tscc2', {}); % September 29
        set(evt4(i), 'Name', ['Second change in fraction hospitalized ' ageCat{i}]);
        set(evt4(i), 'EventFcns', {['iRecovery' num2str(i) ' = iRecoveryC' num2str(i)], ...
            ['iSus' num2str(i) ' = ((1 - iMult)/iMult)*iRecoveryC' num2str(i)], ...
            ['hHospitalized' num2str(i) ' = hHospitalizedC' num2str(i)]});
    end
    
    % Change in hospitalization rate parameters for variants
    eCounter = 1;
    for vCounter = 1:4
        evt7(eCounter) = addevent(m, ['time >= tVariant' num2str(vCounter)], {});
        if vCounter == 3 % Delta
            for i = 1:Na
                set(evt7(eCounter), 'Name', ['Change in fraction hospitalized due to variants ' num2str(vCounter) ' ' ageCat{i}]);
                set(evt7(eCounter), 'EventFcns', {'vDuration = vDurationDelta', ...
                    ['hospMult = hospMult' num2str(vCounter)], ...
                    ['iRecovery' num2str(i) ' = gamma*(1 - hospMult' num2str(vCounter) '*chosp' num2str(i) ' - dnh' num2str(i) ')'], ...
                    ['hHospitalized' num2str(i) ' = gamma*hospMult' num2str(vCounter) '*chosp' num2str(i)], ...
                    ['hvHospitalized' num2str(i) ' = gamma*hospMult' num2str(vCounter) '*chosp' num2str(i)]});
            end
        elseif vCounter == 4 % Omicron
            for i = 1:Na
                set(evt7(eCounter), 'Name', ['Change in fraction hospitalized due to variants ' num2str(vCounter) ' ' ageCat{i}]);
                set(evt7(eCounter), 'EventFcns', {'vDuration = vDuration', ...
                    ['hospMult = hospMult' num2str(vCounter)], ...
                    ['iRecovery' num2str(i) ' = gamma*(1 - hospMult' num2str(vCounter) '*chosp' num2str(i) ' - dnh' num2str(i) ')'], ...
                    ['hHospitalized' num2str(i) ' = gamma*hospMult' num2str(vCounter) '*chosp' num2str(i)], ...
                    ['hvHospitalized' num2str(i) ' = gamma*hospMult' num2str(vCounter) '*chosp' num2str(i)]});
            end
        else
            for i = 1:Na
                set(evt7(eCounter), 'Name', ['Change in fraction hospitalized due to variants ' num2str(vCounter) ' ' ageCat{i}]);
                set(evt7(eCounter), 'EventFcns', {['hospMult = hospMult' num2str(vCounter)], ...
                    ['iRecovery' num2str(i) ' = gamma*(1 - hospMult' num2str(vCounter) '*chosp' num2str(i) ' - dnh' num2str(i) ')'], ...
                    ['hHospitalized' num2str(i) ' = gamma*hospMult' num2str(vCounter) '*chosp' num2str(i)], ...
                    ['hvHospitalized' num2str(i) ' = gamma*hospMult' num2str(vCounter) '*chosp' num2str(i)]});
            end
        end
    end
    
    eCounter = 1;
    for i = 1:size(vaccMult, 1)
        evt5(eCounter) = addevent(m, ['(time >= (348 + ' num2str(i) ')) && (time < (349 + ' num2str(i) '))'], {});
        set(evt5(eCounter), 'Name', ['Vaccination rate ' num2str(i)]);
        set(evt5(eCounter), 'EventFcns', {['vaccNumA1 = vaccNum' num2str(i) 'A1'], ...
            ['vaccNumA2 = vaccNum' num2str(i) 'A2'], ...
            ['vaccNumA3 = vaccNum' num2str(i) 'A3'], ...
            ['vaccNumA4 = vaccNum' num2str(i) 'A4'], ...
            ['boostNumA1 = boostNum' num2str(i) 'A1'], ...
            ['boostNumA2 = boostNum' num2str(i) 'A2'], ...
            ['boostNumA3 = boostNum' num2str(i) 'A3'], ...
            ['boostNumA4 = boostNum' num2str(i) 'A4'], ...
            'betaAsviv = b*(((1 - epsilon)/NV) + epsilon/N)', ...
            'betaIsviv = lambda*betaAsviv', ...
            'betaAsi = b*(((1 - epsilon)/(N - NV)) + epsilon/N)', ...
            'betaIsi = lambda*betaAsi'});
        eCounter = eCounter + 1;
    end
    
    % Bi-weekly change in beta
    eCounter = 1;
    for i = 1:length(tVec)
        evt6(eCounter) = addevent(m, ['(time >= tOptDay' num2str(tVec(i) + 1) ') && (time < tOptDay' num2str(tVec(i) + 1) 'bd)'], {});
        set(evt6(eCounter), 'Name', ['Change in beta, day ' num2str(tVec(i) + 1)]);
        set(evt6(eCounter), 'EventFcns', {['b = bOptDay' num2str(tVec(i) + 1)], ...
            ['betaAuv = bOptDay' num2str(tVec(i) + 1) '*epsilon/N'], ...
            'betaIuv = lambda*betaAuv', ...
            ['betaAsviv = bOptDay' num2str(tVec(i) + 1) '*(((1 - epsilon)/NV) + epsilon/N)'], ...
            'betaIsviv = lambda*betaAsviv', ...
            ['betaAsi = bOptDay' num2str(tVec(i) + 1) '*(((1 - epsilon)/(N - NV)) + epsilon/N)'], ...
            'betaIsi = lambda*betaAsi'});
        eCounter = eCounter + 1;
    end
    
    %% Fitting and Optimization
        
        % Total Hospitalizations
        strOpt1 = ['[' ageCat{1} '].H + [' ageCat{1} '].HV'];
        strOpt2 = ['1 + [' ageCat{1} '].V + [' ageCat{1} '].RV + [' ageCat{1} '].SV + [' ageCat{1} '].EV + [' ageCat{1} '].IV + [' ageCat{1} '].AV + [' ageCat{1} '].B+ [' ageCat{1} '].HV'];
        for i = 2:Na
            strOpt1 = [strOpt1 ' + [' ageCat{i} '].H + [' ageCat{i} '].HV'];
            strOpt2 = [strOpt2 ' + [' ageCat{i} '].V + [' ageCat{i} '].RV + [' ageCat{i} '].SV + [' ageCat{i} '].EV + [' ageCat{i} '].IV + [' ageCat{i} '].AV + [' ageCat{i} '].B+ [' ageCat{i} '].HV'];
        end
        strOpt1 = ['[' ageCat{1} '].hosptot = ' strOpt1];
        strOpt2 = ['NV = ' strOpt2];
        addrule(m, strOpt1, 'RuleType', 'repeatedAssignment');
        addrule(m, strOpt2, 'RuleType', 'repeatedAssignment');
        
        % Vaccination rates and Nbar (population eligible for vaccination)
        pEffective = [0.96 0.96 0.96 0.92];
        for i = 1:Na
            strOpt2 = ['[' ageCat{i} '].H + [' ageCat{i} '].D + [' ageCat{i} '].V + [' ageCat{i} '].RV + [' ageCat{i} '].SV + [' ageCat{i} '].EV + [' ageCat{i} '].IV + [' ageCat{i} '].AV + [' ageCat{i} '].B+ [' ageCat{i} '].HV'];
            strOpt2 = ['NbarA' num2str(i) ' = N' num2str(i) ' - (' strOpt2 ')'];
            addrule(m, strOpt2, 'RuleType', 'repeatedAssignment');
            strOpt3 = ['vaccRateA' num2str(i) ' = (' num2str(pEffective(i)) '*vaccNumA' num2str(i) ')/max(1, NbarA' num2str(i) ')'];
            addrule(m, strOpt3, 'RuleType', 'repeatedAssignment');
        end
        
        % Boosters
        for i = 1:Na
            strOpt2 = ['NboostA' num2str(i) ' = [' ageCat{i} '].SV + [' ageCat{i} '].V'];
            addrule(m, strOpt2, 'RuleType', 'repeatedAssignment');
            strOpt3 = ['boostRateA' num2str(i) ' = (boostNumA' num2str(i) ')/max(1, (NboostA' num2str(i) '))'];
            addrule(m, strOpt3, 'RuleType', 'repeatedAssignment');
        end
end