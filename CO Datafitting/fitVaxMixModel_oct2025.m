% RUNS OPTIMIZATION ON COVID MOB MODEL -- CHECK WHICH MODEL YOU ARE USING!
function [bMatrix, fitConstMAT] = fitVaxMixModel_oct(modelType, eps, indx) %, endOfData)
   
    endOfData = 1095;%825;
    metaStart = 58;
    epsilon = eps;%fliplr(0:0.05:1); %[1 0.9999 0.999 0.99]; %
    load('T_extend_oct9.mat')%data to end of 2022
    T=T_extend_oct9;
    T=groupedData(T);
    T.Properties.DimensionNames={'Row','Variables'};
    T.Properties.VariableNames={'days','hTotalCO'};
    T.Properties.IndependentVariableName='days';

    load('popInfo.mat', 'lpha1', 'metro')
    typeCell = {'metro', 'lpha', 'allCO', 'waState'; 12, 11, 1, 1; metro, lpha1, {'hTotalCO'}, {'washingtonState'}}; 
    regions = typeCell{3, modelType};

    % tVec = 24:14:endOfData;
    t_indices =[24:14:822 830:4:868];%switch data frequency from daily to weekly occurs at entry 830
    tVec=T.days(t_indices)'+1;%switch data frequency from daily to weekly occurs at entry 830
    tL = length(tVec);


    load('bMatrixOCT1_IG.mat') %21x124
    bMatrix = bMatrixOCT1_IG;
        bMatrix = bMatrix(indx,:);
    
    while size(bMatrix, 2) < tL
        bMatrix = [bMatrix, bMatrix(:, size(bMatrix, 2))];
    end
   fitConstMAT=zeros(length(bMatrix(1,:)));
    %% set up optimization
    rNum = 1; rCounter = 1;
    for eInd = 1:length(epsilon) %1:Nr
        
        startOpt = metaStart;
        % disp(regions{rNum})
        responseMap = {['[0to19].hosptot = ' regions{rNum}]};
        paramsToEstimate = {1, tL};
        for j = 1:tL
            paramsToEstimate{j} = ['bOptDay' num2str(tVec(j))];
        end
    
        method = 'particleswarm';
        options = optimoptions(method);
        options = optimoptions(options,'SwarmSize',200, 'HybridFcn', 'fmincon', 'UseParallel', true);%,'StepTolerance', 1e-10, 'FunctionTolerance', 1e-8,'TolCol', 1e-10);
    
        for paramIter = 1:length(metaStart:3:tL) %:tL
           
            endOpt = startOpt  + 4;

            if endOpt >= tL
                endOpt = tL;
            end 
            tend = tVec(endOpt)+42;
            % else
            %     tend = tVec(endOpt)+42;
            % end
                
            if length(1:(tend + 1)) <= length(T.days)
                gDataTemp = T(1:(tend + 1), :);
                % gDataTemp = groupedData(T2);
                % gDataTemp.Properties.IndependentVariableName = 'days';
            else
                gDataTemp = T(1:end,:);
            end
            
            % tic
            % disp(eInd), disp(rCounter), disp(eInd - (eInd - rCounter)), disp(epsilon(eInd))
            m = covidVaxMixModel(tVec(1:endOpt), bMatrix((eInd - (eInd - rCounter)), 1:endOpt), epsilon(eInd));
            csObj = getconfigset(m, 'active');
            csObj.SolverOptions.MaxStep = 1;
            set(csObj, 'StopTime', tend, 'SolverType', 'Sundials');
%             verify(m);
            sbioaccelerate(m, csObj);
            % % disp('Time to (verify and) accelerate the model'), toc
            % % 
            if paramIter == 1 % metaStart
                pGap = 1;
                for i = 1:numel(m.parameters)
                    if strcmp(m.Parameter(i).Name, 'bOptDay24')
                        pStart = i;
                    end

                    if strcmp(m.Parameter(i).Name, 'bOptDay38')
                        pGap = i - pStart;
                    end
                end
            end

            pToE = paramsToEstimate(startOpt:endOpt);
            Bnds = repmat([0.01 1], length(pToE), 1); % 0.4793
            % disp(['Parameters: [' num2str(startOpt) ' ' num2str(endOpt) ']'])
            % 
            %% Optimization 
            for nOpt = 1 %:2

                IV = bMatrix((eInd - (eInd - rCounter)), startOpt:endOpt);
                estimatedParams = estimatedInfo(pToE, 'InitialValue', IV, 'Bounds', Bnds);
                
                % tic
                fitConst = sbiofit(m, gDataTemp, responseMap, estimatedParams, [], method, options, 'UseParallel', true);

                pNum = 1;
                pVec = pStart:pGap:numel(m.parameters);
                for k = startOpt:length(pVec)
                    m.Parameters(pVec(k)).Value = fitConst.ParameterEstimates.Estimate(pNum);
                    if nOpt == 1    
                        if (~strcmp(m.Parameters(pVec(k)).Name, pToE(pNum)))
                            disp('Error! Parameter mismatch: ')
                            disp(m.Parameters(pVec(k)).Name), disp(pToE(pNum))
                            disp(['Index = ' num2str(k) ', pVec(k) = ' num2str(pVec(k)) ', and pNum = ' num2str(pNum)])
                        end
                    end
                    bMatrix((eInd - (eInd - rCounter)), k) = fitConst.ParameterEstimates.Estimate(pNum);
                    fitConstMAT(k) = fitConst.AIC;
                    pNum = pNum + 1;
                end

            end
            
%             disp(['run number: ' num2str(nOpt)]) 
%             disp(bMatrix(eInd, :))
            startOpt = startOpt + 3;
            if startOpt > tL
                startOpt = tL;
            end
            
        end

        %disp(bMatrix((eInd - (eInd - rCounter)), :))
        rCounter = rCounter + 1;

    end
    
    %disp(bMatrix(:, metaStart:tL)')

end