NoOfSteps = 21;

SimResults_Time = cell(NoOfSteps,1);

Config = cell(NoOfSteps,1);

SimResults_Y = cell(NoOfSteps,1);

SimResults_T = cell(NoOfSteps,1);

for j = 1 : 10
    
    parfor i = 1 : 21
        
        [c,m] = ShakerDefaultFuerSim.createTestingConfig(30+(j-1)*5,3,0 + (i-1)*.05,1,.3);
        m.dt = .25;
        m.DefaultMu(2) = 5;
        
        %% The Simulation Part with saving the used time into a vector
        [t,y] = m.simulate;
        
        if t(end) < 100
            m.ODESolver.RelTol = 10^(-6);
            [t,y] = m.simulate;
        end
        
        if t(end) < 100
            m.ODESolver.RelTol = 10^(-10);
            m.setGaussIntegrationRule(5);
            [t,y] = m.simulate;
        end
        
        SimResults_Y{i,1} = y;
        SimResults_T{i,1} = t;
        Config{i} = c;
        
        
        
    end
    save(['FrequenzAktivierungSmall' num2str(j) '.mat'],'SimResults_T','SimResults_Y','Config');
    
end

ErgebnisMatrix = zeros(10,21);
%% Auswertung
for i = 1 : 10
    load(['FrequenzAktivierungSmall' num2str(i)])
    % Jeweils die gewünschte Daten extrahieren und in eine Matrix
    % speichern
    for j = 1 : 21
        
        c = Config{j};
        t = SimResults_T{j};
        y = SimResults_Y{j};
        
        [~,~,~,meanDifference] = c.getOutputOfInterest(t,y);
        
        % In einer Zeile stehen die Ergebnisse zu einer Datei
        ErgebnisMatrix(i,j) = meanDifference;
        
        
    end
    
    
end

save('ErgebnisMatrixSmall3mmDefaultMu','ErgebnisMatrix');