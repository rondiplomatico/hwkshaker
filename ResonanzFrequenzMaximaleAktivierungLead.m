% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden
NoOfSteps = 22;

% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(NoOfSteps,1);

Config = cell(NoOfSteps,1);

Model = cell(NoOfSteps,1);

SimResults_Y = cell(NoOfSteps,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(NoOfSteps,1);

parfor i = 1 : 22
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(160+(i-1)*5,3,1,1,.3);
    m.dt = .25;
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 100
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 100
        m.ODESolver.AbsTol = 10^(-2);
        m.setGaussIntegrationRule(5);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    
    %% Saving the Results into a cell
    SimResults_Time{i,1}=t;
    SimResults_Y{i,1}=y;
    % Speichert Model und Configuration in einer Cell
    Config{i,1} = c;
    Model{i,1} = m;
    disp(['Aktivierung ist 1 und Frequenz ist: ' num2str(c.Fr)]);
end
save('/home/kraschewski/SimulationsErgebnisse/ResonanzFrequenzTest_3mm.mat','SimResults_Time','SimResults_Y','ElapsedTime','Config','Model');
