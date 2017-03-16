% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part
j = 1;

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');

j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');


j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');


j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');

j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');

j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');

j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');

j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1]; 
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');


j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1];
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');

j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1];
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');


j = j+1;

% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!

% Wie viele Steps gemacht werden


% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(40,1);

Config = cell(40,1);



SimResults_Y = cell(40,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(40,1);


viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1];
%% Schleifen fuer den Aktivierungs/Frequenzen Part

parfor i = 1 : 40
    
    
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(i-1)*5,3,0.25,1,.3);
    m.dt = .25;
    m.DefaultMu(1) = viscosity(j); %#ok<PFBNS>
    m.DefaultMu(2) = 5;
    
    %% The Simulation Part with saving the used time into a vector
    [t,y] = m.simulate;
    ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    %
    if t(end) < 90
        m.ODESolver.RelTol = 10^(-1);
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
    end
    %
    if t(end) < 90
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
    
    
end

%% Speichert die meisten Variablen in eine Datei "Results".
save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(j) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config');