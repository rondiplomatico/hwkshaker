%% Vorinitialisierung
% Insgesamt werden die generierten Daten in cells und Vektoren gespeichert
% dabei sind alle Daten, die zusammengehoeren, in der gleichen Zeile der
% Datenstruktur gespeichert!


NoOfSteps = 21; % Wie viele Steps gemacht werden

% Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
SimResults_Time = cell(NoOfSteps,1);

Config = cell(NoOfSteps,1);

Model = cell(NoOfSteps,1);

SimResults_Y = cell(NoOfSteps,1);
% Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% fruehzeitige Abbrueche zu erkennen
ElapsedTime = cell(NoOfSteps,1);




%% Schleifen fuer den Aktivierungs/Frequenzen Part
for j = 1 : 10
    
    
    parfor i = 1 : 21
        [c,m] = ShakerDefaultFuerSim.createTestingConfig(30+(j-1)*5,3,0+(i-1)*(0.05),1);
        m.dt = .25;
        m.DefaultMu(2) = 5;
        %% The Simulation Part with saving the used time into a vector
        [t,y] = m.simulate;
        ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        
        if t(end) < 100
            m.ODESolver.RelTol = 10^(-6);
            [t,y] = m.simulate;
            ElapsedTime{i,1} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
        end
        
        if t(end) < 100
            m.ODESolver.RelTol = 10^(-10);
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
        disp(['Aktivierung ist ' num2str(0+(i-1)*(0.05)) ' und Frequenz: ' num2str(c.Fr)]);
    end
    %% Speichert die meisten Variablen in eine Datei "Results".
    save(['/home/kraschewski/SimulationsErgebnisse/FrequenzAktivierungSmall' num2str(j) '_3mm.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config','Model');
end
% 
% NoOfSteps = 8; % Wie viele Steps gemacht werden
% 
% % Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% % der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
% SimResults_Time = cell(NoOfSteps,1);
% 
% Config = cell(NoOfSteps,1);
% 
% Model = cell(NoOfSteps,1);
% 
% SimResults_Y = cell(NoOfSteps,1);
% % Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% % fruehzeitige Abbrueche zu erkennen
% ElapsedTime = cell(NoOfSteps,1);
% 
% %% Feste Frequenz unterschiedliche Amplituden und Aktivierung
% for i = 1 : 4
%     parfor j = 1 : 8
%         [c,m] = ShakerDefaultFuerSim.createTestingConfig(120,(12-(j-1)),1-(i-1)*.25,1);
%         m.dt = .25
%         m.DefaultMu(2) = 5;
%         %% The Simulation Part with saving the used time into a vector
%         [t,y] = m.simulate;
%         ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
% 
%         if t(end) < 120
%             m.ODESolver.RelTol = 10^(-6);
%             [t,y] = m.simulate;
%             ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         end
%         
%         if t(end) < 120
%             m.ODESolver.RelTol = 10^(-10);
%             [t,y] = m.simulate;
%             ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         end
%         
%         %% Saving the Results into a cell
%         SimResults_Time{j}=t;
%         SimResults_Y{j}=y;
%         % Speichert Model und Configuration in einer Cell
%         Config{j} = c;
%         Model{j} = m;
%     end
% 
%   save(['/home/kraschewski/SimulationsErgebnisse/Amplitude_ActivationChange2' num2str(i) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config','Model');
% end

% NoOfSteps = 9; % Wie viele Steps gemacht werden
% 
% % Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% % der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
% SimResults_Time = cell(NoOfSteps,1);
% 
% Config = cell(NoOfSteps,1);
% 
% Model = cell(NoOfSteps,1);
% 
% SimResults_Y = cell(NoOfSteps,1);
% % Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% % fruehzeitige Abbrueche zu erkennen
% ElapsedTime = cell(NoOfSteps,1);
% 
% 
% %% Mu(13) geaendert mit fester Frequenz und unterschiedlichen Amplituden
% for i = 1 : 2
%     parfor j = 1 : 9
%         [c,m] = ShakerDefaultFuerSim.createTestingConfig(120,(12-(j-1)),.2,1);
%         m.dt = .25
%         m.DefaultMu(2) = 5;
%         switch i
%             case 1
%                 m.DefaultMu(13) = .5
%             case 2
%                 m.DefaultMu(13) = .4
%         end
%         %% The Simulation Part with saving the used time into a vector
%         [t,y] = m.simulate;
%         ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         
%         if t(end) < 120
%             m.ODESolver.RelTol = 10^(-6);
%             [t,y] = m.simulate;
%             ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         end
%         
%         if t(end) < 120
%             m.ODESolver.RelTol = 10^(-10);
%             [t,y] = m.simulate;
%             ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         end        
% 
%         %% Saving the Results into a cell
%         SimResults_Time{j}=t;
%         SimResults_Y{j}=y;
%         % Speichert Model und Configuration in einer Cell
%         Config{j} = c;
%         Model{j} = m;
%     end
% 
%   save(['/home/kraschewski/SimulationsErgebnisse/Mu13_AmplitudeChange2' num2str(i) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config','Model');
% end
% 
% NoOfSteps = 5; % Wie viele Steps gemacht werden
% 
% % Speichert im Eintrag i des Vektors die Dauer der Simulation fuer die i.
% % Iteration
% SimulationTime = cell(NoOfSteps,1);
% 
% % Speichert in cell Spalte 1 die Zeiten und in der Spalte 2 die Ergebnisse
% % der Rechnung. Zeilenweise sind die untersch. Simulationen gespeichert
% SimResults_Time = cell(NoOfSteps,1);
% 
% Config = cell(NoOfSteps,1);
% 
% Model = cell(NoOfSteps,1);
% 
% SimResults_Y = cell(NoOfSteps,1);
% % Der letzte Zeitpunkt der Rechnung wird gespeichert um eventuelle
% % fruehzeitige Abbrueche zu erkennen
% ElapsedTime = cell(NoOfSteps,1);
% 
% 
% %% Mu(1) geaendert mit festen Amplituden und variabler Frequenz
% for i = 1 : 7
%     parfor j = 1 : 18
%         [c,m] = ShakerDefaultFuerSim.createTestingConfig(200-(j-1)*(5),12,.2,1);
%         m.dt = .25
%         m.DefaultMu(2) = 5;
%         switch i
%             case 1
%                 m.DefaultMu(1) = .01;
%             case 2
%                 m.DefaultMu(1) = .05;
%             case 3
%                 m.DefaultMu(1) = .1;
%             case 4
%                 m.DefaultMu(1) = .5;
%             case 5
%                 m.DefaultMu(1) = 1;
%             case 6
%                 m.DefaultMu(1) = 1.5;
%         end
%         %% The Simulation Part with saving the used time into a vector
%         [t,y] = m.simulate;
%         ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         
%         if t(end) < 120
%             m.ODESolver.RelTol = 10^(-6);
%             [t,y] = m.simulate;
%             ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         end
%         
%         if t(end) < 120
%             m.ODESolver.RelTol = 10^(-10);
%             [t,y] = m.simulate;
%             ElapsedTime{j} = t(size(t,2)); % Um vorzeitigen Abbruch zu erkennen
%         end        
% 
%         %% Saving the Results into a cell
%         SimResults_Time{j}=t;
%         SimResults_Y{j}=y;
%         % Speichert Model und Configuration in einer Cell
%         Config{j} = c;
%         Model{j} = m;
%     end
% 
%   save(['/home/kraschewski/SimulationsErgebnisse/Mu1_AmplitudeChange2' num2str(i) '.mat'],'SimResults_Time','SimResults_Y','ElapsedTime','Config','Model');
% end