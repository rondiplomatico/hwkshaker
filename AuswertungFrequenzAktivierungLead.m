%% Variablen festsetzen

% Variable fuer die Aktivierungen
% 
activations = 0 : .05 : 1;
% 
%Fuer die Frequenzen, wird sp�ter automatisch mit Werten bef�llt, NaN
%werden nicht geplottet
frequencies = NaN(39,1);

%Die Ergebnismatrix initialisieren
ErgebnisMatrix = NaN(length(activations),length(frequencies));




%% Schleife ueber die Daten
for i = 1 : 39
    % Dateien laden
    load(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienDefaultMu/FrequenzAktivierung' num2str(i) '_3mm.mat'])
    
    % Jeweils die gew�nschte Daten extrahieren und in eine Matrix
    % speichern
    for j = 1 : 21
        
        c = Config{j};
        t = SimResults_Time{j};
        y = SimResults_Y{j};
        
        [~,~,~,meanDifference] = c.getOutputOfInterest(t,y);
        
        % In einer Zeile stehen die Ergebnisse zu einer Datei
        ErgebnisMatrix(i,j) = meanDifference;
        
        
    end
    
end

save('/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/ErgebnisMatrizen/ErgebnisMatrix3mmComplete','ErgebnisMatrix');