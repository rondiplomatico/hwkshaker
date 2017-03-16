%% Variablen festsetzen

% Variable fuer die Aktivierungen
% 
activations = 0 : .05 : 1;
% 
%Fuer die Frequenzen, wird später automatisch mit Werten befüllt, NaN
%werden nicht geplottet
frequencies = NaN(23,1);

%Die Ergebnismatrix initialisieren
ErgebnisMatrix = NaN(length(activations),length(frequencies));




%% Schleife ueber die Daten
for i = 1 : 23
    % Dateien laden
    load(['/home/kraschewski/SimulationsErgebnisse/FrequenzAktivierungMu13_' num2str(i) '_3mm2.mat'])
    
    % Jeweils die gewünschte Daten extrahieren und in eine Matrix
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

save('/home/kraschewski/SimulationsErgebnisse/ErgebnisMatrixFrequenzAktivierung3mm_Mu13_5','ErgebnisMatrix');