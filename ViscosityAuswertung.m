viscosity = [0.001, 0.0025, 0.005, 0.0075,0.01,0.0125,0.025,0.0375,0.05, 0.075 0.1];


% In die Zeilen kommen die Werte für eine Viskosität rein!
ErgebnisMatrix2 = NaN(11,40);

for i = 1 : 11
    
    load(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/TrajektorienMu/ViscosityResonanz/ViscosityResonanz' num2str(i) '.mat']);
    
    parfor j = 1 : 40
        
        [c,m] = ShakerDefaultFuerSim.createTestingConfig(15+(j-1)*5,3,0.25,1,.3);
        
        t = SimResults_Time{j};
        y = SimResults_Y{j};
        
        
        [~,MeanAmp,~,~] = c.getOutputOfInterest(t,y);
        if numel(MeanAmp) == 1
        ErgebnisMatrix2(i,j) = MeanAmp;
            
        else
        a = MeanAmp(2);
        ErgebnisMatrix2(i,j) = a;
        end
    end
    
    save(['/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/ErgebnisMatrizen/ErgebnisMatrix2Viscosity.mat'],'ErgebnisMatrix2');
    
end