ErgebnisMatrixNeu = NaN(10,21);
load('/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/ErgebnisMatrizen/ErgebnisMatrix3mmComplete');
for i = 1 : 10
    
    save(['/home/kraschewski/SimulationsErgebnisse/FrequenzAktivierungSmall' num2str(i) '_3mm.mat'])
    
    for j = 1 : 21
        [c,m] = ShakerDefaultFuerSim.createTestingConfig(30+(i-1)*5,3,0+(j-1)*(0.05),1);
        t = SimResults_Time{j,1};
        y = SimResults_Y{j,1};
        [~,~,~,MeanDifference] = c.getOutputOfInterest(t,y);
        
        ErgebnisMatrixNeu(i,j) = MeanDifference;
    end
    
end

ErgebnisMatrixNeu = ErgebnisMatrixNeu(1:10,1:21);
ErgebnisMatrix3mmFinal = [ErgebnisMatrixNeu;Ergebnis];
save('/home/kraschewski/SimulationsErgebnisse/Amplitude3mm/ErgebnisMatrizen/ErgebnisMatrix3mmFinal','ErgebnisMatrix3mmFinal');