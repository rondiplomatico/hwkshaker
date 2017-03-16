%%Testet die Auswirkung der Ortsdiskretisierung
model = cell(10,1);
config = cell(10,1);
SimResults_Time = cell(10,1);
SimResults_Y = cell(10,1);
parfor i = 1 : 10
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(125,12,.2,1,.25+(i-1)*(.025));
    m.dt = .25;
    m.DefaultMu(2) = 5;
    model{i} = m;
    config{i} = c;
    [t,y] = m.simulate;
    
    
    if t(end) < 100
        m.ODESolver.RelTol = 10^(-6);
        [t,y] = m.simulate;
    end
    a = length(t);
    SimResults_Time{i} = t;
    SimResults_Y{i} = y;
    
    
end
save('/home/kraschewski/SimulationsErgebnisse/FunktionierendeConfig1','SimResults_Time','SimResults_Y','config','model');
SimResults_T = zeros(300,1);
for i = 1 : 10
    t = SimResults_Time{i};
    a = length(t);
    SimResults_T(1:a,i) = t;
end
[~,index] = find(max(max(SimResults_T))==SimResults_T);
c = config{index};
lastTolerance = c.TOL;

%%Testet die Auswirkung der Relativen Solver Toleranz, nimmt dabei die
%%letzte/beste Diskretisierung der vorangegangenen Reihe
parfor i = 1 : 5
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(125,12,.2,1,lastTolerance);
    m.dt = .25;
    m.DefaultMu(2) = 5;
    
    model{i} = m;
    config{i} = c;
    [t,y] = m.simulate;
    
       
    if t(end) < 100
        m.ODESolver.RelTol = 3*10^(-i);
        [t,y] = m.simulate;
    end
    
    SimResults_Time{i} = t;
    SimResults_Y{i} = y;      
    
end

save('/home/kraschewski/SimulationsErgebnisse/FunktionierendeConfig2','SimResults_Time','SimResults_Y','config','model');

%% Analog fuer die absoulte Toleranz
parfor i = 1 : 5
    [c,m] = ShakerDefaultFuerSim.createTestingConfig(125,12,.2,1,lastTolerance);
    m.dt = .25;
    m.DefaultMu(2) = 5;
    [t,y] = m.simulate;
    
    if t(end) < 100
        m.ODESolver.AbsTol = 3*10^(-i);
        [t,y] = m.simulate;
    end
    
    SimResults_Time{i} = t;
    SimResults_Y{i} = y;     
    
    
end
save('/home/kraschewski/SimulationsErgebnisse/FunktionierendeConfig3','SimResults_Time','SimResults_Y','config','model');
