clear all
close all
%% Parameter hier anpassen
timedotsActivation = [0 10 90 100 180 190 270 280 360];
% Zeitpunkt an dem Aktivierung erreicht
timedotsActivationReached = timedotsActivation(2:2:end);
% Zeitpunkt an dem lineare Aktivierung anfaengt
timedotsActivationStarted = timedotsActivation(1:2:end);
activationDots = [0 .3 .3 .5 .5 1 1 .5 .5];
% Unterschiedliche Aktivierungen
activationDotsRelevant = activationDots(2:2:end);
% Erster Durchlauf: 45Hz 3mm Amp 25% Aktivierung nach 1ms Toleranz 0.3

% [c,m] = ShakerDefaultFuerSim.createTestingConfig(45,3,0.25,1,.3);

% Setup wie in Paper Figure
[c,m] = ShakerDefaultFuerSim.createTestingConfig(95,3,timedotsActivation,activationDots,.3);
m.T = 360;

%% Simulation und Kraefte berechnen
[t,y] = m.simulate;
%
% [df,nf] = m.getResidualForces(t,y);

% m.plotGeometrySetup;
% hold on
% load('PaperFigure2Daten.mat');
%% Jede Zeile enthaelt Kraefte eines Knotens ueber die Zeit t
xforces = zeros(50,length(t));
yforces = zeros(50,length(t));
zforces = zeros(50,length(t));

% Jeder Rand besitzt 25 DOFs mit Knotenkraeften
% Ersten 25 gehoeren an den Linken Rand Zweiten 25 an den rechten
% (hoffentlich)
for i = 1 : length(t)
    xforces(1:50,i) = df(1:3:end,i); % Alle DirichletKomponenten in X Richtung
    yforces(1:50,i) = df(2:3:end,i); % Alle DirichletKomponenten in Y Richtung
    zforces(1:50,i) = df(3:3:end,i); % Alle DirichletKomponenten in Z Richtung
    
end
%% Vergleich welcher Rand groessere Kraefte aufweist offenbar immer abwechselnd
% Erst der linke Rand exemplarisch
xforces_left = xforces(1:25,:);
yforces_left = yforces(1:25,:);
zforces_left = zforces(1:25,:);

% Rechter Rand
xforces_right = xforces(26:50,:);
yforces_right = yforces(26:50,:);
zforces_right = zforces(26:50,:);

%Resultierender Kraftvektor 3 Koordinaten in Zeilen ueber die Zeit
F_res_left = zeros(3,length(t));
F_res_left(1,:) = sum(xforces_left);
F_res_left(2,:) = sum(yforces_left);
F_res_left(3,:) = sum(zforces_left);

F_res_right = zeros(3,length(t));
F_res_right(1,:) = sum(xforces_right);
F_res_right(2,:) = sum(yforces_right);
F_res_right(3,:) = sum(zforces_right);

F_res_gesamt = F_res_left + F_res_right;


%% Der Plot der Kraftvektoren an den Rändern über die Zeit
F_res_abs_left = zeros(1,length(t));
F_res_abs_right = zeros(1,length(t));
F_res_abs_gesamt = zeros(1,length(t));

% Vorzeichen in y-Richtung, da klar dominante Richtung
vorzeichenLinks = sign(F_res_left(2,:));
vorzeichenRechts = sign(F_res_right(2,:));
vorzeichenGesamt = sign(F_res_gesamt(2,:));
for i = 1 : length(t)
    %     quiver3([0,0],[0,100],[0,0],[F_res_left(1,i),F_res_right(1,i)],...
    %         [F_res_left(2,i),F_res_right(2,i)],[F_res_left(3,i),F_res_right(3,i)],1,'b', 'MarkerSize',14);
    %     pause(0.1);
    F_res_abs_left(i) = norm(F_res_left(:,i));
    F_res_abs_right(i) = norm(F_res_right(:,i));
    F_res_abs_gesamt(i) = norm(F_res_gesamt(:,i));
    
end
%
%plot(t,(y(3,:)));
hold on
% Die Resultierende Kraft wird mit einem Vorzeichen versehen welches der
% y-Richtung entspricht (positive Kraft -> positive y-Richtung. negativ
% entsprechend.
plot(t,vorzeichenLinks.*F_res_abs_left);
plot(t,vorzeichenRechts.*F_res_abs_right);

% 1. Zeile linke resultierende 2. Zeile rechte resultierende
F_res_abs = [vorzeichenLinks.*F_res_abs_left;vorzeichenRechts.*F_res_abs_right];

% gesamtkraft der Resultierenden an beiden Rändern zusammen
F_res_abs_gesamt = vorzeichenGesamt.*F_res_abs_gesamt;

plot(t,F_res_abs_gesamt);
% [Axes,h(1),h(2)] = plotyy(t,(y(1005,:)-50)/3,t,F_res_abs(1,:));
% h(3) = line(t,F_res_abs(2,:),'Parent',Axes(2));
% set(h(3),'color','r');

%% Mittelwerte der Amplituden und der Kraefte
hold off
proc = models.musclefibre.experiments.Processor;
% Minimale Amplitude abhängig von der Zeit machen, da sonst auch minima
% berücksichtig würden.
kraftMittelwerte = zeros(4,1);
timeTempIndices = cell(8,1);

t=t(t<=350);

%% Herausfiltern der Maxima
for i = 1 : length(activationDotsRelevant)
    % Zeitpunkte in denen das Extrema gesucht werden soll (bereits voll
    % aktiviert und bis die Aktivierung sich wieder aendert.
    timeTempIndices{i} = (t>= timedotsActivationReached(i) & t <= timedotsActivationStarted(i+1));
    proc.minV = max(F_res_abs_gesamt(timeTempIndices{i}))-3;
    
    % Gibt die lokalen Indices von dem Eingangssignal zurück. Das heißt ich
    % bekomme die Indices auf den Intervallen zurück, die aber nicht die
    % globalen Indices sind!!
    maximaTempIndices = proc.getPeakIdx(t(timeTempIndices{i}),...
        F_res_abs_gesamt(timeTempIndices{i}));
    
    % Index i kennzeichnet hier den einer Aktivierung
    % Der additive Term in der 2. Zeile ist die Korrektur von lokal auf
    % globalen Index. Stimmt. Teilweise faellt ein maximum in die
    % Aktivierungsphase!
    kraftMittelwerte(i) = mean(F_res_abs_gesamt(maximaTempIndices+...
        find(t==timedotsActivationReached(i))-1));
    
end

kraftMittelwerte(2) = 1/2*(kraftMittelwerte(2)+kraftMittelwerte(4));
kraftMittelwerte = kraftMittelwerte(1:3);

%% Berechnung der Amplituden Mittelwerte

amplitudenMittelwerte = zeros(4,1);
y_bar = y(1005,:)-50-y(3,:);

%% Herausfiltern der Maxima
% Wie eben nur anderes Signal TODO!!!!!!!!!!!

for i = 1 : length(activationDotsRelevant)
    % Zeitpunkte in denen das Extrema gesucht werden soll (bereits voll
    % aktiviert und bis die Aktivierung sich wieder aendert.
    timeTempIndices{i} = (t>= timedotsActivationReached(i) & t <= timedotsActivationStarted(i+1));
    proc.minV = max(y_bar(timeTempIndices{i}))-3;
    
    % Gibt die lokalen Indices von dem Eingangssignal zurück. Das heißt ich
    % bekomme die Indices auf den Intervallen zurück, die aber nicht die
    % globalen Indices sind!!
    maximaTempIndices = proc.getPeakIdx(t(timeTempIndices{i}),...
        y_bar(timeTempIndices{i}));
    
    % Index i kennzeichnet hier den einer Aktivierung
    % Der additive Term in der 2. Zeile ist die Korrektur von lokal auf
    % globalen Index. Stimmt. Teilweise faellt ein maximum in die
    % Aktivierungsphase!
    amplitudenMittelwerte(i) = mean(y_bar(maximaTempIndices+...
        find(t==timedotsActivationReached(i))-1));
    
end

amplitudenMittelwerte(2) = 1/2*(amplitudenMittelwerte(2) +amplitudenMittelwerte(4));
amplitudenMittelwerte = amplitudenMittelwerte(1:3);
activations = activationDotsRelevant(1:3);

%% Plot der Kräfte als Funktion der relativen Amplitude (Hysteresen)
close all
% Jeder Abschnitt eine Cell
legendliste = cell(8,1);
plotliste = cell(8,1);
yRegressionListe = cell(8,1);
xRegressionListe = cell(8,1);


for i = 1 : 8
    
    % Plot gestrichelt Aktivierungsphase
    if mod(i,2)==1
        timeTempIndices{i} = (t>= timedotsActivation(i) & t<= timedotsActivation(i+1));
        plotliste{i} = plot(y_bar(timeTempIndices{i}),F_res_abs_gesamt(timeTempIndices{i}),'black--');
        legendliste{i} = '{\boldmath$$a_{ramp}$$}';
        
        
        
        
        % Plot der Aktiviertenkurve 50% aktiviert
    elseif i==4 || i==8
        timeTempIndices{i} = (t>= timedotsActivation(i) & t <= timedotsActivation(i+1));
        plotliste{i} = plot(y_bar(timeTempIndices{i}),F_res_abs_gesamt(timeTempIndices{i}),'g');
        legendliste{i} = '{\boldmath$$a = 0.5 \ \ \ k = 9.4 \frac{N}{mm} \ \ \ E_{diss} = 9 mJ$$}';
        yRegressionListe{i} = F_res_abs_gesamt(timeTempIndices{i});
        xRegressionListe{i} = y_bar(timeTempIndices{i});
        
        
        % Plot der anderen Aktivierung
    else
        timeTempIndices{i} = (t>= timedotsActivation(i) & t<= timedotsActivation(i+1));
        plotliste{i} = plot(y_bar(timeTempIndices{i}),F_res_abs_gesamt(timeTempIndices{i}));
        
        yRegressionListe{i} = F_res_abs_gesamt(timeTempIndices{i});
        xRegressionListe{i} = y_bar(timeTempIndices{i});
    end
    
    % Plot der Gesamtkraft!
    %         F_res_abs(1,(i-1)*180+1:i*180):%...
    %         +F_res_abs(2,(i-1)*180+1:i*180));
    hold on
    
%     waitforbuttonpress;
end
xlabel('{\boldmath$$\bar{y}(t)\ [mm]$$}','Interpreter','latex');
ylabel('{\boldmath$$F_R\ [N]$$}','Interpreter','latex');
legendliste{2} = '{\boldmath$$a = 0.3 \ \ \ k = 7.3 \frac{N}{mm} \ \ \ E_{diss} = 18 mJ$$}';
legendliste{6} = '{\boldmath$$a = 1.0 \ \ \ k = 12.9 \frac{N}{mm} \  E_{diss} = 5 mJ$$}';

h = legend([plotliste{1} plotliste{2:2:6}],{legendliste{1},legendliste{2:2:6}},'location','southeast','fontsize',16);
%     '{\boldmath$$0.5 < \alpha \leq 1$$}','{\boldmath$$1 > \alpha \ge 0.5$$}',...
%     'location','best');
set(h,'Interpreter','latex');
%
set(gca,'FontSize',25);

% h = legend('$$\bar{y}(t)$$','$$F_RLinks$$','$$F_{R}Rechts$$');
% set(h,'Interpreter','latex');


    %% Regressionen machen.. Werte so kleiner als die vom plot tool linear fit.
    
    % Hier wird quasi polyfit gemacht...
    
    % Zusammenführen der gleichen Aktivierungen Index 4 und 8
    
    xRegressionListe{4} = [xRegressionListe{4}';xRegressionListe{8}'];
    yRegressionListe{4} = [yRegressionListe{4}';yRegressionListe{8}'];
    regressionsGeradenListe = cell(3,1);
    
    
    % Nur für die bereits voll aktivierten Abschnitte
    koeffizientenListe = cell(3,1);
    j = 1;
    for k = 2 :2 : 6
        
        if k == 4
            X = [ones(length(xRegressionListe{k}),1), xRegressionListe{k}];
            Y = yRegressionListe{k};
            
            koeffizientenListe{k-j} = X\Y;
        else
            % Für den X Vektor noch die einser Spalte hinzufügen
            X = [ones(length(xRegressionListe{k}'),1), xRegressionListe{k}'];
            Y = yRegressionListe{k}';
            
            koeffizientenListe{k-j} = X\Y;
        end
        % ax+b
        b = koeffizientenListe{k-j}(1);
        a = koeffizientenListe{k-j}(2);
        
        regressionsGeradenListe{k-j} =@(x) a.*x+b;
        plot(y_bar(timeTempIndices{k}),regressionsGeradenListe{k-j}(y_bar(timeTempIndices{k})));
        
        j = j+1;
    end