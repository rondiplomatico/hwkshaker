classdef ShakerDefaultFuerSim < models.muscle.AMuscleConfig
    
    properties
        Fr = 50; % Shaker frequency [Hz]
        Amp = 2; % Shaker amplitude [mm]
        % New Parameters for frequency, a tolerance for the maxDy of
        % TMRFunction (for the discretisation part) and the maximal amount
        % of points on y-axis(not the total amount of nodes) that are usable
        TOL;
        maxN;
        maxYLength;
        
        % New Property to define "other" stretchfunctions
        % is a vector where the first argument sets the end of the constant tendon part and
        % the second argument sets the start of the const muscle part.
        % Parameters must be between 0 and ylen
        constFromTo;
    end
    %GEAENDERT! Properties von Private auf Public (Fuer das Skript zum Testen)
    properties(SetAccess=private)
        ylen;
        radfun;
        stretchfun;
        
        % Neues Attribut, das nur fuer die TMR Berechnung verwendet wird
        radForTMR
    end
    
    methods
        function this = ShakerDefaultFuerSim(varargin)
            this = this@models.muscle.AMuscleConfig(varargin{:});
            
            % Default Stretch is exp
            this.addOption('Stretch','exp');
            % Default tolerance 0.3
            this.addOption('TOL',0.3);
            % Default maxYLength to 20
            this.addOption('maxYLength',20);
            % Default N points on yAxis (not the total amount of nodes!)
            this.addOption('maxN',500);
            
            this.addOption('Frequency',50);
            
            this.addOption('Amplitude',2);
            
            
            
            % Default constFromTo is zero
            this.addOption('constFromTo',[0,0]);
            
            % Compute outer shape
            this.ylen = 100;
            k = kernels.GaussKernel(20);
            
            % Wichtig!! Sollte diese Funktion hier geaendert werden ist
            % diese unter Umstaenden in der
            % getTendonMuscleRatio(this,Points) Methode ebenfalls
            % anzupassen!
            this.radfun = @(x)10*k.evaluate(this.ylen/2,x)+2;
            this.radForTMR = @(x)10*k.evaluate(this.ylen/2,x)+2;
            
            % Also invokes geo setup
            this.init;
            
            this.maxYLength=this.Options.maxYLength;
        end
        
        function prepareSimulation(this, ~, ~)
            % Sets the shaking boundary conditions according to currently
            % set frequency (Fr) and amplitude (Amp)
            this.VelocityBCTimeFun = general.functions.Sinus('Frequency',this.Fr,'Amplitude',this.Amp);
        end
        
        function plotTMR(this)
            y = linspace(0,this.ylen,1000);
            tmr = this.getTendonMuscleRatio([zeros(size(y)); y]);
            pm = PlotManager;
            pm.LeaveOpen = true;
            ax = pm.nextPlot('tmr','Prescribed stretches and resulting TMR','y-axis [mm]','Value');
            plot(ax, y,this.stretchfun(y),'r',y,tmr,'b');
            legend('Stretch','TMR');
            pm.done;
        end
        
        function configureModel(this, m)
            configureModel@models.muscle.AMuscleConfig(this, m);
            m.T = 75;
            m.dt = .5;
            
            mu = m.DefaultMu;
            % Small viscosity
            mu(1) = .005;
            m.DefaultMu = mu;
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only.
            
            f = this.Model.System.f;
            % This is empty if no prepareSimulation has been called yet
            if isempty(f.AnisoPassiveTendon)
                tmr = zeros(1,size(points,2));
                return;
            end
            
            % Berechne zunaechst die Kraft die an der Sehne aufgebracht
            % werden muss, um die bei x=0 vorgegebene Sehnenverzerrung zu erreichen
            
            F = f.AnisoPassiveTendon(this.stretchfun(0))*pi*this.radForTMR(0)^2;
            
            % Flaeche als Funktion des Ortes
            A = @(y)this.radForTMR(y).^2*pi;
            
            % Extrahiere y-Koordinaten aus den gegebenen Punkten
            ycoord = points(2,:);
            
            % Berechne Stress der einzelkomponenten bei gegebenem Stretch
            % an den gegebenen y-Positionen
            musclestress = f.AnisoPassiveMuscle(this.stretchfun(ycoord));
            tendonstress = f.AnisoPassiveTendon(this.stretchfun(ycoord));
            
            % Umstellen von
            % F = \sigma*A = (1-r)\sigma_m + r\sigma_t
            % nach r ergibt
            % r = (F/A - }\sigma_m) / (sigma_t - sigma_m)
            tmr = (F./A(ycoord) - musclestress)./(tendonstress-musclestress);
            
        end
        
        
        %%  Discr due to function values, TMRFunc must be a function handle
        function [ypoints] = getDiscr(this,TMRFunc)
            % Kurze Fehlerbehandlung
            if ~isa(TMRFunc,'function_handle')
                error('The function argument must be of type function_handle');
            end
            if this.maxYLength*this.maxN*2 < this.ylen
                error('Can not reach maxYLength with this amount of points');
            end
            
            
            %% Vorinitialisierung
            N = this.maxN;
            N = N/2; % Nur bis zur haelfte
            ypoints = zeros(N,1);
            distance = zeros(N,1); % Speichere die Abstaende um spaeter an der Mitte zu spiegeln
            j = 2; % Erster Punkt ist die 0, fuer den Index in distance und ypoints
            yDefault = linspace(0,this.ylen/2,1000);
            TMRFuncComp = TMRFunc(yDefault); % Vergleichswerte fuer spaeter
            refValue = TMRFunc(0); % Referenzwert um die Toleranz zu checken
            posChange = 0; % Wird in der Hauptlogik benoetigt um die variable Laenge der YPunkte zu beruecksichtigen
            
            %% Abschnitt um Genauigkeit in den Vergleichsfunktionswerten zu erreichen
            compare = abs(TMRFuncComp(2:length(TMRFuncComp))-TMRFuncComp(1:length(TMRFuncComp)-1))>this.TOL;
            compare = find(compare);
            while ~isempty(compare)
                for k = 1 : length(compare)
                    yDefault = [yDefault(1:compare(k)) (yDefault(compare(k))+yDefault(compare(k)+1))/2 yDefault(compare(k)+1:length(yDefault))];
                end
                TMRFuncComp = TMRFunc(yDefault);
                compare = abs(TMRFuncComp(2:length(TMRFuncComp))-TMRFuncComp(1:length(TMRFuncComp)-1))>this.TOL;
                compare = find(compare);
            end
            
            
            %% Hauptlogik
            for i = 1 : length(yDefault);
                % Erster Punkt der TOL ueberschreitet -> Vorheriger in ypoints
                if(abs(TMRFuncComp(i)-refValue)>this.TOL && i<=length(yDefault))
                    refValue = TMRFuncComp(i-1);
                    ypoints(j) = yDefault(i-1);
                    distance(j-1) = ypoints(j)-ypoints(j-1);
                    j = j+1;
                    if j>N % Alle Punkte verwendet?
                        
                        if(abs(refValue-TMRFunc(this.ylen/2))>this.TOL)
                            error('Tolerance can not be reached with this amount of points')
                        end
                        if(abs(max(distance))>this.maxYLength)
                            error('Can not fall below maxYLength with this amount of points');
                        end
                        if (2*abs(ypoints(N)-this.ylen)>this.maxYLength)
                            error('can not fall below maxYLength with this amount of points');
                        end
                    end
                    
                    % Toleranz wird am Ende nicht mehr ueberschritten -> verfuegbarer Punkt wird Mittelpunkt
                elseif i+1==length(yDefault)
                    ypoints(j) = this.ylen/2;
                    ypoints = ypoints(1:j);
                    distance(j-1) = ypoints(j)-ypoints(j-1);
                    distance = distance(1:j-1);
                end
            end
            %% Bis hier sind die Punkte gesetzt um die Toleranz einzuhalten!
            %%weiter sind hier noch Punkte uebrig um zu Verfeinern ->
            %%maxYLength einhalten
            k = find(distance>this.maxYLength);
            usablePoints = 2*N-2*j+1;
            while max(distance>this.maxYLength)==1 && length(k)*2 <= usablePoints
                
                % Iteration ueber die Intervalle in denen maxX ueberschritten wird, zwischen die Punkte einene weiteren setzen.
                for i=1:length(k)
                    ypoints = [ypoints(1:k(i)+posChange); (ypoints(k(i)+posChange)+ypoints(k(i)+1+posChange))/2; ypoints(k(i)+1+posChange:length(ypoints))]; %y veraendert sich nach einem Schleifendurchlauf -> PosChange+1
                    distance = [distance(1:k(i)+posChange-1);distance(k(i)+posChange)/2; distance(k(i)+posChange)/2 ; distance(k(i)+1+posChange:length(distance))]; %Abstand entsprechend angepasst -> halbiert dafuer 2mal und posChange+2
                    posChange = posChange+1;
                    j = j+1; % Fuer das Spiegeln an der Mitte
                end
                usablePoints = usablePoints-2*length(k); % length(k) Punkte auf einer Haelfte gesetzt
                k = find(distance>this.maxYLength);
                posChange = 0;
            end
            if(max(k)>=1)
                error('Can not fall below maxYLength with this amount of points');
            end
            
            % Zum "Spiegeln" an der Mitte
            ypointsTemp = zeros(length(ypoints),1);
            distance = flipud(distance);
            ypointsTemp(1) = ypoints(j);
            ypointsTemp(length(ypoints)) = ypoints(j)+this.ylen/2;
            for i=2 : length(ypoints)-1
                % Mit distance wird dafuer gesorgt, dass die Abstaende zwischen den Punkten rechts erhalten bleibt
                ypointsTemp(i) = ypointsTemp(i-1)+distance(i-1);
            end
            
            ypoints = [ypoints(1:j-1) ; ypointsTemp];
            
            %% Distanz wird maximal kleiner, falls Stueckweise Konstant,
            %%setze an die entsprechenden Stellen noch ein Element:
            % Scheint ein Problem mit der plotGeometrySetup hervorzurufen
            
            
            % sollte an den Stellen nicht zufaellig ein Punkt gesetzt
            % worden sein, setze einen.
            if ~any(ypoints == this.constFromTo(1))
                
                [~,index] = min(abs(ypoints-this.constFromTo(1)));
                ypoints = [ypoints(1:index); this.constFromTo(1); ypoints(index+1:end)];
                
            end
            
            
            if ~any(ypoints == this.constFromTo(2))
                
                [~,index] = min(abs(ypoints-this.constFromTo(2)));
                ypoints = [ypoints(1:index); this.constFromTo(2); ypoints(index+1:end)];
                
            end
            
            % Um die Symmetrie zu erhalten. Es wird aktuell noch ein
            if ~any(ypoints == this.constFromTo(2)+2*abs(this.ylen/2 - this.constFromTo(2))-.1)
                
                [~,index] = min(abs(ypoints - (this.constFromTo(2) + 2*abs(this.ylen/2 - this.constFromTo(2)))-.1));
                ypoints = [ypoints(1:index); this.constFromTo(2)+2*abs(this.ylen/2 - this.constFromTo(2))-.1...
                    ;ypoints(index+1:end)];
                
            end
            
            if ~any(ypoints == this.ylen-this.constFromTo(1))
                
                [~,index] = min(abs(ypoints - (this.ylen - this.constFromTo(1))));
                ypoints = [ypoints(1:index); this.ylen - this.constFromTo(1); ypoints(index+1:end)];
                
            end
            % Da ich nur "links" einfuege, wird am Ende der Vektor sortiert
            ypoints = sort(ypoints);
            
            %% Kleine Ausgaben fuer den Benutzer
            disp(['Maximal distance between yPoints is ',num2str(max(distance)) ,' mm']);
            % Da zwischen den Punkten immer noch ein Punkt gesetzt wird
            disp(['Tolerance of ' num2str(this.TOL) ' is reached and ' num2str(length(ypoints)*2+1) ' points are used']);
        end
        
        
        
        %% Plots the Position of chosen y-points over the time
        % here we've chosen the first, the mid, and the last point
        % Plottet die entsprechenden Trajektorien und bestimmt die
        % maximalen Amplituden, die Zeitpunkte an denen diese Auftreten und
        % den time-shift zwischen dem mittleren Punkt und den R�ndern
        
        function [timeShift,meanAmp_Mid,maxAmpTime,meanDifference] = getOutputOfInterest(this,t,y)
%             close all
            %% Processor fuer die Auswertung der Amplituden/Peaks
            proc = models.musclefibre.experiments.Processor;
            proc.minV = this.Amp/2;
            
            
            
            %% Bestimme die gewuenschten Punkte aus den Daten
            
            geo = this.Geometry;
            midIndex = (geo.NumElements/2)+1;
            midPoint = 3*geo.Elements(midIndex);
            
            firstPoint = 3;
            
            [last1,last2] = size(geo.Elements);
            lastPoint = 3*geo.Elements(last1,last2-3);
            
            Traj_Mid = y(midPoint,:)-y(midPoint,1);
            Traj_First = y(firstPoint,:);
            Traj_Last = y(lastPoint,:)-this.ylen;
            
%             % Die Kurven uebereinander gelegt
%             title('Alle Kurven uebereinandergelegt und Aktivierung');
%             hold on
%             plot(t,Traj_Mid,'g');
%             plot(t,Traj_Last,'b');
%             plot(t,Traj_First,'k');
%             xlabel('time [s]');
%             ylabel('Position on y-axe');
% %             hold off
            
            
            
            if t(end) > 20
                %% Extract the Peaks and Amplitudes
                % Vectors with all the Peak indices
                PeaksMid = proc.getPeakIdx(t,Traj_Mid);
                PeaksFirst = proc.getPeakIdx(t,Traj_First);
                PeaksLast = proc.getPeakIdx(t,Traj_Last);
                
                %% Um nur die voll aktivierten in die Daten aufzunehmen
                
                if length(PeaksFirst) > 4
                    % Falls die Laengen nicht passen -> "Ueberfluessiges" weglassen sd. keine
                    % Auswirkungen auf timeShift
                    if length(PeaksMid) ~= length(PeaksFirst)
                        if length(PeaksMid) > length(PeaksFirst)
                            PeaksMid = PeaksMid(1:length(PeaksFirst));
                        else
                            PeaksFirst = PeaksFirst(1:length(PeaksMid));
                        end
                    end
                    
                    % Pos Values -> Peak mid after Peak of first
                    timeShiftAbs = t(PeaksMid)-t(PeaksFirst);
                    % Mean of the time differences negative value means that
                    % the Mid "walks" behind the first/last Points
                    timeShift = mean(timeShiftAbs);
                    
%                     % Bereinige die Peaks um die Verschiebung alle um den Wert der Randverschiebung..?!
%                     for i = 1 : length(PeaksMid)-1
%                         % Ich teile die Trajektoren in Abschnitte Peak to
%                         % Peak auf und Verschiebe den jeweils rechten Peak
%                         % um die Differenz zwischen lokales Minimum und
%                         % Null der Raender
%                         TrajSplitter = Traj_First(PeaksFirst(i):PeaksFirst(i+1));
%                         Traj_First(PeaksFirst(i+1)) = Traj_First(PeaksFirst(i+1))-min(TrajSplitter);
%                         Traj_Mid(PeaksMid(i+1)) = Traj_Mid(PeaksMid(i+1))-min(TrajSplitter);
%                         Traj_Last(PeaksLast(i+1)) = Traj_Last(PeaksLast(i+1))-min(TrajSplitter);
%                     end
                    
                    % Die ersten Rauswerfen
                    if length(PeaksMid) < 4
                    else
                        
                        PeaksMid = PeaksMid(4:end);
                        PeaksFirst = PeaksFirst(4:end);
                        PeaksLast = PeaksLast(4:end);
                        
                    end
                    
                    % Jetzt ueberall die MW genommen, sowie auch negative
                    % Werte zugelassen...
                    
                    meanAmp_Mid = mean(Traj_Mid(PeaksMid));
                    
                    meanAmp_First = mean(Traj_First(PeaksFirst));
                    
                    meanAmp_Last = mean(Traj_Last(PeaksLast));
                    
                    %Interessant
                    meanDifference = mean(Traj_Mid(PeaksMid)-Traj_First(PeaksFirst));
                    
                    meanAmp = [meanAmp_First,meanAmp_Mid,meanAmp_Last];
                    
                    % Evtl. einfach den 2. Peak nach Aktivierung
                    % Max Amplitudes and time
                    maxAmp_Mid = max(Traj_Mid(PeaksMid));
                    maxAmp_MidTime = t(maxAmp_Mid == Traj_Mid);
                    
                    maxAmp_First = max(Traj_First(PeaksFirst));
                    maxAmp_FirstTime = t(maxAmp_First == Traj_First);
                    
                    maxAmp_Last = max(Traj_Last(PeaksLast));
                    maxAmp_LastTime = t(maxAmp_Last == Traj_Last);
                    
                    maxAmp = [maxAmp_First,maxAmp_Mid,maxAmp_Last];
                    maxAmpTime = [maxAmp_FirstTime,maxAmp_MidTime,maxAmp_LastTime];
                    
                    
                    maxDifference = (maxAmp_Mid-maxAmp_First)./maxAmp_First;
                    
                else
                    timeShift = NaN;
                    meanAmp = NaN;
                    maxAmpTime = NaN;
                    meanDifference = NaN;
                    
                end
                
            else
                timeShift = NaN;
                meanAmp = NaN;
                maxAmpTime = NaN;
                meanDifference = NaN;
                
            end
            
        end
        
        %% Calculates the TMRFunc outside of the getGeo Method
        function tmrFunc = getTMRFunc(this)
            
            %% Einige Initialisierungen
            m  = models.muscle.Model(models.muscle.examples.Belly);
            mu = m.DefaultMu;
            this.TOL = this.Options.TOL;
            this.maxN = this.Options.maxN;
            this.maxYLength = this.Options.maxYLength;
            this.constFromTo = this.Options.constFromTo;
            ref = [1.03, 1.4];
            this.Fr = this.Options.Frequency;
            this.Amp = this.Options.Amplitude;
            
            %% Stretchfun festlegen
            % Man kann nun ueber 'Gauss Gamma' einen GaussKernel als
            % Stretchfunktion verwenden, wichtig ist das Leerzeichen
            % zwischen Gauss und dem gewuenschten Gamma.
            % Bei groesseren Gamma -> Ungenauigkeiten im maximal Wert, fuehrt zur Unbrauchbarkeit!
            if length(this.Options.Stretch) >= 5 && strcmp(this.Options.Stretch(1:5),'Gauss')
                Gamma = str2double(this.Options.Stretch(6:length(this.Options.Stretch)));
                
                k = kernels.GaussKernel(Gamma);
                kexp = kernels.KernelExpansion;
                kexp.Kernel = k;
                % Funktion soll auf [0,1] definiert sein
                kexp.Centers.xi = 1;
                % Damit insgesamt das Maximum wieder bei ref(2) liegt
                kexp.Ma = (ref(2)-ref(1));
                
                % Der erste Term ist uA da um "Ungenauigkeiten" auszugleichen
                stretchfun_0_1 = @(x)(ref(1)-kexp.evaluate(0))+kexp.evaluate(x);
            else
                switch this.Options.Stretch
                    case 'exp'
                        stretchfun_0_1 =@(x) ref(1)+(ref(2)-ref(1)).*x.^4;
                    case 'log'
                        stretchfun_0_1 =@(x)((ref(2)-ref(1))/log(51)).*log(50.*(1/50+x))+ref(1);
                    case 'linear'
                        stretchfun_0_1 =@(x) (ref(2)-ref(1)).*x+ref(1);
                end
            end
            
            len = this.ylen;
            % Transfer the [0,1] argument version into the given domain
            % length (mirror at half)
            this.stretchfun = @(y)(y<=len/2) .* (stretchfun_0_1(2*y/len)) ...
                + (y>len/2) .* (stretchfun_0_1(1-2*(y-len/2)/len));
            
            % Sollten Konstante Abschnitte erwuenscht sein werden hier die
            % Anpassungen gemacht, um im TMR auch Konstante
            % Abschnitte zu haben.
            % Die Funktionen werden hierfuer auf das Intervall
            % [constFromTo(1),constFromTo(2)] transformiert. Von
            % [0,constFromTo(1)] haben wir nur Sehne, auf [constFromTo(2),ylen/2] nur
            % Muskel. Aufgrund der Symmetrie gilt entsprechendes ab ylen/2.
            
            if (norm(zeros(1,2) - this.constFromTo) ~= 0)
                
                k = kernels.GaussKernel(20);
                helpRadFun = @(x)10*k.evaluate(this.ylen/2,x)+2;
                helpRadFunLeft = @(x) (x <= this.ylen/2).*helpRadFun(x);
                helpRadFunRight = @(x) (x > this.ylen/2).*helpRadFun(x);
                
                this.radForTMR = @(x)  (x < this.constFromTo(1)).*helpRadFun(0)...
                    ... Linke Seite
                    +helpRadFunLeft((this.ylen/2 - 0)/(this.constFromTo(2) - this.constFromTo(1))*x +...
                    (0 - (this.ylen/2 - 0)/(this.constFromTo(2) - this.constFromTo(1)) * this.constFromTo(1)))...
                    .* ((x >= this.constFromTo(1)) .* (x < this.constFromTo(2)))...
                    + helpRadFun(this.ylen/2) .* (x >= this.constFromTo(2)).*...
                    (x < (this.constFromTo(2)+2 * abs(this.ylen/2 - this.constFromTo(2))))...
                    ...  Rechte Seite
                    +(x >= (this.constFromTo(2) + 2*abs(this.ylen/2 - this.constFromTo(2)))) .*...
                    helpRadFunRight((this.ylen-this.ylen/2) /...
                    (this.ylen - this.constFromTo(1) - (this.constFromTo(2) +...
                    2 * abs(this.ylen/2 - this.constFromTo(2)))).*x +...
                    this.ylen/2 - (this.ylen-this.ylen/2) /...
                    (this.ylen - this.constFromTo(1) - (this.constFromTo(2) +...
                    2 * abs(this.ylen/2 - this.constFromTo(2)))) * ...
                    (this.constFromTo(2) + 2*abs(this.ylen/2 - this.constFromTo(2))))...
                    .*(x <= this.ylen - this.constFromTo(1))...
                    +(x > this.ylen - this.constFromTo(1)).*helpRadFun(0);
                
                % Erleichtert die Transformationen fuer Konstante Abschnitte
                StretchLeft = @(y) (y<=len/2) .* (stretchfun_0_1(2*y/len));
                StretchRight = @(y) (y>len/2) .* (stretchfun_0_1(1-2*(y-len/2)/len));
                
                this.stretchfun = @(x) (x < this.constFromTo(1)).*ref(1)...
                    ... Linke Seite
                    +StretchLeft((this.ylen/2 - 0)/(this.constFromTo(2) - this.constFromTo(1))*x +...
                    (0 - (this.ylen/2 - 0)/(this.constFromTo(2) - this.constFromTo(1)) * this.constFromTo(1)))...
                    .* ((x >= this.constFromTo(1)) .* (x < this.constFromTo(2)))...
                    + ref(2) .* (x >= this.constFromTo(2)).*...
                    (x < (this.constFromTo(2)+2 * abs(this.ylen/2 - this.constFromTo(2))))...
                    ...  Rechte Seite
                    +(x >= (this.constFromTo(2) + 2*abs(this.ylen/2 - this.constFromTo(2)))) .*...
                    StretchRight((this.ylen-this.ylen/2) /...
                    (this.ylen - this.constFromTo(1) - (this.constFromTo(2) +...
                    2 * abs(this.ylen/2 - this.constFromTo(2)))).*x +...
                    this.ylen/2 - (this.ylen-this.ylen/2) /...
                    (this.ylen - this.constFromTo(1) - (this.constFromTo(2) +...
                    2 * abs(this.ylen/2 - this.constFromTo(2)))) * ...
                    (this.constFromTo(2) + 2*abs(this.ylen/2 - this.constFromTo(2))))...
                    .*(x <= this.ylen - this.constFromTo(1))...
                    +(x > this.ylen - this.constFromTo(1)).*ref(1);
                
            end
            
            
            %% TMRFunction berechnen
            S_m = models.muscle.functions.MarkertLawOriginal(mu(5),mu(6));
            AnisoPassiveMuscle = S_m.getFunction;
            S_t = general.functions.CubicToLinear(mu(7),mu(8));
            AnisoPassiveTendon = S_t.getFunction;
            F = AnisoPassiveTendon(this.stretchfun(0))*pi*this.radForTMR(0)^2;
            
            % Flaeche als Funktion des Ortes
            A = @(y)this.radForTMR(y).^2*pi;
            
            tmrFunc=@(y) (F./A(y)-AnisoPassiveMuscle(this.stretchfun(y)))./...
                (AnisoPassiveTendon(this.stretchfun(y))-AnisoPassiveMuscle(this.stretchfun(y)));
            
        end
       
        
        function alpha = getAlphaRamp(this, timeDots, activations)
            % Creates a linearly increasing scalar function starting at
            % starttime milliseconds ranging from zero to alphamax over
            % ramptime.
            %
            % Parameters:
            % ramptime: The time over which to increase to alphamax. If
            % less or equal to zero, an all zero function is returned.
            % alphamax: The maximum value to achieve. @type double @default
            % AMuscleConfig.ActivationRampMax
            % starttime: The offset time (in milliseconds) to wait before
            % increasing the signal. @type double 
            % @default AMuscleConfig.ActivationRampOffset
%             alpha = @(t)1;
            if nargin < 4
                activations = this.ActivationRampOffset;
                if nargin < 3
                    timeDots = this.ActivationRampMax;
                end
            end
            ramp = general.functions.PiecewiseLinear(timeDots, activations);
            alpha = ramp.getFunction;
        end        
        
        
    end
    methods(Access=protected)
        
        function geo = getGeometry(this)
            
            tmrFunc = getTMRFunc(this);
            [ypoints] = getDiscr(this,tmrFunc);
            geo = fem.geometry.Belly(ypoints,'Radius',this.radfun,'InnerRadius',0);
            
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix ends in xz direction
            displ_dir([1 3],geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            displ_dir([1 3],geo.Elements(geo.NumElements-3:geo.NumElements,...
                geo.MasterFaces(4,:))) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix ends in xz direction
            velo_dir(2,geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            velo_dir(2,geo.Elements(geo.NumElements-3:geo.NumElements,...
                geo.MasterFaces(4,:))) = true;
            % Die Funktion wurde aus Simulationsdaten generiert, passt nicht
            % zu 100% aber verbessert das Problem deutlich!
%             AmplitudenKorrektor =@(x) 317.6237./x + 78.3460./x.^2 - 2.3049.*10^3./x.^3 +1.2644.*10^4./x.^4;
            format long
            AmplitudenKorrektor =@(x) 305 * 1./x;
            velo_dir_val(velo_dir) = 1/(AmplitudenKorrektor(this.Fr));
        end
        
        function anull = seta0(~, anull)
            % Direction is y
            anull(2,:,:) = 1;
        end
                
        
        
    end
    
    methods(Static)
        
        function test_ShakerDefaultFuerSimTMR
            
            c = ShakerDefaultFuerSim('Stretch','log');
            m = c.createModel;
            c.plotTMR;
            m.plotGeometrySetup;
            
            
            c = ShakerDefaultFuerSim('Stretch','exp');
            m = c.createModel;
            c.plotTMR;
            m.plotGeometrySetup;
            
            c = ShakerDefaultFuerSim('Stretch','linear');
            m = c.createModel;
            c.plotTMR;
            m.plotGeometrySetup;
        end
        
        % Erzeugt die Config und das Model fuer das 10^-6_50Hz Video
        % in varargin wird jetzt die Toleranz fuer die Diskretisierung
        % gesetzt
        function [c,m] = createTestingConfig(frequency,amplitude,RampMax,RampOffset,varargin)
            if(isempty(varargin))
                
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',.25,...
                    'constFromTo',[10,40],'maxYLength',20,'Frequency',frequency,'Amplitude',amplitude);
                
            else
                
                c = ShakerDefaultFuerSim('Stretch','Gauss 0.3','TOL',varargin{1},...
                    'constFromTo',[10,40],'maxYLength',20,'Frequency',frequency,'Amplitude',amplitude);
            end
            
            c.ActivationRampMax = RampMax;
            c.ActivationRampOffset = RampOffset;
            m = c.createModel;
            m.T = 150;
            m.ODESolver.RelTol = .01;
            m.ODESolver.AbsTol = .1;
            
            % Ersetzt die TMR's durch die Mittelwerte der bisherigen
            % Wichtig, das Attribut wurde in System auf public gesetzt!!
            m.System.MuscleTendonRatioGP(1:size(m.System.MuscleTendonRatioGP,1),:)...
                = repmat(mean(m.System.MuscleTendonRatioGP),...
                size(m.System.MuscleTendonRatioGP,1),1);
            
            % Kleine Anpassung um an Stellen wirklich nur Tendon/Muscle zu
            % haben.
            m.System.MuscleTendonRatioGP = (m.System.MuscleTendonRatioGP - ...
                min(m.System.MuscleTendonRatioGP(:)))/(max(m.System.MuscleTendonRatioGP(:)...
                -min(m.System.MuscleTendonRatioGP(:))))./(10^(6));
            
            
            
        end
        
    end
end