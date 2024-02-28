%% Code Documentation
% Numerical Evaluation of Analytical Procedure for elastic-plastic
% waves in 1D longitudinal rods with variable cross-sectional
% area. Inspired by Hettche (1968) from paper entitled "Theoretical and
% experimental study on longitudinal impact of tapered rods".

% Constitutive behavior is based on a power-law representation such that
% stress = K*strain^n when stress is greater than the yield stress
% where n is the strain hardening exponent and K is the strength coefficient.
% The accessory function Wavespeed goes into further detail

% Procedure is laid out as a 6 step iterative loop for each spatial and
% time node. 
% Steps are as follows:
%   0) Assume the wave speed at the target point is equal to the elastic wavespeed 
%   1) Calculate wavespeed at relevant (nearby) spatial nodes at the previous time step
%      from constitutive relationship (function Wavespeed)
%   2) Approximate position (alpha) of inception points 1 and/or 2 (denoted as B and
%      D in paper) at previous time step such that a rightward characteristic
%      will propagate from 1 to the target point and a leftward characteristic
%      will propagate from 2 to the target point.
%   3) Interpolate state and wave speed at inception points (stress and velocity)
%      from adjacent mesh points
%   4) Calculate target point state with inception point states
%   5) Calculate target point wave speed based on the target stress found in 
%      (4) using the constitutive relationship (function Wavespeed) 
%   6) Check if wave speed from Step 5 matches assumed value. If they
%      match, then we have the correct state at the target point and we can
%      move on to the next node. If they do not match, then we adopt the 
%      wave speed from Step 6 and restart the loop. Match is accepted if it
%      is within user-defined Epsilon

% Nomenclature: 
% Target Point = Position (spatial node) and time (time increment) at which
%                we wish to find the state. Specified by indices i and j 
%                for spatial node and time increment respectively.
% Inception Point = Position at which the defining characteristic will
%                      propagate from the intermediate point to the target
%                      point. The two intermediate points uniquely define 
%                      the target point and are denoted by 1 and 2 
%                      (i.e. C1, C2, s1, s2, etc.) 
% alpha = Position of inception point (tied to subscripts 1 and 2)
% C = Wave speed
% S (or s) = Stress. Capital letter form denotes matrix while lower case
%            denotes single value (stress at a given point)
% V (or v) = Element velocity. Capital letter form denotes matrix while 
%            lower case denotes single value (velocity at a given point)
% Sy = Yield strength. (Sy0=Yield strength of virgin material)
% CP, CQ, CR = Wave speeds at mesh points (i-1,i,i+1) at previous time step
%              (j-1) where the target point is at the mesh point (j,i)
% g = Shorthand notation for (1/A)(d(ln(A))/dX) whereby g1 = g at inception
%     point 1, g2 = g at inception point 2, gT = g at target point


%Units herein are assumed in metric (meters, Pascals, etc)

%Questions on the code can be addressed to Greg Dorgant at
%gdorgant3@gatech.edu 

%Note: The wave speed assumed at the beginning of the iterative loop does
%not need to be the elastic wave speed (as suggested by Hettche). Similar
%results were found assuming the wave speed was equal to the wave speed at
%the previous time step. If one is able to converge on a correct solution
%over repeated iterations, then the initial guess of wave speed shouldn't
%matter


%% Code Executable
clear all %Boilerplate for me to clear any residual variables just in case

t2 = tic; %Timer for whole code

AmpIt = [5,350]; %Input amplitudes to sweep through
fIt = 0.5:0.5:25; %Input frequencies to sweep through
%Create iteration matrix containing all input variables to try
p = 1;
for j = 1:length(fIt)
    for i = 1:length(AmpIt)
        Iterate(p,:) = {-AmpIt(i)*1e6,fIt(j)*1e3};
        p = p + 1;
    end
end

[Iterations,~] = size(Iterate);
for ij = 1:Iterations %For each case we wish to simulate
    Amp = Iterate{ij,1};
    f_Input = Iterate{ij,2};

    %Parameters for Constitutive behavior
    E = 68.9e9; %Elastic Modulus, Pa
    rho = 2700; %Density, kg/m^3
    n = 0.2; %Strain Hardening exponent 
    C0 = sqrt(E/rho); % Elastic wave speed, m/s
    Sy0 = 276e6; %Initial Yield Stress
    K = (Sy0)^(1-n)*E^n; %Strength coefficient
    
    %Initialize space and time variables
    Luc = 0.25; %Unit Cell Length   
    rA = 0.5; %Ratio in area between the largest and smallest cross-sections: A_min/A_max (A3/A0)
    rL = 0.7; %Ratio in length between the conical sections and the unit cell length: LCone/Luc
    rC = 0.5; %Ratio in length between the first conical section and the total length for the conical sections: L2/LCone
    rR = 0.5; %Ratio in length between the first uniform section and the total length for the uniform sections: L1/LCone
    d0 = 0.167; %Diameter at start/Diameter for the first section/Length Scaling Parameter
    num_uc = 30; %Number of unit cells
    num_nodes_uc = 80; %Number of nodes per unit cell
    dX = Luc/num_nodes_uc; %Spatial discretization as specified from unit cell length and number of nodes
    dt = dX/C0; %Time discretization as follows from Courant-Friedrichs-Levy condition
    n1 = floor(rR*(1-rL)*num_nodes_uc); %Nodes in Section 1, larger cross-section uniform section
    n2 = floor(rL*rC*num_nodes_uc); %Nodes in Section 2, tapering cone
    n4 = floor(rL*(1-rC)*num_nodes_uc); %Nodes in Section 4, expanding cone
    n3 = num_nodes_uc - (n1 + n2 + n4); %Nodes in Section 3, smaller cross-section uniform section
    L = Luc*num_uc; %Length of total pulse shaper
    num_nodes = L/dX; %Number of nodes in system
    L1 = rR*(1-rL)*Luc; %Section 1 length
    L2 = rL*rC*Luc; %Section 2 length
    L3 = (1-rR)*(1-rL)*Luc; %Section 3 length
    L4 = rL*(1-rC)*Luc; %Section 4 length
    X = linspace(0,L,num_nodes); %Set spatial mesh, vector of size num_nodes
    T = 2*L/C0; %Simulation duration, seconds
    num_time_inc = floor(T/dt); %Number of time increments
    t = linspace(0,T,num_time_inc); %Create time vector as linearly spaced vector from 0 to simulation duration space by time step dt
    
    %Specify the variation in cross-sectional area over the length of the rod
    dL = sqrt(rA)*d0; %Smallest cross-sectional diameter
    A0 = (pi*d0^2)/4; %Cross-sectional area in section 1 (large)
    AL = (pi*dL^2)/4; %Cross-sectional area in section 3 (small)
    AreaProfile = 'Linear'; %Area profile in cone: Linear=Linearly tapering cone, Exponential=Exponential curve horn
    if strcmp(AreaProfile,'Linear')
        a = (AL-A0)/L2;
        b = A0;
        ACone = a*(linspace(0,L2,n2)) + b; %Cross-sectional area for the 1st cone section
    elseif strcmp(AreaProfile,'Exponential')
        Coef = (1/L2)*log(AL/A0);
        ACone = A0*exp(Coef*linspace(0,L2,n2)); %Cross-sectional area for the 1st cone section
        dCone = sqrt(4*ACone/pi);
    end

    A_UC = [A0*ones(1,n1/2),ACone,AL*ones(1,n3),fliplr(ACone),A0*ones(1,n1/2)]; %Combine all the cross-sectional area vectors to make the cross-sectional area vector over the unit cell
    A = repmat(A_UC,1,num_uc); % Total area profile over phononic crystal
    d = sqrt(4*A/pi); % Diameter of rod over space
    dA = diff(A)./diff(X); %dA/dx, the area slope at each node
    dA(length(dA)+1) = dA(end); %Assume area slope at end is the same as the node previous. Arises since diff() reduces vector length by 1
    dlnA = diff(log(A))./diff(X);
    dlnA(length(dlnA)+1) = dlnA(end);
    
    %Initialize state variable matrices
    S = zeros(num_time_inc,num_nodes); %Stress
    Sy = Sy0*ones(size(S)); %Yield stress
    V = S; %Element velocity
    Psi = S; %Intermediate stress variable
    R1 = S; %Rightward characteristics
    R2 = S; %Leftward characteristics
    C = C0*ones(num_time_inc,num_nodes); %Wave speed
    Cstar = C; %Wave speed guess
    
    %Specify wave form of input
    InputCurveType = 'Modulated Harmonic Wave';
    if strcmp(InputCurveType, 'Square Wave') 
        TDwell = 1e-3; %Duration of dwell
        s0 = [ones(1,floor(TDwell/dt)),zeros(1,num_time_inc-floor(TDwell/dt))]; %Square Load
    elseif strcmp(InputCurveType, 'Step Wave')
        s0 = ones(1,num_time_inc); %Step Load
    elseif strcmp(InputCurveType, 'Sine Wave')
        PhaseShift = 90; %In degrees
        s0 = sin(2*pi*f_Input*t + deg2rad(PhaseShift));
    elseif strcmp(InputCurveType,'Modulated Harmonic Wave')
        SignalPeriod = 5e-4;
        Window = (sin(pi*t/(2*SignalPeriod))).^2.*(t<2*SignalPeriod);
        s0 = [sin(2*pi*f_Input*t(1:floor(SignalPeriod/dt))).*Window(1:floor(SignalPeriod/dt)),sin(2*pi*f_Input*t(floor(SignalPeriod/dt)+1:end))];
    end

    %Scale input wave form to desired amplitude
    s0 = Amp*s0;
    S(:,1) = s0; %Incident Load such that node 1 always has stress specified by s0

    CW = Wavespeed(abs(S(1,1)),Sy0,C0,n,rho,K); %Wave speed at initial time and space node to calculate velocity at that point
    g2 = dlnA(2);
    V(1,1) = -S(1,1)*((3*C0-CW)/(C0^2) - dt*g2)/(2*rho); 

    %Values of characteristics at first time increment
    for i = 1:num_nodes
        Psi(1,i) = Sigma2Psi(S(1,i),Sy0,C0,n,K,rho);
        R1(1,i) = Psi(1,i) - rho*V(1,i) + (dt/2)*S(1,i)*dlnA(i);
        R2(1,i) = -Psi(1,i) - rho*V(1,i) + (dt/2)*S(1,i)*dlnA(i);
    end

    %Initialize conditional variables
    epsilon = 1e-3; %Close Enough condition value for wave speed
    np = 5; %Allowable iterations past the first for wave speed. Only need 1 for linearly hardening curves since there are only two possible wave speeds
    
    for j = 2:num_time_inc %Iterate over time nodes excluding the initial one since that state is pre-defined
        timer1 = tic; %Time iteration timer
        for i = 1:num_nodes %Iterate over all spatial nodes
            tic %Spatial iteration timer
            % CW = Cstar(j-1,i); %Initial guess for wave speed at target node
            CW = C0;
            Cond = false; %Initialize while loop condition
            p = 1; %Initialize wave speed iteration counter
            q=0;
            while Cond == false %While we have not confirmed wave speed at target point
                
                %Check if target's position begets special treament
                if i == num_nodes %If last node
                    %Step 1: Calculate wave speed at second to last point
                    %CS is as assumed since stress is 0 at the free end
                    CP = Wavespeed(abs(S(j-1,i-1)),Sy(j-1,i-1),C0,n,rho,K);
   
                    %Step 2: Calculate position of inception point
                    alpha1(p) = 2*C0*dX/(3*C0 - CP);
    
                    %Step 3: Calculate state at inception point   
                    s1 = S(j-1,i) + (S(j-1,i-1) - S(j-1,i))*alpha1(p)/dX;
                    v1 = V(j-1,i) + (V(j-1,i-1) - V(j-1,i))*alpha1(p)/dX;
                    g1 = dlnA(i) + (dlnA(i-1) - dlnA(i))*alpha1(p)/dX;
                    C1(p) = alpha1(p)/dt;
                    
                    %%% Step 5: Calculate state variables at target point
                    gT = dlnA(i);
                    S(j,i) = 0; %Free boundary condition at right end/last node
                    %Velocity at free end is calculated as such when one
                    %variable is specified and only nodes to the left exist
                    V(j,i) = (  (((3*C1(p)-CW(p))/C1(p)^2) + dt*gT)*S(j,i) + (-((3*C1(p)-CW(p))/C1(p)^2)*s1 + dt*gT) + 2*rho*v1  )/(2*rho);

                    %%%
                    Cstar(j,i) = C0;  
                    break %Exit while loop because we know stress at node definitively
                elseif i == 1 %If target point is first node or input face
                    
                    %Step 1: Calculate wave speed at nearby points
                    CW(p) = Wavespeed(abs(s0(j)),Sy(j,i),C0,n,rho,K); %Known definitively as stress is defined
                    CR = Wavespeed(abs(S(j-1,i+1)),Sy(j-1,i+1),C0,n,rho,K); 
                    CQ = Wavespeed(abs(S(j-1,i)),Sy(j-1,i),C0,n,rho,K); 
    
                    %Step 2: Calculate position of inception point
                    alpha2(p) = 2*CW(p)*dX/(2*C0 - CR + CW(p));                
                    
                    %Step 3: Calculate dependent variables at inception point                
                    s2 = S(j-1,i) + (S(j-1,i+1) - S(j-1,i))*alpha2(p)/dX;
                    v2 = V(j-1,i) + (V(j-1,i+1) - V(j-1,i))*alpha2(p)/dX;
                    g2 = dlnA(i) + (dlnA(i+1) - dlnA(i))*alpha2(p)/dX;
                    C2(p) = alpha2(p)/dt;
                   
                    %Step 5: Calculate state variables at target point
                    gT = dlnA(i);
                    S(j,i) = s0(j);    
                    V(j,i) = (  (-((3*C2(p)-CW(p))/C2(p)^2) + dt*gT)*S(j,i) + (((3*C2(p)-CW(p))/C2(p)^2) + dt*gT)*s2 + 2*rho*v2  )/(2*rho);

                    Cstar(j,i) = CW(p);

                   break %Exit while loop because we know the stress definitively at this point
                else %Target point is not on edges of rod
    
                    %Step 1: Calculate wave speed at nearby points
                    CQ = Wavespeed(abs(S(j-1,i)),Sy(j-1,i),C0,n,rho,K);
                    CP = Wavespeed(abs(S(j-1,i-1)),Sy(j-1,i-1),C0,n,rho,K);
                    CR = Wavespeed(abs(S(j-1,i+1)),Sy(j-1,i+1),C0,n,rho,K);
    
                    %Step 2: Calculate position of inception points
                    alpha1(p) = (CW(p) + CQ)*dX/(2*C0 - CP + CQ);
                    alpha2(p) = (CW(p) + CQ)*dX/(2*C0 - CR + CQ);
                    
                    %Step 3: Calculate state at inception points
                    s1 = S(j-1,i) + (S(j-1,i-1) - S(j-1,i))*alpha1(p)/dX;
                    v1 = V(j-1,i) + (V(j-1,i-1) - V(j-1,i))*alpha1(p)/dX;
                    g1 = dlnA(i) + (dlnA(i-1) - dlnA(i))*alpha1(p)/dX;
                    s2 = S(j-1,i) + (S(j-1,i+1) - S(j-1,i))*alpha2(p)/dX;
                    v2 = V(j-1,i) + (V(j-1,i+1) - V(j-1,i))*alpha2(p)/dX;
                    g2 = dlnA(i) + (dlnA(i+1) - dlnA(i))*alpha2(p)/dX;
                    C2(p) = alpha2(p)/dt;
                    C1(p) = alpha1(p)/dt;
                    
                    %Step 5: Calculate state variables at target point
                    gT = dlnA(i);
                    S(j,i) = (  (((3*C2(p)-CW(p))/C2(p)^2) + dt*gT)*s2 + 2*rho*v2 - (-((3*C1(p)-CW(p))/C1(p)^2) + dt*gT)*s1-2*rho*v1  )/(((3*C2(p)-CW(p))/C2(p)^2) + ((3*C1(p)-CW(p))/C1(p)^2) + dt*(g1-g2));
                    V(j,i) = (  (((3*C1(p)-CW(p))/C1(p)^2) + dt*gT)*S(j,i) + (-((3*C1(p)-CW(p))/C1(p)^2) + dt*gT)*s1 + 2*rho*v1  )/(2*rho);

                end
            
                %Step 6: Calculate wave speed based on constitutive behavior
                Cstar(j,i) = Wavespeed(abs(S(j,i)),Sy(j,i),C0,n,rho,K);
    
                pCount(j,i) = p;
                C(j,i) = CW(p);
                %%% Step 7: Is our solution close enough?
                if abs((Cstar(j,i)-CW(p))/Cstar(j,i)) < epsilon %If the wave speed calculated from the constitutive relation
                %using the calculated stress is close enough to our guessed wave speed
                    Cond = true; %Escape loop, we've found the correct stress at target point
                    break
                else
                    CW(p+1) = Cstar(j,i); %Update guess on wave speed at target point
                    
                end
                %%%

                %%% Update iteration counter. Doesn't matter if Step 7 was successful and escape condition is now true
                if p == np+1 %If we exceed number of free iterations, revert to specified wave speed
%                     CW(p+1) = min(CW);
%                     CW(p+1) = C(j-1,i);
                    CW(p+1) = C0;
                    p = p + 1;
                elseif p == np+2  %After allowing procedure to run once more, escape loop
                    Cond = true;
                else
                    p = p + 1; %Increment counter
                end
            end

            %Post process stress and yield stress fields to find field of
            %intermediate stress variables
            Psi(j,i) = Sigma2Psi(S(j,i),Sy(j-1,i),C0,n,K,rho);
            Phi = dt*S(j,i)*dlnA(i);
            %Find rightward and leftward characteristics
            R1(j,i) = Psi(j,i) - rho*V(j,i) + Phi;
            R2(j,i) = -Psi(j,i) - rho*V(j,i) + Phi;

            %Update yield stress at target point if plastic yielding occurred
            if abs(S(j,i)) > Sy(j,i)
                Sy(j+1:end,i) = abs(S(j,i));
            end
            tNode(j,i) = toc; %Record time spent at this target point
        end

        tTime(ij,j) = toc(timer1); %Record time spent during this time increment

    end

    SA(ij).R1 = R1;
    SA(ij).R2 = R2;
    SA(ij).Stress = S;
    SA(ij).Psi = Psi;
    SA(ij).Time = t;
    SA(ij).Position = X;

    if n > 1
        Description_Constitutive = 'Increasingly Hardening';
    elseif n < 1
        Description_Constitutive = 'Decreasingly Hardening';
    else
        Description_Constitutive = 'Linearly Hardening';
    end
    if strcmp(InputCurveType,'Square Wave') || strcmp(InputCurveType,'Step Wave')
        SA(ij).Description = sprintf(strjoin({'%d MPa with' InputCurveType,'Input',Description_Constitutive},' '),[-Iterate{ij,1}/1e6]);
    else
        SA(ij).Description = sprintf(strjoin({'%0.0f MPa at %0.1f kHz with ' InputCurveType,' Input PnC Explicit Mat1'},''),[-Iterate{ij,1}/1e6,Iterate{ij,2}/1e3]);
    end
        
    [Strain{ij},Residual_Strain{ij},YieldStrength{ij},Temp] = StrainFcn(S,E,n,Sy0,K);
%     [Strain{ij},Residual_Strain{ij}] = StrainFcn(Abaqus(ij).Stress,E,n,Sy0,K);

    SA(ij).ResidualStrain = Residual_Strain{ij};
    SA(ij).Strain = Strain{ij};
    SA(ij).YieldStrength = YieldStrength{ij};

end

%% Load data from Abaqus
%Data in Abaqus stored under same file nomenclature as SA description field
%as text file
TypeOpts = delimitedTextImportOptions('Delimiter',{' ',':','E:','N:'},...
    'DataLines',[1,2],'ConsecutiveDelimitersRule','join',...
    'LeadingDelimitersRule','ignore');
File_Location = 'C:\Users\gdorgant3\Documents\';
for i = 1:length(SA)
    StrainFile = strjoin({File_Location,'Residual Strain ',SA(i).Description,'.txt'},'');
    File = strjoin({File_Location,SA(i).Description,'.txt'},'');

    VarTypes = rows2vars(readtable(File,TypeOpts));
    VarTypes = VarTypes(:,2);
    temp = VarTypes(2:2:end,1);
    IDs = table2cell(unique(temp,'stable'));
    num_Measurement_Locations = height(temp)/length(IDs);
    temp = table2array(varfun(@(x)str2double(x),VarTypes(3:2:end,1)));
    DataOpts = delimitedTextImportOptions('DataLines',3,'Delimiter',{' '},...
        'VariableNames',repmat({'Var'},1,length(IDs)*num_Measurement_Locations+1),...
        'VariableTypes',repmat({'double'},1,length(IDs)*num_Measurement_Locations+1),...
        'ConsecutiveDelimitersRule','join','LeadingDelimitersRule','ignore','ExtraColumnsRule','ignore');
    Data = table2array(readtable(File,DataOpts));
    Abaqus(i).Time = Data(:,1);
    FE_dt = Abaqus(i).Time(2)-Abaqus(i).Time(1);
    Data(:,1) = [];

    for j = 1:length(IDs)
        if strcmp(IDs{j},'S')
            Name = 'Stress';
        elseif strcmp(IDs{j},'PE')
            Name = 'ResidualStrain';
        elseif strcmp(IDs{j},'A')
            Name = 'Acceleration';
        elseif strcmp(IDs{j},'V')
            Name = 'Velocity';
        elseif strcmp(IDs{j},'U')
            Name = 'Displacement';
        elseif strcmp(IDs{j},'KE')
            Name = 'KineticEnergy';
        elseif strcmp(IDs{j},'PD')
            Name = 'PlasticDissipationEnergy';
        elseif strcmp(IDs{j},'SE')
            Name = 'ElasticStrainEnergy';
        elseif strcmp(IDs{j},'LE')
            Name = 'Strain';
        end
        Abaqus(i).(Name) = Data(:,1+(j-1)*num_Measurement_Locations:(j)*num_Measurement_Locations);
    end
    Abaqus(i).Description = SA(i).Description;

end

%% Plastic Strain Surface (Fig. 22)
clear Q
p=0;
I_Measure=ceil((L/C0)/dt); %Take measurement at time required for elastic wave to travel one length of system
for i = 1:length(SA)
    if contains(SA(i).Description,'350 MPa')
        Q(p,:) = SA(i).ResidualStrain(:,I_Measure);
        p = p + 1;
    end
end

ds_X = 10;
figure
surface(X(1:ds_X:end)/Luc,fIt,Q(:,1:ds_X:end))
shading interp
xlabel('Unit Cell Number')
ylabel('Input Frequency (kHz)')
zlabel('Plastic Strain')
ax = gca;
colo = colorbar('northoutside');
clim([-2,2])
zlim([-10,10])
cmap = colormap;
cmap(129,:) = [0.75,0,0];
tmp1 = 2;
tmp2 = 0.1;
yline(8.36,'--','Color',[tmp2,tmp2,tmp2],'LineWidth',tmp1)
yline(12,'--','Color',[tmp2,tmp2,tmp2],'LineWidth',tmp1)
colormap(cmap);
sizeFigureFonts('Times New Roman',14)



%% PnC Validation Surface (Fig. 21)
clear Q 

%Take stress during runs where input was elastic 
p=0;
I_Measure=(10-1)*num_nodes_uc+1; %Take measurement at 10th unit cell
for i = 1:length(SA)
    if contains(SA(i).Description,'5 MPa')
        Q(p,:) = SA(i).Stress(:,I_Measure);
        p = p + 1;
    end
end


ds_t = 10;
figure
surface(t(1:ds_t:end)*1e3,fIt,Q(:,1:ds_t:end)/1e6)
shading interp
xlabel('Time (ms)')
ylabel('Input Frequency (kHz)')
zlabel('Stress (MPa)')
ax = gca;
colo = colorbar('northoutside');
colo.Label.String = 'Stress (MPa)';
colo.Label.FontSize = 14;
colo.Label.FontName = 'Times New Roman';
tmp1 = 2;
tmp2 = 0.1;
yline(8.36,'--','Color',[tmp2,tmp2,tmp2],'LineWidth',tmp1)
yline(12,'--','Color',[tmp2,tmp2,tmp2],'LineWidth',tmp1)
clim([-10,10])
cmap = colormap;
sizeFigureFonts('Times New Roman',14)


%% Stress Surface at 3 frequencies (Fig. 23)

ds_X = 10;
ds_t = 10;
for i = 1:length(SA)
    if contains(SA(i).Description,'350 MPa')
        if contains(SA(i).Description,'5 kHz') || contains(SA(i).Description,'10 kHz') || contains(SA(i).Description,'15 kHz')
            figure
            surface(X(1:ds_X:end)/Luc,t(1:ds_t:end)*1e3,SA(i).Stress(1:ds_t:end,1:ds_X:end)/1e6)
            shading interp
            ylabel('Time (ms)')
            xlabel('Unit Cell Number')
            zlabel('Stress (MPa)')
            ylim([0,3])
            sizeFigureFonts('Times New Roman',16)
            zlim([Sy0,1e9]/1e6)
        end
    end
end

%% Characteristic Surfaces Fig.24
ds_X = 10;
ds_t = 10;

% Rightward characteristics
figure
surface(X(1:ds_X:end)/Luc,t(1:ds_t:end)*1e3,1e-5*SA(1).R1(1:ds_t:end,1:ds_X:end))
shading interp
ylabel('Time (ms)')
xlabel('Unit Cell Number')
c1 = colorbar('northoutside');
sizeFigureFonts('Times New Roman',16)
clim([-2,2])
g1 = gca;
g1.Position = [0.125,0.15,0.775,0.65];
c1.Position = [0.125,0.825,0.775,0.05];
c1.Ticks = [-2,-1,0,1,2];
sizeFigureFonts('Times New Roman',16)

%Leftward characteristics
figure
surface(X(1:ds_X:end)/Luc,t(1:ds_t:end)*1e3,1e-5*(SA(1).R2(1:ds_t:end,1:ds_X:end)))
shading interp
ylabel('Time (ms)')
xlabel('Unit Cell Number')
c2 = colorbar('northoutside');
clim([-1,1])
g2 = gca;
g2.Position = [0.125,0.15,0.775,0.65];
c2.Position = [0.125,0.825,0.775,0.05];
c2.Ticks = [-1,-0.5,0,0.5,1];
sizeFigureFonts('Times New Roman',16)



%% Break for end
tTotal = toc(t2); %Record total time spent during procedure

