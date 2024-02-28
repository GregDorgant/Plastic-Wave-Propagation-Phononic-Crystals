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
InputCurveType = {'Step Wave','Sine Wave'}; % Wave profiles to sweep through
AmpIt = [5,350]; %Amplitudes to sweep through
fIt = [5]; %Input frequencies to sweep through
p = 1;

% Create iterate through which the code will loop
for i = 1:length(AmpIt)
    for j = 1:length(InputCurveType)
        if contains(InputCurveType{j},'ine') || contains(InputCurveType{j},'Harmonic')
            for k = 1:length(fIt)
    
                Iterate(p,:) = {-AmpIt(i)*1e6,fIt(k)*1e3,InputCurveType{j}};
                p = p + 1;
            end
        else
                Iterate(p,:) = {-AmpIt(i)*1e6,0,InputCurveType{j}};
                p = p + 1;
        end
    end
end

[Iterations,~] = size(Iterate);
for ij = 1:Iterations %for each specified configuration
    timer2 = tic;
    Amp = Iterate{ij,1}; %Adopt specified amplitude
    f_Input = Iterate{ij,2}; %Adopt specified frequency input

    %Parameters for Constitutive behavior
    E = 68.9e9; %Elastic Modulus, Pa
    rho = 2700; %Density, kg/m^3
    n = 0.2; %Strain Hardening exponent 
    Ce = sqrt(E/rho); % Elastic wave speed, m/sLCone
    Sy0 = 276e6; %Initial Yield Stress
    K = (Sy0)^(1-n)*E^n; %Strength Coefficient
    
    %Initialize space and time variables
    L = 0.25*30; %Unit Cell Length   
    rA = 0.5; %Ratio in area between the largest and smallest cross-sections: A_min/A_max (A3/A0)
    rL = 0.7; %Ratio in length between the conical sections and the unit cell length: LCone/Luc
    rC = 0.5; %Ratio in length between the first conical section and the total length for the conical sections: L2/LCone
    rR = 0.5; %Ratio in length between the first uniform section and the total length for the uniform sections: L1/LCone
    d0 = 0.167; %Diameter at start/Diameter for the first section/Length Scaling Parameter
    num_uc = 30; %Number of unit cells
    dt = 5e-7; %Set time increment to approximately match FE simulations, seconds
    dX = Ce*dt; %Set spatial mesh distance by time increment and elastic wave speed based on Courant Condition
    num_nodes = floor(L/dX); %Find number of nodes as integer from total length and spatial discretization
    dX = floor(L/num_nodes); %Ensure spatial discretization leads to integer number of nodes
    dt = dX/Ce; %Refinforce CFL condition
    n1 = floor((1-rL)*num_nodes); %Number of nodes in uniform section
    n2 = num_nodes-n1; %Number of nodes in cone section
    L1 = (1-rL)*L; %Length of uniform section
    L2 = rL*L; %Length of cone section
    X = linspace(0,L,num_nodes); %Set spatial mesh, vector of size num_nodes
    T = 2*L/Ce; %Simulation duration, seconds
    num_time_inc = floor(T/dt); %Number of time increments
    t = linspace(0,T,num_time_inc); %Create time vector as linearly spaced vector from 0 to simulation duration space by time step dt
    
    %Specify the variation in cross-sectional area over the length of the rod
    dL = sqrt(rA)*d0; %Smallest cross-sectional diameter
    A0 = (pi*d0^2)/4; %Larger cross-sectional area (uniform section)
    AL = (pi*dL^2)/4; %Smaller cross-sectional area (cone's final cross-section)
    AreaProfile = 'Linear'; %Cone's taper profile (linear=normal cone, exponential=horn)
    if strcmp(AreaProfile,'Linear')
        ACone = (AL-A0)/L2*(linspace(0,L2,n2)) + A0; %Cross-sectional area for the 1st cone section
    elseif strcmp(AreaProfile,'Exponential')
        Coef = (1/L2)*log(AL/A0);
        ACone = A0*exp(Coef*linspace(0,L2,n2));
        dCone = sqrt(4*ACone/pi);
    end
    A = [A0*ones(1,n1),ACone]; %Cross-sectional area vector for unit cell is only comprised of the tapering section and larger cross-section uniform rod
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
    C = Ce*ones(num_time_inc,num_nodes); %Wave speed
    Cstar = C; %Wave speed guesses
    
    %Specify wave form of input
    InputType = Iterate{ij,3};
    if strcmp(InputType, 'Step Wave')
        s0 = ones(1,num_time_inc); %Step Load
    elseif strcmp(InputType, 'Sine Wave')
        PhaseShift = 0; %In degrees
        s0 = sin(2*pi*f_Input*t + deg2rad(PhaseShift));
    elseif strcmp(InputType,'Modulated Harmonic Wave')
        SignalPeriod = 5e-4;
        Window = (sin(pi*t/(2*SignalPeriod))).^2.*(t<2*SignalPeriod);
        s0 = [sin(2*pi*f_Input*t(1:floor(SignalPeriod/dt))).*Window(1:floor(SignalPeriod/dt)),sin(2*pi*f_Input*t(floor(SignalPeriod/dt)+1:end))];
    end
    %Scale input wave form to desired amplitude
    s0 = Amp*s0;
    S(:,1) = s0; %Incident Load such that node 1 always has stress s0

    CT = Wavespeed(abs(S(1,1)),Sy0,Ce,n,rho,K); %Wave speed at initial time and space node to calculate velocity at that point
    g2 = dlnA(2);
    V(1,1) = -S(1,1)*((3*Ce-CT)/(Ce^2) - dt*g2)/(2*rho);

    %Values of characteristics at first time increment
    for i = 1:num_nodes
        Psi(1,i) = Sigma2Psi(S(1,i),Sy0,Ce,n,K,rho);
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
%             CW = Cstar(j-1,i); %Initial guess for wave speed at target node
            CT = Ce;  %Initial guess for wave speed at target node
            Cond = false; %Initialize while loop condition
            p = 1; %Initialize wave speed iteration counter
            while Cond == false %While we have not confirmed wave speed at target point
                
                %Check if target's position begets special treament
                if i == num_nodes %If last node
                    %Step 1: Calculate wave speed at second to last point
                    %CS is as assumed since stress is 0 at the free end
                    CP = Wavespeed(abs(S(j-1,i-1)),Sy(j-1,i-1),Ce,n,rho,K);
   
                    %Step 2: Calculate position of inception point
                    alpha1(p) = 2*Ce*dX/(3*Ce - CP);
    
                    %Step 3: Calculate state at inception point   
                    s1 = S(j-1,i) + (S(j-1,i-1) - S(j-1,i))*alpha1(p)/dX;
                    v1 = V(j-1,i) + (V(j-1,i-1) - V(j-1,i))*alpha1(p)/dX;
                    g1 = dlnA(i) + (dlnA(i-1) - dlnA(i))*alpha1(p)/dX;
                    gT = dlnA(i);

                    C1(p) = alpha1(p)/dt; %Calculate wave speed from inception point position

                    %%% Step 5: Calculate state variables at target point
                    S(j,i) = 0; %Free boundary condition at right end/last node
                    %Velocity at free end is calculated as such when one
                    %variable is specified and only nodes to the left exist                   
                    V(j,i) = (  (((3*C1(p)-CT(p))/C1(p)^2) + dt*gT)*S(j,i) + (-((3*C1(p)-CT(p))/C1(p)^2)*s1 + dt*gT) + 2*rho*v1  )/(2*rho);
                    %%%
                    Cstar(j,i) = Ce;    
                    break %Exit while loop because we know stress at node definitively
                elseif i == 1 %If target point is first node or input face
                    
                    %Step 1: Calculate wave speed at nearby points
                    CT(p) = Wavespeed(abs(s0(j)),Sy(j,i),Ce,n,rho,K); %Known definitively as stress is defined
                    CR = Wavespeed(abs(S(j-1,i+1)),Sy(j-1,i+1),Ce,n,rho,K); 
                    CQ = Wavespeed(abs(S(j-1,i)),Sy(j-1,i),Ce,n,rho,K); 
    
                    %Step 2: Calculate position of inception point
                    alpha2(p) = (CT(p)+CQ)*dX/(2*Ce - CR + CT(p));                
                    
                    %Step 3: Calculate dependent variables at inception point                
                    s2 = S(j-1,i) + (S(j-1,i+1) - S(j-1,i))*alpha2(p)/dX;
                    v2 = V(j-1,i) + (V(j-1,i+1) - V(j-1,i))*alpha2(p)/dX;
                    g2 = dlnA(i) + (dlnA(i+1) - dlnA(i))*alpha2(p)/dX;
                    gT = dlnA(i);
            
                    C2(p) = alpha2(p)/dt; %Calculate wave speed from inception point position
                    %Step 5: Calculate state variables at target point
                    S(j,i) = s0(j);  
                    V(j,i) = (  (-((3*C2(p)-CT(p))/C2(p)^2) + dt*gT)*S(j,i) + (((3*C2(p)-CT(p))/C2(p)^2) + dt*gT)*s2 + 2*rho*v2  )/(2*rho);

                    Cstar(j,i) = Wavespeed(abs(s0(j)),Sy(j,i),Ce,n,rho,K);  
                   break %Exit while loop because we know the stress definitively at this point
                else %Target point is not on edges of rod
    
                    %Step 1: Calculate wave speed at nearby points
                    CQ = Wavespeed(abs(S(j-1,i)),Sy(j-1,i),Ce,n,rho,K);
                    CP = Wavespeed(abs(S(j-1,i-1)),Sy(j-1,i-1),Ce,n,rho,K);
                    CR = Wavespeed(abs(S(j-1,i+1)),Sy(j-1,i+1),Ce,n,rho,K);

                    %Step 2: Calculate position of inception points
                    alpha1(p) = (CT(p) + CQ)*dX/(2*Ce - CP + CQ);
                    alpha2(p) = (CT(p) + CQ)*dX/(2*Ce - CR + CQ);
       
                    %Step 3: Calculate state at inception points
                    s1 = S(j-1,i) + (S(j-1,i-1) - S(j-1,i))*alpha1(p)/dX;
                    v1 = V(j-1,i) + (V(j-1,i-1) - V(j-1,i))*alpha1(p)/dX;
                    g1 = dlnA(i) + (dlnA(i-1) - dlnA(i))*alpha1(p)/dX;
                    s2 = S(j-1,i) + (S(j-1,i+1) - S(j-1,i))*alpha2(p)/dX;
                    v2 = V(j-1,i) + (V(j-1,i+1) - V(j-1,i))*alpha2(p)/dX;
                    g2 = dlnA(i) + (dlnA(i+1) - dlnA(i))*alpha2(p)/dX;

                    %Calculate wave speed from inception point position
                    C2(p) = alpha2(p)/dt;
                    C1(p) = alpha1(p)/dt;

                    gT = dlnA(i);
                    %Step 5: Calculate state variables at target point
                    S(j,i) = (  (((3*C2(p)-CT(p))/C2(p)^2) + dt*gT)*s2 + 2*rho*v2 - (-((3*C1(p)-CT(p))/C1(p)^2) + dt*gT)*s1-2*rho*v1  )/(((3*C2(p)-CT(p))/C2(p)^2) + ((3*C1(p)-CT(p))/C1(p)^2) + dt*(g1-g2));
                    V(j,i) = (  (((3*C1(p)-CT(p))/C1(p)^2) + dt*gT)*S(j,i) + (-((3*C1(p)-CT(p))/C1(p)^2) + dt*gT)*s1 + 2*rho*v1  )/(2*rho);

                end
            
                %Step 6: Calculate wave speed based on constitutive behavior
                Cstar(j,i) = Wavespeed(abs(S(j,i)),Sy(j,i),Ce,n,rho,K);
                
                %%% Step 7: Is our solution close enough?
                if abs((Cstar(j,i)-CT(p))/Cstar(j,i)) < epsilon %If the wave speed calculated from the constitutive relation
                %using the calculated stress is close enough to our guessed wave speed
                    Cond = true; %Escape loop, we've found the correct stress at target point
                else
                    CT(p+1) = Cstar(j,i); %Update guess on wave speed at target point
                end
                %%%

                pCount(j,i) = p;
                C(j,i) = CT(p);
                %%% Update iteration counter. Doesn't matter if Step 7 was successful and escape condition is now true
                if p == np+1 %If we exceed number of free iterations
                    CT(p+1) = C(j-1,i); %If we cannot converge on solution, assume wave speed does not change from previous time step
                    % CW(p+1) = C0; %If we cannot converge on solution, assume elastic wave speed
                    p = p + 1; %Allow one more iteration to calculate state with final assumed wave speed
                elseif p == np+2 %Escape loop
                    Cond = true;
                else
                    p = p + 1; %Increment counter
                end
            end
            %Post process stress and yield stress fields to find field of
            %intermediate stress variables
            Psi(j,i) = Sigma2Psi(S(j,i),Sy(j-1,i),Ce,n,K,rho);
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
    if strcmp(InputType,'Square Wave') || strcmp(InputType,'Step Wave')
        SA(ij).Description = sprintf(strjoin({'%d MPa with' InputType,'Input Rod Cone Explicit Mat1'},' '),[-Iterate{ij,1}/1e6]);
    else
        SA(ij).Description = sprintf(strjoin({'%0.0f MPa at %0.1f kHz with ' InputType,' Input Rod Cone Explicit Mat1'},''),[-Iterate{ij,1}/1e6,Iterate{ij,2}/1e3]);
    end
        
    [Strain{ij},Residual_Strain{ij}] = StrainFcn(S,E,n,Sy0,K);
%     [Strain{ij},Residual_Strain{ij}] = StrainFcn(Abaqus(ij).Stress,E,n,Sy0,K);

    SA(ij).ResidualStrain = Residual_Strain{ij};
    SA(ij).Strain = Strain{ij};

    ItTimer(ij) = toc(timer2); %Record time spent on this iteration case
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
    ResStrainDataAtEnd = readtable(StrainFile);
    ResStrainDataAtEnd = table2array(ResStrainDataAtEnd(:,2));
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
%% SA vs FE Validation in Residual Strain (Fig. 17b, 18c)

for i = 1:length(SA)
    if contains(SA(i).Description,'350 MPa')
        dt_fe = Abaqus(i).Time(2)-Abaqus(i).Time(1);
        [~,tmp] = min(abs(Abaqus(i).Time-3e-3));
        figure
        plot(SA(i).Time(1:ceil(3e-3/dt))*1e3,SA(i).ResidualStrain(1:ceil(3e-3/dt),ceil(3/dX))/(Sy0/E),'Color',[1,0,0], 'LineWidth',1.5)
        hold on
        plot(Abaqus(i).Time(1:tmp)*1e3,Abaqus(i).ResidualStrain(1:tmp,25)/(Sy0/E),'Color',[0,0,1], 'LineWidth',1.5,'LineStyle','--')
        hold off
        xlabel('Time (ms)')
        xticks(0:3)
        xlim([0,3])
        ylabel('$\varepsilon^{p}/\varepsilon_{Y}$','Interpreter','latex')
        sizeFigureFonts('Times New Roman',18)
    end
end

%% Validating Residual Strain Calculation (Fig. 18d)

for i = 1:length(SA)
    if contains(SA(i).Description,'350 MPa') && contains(SA(i).Description,'Sine Wave')
        j = i;
    end
end
[Strain,Residual_Strain] = StrainFcn(Abaqus(j).Stress,E,n,Sy0,K);
figure
plot(Abaqus(j).Time*1e3,Residual_Strain(:,25)/(Sy0/E),'Color',[1,0,0], 'LineWidth',1.5)
hold on
plot(Abaqus(j).Time*1e3,Abaqus(j).ResidualStrain(:,25)/(Sy0/E),'Color',[0,0,1], 'LineWidth',1.5,'LineStyle','--')
hold off
xlabel('Time (ms)')
xticks(0:3)
xlim([0,3])
% ylim([1.1*min(Abaqus(j).ResidualStrain(:,25)/(Sy0/E)),1.1*max(Abaqus(j).ResidualStrain(:,25)/(Sy0/E))])
yticks([-0.02,-0.01,0,0.01])
sizeFigureFonts('Times New Roman',18)
ylabel('$\varepsilon^{p}/\varepsilon_{Y}$','Interpreter','latex')

%% SA vs FE Stress Validation (Figs. 16, 17a, 18a, 18b)

for i = 1:length(SA)
    dt_fe = Abaqus(i).Time(2)-Abaqus(i).Time(1);
    [~,tmp] = min(abs(Abaqus(i).Time-3e-3));
    figure
    plot(SA(i).Time(1:ceil(3e-3/dt))*1e3,SA(i).Stress(1:ceil(3e-3/dt),ceil(3/dX))/1e6,'Color',[1,0,0], 'LineWidth',1.5)
    hold on
    plot(Abaqus(i).Time(1:tmp)*1e3,Abaqus(i).Stress(1:tmp,25)/1e6,'Color',[0,0,1], 'LineWidth',1.5,'LineStyle','--')
    hold off
    xlabel('Time (ms)')
    ylabel('Stress (MPa)') 
    xticks(0:3)
    xlim([0,3])
    if contains(SA(i).Description,'350 MPa') && contains(SA(i).Description,'Step')
        ylim([-400,0])
        sizeFigureFonts('Times New Roman',18)
    elseif contains(SA(i).Description,'350 MPa') && contains(SA(i).Description,'Sine')
        ylim([-400,400])
        sizeFigureFonts('Times New Roman',18)

        figure
        plot(SA(i).Time(1:ceil(3e-3/dt))*1e3,SA(i).Stress(1:ceil(3e-3/dt),ceil(3/dX))/1e6,'Color',[1,0,0], 'LineWidth',1.5)
        hold on
        plot(Abaqus(i).Time(1:tmp)*1e3,Abaqus(i).Stress(1:tmp,25)/1e6,'Color',[0,0,1], 'LineWidth',1.5,'LineStyle','--')
        hold off
        xlabel('Time (ms)')
        ylabel('Stress (MPa)') 
        xlim([2.7,2.8])
        xticks([2.7,2.75,2.8])
        ylim([280,290])
        yticks([280,285,290])
        sizeFigureFonts('Times New Roman',18)
    elseif contains(SA(i).Description,'5 MPa') && contains(SA(i).Description,'Sine')
        ylim([-6,6])
        yticks(-5:2.5:5)
        sizeFigureFonts('Times New Roman',18)
    elseif contains(SA(i).Description,'5 MPa') && contains(SA(i).Description,'Step')
        ylim([-6,1.5])
        yticks(-5:2.5:0)
        sizeFigureFonts('Times New Roman',18)
    end

end

%% Break for end
tTotal = toc(t2); %Record total time spent during procedure

