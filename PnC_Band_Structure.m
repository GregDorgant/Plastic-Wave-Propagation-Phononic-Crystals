%% Documentation
% 1D adaptation of finite element procedure presented by Phani et al. (2006)
% to derive the band structure for a 2D lattice. 
% The real part of the wavenumber is found by solving the eigenvalue
% problem set up by Phani et al. in which we assume a value for the wave
% number and solve for the various frequencies corresponding to the
% wavenumber. The number of frequencies for each wavenumber should
% correspond to the number of elements used
% The imaginary part of the wavenumber is found by assuming a frequency and
% a wavenumber of real part=pi and an assumed imaginary part. The correct
% imaginary part of the wavenumber is that which minimizes the
% determinant of the problem for a given frequency. This was found to be
% very sensitive so be careful in recreating it.


%% Executable
clear all
Luc = 0.25; %Unit Cell Length   
rA = 0.5; %Ratio in area between the largest and smallest cross-sections: A_min/A_max (A3/A0)
rL = 0.7; %Ratio in length between the conical sections and the unit cell length: LCone/Luc
rC = 0.5; %Ratio in length between the first conical section and the total length for the conical sections: L2/LCone
rR = 0.5; %Ratio in length between the first uniform section and the total length for the uniform sections: L1/LCone
d0 = 0.167; %Diameter at start/Diameter for the first section/Length Scaling Parameter
E = 68.9e9; %Elastic modulus
rho = 2700; %Mass density
C0 = sqrt(E/rho); %Elastic wave speed
L = Luc; %Set system length as unit cell length
n = 80; %Number of nodes in system
dX = L/n; %Set spatial discretization using system length and nodes
n1 = floor(rR*(1-rL)*n); %Nodes in Section 1, larger cross-section rod
n2 = floor(rL*rC*n); %Nodes in Section 2, tapering cone
n3 = floor((1-rR)*(1-rL)*n); %Nodes in Section 3, smaller cross-section rod
n4 = floor(rL*(1-rC)*n); %Nodes in Section 4, expanding cone
L1 = rR*(1-rL)*Luc; %Length of Section 1
L2 = rL*rC*Luc; %Length of Section 2
L3 = (1-rR)*(1-rL)*Luc; %Length of Section 3
L4 = rL*(1-rC)*Luc; %Length of Section 4

dL = sqrt(rA)*d0;
A0 = (pi*d0^2)/4;
AL = (pi*dL^2)/4;
ACone = ((AL-A0)/L2)*linspace(0,L2,n2) + A0; %Cross-sectional area for the tapering section of the rod
A1_Short = A0*ones(1,n1/2); %Larger cross-sectional uniform rod
A1 = A0*ones(1,n1); %Larger cross-sectional uniform rod
A3 = AL*ones(1,n3); %Smaller cross-sectional uniform rod
A = [A1_Short,ACone,A3,fliplr(ACone),A1_Short]; %Combine all the cross-sectional area vectors to make the cross-sectional area vector over the unit cell

K = E/(dX)*(diag(-A,1) + diag(-A,-1) + diag([0,A]) + diag([A,0])); %Form stiffness matrix of system
M = (rho*dX/4)*(diag(A,1) + diag(A,-1) + diag([0,A]) + diag([A,0])); %Form mass matrix of system
k = []; %Initialize wave number vector
omega = []; %Initialize frequency vector 
k1 = 0:pi/100:pi; %Trial wave numbers
T = eye(n+1); %Floquet Bloch matrix
T(:,end) = []; %Empty out the last column for Floquet Bloch condition

for i = 1:length(k1) %For each of the trial wave numbers, solve eigenvalue problem
    T(end,1) = exp(1i*k1(i)); 
    AA = T'*K*T;
    BB = T'*M*T;
    [~,temp] = eig(AA,BB);
    temp = sqrt(diag(temp));
    k = [k,k1(i)*ones(1,length(temp))];
    omega = [omega,temp'];
end

[omega1,I] = sort(omega); %Sort frequencies in ascending order
kR = k(I); %Sort wavenumbers by same order as frequencies

%Redo process for imaginary wavenumber using a longer and more sensitive
%method, so restrict the range
T = eye(n+1);
T(:,end) = [];
kI = 0:pi/1e4:pi/4;
omega2 = 2*pi*(8350:10:12010);
for i = 1:length(omega2)
    for j = 1:length(kI)
        T(end,1) = exp(1i*(pi + 1i*kI(j)));
        D = T'*(K-(omega2(i)^2)*M)*T;
        detD(i,j) = det(D/(3e11));
    end
    detD(i,:) = abs(detD(i,:).*conj(detD(i,:)));
end
[~,ktemp] = min(detD,[],2); %Find appropriate wavenumber for the frequency as the value which minimizes determinant detD

%% Plot Band Structure (Fig. 20)

figure
plot(kR(1:300),real(omega1(1:300))/(1e3*2*pi),'Color',[0,0,1],'LineWidth',1.5)
hold on
plot(kI(ktemp),(omega2)/(1e3*2*pi),'Color',[1,0,0],'LineWidth',1.5)
% plot(kappa3(:,1)*L,kappa3(:,2)/(1e3*2*pi))
hold off
ylim([0,25])
xlim([0,pi])
xlabel('Wavenumber $(\frac{\kappa}{a})$','Interpreter','latex','FontSize',16)
ylabel('Frequency (kHz)','Interpreter','latex','FontSize',16)
xticks([0, pi/4, pi/2, 3*pi/4, pi])
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
tmp1 = 2;
tmp2 = 0.25;
yline(8.36,'--','Color',[tmp2,tmp2,tmp2],'LineWidth',tmp1)
yline(12,'--','Color',[tmp2,tmp2,tmp2],'LineWidth',tmp1)
sizeFigureFonts('Times New Roman',16)
text(pi/2-0.2,10.25,'1st Bandgap','FontSize',16,'FontName','Times New Roman')
