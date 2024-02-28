%% Load data and initialize
DataEx1 = importdata("Comparison Data\Abaqus Data Ex1 ElNum1.txt");
DataEx1_2 = importdata("Comparison Data\Abaqus Data Ex1 ElNum2.txt");
DataEx2_Pt1 = importdata("Comparison Data\Abaqus Data Ex 2 Try 2.txt");
DataEx2_Pt2 = importdata("Comparison Data\Abaqus Data Ex 2 Second Pt.txt");
DataEx3_Pt1 = importdata("Comparison Data\Abaqus Data Ex 3 Try 2.txt");
DataEx3_Pt2 = importdata("Comparison Data\Abaqus Data Ex 3 Second Pt.txt");

E = 68.9e9;
sY = 276e6;
eY = (sY/E);
rho = 2700;
Ce = sqrt(E/rho);
L=1;
%% Elastic Step Wave Example (Fig.7e)

eA = (1e8/E);
t1 = linspace(0,L/(3*Ce),100);
t2 = linspace(L/(3*Ce),2*(L/(3*Ce)),100);
t3 = linspace((2*L/(3*Ce)),4*L/(3*Ce),100);
t4 = linspace(4*L/(3*Ce),3*L/(2*Ce),100);

e1 = 0*ones(size(t1));
e2 = -eA*ones(size(t2));
e3 = -2*eA*ones(size(t3));
e4 = -eA*ones(size(t4));
figure
plot(1e3*[t1,t2,t3,t4],[e1,e2,e3,e4]/eY,'LineWidth',1.5)
hold on
plot(1e3*DataEx1(:,1),DataEx1(:,2)/eY,'LineWidth',1.5)
plot(1e3*DataEx1_2(:,1),DataEx1_2(:,2)/eY,'LineWidth',1.5)
hold off
ylabel('Strain $(\frac{\varepsilon}{\varepsilon_{Y}})$','Interpreter','latex')
xlabel('Time (ms)','Interpreter','latex')
title('Strain Elastic Example')
xlim(1e3*[0,3*L/(2*Ce)])
yticks([-0.8,-0.6,-0.4,-0.2,0])
xticks([0,0.05,0.1,0.15,0.2,0.25])
ylim([-0.8,0.05])
sizeFigureFonts('Times New Roman',14)

%% Plastic Step Wave Example (Fig. 8b)

Cp = Ce/3;
sA = 3e8;
K = E/9;
eA = sY/E + (sA-sY)/K;

t1 = linspace(0,L/(3*Ce),100);
t2 = linspace(L/(3*Ce),(L/(3*Cp)),100);
t3 = linspace((L/(3*Cp)),5*L/(3*Ce),100);
t4 = linspace(5*L/(3*Ce),2*L/(Ce),100);

e1 = 0*ones(size(t1));
e2 = -eY*ones(size(t2));
e3 = -eA*ones(size(t3));
s4 = -0.5*(sA-sY)*(1+Ce/Cp);
er3 = -eA+sA/E;
e4 = (s4/E + er3)*ones(size(t4));
figure
plot(1e3*[t1,t2,t3,t4],[e1,e2,e3,e4]/eY,'LineWidth',1.5)
hold on
plot(1e3*DataEx2_Pt1(:,1),DataEx2_Pt1(:,2)/eY,'LineWidth',1.5)
hold off
ylabel('Strain $(\frac{\varepsilon}{\varepsilon_{Y}})$','Interpreter','latex')
xlabel('Time (ms)','Interpreter','latex')
title('Strain L/3 Plastic Example')
ylim([-1.85,0.05])
yticks([-1.5,-1,-0.5,0])
xticks([0,0.1,0.2,0.3,0.4])
sizeFigureFonts('Times New Roman',14)

%Figure 11
xstar = 2*L*Cp/(Ce+Cp);
tstar = 2*L/(Ce+Cp);
t = 1.15*tstar;
X53 = xstar + (t-tstar)*Ce;
x24 = xstar - (t-tstar)*Ce;
s2 = sA;
s4 = (sA-sY)*(Ce+Cp)/(2*Cp);
s5 = s4;
s3 = 0;
e2 = eA;
er2 = eA-s2/E;
e4 = er2 + s4/E;
er4 = er2;
e5 = s5/E;
er5 = 0;
e3 = 0;
er3 = 0;
dx = L/1000;
X2 = 0:dx:x24-dx;
X4 = x24:dx:xstar-dx;
X5 = xstar:dx:X53-dx;
X3 = X53:dx:L;
X = [X2,X4,X5,X3];
Strain = [e2*ones(1,length(X2)),e4*ones(1,length(X4)),e5*ones(1,length(X5)),e3*ones(1,length(X3))];
ResidualStrain = [er2*ones(1,length(X2)),er4*ones(1,length(X4)),er5*ones(1,length(X5)),er3*ones(1,length(X3))];
figure
plot(X,Strain/eY,'LineWidth',1.5)
hold on
plot(X,ResidualStrain/eY,'LineWidth',1.5)
hold off
xlabel('Position $(X)$','Interpreter','latex')
ylabel('Strain','Interpreter','latex')

legend({'$\varepsilon/\varepsilon_{Y}$','$\varepsilon^{p}/\varepsilon_{Y}$'},'Interpreter','latex')
sizeFigureFonts('Times New Roman',14)
%% Plastic Unloading Example (Fig.12b, 11)
%Figure 12b
Cp = Ce/2;
sA = 3e8;
K = E/4;
eA = sY/E + (sA-sY)/K;
T = 1e-4;
t1 = linspace(0,L/(3*Ce),100);
t2 = linspace(L/(3*Ce),(L/(3*Cp)),100);
Xstar = ((Cp*Ce)/(Ce-Cp))*T;
tstar = Ce/(Ce-Cp)*T;
t4s = tstar + (Xstar-L/3)/(Ce);
t3 = linspace(L/(3*Cp),L/(3*Ce)+T,100);
t4 = linspace(L/(3*Ce)+T,t4s,100);
t5 = linspace(t4s,t4s+2*L/(3*Ce),100);
s5 = 0.5*(sA-sY)*(1-Ce/Cp);
e1 = 0*ones(size(t1));
e2 = -eY*ones(size(t2));
er3 = -eA+sA/E;
e3 = -eA*ones(size(t3));
e4 = (er3)*ones(size(t4));
e5 = (s5/E+er3)*ones(size(t5));
tmp = ceil(1/(3*dX));
figure
plot(1e3*[t1,t2,t3,t4,t5],[e1,e2,e3,e4,e5]/eY,'LineWidth',1.5)
hold on
plot(1e3*DataEx3_Pt1(:,1),DataEx3_Pt1(:,2)/eY,'LineWidth',1.5)
plot(1e3*t,SA(1).Strain(:,tmp)/eY,'--','LineWidth',1.5)
hold off
ylabel('Strain $(\frac{\varepsilon}{\varepsilon_{Y}})$','Interpreter','latex')
xlabel('Time (ms)','Interpreter','latex')
xlim([0,0.3])
ylim([-1.4,0.05])
yticks([-1.25,-1,-0.75,-0.5,-0.25,0])
xticks([0,0.1,0.2,0.3])
sizeFigureFonts('Times New Roman',14)

