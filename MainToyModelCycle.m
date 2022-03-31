%% Main script to generate Fig. 1
%%%Compute phenotypes of the cycle model
%
%   Copyright (c) 2022, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree.
clear
clc
close all

% Set figure size
Wi = 500;
Le = 330;

% Set number of states in the cycle
N = 3;

% Set active states, here only the first state leads to expression
% especially relevant when computing expression levels and propagated noise
Ip = false(N,1);
Ip(1) = true;
Np = sum(Ip);

% Set overall state occupancies
pa = 1/N;
P = [pa*ones(Np,1)/Np;ones(N-Np,1)*(1-pa)/(N-Np)];

% Set time scales
% completion time of the cycle/period (when alpha=1)
% define the max current 1/Tcyc
Tcycle = N;
% residence time of first state
T = pa*Tcycle/Np;

% Revesibility parameter for the cycle
% alpha = 0 fully reversible
% alpha = 1 fully irreversible
Na = 200;
la = logspace(-5,log10(0.5),Na/2-1);
alpha = [0,la,1-la(end:-1:1),1];

% Custom colormap
CMP = importdata('mycmap2.mat');
cmap = CMP(round(255*alpha+1),:);

% W acts as either a time or frequency vector
Nw = 200;
W = logspace(-2,2,Nw);

% Storing phenotypes
Pa = zeros(Na,1);
Phi = zeros(Na,Nw);
S = zeros(Na,1);
J = zeros(Na,1);

Pws = zeros(Na,Nw);
Ac = zeros(Na,Nw);
Rx = zeros(Na,Nw);
Ssys = zeros(Na,Nw);
Smed = zeros(Na,Nw);

m1Ta = zeros(Na,1);
s1Ta = zeros(Na,1);
Fa = zeros(Na,Nw);
Ta = zeros(Na,Nw);

m1Ti = zeros(Na,1);
s1Ti = zeros(Na,1);
Fi = zeros(Na,Nw);
Ti = zeros(Na,Nw);

for i=1:Na
    % Compute transition rate as a function of alpha
    % forward wf and backward wb rate are defined such that overall
    % occupancies and residence times remain the same.
    a = alpha(i);
    w = 1/(T*(2-a));
    Z = P(1)./P;
    wf = w*Z;
    wb = (1-a)*w*Z;
    
    % Build state rate matrix
    M = makeRateMatrixCycle(wf,wb);
    
    % Compute various regulatory phenotypes
    [p,dph,s,k,Resid] = computePheno(M,Ip,W);
    
    % Occupancy of the first state, should be equal to pa above
    Pa(i) = sum(p(Ip)); 
    
    % Dynamical phenotypes
    Phi(i,:) = dph(1,:);
    Pws(i,:) = dph(2,:);
    Ac(i,:) = dph(3,:);
    Rx(i,:) = dph(4,:);
    Ssys(i,:) = dph(5,:);
    Smed(i,:) = dph(6,:);
    
    % Entropy and current at steady state
    % since the flux is constant through the cycle at steady state,
    % one can pick any current between adjacent states to get the current,
    % such as K = K(2,1);
    S(i) = s;
    J(i) = k(2,1);
    
    % Residence time
    m1Ta(i) = Resid.m1Ta;
    s1Ta(i) = Resid.s1Ta;
    Fa(i,:) = Resid.FTa;
    Ta(i,:) = Resid.tTa;
        
    m1Ti(i) = Resid.m1Ti;
    s1Ti(i) = Resid.s1Ti;
    Fi(i,:) = Resid.FTi;
    Ti(i,:) = Resid.tTi;
end

% 2-state model as reference
pa = sum(P(Ip));
kf = P(1)/(pa*T);
kb = kf*pa/(1-pa);
tc = 1/(kf+kb);
M = makeRateMatrixCycle(kf,kb);
ip = false(2,1);
ip(1) = true;
[p,dph,s,k] = computePheno(M,ip,W);
phi = dph(1,:);
pws = dph(2,:);
ssys = dph(5,:);
smed = dph(6,:);
j = k(2,1);

%% Make reaction network graph (make Fig. 1B)
clc
close all

% Pick alpha
% alpha = 0 fully reversible
a = 0;
% alpha = 1 fully irreversible
%a = 1;

% Build state rate matrix
w = 1/(T*(2-a));
Z = P(1)./P;
wf = w*Z;
wb = (1-a)*w*Z;

M = makeRateMatrixCycle(wf,wb);
M(1:(N+1):end) = 0;

%%% Fig. 1B
H1=figure(1);
set(H1,'position',[50 50 Wi 0.9*Wi],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);

cedges = lines(2);
G = digraph(M');
LWidths = 5*abs(G.Edges.Weight/max(G.Edges.Weight));
EColor = zeros(size(G.Edges,1),3);

vv = G.Edges.EndNodes;
Ib = mod(vv(:,1),N) == mod(vv(:,2)+1,N);
If = ~Ib;
EColor(If,:) = repmat(cedges(2,:),sum(If),1);
EColor(Ib,:) = repmat(cedges(1,:),sum(Ib),1);
if ~isempty(LWidths)
plot(h1,G,'LineWidth',LWidths,'EdgeColor',EColor,'ArrowSize',18,'MarkerSize',10,'NodeColor','k','NodeFontSize',28)
end
set(h1,'Xtick',[],'Ytick',[],'linewidth',1.5,'fontsize',24)
set(h1,'XColor','none','YColor','none')

%% Plot dynamical phenotypes (make Fig. 1D-E)
clc
close all

%%% Occupancy
H1=figure(1);
set(H1,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);
hold(h1,'on')

plot(h1,[0,1],pa*[1,1],'-k','linewidth',2)
for i=1:Na
    plot(h1,alpha(i:end),Pa(i:end),'color',cmap(i,:),'linewidth',2)
end
set(h1,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h1,'\alpha')
ylabel(h1,'P_a')
ylim(h1,[0,1])

%%% Entropy production
H2=figure(2);
set(H2,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h2 = axes('parent',H2);
hold(h2,'on')

plot(h2,[0,1],s*[1,1],'-k','linewidth',2)
for i=1:Na
    plot(h2,alpha(i:end),S(i:end),'color',cmap(i,:),'linewidth',2)
end
set(h2,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h2,'\alpha')
ylabel(h2,'Entropy prod S [k_B/min]')
Y = ylim(h2);
ylim(h2,[-1,Y(2)])

%%% Current
H3=figure(3);
set(H3,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h3 = axes('parent',H3);
hold(h3,'on')

plot(h3,[0,1],j*[1,1],'-k','linewidth',2)
for i=1:Na
    plot(h3,alpha(i:end),J(i:end),'color',cmap(i,:),'linewidth',2)
end
set(h3,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h3,'\alpha')
ylabel(h3,'Current J [1/min]')

%%% Entropy production vs Current, Fig. 1D
H4=figure(4);
set(H4,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h4 = axes('parent',H4);
hold(h4,'on')

plot(h4,[1,1]/Tcycle,[0,20],'--k','linewidth',2)
for i=1:Na
    plot(h4,J(i:end),S(i:end),'color',cmap(i,:),'linewidth',3)
end
set(h4,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h4,'Current J [1/min]')
ylabel(h4,'Entropy prod S [k_B/min]')
ylim(h4,[0,8])

%%% Mean residence time
H5=figure(5);
set(H5,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h5 = axes('parent',H5);
hold(h5,'on')

plot(h5,alpha,m1Ta,'linewidth',1.5)
plot(h5,alpha,m1Ti,'linewidth',1.5)
set(h5,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h5,'\alpha')
ylabel(h5,'<T> [min]')
Y=ylim(h5);
ylim(h5,[0,Y(2)])

%%% Std residence time
H6=figure(6);
set(H6,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h6 = axes('parent',H6);
hold(h6,'on')

plot(h6,alpha,s1Ta,'linewidth',1.5)
plot(h6,alpha,s1Ti,'linewidth',1.5)
set(h6,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h6,'\alpha')
ylabel(h6,'\sigma_T [min]')
Y=ylim(h6);
ylim(h6,[0,Y(2)])

%%% Active residence time distribution, Fig. 1E
H7=figure(7);
set(H7,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h7 = axes('parent',H7);
hold(h7,'on')

for i=1:Na
    plot(h7,Ta(i,:),Fa(i,:),'color',cmap(i,:),'linewidth',1.5)
end
plot(h7,Ta(1,:),kf*exp(-kf*Ta(1,:)),'-k','linewidth',2)
set(h7,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h7,'Residence time T_1 [min]')
ylabel(h7,'P(T_1)')
xlim(h7,[0,3])

%%% Inactive residence time distribution, Fig. 1E
H8=figure(8);
set(H8,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h8 = axes('parent',H8);
hold(h8,'on')

for i=1:Na
    plot(h8,Ti(i,:),Fi(i,:),'color',cmap(i,:),'linewidth',1.5)
end
plot(h8,Ti(1,:),kb*exp(-kb*Ti(1,:)),'-k','linewidth',2)
set(h8,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h8,'Residence time T_{2\rightarrow3} [min]')
ylabel(h8,'P(T_{2\rightarrow3})')
xlim(h8,[0,6])

%%% Noise
H9=figure(9);
set(H9,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h9 = axes('parent',H9);
hold(h9,'on')

for i=1:Na
    plot(h9,W,Phi(i,:),'color',cmap(i,:),'linewidth',1.5)
end
plot(h9,W,phi,'-k','linewidth',2)
set(h9,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h9,'\tau [min]')
ylabel(h9,'Propagated noise \Phi(\tau)')

%%% Noise reduction (3-state vs 2-state)
H10=figure(10);
set(H10,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h10 = axes('parent',H10);
hold(h10,'on')

for i=1:Na
    plot(h10,W,Phi(i,:)./phi,'color',cmap(i,:),'linewidth',1.5)
end
set(h10,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h10,'\tau [min]')
ylabel(h10,'\Phi_N(\tau)/\Phi_2(\tau)')

%%% Noise reduction (3-state vs 2-state)
H11=figure(11);
set(H11,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h11 = axes('parent',H11);
hold(h11,'on')

plot(h11,[0,1],[1,1],'-k','linewidth',2)
for i=1:Na
    plot(h11,alpha(i:end),Phi(i:end,end)./phi(end),'color',cmap(i,:),'linewidth',2)
end
set(h11,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h11,'\alpha')
ylabel(h11,'\Phi_N(\tau\rightarrow\infty)/\Phi_2(\tau\rightarrow\infty)')

%%% Power spectrum
H12=figure(12);
set(H12,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h12 = axes('parent',H12);
hold(h12,'on')

for i=1:Na
    plot(h12,W,Pws(i,:),'color',cmap(i,:),'linewidth',1.5)
end
set(h12,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h12,'\omega [1/min]')
ylabel(h12,'Spec dens S(\omega) [1/min]')
plot(h12,W,pws,'-k','linewidth',1.5)

%%% Auto-correlation
H13=figure(13);
set(H13,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h13 = axes('parent',H13);
hold(h13,'on')

for i=1:Na
    plot(h13,W,Ac(i,:),'color',cmap(i,:),'linewidth',1.5)
end
set(h13,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h13,'\tau [min]')
ylabel(h13,'AC(\tau)')
plot(h13,W(1,:),pa*(1-pa)*exp(-W(1,:)/tc),'-k','linewidth',1.5)
xlim(h13,[0,2*Tcycle])

%%% Transient relaxation
H14=figure(14);
set(H14,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h14 = axes('parent',H14);
hold(h14,'on')

for i=1:Na
    plot(h14,W,Rx(i,:),'color',cmap(i,:),'linewidth',1.5)
end
plot(h14,W(1,:),1-exp(-W(1,:)/tc),'-k','linewidth',1.5)
set(h14,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h14,'t [min]')
ylabel(h14,'Relaxation P_1(t)/P_1(t\rightarrow\infty)')
xlim(h14,[0,2*Tcycle])

%%% Entropy production during transient (medium)
H15=figure(15);
set(H15,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h15 = axes('parent',H15);
hold(h15,'on')

for i=1:Na
    plot(h15,W,Ssys(i,:),'color',cmap(i,:),'linewidth',1.5)
end
plot(h15,W,ssys,'-k','linewidth',1.5)
set(h15,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h15,'t [min]')
ylabel(h15,'S_{sys} [k_B/min]')
xlim(h15,[0,2*Tcycle])

%%% Entropy production during transient (system)
H16=figure(16);
set(H16,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h16 = axes('parent',H16);
hold(h16,'on')

for i=1:Na
    plot(h16,W,Smed(i,:),'color',cmap(i,:),'linewidth',1.5)
end
plot(h16,W,smed,'-k','linewidth',1.5)
set(h16,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h16,'t [min]')
ylabel(h16,'S_{med} [k_B/min]')
xlim(h16,[0,2*Tcycle])

%% Generate stochastic trajectories (make Fig. 1C)
clc
close all

% Production rate of molecules in the active states kr
kr = 5;

% Length of time trace
tend = Tcycle*10;
% Temporal resolution for display
dt = T/100;
% Initial condition, i.e. the initial state
k0 = find(~Ip,1,'first');

% Reversible case alpaha=0
a = 0;
w = 1/(T*(2-a));
Z = P(1)./P;
wf = w*Z;
wb = (1-a)*w*Z;

% Make state rate matrix
M = makeRateMatrixCycle(wf,wb);
% Generate trajectories
[S0,T0] = genStoTrajectories(M,Ip,kr,dt,tend,k0);

% Irreversible case alpaha=1
a = 1;
w = 1/(T*(2-a));
Z = P(1)./P;
wf = w*Z;
wb = (1-a)*w*Z;

% Make state rate matrix
M = makeRateMatrixCycle(wf,wb);
% Generate trajectories
[S1,T1] = genStoTrajectories(M,Ip,kr,dt,tend,k0);

%%% State trajectory & winding number, alpha=0, Fig. 1C left
H1=figure(1);
set(H1,'position',[50 50 1.2*Wi 1.2*Le],'paperpositionmode','auto','color','w');
sgtitle(H1,'\alpha=0','fontsize',22,'fontweight','bold')

h1a = subplot(2,1,1,'parent',H1);
hold(h1a,'on')
stairs(h1a,T0,S0(:,1),'linewidth',2,'color',cmap(1,:))
set(h1a,'fontsize',20,'linewidth',2,'ytick',1:N,'tickdir','out')
ylim(h1a,[0.5,N+0.5])
ylabel(h1a,'States')

h1b = subplot(2,1,2,'parent',H1);
hold(h1b,'on')
plot(h1b,[0,tend],[0,0],'--k','linewidth',2)
stairs(h1b,T0,S0(:,3)/N,'linewidth',2,'color',cmap(1,:))
set(h1b,'fontsize',20,'linewidth',2,'tickdir','out')
xlabel(h1b,'Time [min]')
ylabel(h1b,'Winding #')

%%% State trajectory & winding number, alpha=1, Fig. 1C right
H2=figure(2);
set(H2,'position',[50 50 1.2*Wi 1.2*Le],'paperpositionmode','auto','color','w');
sgtitle(H2,'\alpha=1','fontsize',22,'fontweight','bold')

h2a = subplot(2,1,1,'parent',H2);
hold(h2a,'on')
stairs(h2a,T1,S1(:,1),'linewidth',2,'color',cmap(end,:))
set(h2a,'fontsize',20,'linewidth',2,'ytick',1:N,'tickdir','out')
ylim(h2a,[0.5,N+0.5])
ylabel(h2a,'States')

h2b = subplot(2,1,2,'parent',H2);
hold(h2b,'on')
plot(h2b,[0,tend],[0,tend/Tcycle],'--k','linewidth',2)
stairs(h2b,T1,S1(:,3)/N,'linewidth',2,'color',cmap(end,:))
set(h2b,'fontsize',20,'linewidth',2,'tickdir','out')
xlabel(h2b,'Time [min]')
ylabel(h2b,'Winding #')

rmax = max([S0(:,2);S1(:,2)]);

%%% Occupancy of active state & production events, alpha=0
H3=figure(3);
set(H3,'position',[50 50 1.2*Wi 1.2*Le],'paperpositionmode','auto','color','w');
sgtitle(H3,'\alpha=0','fontsize',22,'fontweight','bold')

h3a = subplot(2,1,1,'parent',H3);
hold(h3a,'on')
stairs(h3a,T0,double(S0(:,1)==1),'linewidth',2,'color',cmap(1,:))
set(h3a,'fontsize',20,'linewidth',2,'ytick',0:1,'tickdir','out')
ylim(h3a,[-0.1,1+0.1])
ylabel(h3a,'Occupancy')

h3b = subplot(2,1,2,'parent',H3);
hold(h3b,'on')
bar(h3b,T1,S0(:,2),1,'FaceColor',cmap(1,:))
set(h3b,'fontsize',20,'linewidth',2,'ytick',0:rmax,'tickdir','out')
ylim(h3b,[0,rmax+0.1])
xlabel(h3b,'Time [min]')
ylabel(h3b,'Production')

%%% Occupancy of active state & production events, alpha=1
H4=figure(4);
set(H4,'position',[50 50 1.2*Wi 1.2*Le],'paperpositionmode','auto','color','w');
sgtitle(H4,'\alpha=1','fontsize',22,'fontweight','bold')

h4a = subplot(2,1,1,'parent',H4);
hold(h4a,'on')
stairs(h4a,T1,double(S1(:,1)==1),'linewidth',2,'color',cmap(end,:))
set(h4a,'fontsize',20,'linewidth',2,'ytick',0:1,'tickdir','out')
ylim(h4a,[-0.1,1+0.1])
ylabel(h4a,'Occupancy')

h4b = subplot(2,1,2,'parent',H4);
hold(h4b,'on')
bar(h4b,T1,S1(:,2),1,'FaceColor',cmap(end,:))
set(h4b,'fontsize',20,'linewidth',2,'ytick',0:rmax,'tickdir','out')
ylim(h4b,[0,rmax+0.1])
xlabel(h4b,'Time [min]')
ylabel(h4b,'Production')
