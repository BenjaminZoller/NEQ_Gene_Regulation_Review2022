%% Main script to generate Fig. 2 left
%%% Compute phenotypes of the proof-reading model
clear
clc
close all

% Set figure size
Wi = 500;
Le = 330;

% Custom colormap
cmap = importdata('mycmap2.mat');

% TF binding rate, proportional to TF concentration
Nc = 200;
kb = logspace(-3,5,Nc);

% Sampling the Proof-reading rate parameter
Nq = 200;
kq = logspace(-5,5,Nq);

% TF unbinding rate, related to affinity
% for specific (strong) binding site
ku = 1;
% for non-specific binding site
ku_ns = 1e2;

% Expression levels
ES = zeros(Nc,Nq);
ENS = zeros(Nc,Nq);

for i=1:Nc
    for j=1:Nq
        % One could use the numerical solution for occupancies based on
        % the state rate matrix of the model, cf. "computePheno.m" function
        %M = makeRateMatrixBinding(kb(i),ku,kq(j));
        %Ip = false(3,1);
        %Ip(1) = true;
        %ES(i,j) = getExp(M,Ip)
        
        % Since the model is easy to solve analyticaly,
        % we for the analytical solution
        % Compute occupancy for specific binding site
        P = zeros(3,1);
        P(3) = ku^2+kq(j)*ku;
        P(2) = kb(i)*ku;
        P(1) = kb(i)*kq(j);
        P = P/sum(P);
        ES(i,j) = P(1);
        
        % Compute occupancy for non-specific
        P = zeros(3,1);
        P(3) = ku_ns^2+kq(j)*ku_ns;
        P(2) = kb(i)*ku_ns;
        P(1) = kb(i)*kq(j);
        P = P/sum(P);
        ENS(i,j) = P(1);
    end
end

% Definition of specificity
% ratio of specific expression level over non-specific
S = ES./ENS;

%% Make reaction network graph (make Fig. 2A left)
clc
close all

% Build state rate matrix
M = makeRateMatrixBinding(1e0,1e0,1e0);
N = size(M,1);
M(1:(N+1):end) = 0;

%%% Fig. 2A left
H1=figure(1);
set(H1,'position',[50 50 Wi 0.9*Wi],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);

cedges = lines(3);
G = digraph(M');
LWidths = 5*abs(G.Edges.Weight/max(G.Edges.Weight));
EColor = zeros(size(G.Edges,1),3);

vv = G.Edges.EndNodes;
Iu = false(4,1);
Ib = false(4,1);
Iq = false(4,1);
Iu([1,3]) = true;
Ib(4) = true;
Iq(2) = true;
EColor(Iu,:) = repmat(cedges(1,:),sum(Iu),1);
EColor(Ib,:) = repmat(cedges(2,:),sum(Ib),1);
EColor(Iq,:) = repmat(cedges(3,:),sum(Iq),1);
if ~isempty(LWidths)
pp=plot(h1,G,'LineWidth',LWidths,'EdgeColor',EColor,'ArrowSize',12,'MarkerSize',10,...
    'NodeColor','k','NodeFontSize',20,'EdgeFontSize',20,'EdgeFontAngle','normal');
end
set(h1,'Xtick',[],'Ytick',[],'linewidth',1.5,'fontsize',24)
set(h1,'XColor','none','YColor','none')
labelnode(pp,1,'Active')
labelnode(pp,2,'Bound')
labelnode(pp,3,'Unbound')
pp.EdgeLabel = {'k_{-}','k_q','k_{-}','k_{+}'};

%% Plot phenotypes (make Fig. 2B left)
clc
close all

% Targeted level of expression
E0 = 0.5;

% Optimal kq at E=0.5, see lines 576-583 below for computation
kq_opti = 2.4392;
[~,k] = min(abs(kq-kq_opti));

% Set proof-reading paramter to optimal value for non-eq model
kq_neq = kq(k);
% Set proof-reading paramter to large value (~infinity) for eq model
kq_eq = 1e5;

% Calculate the concentration required to target E=0.5 for eq and non-eq model
zeq = 1+ku/kq_eq;
zneq = 1+ku/kq_neq;
kbeq = ku*E0*zeq/(1-E0*zeq);
kbneq = ku*E0*zneq/(1-E0*zneq);

myc = lines(2);

%%% Expression vs Concentration, Fig. 2B left
H1=figure(1);
set(H1,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);
hold(h1,'on')

plot(h1,kb,ES(:,end),'-','color',myc(1,:),'linewidth',2)
plot(h1,kb,ENS(:,end),'--','color',myc(1,:),'linewidth',2)
plot(h1,kb,ES(:,k),'-','color',myc(2,:),'linewidth',2)
plot(h1,kb,ENS(:,k),'--','color',myc(2,:),'linewidth',2)
plot(h1,kbeq,E0,'ok','markersize',10,'linewidth',2)
plot(h1,kbneq,E0,'*k','markersize',12,'linewidth',2)
set(h1,'fontsize',22,'linewidth',2,'xscale','log','ytick',0:1,'xtick',[1e-3,1e1,1e5],'tickdir','out')
xlabel(h1,'Concentration k_{+}')
ylabel(h1,'Expression E')
xlim(h1,[kb(1),kb(end)])

%%% Specificity vs Concentration
H2=figure(2);
set(H2,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h2 = axes('parent',H2);
hold(h2,'on')

plot(h2,kb,S(:,end),'-k','linewidth',2)
plot(h2,kb,S(:,k),'-','color',myc(2,:),'linewidth',2)
set(h2,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h2,'Concentration k_{+}')
ylabel(h2,'Specificity S=E^S/E^{NS}')
xlim(h2,[kb(1),kb(end)])

%%% Expression as a function of Concentration and Proof-reading rate
H3=figure(3);
set(H3,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h3 = axes('parent',H3);
hold(h3,'on')

imagesc(h3,kb,kq,ES)
set(h3,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h3,'Concentration k_{+}')
ylabel(h3,'k_q')
cb=colorbar(h3,'LineWidth',2);
set(get(cb,'ylabel'),'String','E','fontsize',22);
xlim(h3,[kb(1),kb(end)])
colormap(h3,cmap)

%%% Specificity as a function of Concentration and Proof-reading rate
H4=figure(4);
set(H4,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h4 = axes('parent',H4);
hold(h4,'on')

imagesc(h4,kb,kq,S)
set(h4,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','colorscale','log','tickdir','out')
xlabel(h4,'Concentration k_{+}')
ylabel(h4,'k_q')
cb=colorbar(h4,'LineWidth',2);
set(get(cb,'ylabel'),'String','S','fontsize',22);
xlim(h4,[kb(1),kb(end)])
colormap(h4,cmap)

%%% Specificity gain & Maximum expression vs Proof-reading rate
H5=figure(5);
set(H5,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h5 = axes('parent',H5);
hold(h5,'on')

Emax = @(kq) kq./(kq+ku);
Sgain = @(kq) (ku_ns+kq)./(ku+kq);
yyaxis(h5,'left')
plot(h5,kq,Sgain(kq),'-','linewidth',2)
set(h5,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
ylabel(h5,'Specificity gain S/S_{EQ}')
yyaxis(h5,'right')
plot(h5,kq,Emax(kq),'-','linewidth',2)
plot(h5,kq(k)*[1,1],[0,1],'--k','linewidth',2)
ylabel(h5,'Maximum expression E_{max}')
xlabel(h5,'k_q')

%% Unbinding identifiability (at fixed concentration kb)

% Binding rate kb, i.e. TF concentration is fixed and the same for both model
kb0 = 1e1;
% Sample specific TF unbinding rate, or equivalently TF specific affinity
kus = logspace(-1,3,Nc);

% Set proof-reading paramter to optimal value for non-eq model
kq_neq = kq_opti;
% Set proof-reading paramter to large value (~infinity) for eq model
kq_eq = 1e5;

% Expression levels
Eneq = zeros(Nc,1);
Eeq = zeros(Nc,1);

% TF occupancy
Pneq = zeros(Nc,1);
Peq = zeros(Nc,1);

for i=1:Nc
    % Analytical solution non-eq
    P = zeros(3,1);
    P(3) = kus(i)^2+kq_neq*kus(i);
    P(2) = kb0*kus(i);
    P(1) = kb0*kq_neq;
    P = P/sum(P);
    Eneq(i) = P(1);
    Pneq(i) = P(1)+P(2);
    
    % Analytical solution eq
    P = zeros(3,1);
    P(3) = kus(i)^2+kq_eq*kus(i);
    P(2) = kb0*kus(i);
    P(1) = kb0*kq_eq;
    P = P/sum(P);
    Eeq(i) = P(1);
    Peq(i) = P(1)+P(2);
end

%%% TF occupancy vs TF binding site affinity
H6=figure(6);
set(H6,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h6 = axes('parent',H6);
hold(h6,'on')

plot(h6,ku_ns./kus,Eeq,'-k','linewidth',2)
plot(h6,ku_ns./kus,Eneq,'-','color',myc(2,:),'linewidth',2)
set(h6,'fontsize',22,'linewidth',2,'xscale','log','ytick',0:1,'tickdir','out')
xlabel(h6,'Affinity k_{-}^{NS}/k_{-}')
ylabel(h6,'Expression E')
xlim(h6,[min(ku_ns./kus),max(ku_ns./kus)])

%%% Expression vs TF binding site affinity
H7=figure(7);
set(H7,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h7 = axes('parent',H7);
hold(h7,'on')

plot(h7,ku_ns./kus,Peq,'-k','linewidth',2)
plot(h7,ku_ns./kus,Pneq,'--','color',myc(2,:),'linewidth',2)
set(h7,'fontsize',22,'linewidth',2,'xscale','log','ytick',0:1,'tickdir','out')
xlabel(h7,'Affinity k_{-}^{NS}/k_{-}')
ylabel(h7,'TF occupancy')
xlim(h7,[min(ku_ns./kus),max(ku_ns./kus)])

%% Unbinding identifiability (with concentration kb set such that E=0.5)

% Targeted level of expression
E0 = 0.5;

% Specific unbinding rate for a typical strong biding site, 
% namely ku = 1e2*ku_ns
ku = 1;

% Sample the TF unbinding rate, i.e. affinity
kus = logspace(-1,3,Nc);

% Sample both neq and eq model
kq = [kq_opti,1e5];

% Expression levels
ES = zeros(length(kq),Nc);
% TF occupancy
PS = zeros(length(kq),Nc);

for i=1:length(kq)
    % Determine concentration to achieve trageted expression level E0
    % with a typical strong binding site
    EE = E0*(1+ku/kq(i));
    kb0 = ku*EE/(1-EE);
    for j=1:Nc        
        %analytical solution specific
        P = zeros(3,1);
        P(3) = kus(j)^2+kq(i)*kus(j);
        P(2) = kb0*kus(j);
        P(1) = kb0*kq(i);
        P = P/sum(P);
        ES(i,j) = P(1);
        PS(i,j) = P(1)+P(2);
    end
end

%%% TF occupancy vs TF binding site affinity
H8=figure(8);
set(H8,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h8 = axes('parent',H8);
hold(h8,'on')

plot(h8,ku_ns./kus,ES(end,:),'-k','linewidth',2)
plot(h8,ku_ns./kus,ES(1,:),'-','color',myc(2,:),'linewidth',2)
set(h8,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h8,'Affinity k_{-}^{NS}/k_{-}')
ylabel(h8,'Expression E')
xlim(h8,[min(ku_ns./kus),max(ku_ns./kus)])

%%% Expression vs TF binding site affinity
H9=figure(9);
set(H9,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h9 = axes('parent',H9);
hold(h9,'on')

plot(h9,ku_ns./kus,PS(end,:),'-k','linewidth',2)
plot(h9,ku_ns./kus,PS(1,:),'-','color',myc(2,:),'linewidth',2)
set(h9,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h9,'Affinity k_{-}^{NS}/k_{-}')
ylabel(h9,'TF occupancy')
xlim(h9,[min(ku_ns./kus),max(ku_ns./kus)])


%% Sampling dynamical phenotypes at fixed expression
clear
clc
close all

% Set figure size
Wi = 500;
Le = 330;

% Custom colormap
cmap = importdata('mycmap2.mat');
ncmap = size(cmap,1)-1;

% Targeted level of expression
E0 = 0.5;

% Define which states are active
Ip = false(3,1);
% for active transcription
Ip(1) = true;
% for TF binding
%Ip(1:2) = true;

% Specific unbinding rate for a typical strong biding site, 
% namely ku = 1e2*ku_ns
ku = 1;

% Sampling the Proof-reading rate parameter
% kqmin is the minimal value that allows E0
Nq = 200;
kqmin = E0*ku/(1-E0);
lkqmax = 5;
kq = logspace(log10(kqmin),lkqmax,Nq);

Nw = 200;
W = logspace(-2,2,Nw);
Phi = zeros(Nq,Nw);
Pws = zeros(Nq,Nw);
Ac = zeros(Nq,Nw);
Rx = zeros(Nq,Nw);

m1Ta = zeros(Nq,1);
s1Ta = zeros(Nq,1);
Fa = zeros(Nq,Nw);
Ta = zeros(Nq,Nw);

m1Ti = zeros(Nq,1);
s1Ti = zeros(Nq,1);
Fi = zeros(Nq,Nw);
Ti = zeros(Nq,Nw);

for i=1:Nq
    % Determine concentration to achieve trageted expression level E0
    % with a typical strong binding site
    EE = E0*(1+ku/kq(i));
    kb0 = ku*EE/(1-EE);
    
    M = makeRateMatrixBinding(kb0,ku,kq(i));
    [p,dph,~,~,Resid] = computePheno(M,Ip,W);
    
    Phi(i,:) = dph(1,:); 
    Pws(i,:) = dph(2,:);
    Ac(i,:) = dph(3,:);
    Rx(i,:) = dph(4,:);
    
    m1Ta(i) = Resid.m1Ta;
    s1Ta(i) = Resid.s1Ta;
    Fa(i,:) = Resid.FTa;
    Ta(i,:) = Resid.tTa;
        
    m1Ti(i) = Resid.m1Ti;
    s1Ti(i) = Resid.s1Ti;
    Fi(i,:) = Resid.FTi;
    Ti(i,:) = Resid.tTi;  
end

%% Plot dynamical phenotypes
clc
close all

lkqi = -log10(kq);

%%% Noise
H1=figure(1);
set(H1,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);
hold(h1,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h1,W,Phi(i,:),'color',cmap(cr,:),'linewidth',1.5)
end
plot(W,Phi(end,:),'-k','linewidth',2)
set(h1,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h1,'\tau [min]')
ylabel(h1,'Propagated noise \Phi(\tau)')

%%% Active residence time distribution
H2=figure(2);
set(H2,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h2 = axes('parent',H2);
hold(h2,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h2,Ta(i,:),Fa(i,:),'color',cmap(cr,:),'linewidth',1.5)
end
plot(h2,Ta(end,:),Fa(end,:),'-k','linewidth',2)
set(h2,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h2,'T_a [min]')
ylabel(h2,'P(T_a)')

%%% Inactive residence time distribution
H3=figure(3);
set(H3,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h3 = axes('parent',H3);
hold(h3,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h3,Ti(i,:),Fi(i,:),'color',cmap(cr,:),'linewidth',1.5)
end
plot(h3,Ti(end,:),Fi(end,:),'-k','linewidth',2)
set(h3,'fontsize',22,'linewidth',2,'tickdir','out')
xlabel(h3,'T_i [min]')
ylabel(h3,'P(T_i)')
ylim(h3,[0,1])

%%% Power spectrum
H4=figure(4);
set(H4,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h4 = axes('parent',H4);
hold(h4,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h4,W,Pws(i,:),'color',cmap(cr,:),'linewidth',1.5)
end
plot(h4,W,Pws(end,:),'-k','linewidth',1.5)
set(h4,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h4,'\omega [1/min]')
ylabel(h4,'Spec dens S(\omega) [1/min]')

%%% Auto-correlation
H5=figure(5);
set(H5,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h5 = axes('parent',H5);
hold(h5,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h5,W,Ac(i,:),'color',cmap(cr,:),'linewidth',1.5)
end
set(h5,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h5,'\tau [min]')
ylabel(h5,'AC(\tau)')
plot(h5,W,Ac(end,:),'-k','linewidth',1.5)

%%% Transient relaxation
H6=figure(6);
set(H6,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h6 = axes('parent',H6);
hold(h6,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h6,W,Rx(i,:),'color',cmap(cr,:),'linewidth',1.5)
end
plot(h6,W,Rx(end,:),'-k','linewidth',1.5)
set(h6,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h6,'t [min]')
ylabel(h6,'Relaxation P_a(t)/P_a')

%% Sampling specificity at fixed E=0.5 (Fig. 2C-D left)
clear
clc
close all

% Set figure size
Wi = 500;
Le = 330;

% Custom colormap
cmap = importdata('mycmap2.mat');
ncmap = size(cmap,1)-1;

% Targeted level of expression
E0 = 0.5;

% TF unbinding rate, related to affinity
% for non-specific binding site
ku_ns = 1e2;

% Sampling the specific unbiding rate
Nu = 200;
ku = logspace(log10(0.5),2,Nu);

% Sampling the proof-reading rate
Nq = 200;
lkqmax = 5;

% Setting the active state leading to expression
Ip = false(3,1);
Ip(1) = true;

% Setting the lifetime of expressed molecules (proteins)
% this time scale defines the typical duration over which expression noise
% is averaged
tp = 6e4;

% Storing the phenotypes
ES = zeros(Nu,Nq);
ENS = zeros(Nu,Nq);
PHI = zeros(Nu,Nq);
KQ = zeros(Nu,Nq);

for i=1:Nu
    % Sampling the proof-reading rate
    kqmin = E0*ku(i)/(1-E0);
    kqmin = kqmin*(1+1e-4);
    kq = logspace(log10(kqmin),lkqmax,Nq);
    KQ(i,:) = kq;
    for j=1:Nq
        % Determine concentration to achieve trageted expression level E0
        EE = E0*(1+ku(i)/kq(j));
        kb0 = ku(i)*EE/(1-EE);
        
        % Analytical solution for specific expression
        P = zeros(3,1);
        P(3) = ku(i)^2+kq(j)*ku(i);
        P(2) = kb0*ku(i);
        P(1) = kb0*kq(j);
        P = P/sum(P);
        ES(i,j) = P(1);
        
        % Compute noise
        M = makeRateMatrixBinding(kb0,ku(i),kq(j));
        
        PHI(i,j) = getNoise(M,Ip,tp);
        
%         Mr = M;
%         N = size(Mr,1)-1;
%         k0 = Mr(1:N,end);
%         Mr = Mr(1:N,1:N) - k0*ones(1,N);
%         
%         sn = P(1)*(1-P(1));
%         v0 = [1,0];
%         v1 = zeros(N,1);
%         v1(1) = sn;
%         v1(2) = -P(1)*P(2);
% 
%         PHI(i,j) = v0*((eye(N)-tp*Mr)\v1)/sn;

        % Analytical solution for non-specific expression
        P = zeros(3,1);
        P(3) = ku_ns^2+kq(j)*ku_ns;
        P(2) = kb0*ku_ns;
        P(1) = kb0*kq(j);
        P = P/sum(P);
        ENS(i,j) = P(1);
    end
end

% Compute specificity
S = ES./ENS;

%% Plot regulatory phenotypes at fixed E (Fig. 2C-D left)
clc
close all

% Find optimal proof-reading rate kq_opti.
% The optimal value depends on the TF biding site affinity.
% Let's set ku=1 for a typical strong binding site
j = find(ku>=1,1,'first');
% For this specific ku, let's find kq that maximizes the specificity S
[~,k] = max(S(j,:));
% The optimal kq is then
kq_opti = KQ(j,k);

% Determine TF residence time in the bound state (Bound+Active)
KU = repmat(ku(:),1,Nq);
T = ku_ns./KU;

lkqi = -log10(KQ(j,:));

%%% Specificity vs Proof-reading rate
H1=figure(1);
set(H1,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);
hold(h1,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h1,KQ(j,i:end),S(j,i:end),'-','color',cmap(cr,:),'linewidth',3)
end
plot(h1,KQ(j,k),S(j,k),'*k','linewidth',2,'markersize',12)
plot(h1,KQ(j,end),S(j,end),'ok','linewidth',2,'markersize',10)
set(h1,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h1,'k_q')
ylabel(h1,'Specificity S')

%%% Specificity vs Proof-reading ratio, Fig. 2C left
H2=figure(2);
set(H2,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h2 = axes('parent',H2);
hold(h2,'on')

for i=1:Nq
    cr=1+round(ncmap*(lkqi(i)+lkqmax)/lkqmax);
    plot(h2,ku(j)./KQ(j,i:end),S(j,i:end),'-','color',cmap(cr,:),'linewidth',3)
end
plot(h2,ku(j)./KQ(j,k),S(j,k),'*k','linewidth',2,'markersize',12)
plot(h2,ku(j)./KQ(j,end),S(j,end),'ok','linewidth',2,'markersize',10)
set(h2,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h2,'Proof-reading ratio k_{-}/k_q')
ylabel(h2,'Specificity S')

%%% Specificity vs Residence time, Fig. 2D left
H3=figure(3);
set(H3,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h3 = axes('parent',H3);
hold(h3,'on')

[~,Is] = sort(1./KQ(:),'descend');
scatter(h3,T(Is),S(Is),[],1./KQ(Is),'filled')
plot(h3,ku_ns./ku,S(:,end),'-k','linewidth',3)
plot(h3,ku_ns./ku(j),S(j,k),'*k','linewidth',2,'markersize',12)
plot(h3,ku_ns./ku(j),S(j,end),'ok','linewidth',2,'markersize',10)
set(h3,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','colorscale','log','tickdir','out')
xlabel(h3,'Residence time T_{A} [1/k_{-}^{NS}]')
ylabel(h3,'Specificity S')
xlim(h3,[ku_ns./ku(end),ku_ns./ku(1)])
colormap(h3,cmap)

%%% Noise vs Residence time
H4=figure(4);
set(H4,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h4 = axes('parent',H4);
hold(h4,'on')

scatter(h4,T(Is),PHI(Is),[],1./KQ(Is),'filled')
plot(h4,ku_ns./ku,PHI(:,end),'-k','linewidth',3)
plot(h4,ku_ns./ku(j),PHI(j,k),'*k','linewidth',2,'markersize',12)
plot(h4,ku_ns./ku(j),PHI(j,end),'ok','linewidth',2,'markersize',10)
set(h4,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','colorscale','log','tickdir','out')
xlabel(h4,'Residence time T_{A} [1/k_{-}^{NS}]')
ylabel(h4,'Propagated noise \Phi')
xlim(h4,[ku_ns./ku(end),ku_ns./ku(1)])
colormap(h4,cmap)
