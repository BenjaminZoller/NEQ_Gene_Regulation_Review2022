%% Main script to generate Fig. 2 right
%%% Compute phenotypes of the cooperativity model
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

% Custom colormap
cmap = importdata('mycmap2.mat');
ncmap = size(cmap,1)-1;

warning('off')
spmd
    warning('off')
end

% Number of binding sites N
N = 3;
% Unbinding rate ku
ku = 1;
% Time averaging time scale for noise, i.e. lifetime of produced proteins
tp = 6e4;
% Define active states leading to expression
Ip = false(2^N,1);
% All-or-nothing expression scheme
Ip(1) = true;
% One-or-more expression scheme
%Ip(1:(end-1)) = true;

% Expression E is proprotional to the occupancy of the active states
% Here for simplification, we define E as being the occupancy,
% namely E in [0,1]
% Trageted expression level E0
E0 = 0.5;
Emin = 0;
Emax = 1;

% Sampling cooperativities, a1 & a2
Ns = 200;
amax = 3; %3 seems sufficient to sample the intersting region
alpha = logspace(0,amax,Ns);
% (Optional) Use this range for a better sampling of induction curves
% Ns = 300;
% alpha = logspace(1.9,2.3,Ns);

% Storing regulatory phenotypes
% non-eq phenotypes
A1 = nan(Ns,Ns);
A2 = nan(Ns,Ns);
KB = nan(Ns,Ns);
E = nan(Ns,Ns);
T = nan(Ns,Ns);
ST = nan(Ns,Ns);
S = nan(Ns,Ns);
PHI = nan(Ns,Ns);
HILL = nan(Ns,Ns);
% eq phenotypes
CEQ = nan(Ns,Ns);
KEQ = nan(Ns,Ns);
HEQ = nan(Ns,Ns);
PEQ = nan(Ns,Ns);

% Some numerical constants necessary for search & optimization function
B1 = -20*log(10);
B2 = 10*log(10);
options = optimset('TolX',1e-16);
myeps = 1e-3; %for derivative
inveps = 1/myeps;

% Main loop to sample regulatory phenotypes over cooperativities
% Parallelized loop
parfor i=1:Ns
% Standard loop in case parfor cannot be use or for debugging    
%for i=1:Ns 
    disp(i)
    for j=1:Ns
        % Sample a1 & a2 randomely, such that a2>=a1
        if j==1
            co = alpha(i)*ones(1,N);
        else            
            r = rand(1);
            ai = log10(alpha(i));
            co = 10.^(ai+[0,(N*(ai+1))*r*ones(1,N-1)]);
        end
        A1(i,j) = co(1);
        A2(i,j) = co(2);
                
        myfun = @(x) getExp(makeRateMatrixCooperativity(exp(x),ku,co),Ip)-E0;
        if myfun(B1)*myfun(B2) < 0
            % Find kbind (concentration) to meet targeted expression E0
            [x,~,exitflag] = fzero(myfun,[B1,B2],options);
            if exitflag ~= 1
                exitflag
            end
            
            kb = exp(x);
            KB(i,j) = kb;
            
            % Make state rate matrix
            Mneq = makeRateMatrixCooperativity(kb,ku,co);
            
            % Compute expression
            [E(i,j),p] = getExp(Mneq,Ip)
            
            % Compute residence time
            [t,st] = getResid(Mneq,Ip,p);
            T(i,j) = t;
            ST(i,j) = st;
            
            % Compute noise
            PHI(i,j) = getNoise(Mneq,Ip,tp);
            
            % Compute entropy
            S(i,j) = getEntropy(Mneq,p);
            
            % Compute Hill coef
            ekb = myeps*kb;
            kb1 = kb-ekb;
            kb2 = kb+ekb;
            
            E1 = getExp(makeRateMatrixCooperativity(kb1,ku,co),Ip);
            E2 = getExp(makeRateMatrixCooperativity(kb2,ku,co),Ip);
            
            dy = E2 - E1;
            HILL(i,j) = inveps*2*dy/(Emax-Emin);
            
            %%% Build equivalent eq model achieving same expression level E
            % and same residence time T for the active state.
            % For the One-or-more expression scheme, the equ cooperativity
            % leading to same T must be computed numerically.
            % For the All-or-nothing expression scheme, the equ cooperativity
            % leading to same T can be calculated analytically easily.
            ceq = (N*t*ku)^(1/(N-1)); 
            CEQ(i,j) = ceq;
            
            % Find kbind (concentration) to meet targeted expression
            myfun = @(x) getExp(makeRateMatrixCooperativity(exp(x),ku,ceq*ones(1,N)),Ip)-E0;
            [x,~,exitflag] = fzero(myfun,[B1,B2],options);
            
            kb = exp(x);
            KEQ(i,j) = kb;
            
            % Make state rate matrix
            Meq = makeRateMatrixCooperativity(kb,ku,ceq*ones(1,N));
            
            % Compute noise
            PEQ(i,j) = getNoise(Meq,Ip,tp);
            
            % Compute Hill coef
            ekb = myeps*kb;
            kb1 = kb-ekb;
            kb2 = kb+ekb;
            
            E1 = getExp(makeRateMatrixCooperativity(kb1,ku,ceq*ones(1,N)),Ip);
            E2 = getExp(makeRateMatrixCooperativity(kb2,ku,ceq*ones(1,N)),Ip);
            
            dy = E2 - E1;
            HEQ(i,j) = inveps*2*dy/(Emax-Emin);
        end
    end
end

%% Make reaction network graph (make Fig. 2A right)
clc
close all

% Select intersting models (parameters)
% Reasonable/plausible active residence time
Ig = T > 1e4-3e2 & T < 1e4+3e2;
[ki,kj] = find(Ig);
% High Hill coefficient
hh = HILL(Ig);
[~,r] = max(hh);

% Resulting parameters
co = [A1(ki(r),kj(r)),A2(ki(r),kj(r)),A2(ki(r),kj(r))];
kb = KB(ki(r),kj(r));

% Build state rate matrix
[M,X,V,Ed] = makeRateMatrixCooperativity(kb,ku,co);
m = size(M,1);
M(1:(m+1):end) = 0;

%%% Fig. 2A right
H1=figure(1);
set(H1,'position',[50 50 Wi 0.9*Wi],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);

cedges = lines(4);
G = digraph(M');
LWidths = 5*abs( (log10(G.Edges.Weight)-log10(min(G.Edges.Weight))) ) / (log10(max(G.Edges.Weight))-log10(min(G.Edges.Weight)));
LWidths(LWidths==0) = 1e-2;
EColor = zeros(size(G.Edges,1),3);

vv = G.Edges.EndNodes;
for i=1:length(vv)
    di = Ed(vv(i,2),vv(i,1));
    if di==1
        EColor(i,:) = cedges(2,:);
    elseif di==-1
        EColor(i,:) = cedges(1,:);
    elseif di==-2
        EColor(i,:) = cedges(3,:);
    elseif di==-3
        EColor(i,:) = cedges(3,:);
    end
end
if ~isempty(LWidths)
    xx = [1,0,1,0,2,1,2,1];
    yy = [3,2,2,1,2,1,1,0];
    pp=plot(h1,G,'XData',xx,'YData',yy,'LineWidth',LWidths,'EdgeColor',EColor,'ArrowSize',12,'MarkerSize',10,...
        'NodeColor','k','NodeFontSize',20,'EdgeFontSize',20,'EdgeFontAngle','normal');
end
set(h1,'Xtick',[],'Ytick',[],'linewidth',1.5,'fontsize',24)
set(h1,'XColor','none','YColor','none')

for i=1:m
    nodelabel{i} = arrayfun(@(x) num2str(x),X(i,:));
end
pp.NodeLabel = nodelabel;
edgelabel = cell(length(vv),1);

%% Plot induction curves for models leading to same active residence time at expression E0 (make Fig. 2B & C right)
clc
close all

% Select intersting models (parameters)
% Reasonable/plausible constant active residence time
Ig = T > 1e4-3e2 & T < 1e4+3e2;
[ki,kj] = find(Ig);
Nr = length(ki);

co = zeros(Nr,3);
kb = zeros(Nr,1);
hh = zeros(Nr,1);
ph = zeros(Nr,1);
tt = zeros(Nr,1);
co_eq = zeros(Nr,3);
kb_eq = zeros(Nr,1);
hh_eq = zeros(Nr,1);
ph_eq = zeros(Nr,1);

for r=1:Nr
    co(r,:) = [A1(ki(r),kj(r)),A2(ki(r),kj(r)),A2(ki(r),kj(r))];
    kb(r) = KB(ki(r),kj(r));
    hh(r) = HILL(ki(r),kj(r));
    ph(r) = PHI(ki(r),kj(r));
    
    tt(r) = T(ki(r),kj(r));
    
    co_eq(r,:) = CEQ(ki(r),kj(r))*ones(1,N);
    kb_eq(r) = KEQ(ki(r),kj(r));
    hh_eq(r) = HEQ(ki(r),kj(r));
    ph_eq(r) = PEQ(ki(r),kj(r));
end

% Pick models that meet the trageted Hill coefficient
htarget = [1:0.2:3,max(hh)];
Kr = zeros(size(htarget));
for i=1:length(htarget)
    [~,Kr(i)] = min(abs(hh-htarget(i)));
end

% Compute the induction curves for these models
Nr = length(Kr);
kbe = logspace(log10(median(kb))-2.5,log10(median(kb))+2.5,Ns);
E_eq = zeros(Nr,Ns);
E_neq = zeros(Nr,Ns);
for r=1:Nr
    rr = Kr(r);
    for i=1:Ns
        Meq = makeRateMatrixCooperativity(kbe(i),ku,co_eq(rr,:));
        Mneq = makeRateMatrixCooperativity(kbe(i),ku,co(rr,:));
        
        E_eq(r,i) = getExp(Meq,Ip);
        E_neq(r,i) = getExp(Mneq,Ip);
    end
end

a12 = co(:,2)./co(:,1);
la12 = log10(a12);

%%% Induction curves (Fig. 2B right)
H1=figure(1);
set(H1,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);
hold(h1,'on')

plot(h1,[min(kbe),max(kbe)],[0.5,0.5],'--k','linewidth',2)
for r=1:Nr
    cr=1+round(ncmap*(la12(Kr(r)))/max(la12(:)));
    plot(h1,kbe,E_neq(r,:),'-','color',cmap(cr,:),'linewidth',2)
end
plot(h1,kbe,E_eq(1,:),'-k','linewidth',2)
plot(h1,kb_eq(Kr(end)),0.5,'ok','markersize',10,'linewidth',2)
plot(h1,kb(Kr(end)),0.5,'*k','markersize',12,'linewidth',2)
set(h1,'fontsize',22,'linewidth',2,'xscale','log','xtick',[1e-5,1e-2],'ytick',[0,1],'tickdir','out')
xlabel(h1,'Concentration k_{+}')
ylabel(h1,'Expression E')
xlim(h1,[min(kbe),max(kbe)])

%%% Transformed induction curves
H2=figure(2);
set(H2,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h2 = axes('parent',H2);
hold(h2,'on')

plot(h2,[min(kbe),max(kbe)],[1,1],'--k','linewidth',2)
for r=1:Nr
    cr=1+round(ncmap*(la12(Kr(r)))/max(la12(:)));
    plot(h2,kbe,E_neq(r,:)./(1-E_neq(r,:)),'-','color',cmap(cr,:),'linewidth',2)
end
plot(h2,kbe,E_eq(1,:)./(1-E_eq(1,:)),'-k','linewidth',2)
set(h2,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','xtick',[1e-5,1e-2],'tickdir','out')
xlabel(h2,'Concentration k_{+}')
ylabel(h2,'E/(1-E)')
xlim(h2,[min(kbe),max(kbe)])
ylim(h2,[1e-3,1e3])

%%% Hill coefficient versus cooperativity ratio (Fig. 2C right)
H3=figure(3);
set(H3,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h3 = axes('parent',H3);
hold(h3,'on')

[~,Is] = sort(co(:,2)./co(:,1));
a12 = a12(Is);
la12 = la12(Is);
hhs = hh(Is);
for i=1:length(hh)
    cr=1+round(ncmap*(la12(i))/max(la12(:)));
    plot(h3,a12(i:end),hhs(i:end),'-','color',cmap(cr,:),'linewidth',3)
end
plot(h3,1,hhs(1),'ok','markersize',10,'linewidth',2)
plot(h3,co(Kr(end),2)./co(Kr(end),1),hh(Kr(end)),'*k','markersize',12,'linewidth',2)
set(h3,'fontsize',22,'linewidth',2,'xscale','log','ytick',1:3,'xtick',[1,1e3,1e6,1e9],'tickdir','out')
xlabel(h3,'Cooperativity ratio a_2/a_1')
ylabel(h3,'Hill coefficient H')
xlim(h3,[1,1e9])

%%% Concentration at E0 versus cooperativity ratio
H4=figure(4);
set(H4,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h4 = axes('parent',H4);
hold(h4,'on')

kbs = kb(Is);
for i=1:length(hh)
    cr=1+round(ncmap*(la12(i))/max(la12(:)));
    plot(h4,a12(i:end),kbs(i:end),'-','color',cmap(cr,:),'linewidth',3)
end
plot(h4,1,kbs(1),'ok','markersize',10,'linewidth',2)
plot(h4,co(Kr(end),2)./co(Kr(end),1),kb(Kr(end)),'*k','markersize',12,'linewidth',2)
set(h4,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','xtick',[1,1e3,1e6,1e9],'tickdir','out')
xlabel(h4,'Cooperativity ratio a_2/a_1')
ylabel(h4,'Concentration k_+^{1/2}')

%% Plot parameter space & regulatory phenotypes
clc
close all

% Select intersting models (parameters) for display within phase space
% Reasonable/plausible constant active residence time
Ig = T > 1e4-3e2 & T < 1e4+3e2;

c1max = 10^(1*amax);
c2max = 10^(N*(amax+1));

%%% Hill coefficient as function of cooperativities
H1=figure(1);
set(H1,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);
hold(h1,'on')

[~,Is] = sort(HILL(:),'ascend');
scatter(h1,A1(Is),A2(Is)./A1(Is),[],HILL(Is),'filled')
plot(h1,A1(Ig),A2(Ig)./A1(Ig),'ok')
set(h1,'fontsize',22,'linewidth',2,'xscale','log','yscale','log',...
    'Layer','top','tickdir','out')
xlabel(h1,'a_1')
ylabel(h1,'a_2/a_1')
xlim(h1,[1,c1max])
ylim(h1,[1,c2max])
cb=colorbar(h1,'LineWidth',2);
set(get(cb,'ylabel'),'String','Hill coefficient H','fontsize',20);
colormap(h1,cmap)

%%% Entropy production as function of cooperativities
H2=figure(2);
set(H2,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h2 = axes('parent',H2);
hold(h2,'on')

[~,Is] = sort(S(:),'ascend');
scatter(h2,A1(Is),A2(Is)./A1(Is),[],S(Is),'filled')
plot(h2,A1(Ig),A2(Ig)./A1(Ig),'ok')
set(h2,'fontsize',22,'linewidth',2,'xscale','log','yscale','log',...
    'colorscale','log','Layer','top','tickdir','out')
xlabel(h2,'a_1')
ylabel(h2,'a_2/a_1')
xlim(h2,[1,c1max])
ylim(h2,[1,c2max])
cb=colorbar(h2,'LineWidth',2);
set(get(cb,'ylabel'),'String','Entropy prod S','fontsize',20);
colormap(h2,cmap)

%%% Active residence time as function of cooperativities
H3=figure(3);
set(H3,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h3 = axes('parent',H3);
hold(h3,'on')

[~,Is] = sort(T(:),'ascend');
scatter(h3,A1(Is),A2(Is)./A1(Is),[],T(Is),'filled')
plot(h3,A1(Ig),A2(Ig)./A1(Ig),'ok')
set(h3,'fontsize',22,'linewidth',2,'xscale','log','yscale','log',...
    'colorscale','log','Layer','top','tickdir','out')
xlabel(h3,'a_1')
ylabel(h3,'a_2/a_1')
xlim(h3,[1,c1max])
ylim(h3,[1,c2max])
cb=colorbar(h3,'LineWidth',2);
set(get(cb,'ylabel'),'String','Residence time T_A','fontsize',20);
colormap(h3,cmap)

%%% Noise as function of cooperativities
H4=figure(4);
set(H4,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h4 = axes('parent',H4);
hold(h4,'on')

[~,Is] = sort(PHI(:),'ascend');
scatter(h4,A1(Is),A2(Is)./A1(Is),[],PHI(Is),'filled')
plot(h4,A1(Ig),A2(Ig)./A1(Ig),'ok')
set(h4,'fontsize',22,'linewidth',2,'xscale','log','yscale','log',...
    'colorscale','log','Layer','top','tickdir','out')
xlabel(h4,'a_1')
ylabel(h4,'a_2/a_1')
xlim(h4,[1,c1max])
ylim(h4,[1,c2max])
cb=colorbar(h4,'LineWidth',2);
set(get(cb,'ylabel'),'String','Propagated noise \Phi','fontsize',20);
colormap(h4,cmap)


%% Plot regulatory phase space (make Fig. 2D right)
clc
close all

% Reasonable/plausible constant active residence time
Ig = T > 1e4-3e2 & T < 1e4+3e2;
% Among which, highlight one non-eq model with highest Hill coefficient
[ki,kj] = find(Ig);
hh = HILL(Ig);
[~,r] = max(hh);

co = [A1(ki(r),kj(r)),A2(ki(r),kj(r)),A2(ki(r),kj(r))];
kb = KB(ki(r),kj(r));
hh = HILL(ki(r),kj(r));
ph = PHI(ki(r),kj(r));
ss = S(ki(r),kj(r));
tt = T(ki(r),kj(r));
co_eq = CEQ(ki(r),kj(r))*ones(1,N);
kb_eq = KEQ(ki(r),kj(r));
hh_eq = HEQ(ki(r),kj(r));
ph_eq = PEQ(ki(r),kj(r));

% Compute induction curves for the non-eq model above 
% and its equilibrium counterpart
kbe = logspace(log10(kb)-2,log10(kb)+3,Ns);
E_eq = zeros(Ns,1);
E_neq = zeros(Ns,1);
for i=1:Ns
    E_eq(i) = getExp(makeRateMatrixCooperativity(kbe(i),ku,co_eq),Ip);
	E_neq(i) = getExp(makeRateMatrixCooperativity(kbe(i),ku,co),Ip);
end

%%% Induction curves
H1=figure(1);
set(H1,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h1 = axes('parent',H1);
hold(h1,'on')

plot(h1,[min(kbe),max(kbe)],[0.5,0.5],'--k','linewidth',2)
plot(h1,kbe,E_eq,'-k','linewidth',2)
plot(h1,kbe,E_neq,'-r','linewidth',2)
plot(h1,kb_eq,0.5,'ok','linewidth',2,'markersize',10)
plot(h1,kb,0.5,'*r','linewidth',2,'markersize',12)
set(h1,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h1,'k_{+}')
ylabel(h1,'Expression E')
xlim(h1,[min(kbe),max(kbe)])

%%% Transformed induction curves
H2=figure(2);
set(H2,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h2 = axes('parent',H2);
hold(h2,'on')

plot(h2,[min(kbe),max(kbe)],[1,1],'--k','linewidth',2)
plot(h2,kbe,E_eq./(1-E_eq),'-k','linewidth',2)
plot(h2,kbe,E_neq./(1-E_neq),'-r','linewidth',2)
set(h2,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h2,'k_{+}')
ylabel(h2,'E/(1-E)')
xlim(h2,[min(kbe),max(kbe)])

%%% Hill coefficient as log derivative
H3=figure(3);
set(H3,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h3 = axes('parent',H3);
hold(h3,'on')

Zeq = diff(log10(E_eq./(1-E_eq)))./diff(log10(kbe));
Zneq = diff(log10(E_neq./(1-E_neq)))./diff(log10(kbe));
keq = find(kbe>kb_eq,1,'first');
kneq = find(kbe>kb,1,'first');

plot(h3,kbe(1:(end-1)),Zeq,'-k','linewidth',2)
plot(h3,kbe(1:(end-1)),Zneq,'-r','linewidth',2)
plot(h3,kb_eq,Zeq(keq),'ok','linewidth',2,'markersize',10)
plot(h3,kb,Zneq(kneq),'*r','linewidth',2,'markersize',12)
set(h3,'fontsize',22,'linewidth',2,'xscale','log','tickdir','out')
xlabel(h3,'k_{+}')
ylabel(h3,'d/dk_{+} log(E/(1-E))')
xlim(h3,[min(kbe),max(kbe)])

%%% Regulatory phase space
% Highlight models whose Hill coefficient exceed eq limit (H=3 for N=3)
Ib = HILL > 3 & T < 1e5;

%%% Hill coefficient vs active residence time
H4=figure(4);
set(H4,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h4 = axes('parent',H4);
hold(h4,'on')

plot(h4,T(:),HILL(:),'o','MarkerFaceColor',lines(1))
plot(h4,T(Ib),HILL(Ib),'o','color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(h4,T(:,1),HILL(:,1),'-k','linewidth',3)
plot(h4,[min(T(:)),1e5],N*[1,1],'--k','linewidth',2)
plot(h4,tt,hh_eq,'ok','linewidth',2,'markersize',10)
plot(h4,tt,hh,'*r','linewidth',2,'markersize',12)
set(h4,'fontsize',22,'linewidth',2,'xscale','log','Layer','top','ytick',1:3,'tickdir','out')
xlabel(h4,'Residence time T_{A} [1/k_{-}]')
ylabel(h4,'Hill coefficient H')
xlim(h4,[min(T(:)),1e5])

%%% Noise vs active residence time
H5=figure(5);
set(H5,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h5 = axes('parent',H5);
hold(h5,'on')

plot(h5,T(:),PHI(:),'o','MarkerFaceColor',lines(1))
plot(h5,T(Ib),PHI(Ib),'o','color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(h5,T(:,1),PHI(:,1),'-k','linewidth',3)
plot(h5,tt,ph_eq,'ok','linewidth',2,'markersize',10)
plot(h5,tt,ph,'*r','linewidth',2,'markersize',12)
set(h5,'fontsize',22,'linewidth',2,'xscale','log','ytick',0:1,'tickdir','out')
xlabel(h5,'Residence time T_{A} [1/k_{-}]')
ylabel(h5,'Propagated noise \Phi')
xlim(h5,[min(T(:)),1e5])

%%% Entropy production vs active residence time
H6=figure(6);
set(H6,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h6 = axes('parent',H6);
hold(h6,'on')

plot(h6,T(:),S(:),'o','MarkerFaceColor',lines(1))
plot(h6,T(Ib),S(Ib),'o','color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(h6,tt,ss,'*r','linewidth',2,'markersize',12)
set(h6,'fontsize',22,'linewidth',2,'xscale','log','yscale','log','tickdir','out')
xlabel(h6,'Residence time T_{A} [1/k_{-}]')
ylabel(h6,'Entropy Production S')
xlim(h6,[min(T(:)),1e5])
ylim(h6,[1e-8,5])

%%% Hill coefficient vs active residence time with coopertivity ratio
%%% colorcoded (Fig. 2D right)
H7=figure(7);
set(H7,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h7 = axes('parent',H7);
hold(h7,'on')

A21 = A2./A1;
[~,Is] = sort(A21(:),'ascend');
scatter(h7,T(Is),HILL(Is),[],A21(Is),'filled')
plot(h7,T(:,1),HILL(:,1),'-k','linewidth',3)
plot(h7,[min(T(:)),1e5],N*[1,1],'--k','linewidth',2)
plot(h7,tt,hh_eq,'ok','linewidth',2,'markersize',10)
plot(h7,tt,hh,'*r','linewidth',2,'markersize',12)
set(h7,'fontsize',22,'linewidth',2,'xscale','log','Layer','top','colorscale','log','ytick',1:3,'tickdir','out')
xlabel(h7,'Residence time T_{A} [1/k_{-}]')
ylabel(h7,'Hill coefficient H')
xlim(h7,[min(T(:)),1e5])
caxis(h7,[1,c2max])
cb=colorbar(h7,'LineWidth',2);
set(get(cb,'ylabel'),'String','a_2/a_1','fontsize',20);
colormap(h7,cmap)

%%% Hill coefficient vs active residence time with noise colorcoded
H8=figure(8);
set(H8,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h8 = axes('parent',H8);
hold(h8,'on')

[~,Is] = sort(PHI(:),'ascend');
scatter(h8,T(Is),HILL(Is),[],PHI(Is),'filled')
plot(h8,T(:,1),HILL(:,1),'-k','linewidth',3)
plot(h8,[min(T(:)),1e5],N*[1,1],'--k','linewidth',2)
plot(h8,tt,hh_eq,'ok','linewidth',2,'markersize',10)
plot(h8,tt,hh,'*r','linewidth',2,'markersize',12)
set(h8,'fontsize',22,'linewidth',2,'xscale','log','Layer','top','ytick',1:3,'tickdir','out')
xlabel(h8,'Residence time T_{A} [1/k_{-}]')
ylabel(h8,'Hill coefficient H')
xlim(h8,[min(T(:)),1e5])
cb=colorbar(h8,'LineWidth',2);
set(get(cb,'ylabel'),'String','Propagated noise \Phi','fontsize',20);
colormap(h8,cmap)

%%% Hill coefficient vs active residence time with entropy production colorcoded
H9=figure(9);
set(H9,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
h9 = axes('parent',H9);
hold(h9,'on')

[~,Is] = sort(S(:),'ascend');
scatter(h9,T(Is),HILL(Is),[],S(Is),'filled')
plot(h9,T(:,1),HILL(:,1),'-k','linewidth',3)
plot(h9,[min(T(:)),1e5],N*[1,1],'--k','linewidth',2)
plot(h9,tt,hh_eq,'ok','linewidth',2,'markersize',10)
plot(h9,tt,hh,'*r','linewidth',2,'markersize',12)
set(h9,'fontsize',22,'linewidth',2,'xscale','log','Layer','top','colorscale','log','ytick',1:3,'tickdir','out')
xlabel(h9,'Residence time T_{A} [1/k_{-}]')
ylabel(h9,'Hill coefficient H')
xlim(h9,[min(T(:)),1e5])
cb=colorbar(h9,'LineWidth',2);
set(get(cb,'ylabel'),'String','Entropy prod S','fontsize',20);
colormap(h9,cmap)
