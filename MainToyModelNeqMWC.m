%% Main script to generate Fig. 3
%%% Compute phenotypes of the non-equilibrium MWC model
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

% Number of binding sites
N = 3;

% Targeted expression level
E0 = 0.5;

%%% Model parameters
% Lifetime of proteins
tp = 3.6e6;

% TF unbinding rate, related to affinity
% for specific (strong) binding site and non-specific binding site
kb = 1; % will be chosen to reach E0
ku = 0.01;
ku_ns = 1;

% Mediator binding and unbinding rate
kb_M = 1e-2;
ku_M = 1e4;

% Un-linking rate
klink = 10; % placeholder, will be sampled
kunlink = 0;

% Cooperativity parameter
alpha = 4e3; % placeholder, will be sampled

%%% States
[~,configBound,configTF] = makeRateMatrixNeqMWC(kb,ku,kb_M,ku_M,klink,kunlink,alpha,N);

% States corresponding to first biding site occupied by TF used to compute
% TF residence time.
Itf1 = configTF(:,1) == 1;
% States corresponding to Mediator bound used to compute expression and
% Mediator residence time.
Im = configBound(:,3) == 1;

%%% Constants for sampling
% Some useful theoretical values
eL = kb_M/ku_M;
inva = (eL*(1/E0-1))^(1/N);
Emin = eL/(1+eL);
Smax = E0*ku_M/kb_M;

% Some numerical constants necessary for search & optimization function
options = optimset('TolX',1e-16);
myeps = 1e-3;
inveps = 1/myeps;
B1 = -20*log(10);
B2 = 10*log(10);
A1 = -log(inva)+1e-12;
A2 = 10*log(10);


%% Compute phenotypes at fixed expression E0

Np = 200;

% Sampling klink and alpha
KLINK = logspace(-6,8,Np);
ALPHA = logspace(-log10(inva),8,Np);

% Storing phenotypes
KBIND = nan(Np,Np);
TTF = nan(Np,Np);
PHI = nan(Np,Np);
SPEC = nan(Np,Np);
EMAX = nan(Np,Np);
HILL = nan(Np,Np);

KBIND_EQ = nan(Np,Np); %TF concentration EQ
PHI_EQ = nan(Np,Np);
SPEC_EQ = nan(Np,Np);
HILL_EQ = nan(Np,Np);

% Main loop to sample regulatory phenotypes
% Parallelized loop
parfor i=1:Np
% Standard loop in case parfor cannot be use or for debugging 
%for i=1:Np
    disp(i)
    klink = KLINK(i);
    for j=1:Np
        alpha = ALPHA(j);
        
        myfun = @(x) getExp(makeRateMatrixNeqMWC(exp(x),ku,kb_M,ku_M,klink,kunlink,alpha,N),Im)-E0;
        
        if myfun(B1)*myfun(B2) < 0
            % Find kbind (TF concentration) to meet targeted expression E0
            [x,~,exitflag] = fzero(myfun,[B1,B2],options);
            if exitflag ~= 1
                exitflag
            end
            
            kb = exp(x);
            KBIND(i,j) = kb;
            
            % Make state rate matrix
            Mneq = makeRateMatrixNeqMWC(kb,ku,kb_M,ku_M,klink,kunlink,alpha,N);
            
            % Compute expression
            [ES,P] = getExp(Mneq,Im);
            
            % Compute residence time
            T = getResid(Mneq,Itf1,P);
            TTF(i,j) = T;
            
            % Compute noise
            PHI(i,j) = getNoise(Mneq,Im,tp);
            
            % Compute specificity
            Mneq_ns = makeRateMatrixNeqMWC(kb,ku_ns,kb_M,ku_M,klink,kunlink,alpha,N);
            EU = getExp(Mneq_ns,Im);
            SPEC(i,j) = ES/EU;
                                  
            % Compute Hill coeff
            % We must find the half maximal expression level EH, which is
            % usually different than E0 for this class of model.
            % First compute maximal expression
            Emax = getExp(makeRateMatrixNeqMWC(exp(B2),ku,kb_M,ku_M,klink,kunlink,alpha,N),Im);
            EMAX(i,j) = Emax;
            % Compute EH based on Emax and Emin
            EH = 0.5*(Emax+Emin);
            % Find kb0 (concentration) to meet EH
            myfun2 = @(x) getExp(makeRateMatrixNeqMWC(exp(x),ku,kb_M,ku_M,klink,kunlink,alpha,N),Im)-EH;
            
            if myfun2(B1)*myfun2(B2) < 0 && Emax > 1.1*Emin
                [x,~,exitflag] = fzero(myfun2,[B1,B2],options);
                if exitflag ~= 1
                    exitflag
                end
                
                % Once kb0 is known, we can compute the Hill coef
                kb0 = exp(x);
                ekb = myeps*kb0;
                kb1 = kb0-ekb;
                kb2 = kb0+ekb;
                
                E1 = getExp(makeRateMatrixNeqMWC(kb1,ku,kb_M,ku_M,klink,kunlink,alpha,N),Im);
                E2 = getExp(makeRateMatrixNeqMWC(kb2,ku,kb_M,ku_M,klink,kunlink,alpha,N),Im);
                
                dy = E2 - E1;
                HILL(i,j) = inveps*2*dy/(Emax-Emin);                
            end
            
            %%% Build equivalent equ model achieving same expression 
            % level E0 and same TF residence time T.
            % Below we use analytical solutions for most of the phenotypes
            % since the equ-MWC model is tractable. However, one could
            % compute the phenotypes numericaly using the provided
            % functions.
            % For more details on the derivation of these solutions, see
            % Grah et al. 2020, DOI:10.1073/pnas.2006731117
            myfun3 = @(x) log(getResidMWC(exp(x),ku,eL,E0,N))-log(T);
            
            if myfun3(A1)*myfun3(A2) < 0
                % Find alpha to meet targeted expression & residence time
                [x,~,exitflag] = fzero(myfun3,[A1,A2],options);
                if exitflag ~= 1
                    exitflag
                end
                
                alpha_eq = exp(x);
                
                % Compute TF concentration
                x = (eL*(1/E0-1))^(1/N);
                eE = (x-1)./(1-x.*alpha_eq);
                kb_eq = ku*eE;
                KBIND_EQ(i,j) = kb_eq;
                
                % Compute Noise
                [Meq,configBound] = makeRateMatrixMWC(kb_eq,ku,kb_M,ku_M,alpha_eq,N);
                Im_eq = configBound(:,3) == 1;
                PHI_EQ(i,j) = getNoise(Meq,Im_eq,tp);
                
                % Compute specificity
                eNS = kb_eq/ku_ns;
                A =  ((1+eE)./(1+eE.*alpha_eq)).^N;
                B = ((1+eNS)./(1+eNS.*alpha_eq)).^N;
                SPEC_EQ(i,j) = (eL+B)./(eL+A);
                
                % Compute Hill coeff
                Emax_eq = alpha_eq.^N*eL./(alpha_eq.^N*eL+1);
                Emin_eq = eL./(eL+1);
                EH = 0.5*(Emax_eq+Emin_eq);
                
                x = (eL*(1./EH-1)).^(1/N);
                HILL_EQ(i,j) = 4*N./(Emax_eq-Emin_eq) .* (1-x).*(1-x.*alpha_eq)./(1-alpha_eq) .* eL.*x.^(N-1)./(eL+x.^N).^2;
            end
        end
    end
end

PhenoStruct.klink = KLINK;
PhenoStruct.alpha = ALPHA;
PhenoStruct.kbind = KBIND;
PhenoStruct.Ttf = TTF;
PhenoStruct.Phi = PHI;
PhenoStruct.Spec = SPEC;
PhenoStruct.Emax = EMAX;
PhenoStruct.Hill = HILL;

PhenoStruct.kbind_eq = KBIND_EQ;
PhenoStruct.Spec_eq = SPEC_EQ;
PhenoStruct.Phi_eq = PHI_EQ;
PhenoStruct.Hill_eq = HILL_EQ;

%%% save PhenoStruct
save(['PhenoStruct_N',num2str(N),'_E',strrep(num2str(0.5),'.',''),'.mat'],'PhenoStruct');

DS = SPEC./SPEC_EQ;
TT0 = TTF/T0;

%% Or directly load phenotypes if already computed
clc

load(['PhenoStruct_N',num2str(N),'_E',strrep(num2str(0.5),'.',''),'.mat'],'PhenoStruct');

KLINK = PhenoStruct.klink;
ALPHA = PhenoStruct.alpha;
KBIND = PhenoStruct.kbind;

TTF = PhenoStruct.Ttf;
PHI = PhenoStruct.Phi;
SPEC = PhenoStruct.Spec;
EMAX = PhenoStruct.Emax;
HILL = PhenoStruct.Hill;

KBIND_EQ = PhenoStruct.kbind_eq;
SPEC_EQ = PhenoStruct.Spec_eq;
PHI_EQ = PhenoStruct.Phi_eq;
HILL_EQ = PhenoStruct.Hill_eq;

DS = SPEC./SPEC_EQ;
TT0 = TTF/T0;

%% Make reaction network graph
clc
close all

% Pick a 'good' model 
UF = SPEC/Smax .* (1-PHI).^2 .* (HILL-min(HILL(:)))./(N-min(HILL(:))) .* (1./(5+TT0));
[~,I] = max(UF(:));
[kl,ka] = ind2sub(size(UF),I);
klink = KLINK(kl);
alpha = ALPHA(ka);
kb = KBIND(kl,ka);

% Build state rate matrix
[Mneq,configBound,configTF,configLink,Ed] = makeRateMatrixNeqMWC(kb,ku,kb_M,ku_M,klink,kunlink,alpha,N);

Hi=figure(1);
set(Hi,'position',[50 50 700 500],'paperpositionmode','auto','color','w');
hi = axes('parent',Hi);
hold(hi,'on')

m = size(Mneq,1);
Mg = Mneq;
Mg(1:(m+1):end) = 0;

cedges = lines(6);
G = digraph(Mg');
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
    elseif di==-0.5
        EColor(i,:) = cedges(4,:);
    elseif di==0.5
        EColor(i,:) = cedges(3,:);
    elseif di==-1.5
        EColor(i,:) = cedges(6,:);
    elseif di==1.5
        EColor(i,:) = cedges(5,:);
    end
end
if ~isempty(LWidths)
    if N==1
        xx = [1,1,2,0,1];
        yy = [0,1,2,2,3];
        pp=plot(hi,G,'XData',xx,'YData',yy,'LineWidth',LWidths,'EdgeColor',EColor,'ArrowSize',12,'MarkerSize',10,...
            'NodeColor','k','NodeFontSize',20,'EdgeFontSize',20,'EdgeFontAngle','normal');
    else     
        pp=plot(hi,G,'LineWidth',LWidths,'EdgeColor',EColor,'ArrowSize',12,'MarkerSize',10,...
            'NodeColor','k','NodeFontSize',20,'EdgeFontSize',20,'EdgeFontAngle','normal');
    end
end
set(hi,'Xtick',[],'Ytick',[],'linewidth',1.5,'fontsize',24)
set(hi,'XColor','none','YColor','none')

for i=1:m
    nodelabel{i} = [num2str(configBound(i,end)),';',...
        arrayfun(@(x) num2str(x),configTF(i,:)),';',...
        arrayfun(@(x) num2str(x),configLink(i,1:end))];
end
pp.NodeLabel = nodelabel;

%% Generate stochastic trajectories (make Fig. 3B)
clc
close all

% Constants for raster trajectories
dt = 1e2;
Tend = 3e7;
Tend2 = Tend+2e7;
T = 0:dt:Tend;

% Protein production rate kr and lifetime of proteins tau
kr = 1e3/3.6e5; %1/h
tau = 3.6e6; %10h

%%% Non-equ model raster
% Pick a 'good' model 
UF = SPEC/Smax .* (1-PHI).^2 .* (HILL-min(HILL(:)))./(N-min(HILL(:))) .* (1./(5+TT0));
[~,I] = max(UF(:));
[kl,ka] = ind2sub(size(UF),I);
klink = KLINK(kl);
alpha = ALPHA(ka);
kb = KBIND(kl,ka);

[Mneq,configBound] = makeRateMatrixNeqMWC(kb,ku,kb_M,ku_M,klink,kunlink,alpha,N);
Im = configBound(:,3) == 1;

% Generate raster
Sneq = genRastTrajectories(Mneq,Im,kr,tau,dt,Tend2);
Sneq = Sneq((2e7/dt+1):end,:);

%%% Equ model raster
% Compute equ parameters leading to same expression E and TF residence time 
kb_eq = KBIND_EQ(kl,ka);
x = (eL*(1/E0-1))^(1/N);
eE = kb_eq/ku;
alpha_eq = (1-((x-1)/eE))/x;

[Meq,configBound] = makeRateMatrixMWC(kb_eq,ku,kb_M,ku_M,alpha_eq,N);
Im_eq = configBound(:,3) == 1;

% Generate raster
Seq = genRastTrajectories(Meq,Im_eq,kr,tau,dt,Tend2);
Seq = Seq((2e7/dt+1):end,:);

% Constants for plot
Ts = 1e7;
myc = lines(2);

%%% Mediator time trace (Fig. 3B top)
Hi=figure(1);
set(Hi,'position',[50 50 Wi Le/2],'paperpositionmode','auto','color','w');
hi = axes('parent',Hi);
hold(hi,'on')

x = T/Ts;
y = Sneq(:,2);
stairs(hi,x(1:1e3:end),y(1:1e3:end),'linewidth',2,'color',myc(2,:));
set(hi,'linewidth',2,'fontsize',22,'tickdir','out','ytick',[0,1],'xtick',[0,Tend/Ts])
xlabel('Time [1/k_{-}]')
ylabel('Mediator')
ylim([0,1.1])
xlim(hi,[T(1),T(end)]/Ts)

%%% Protein time trace (Fig. 3B bottom)
Hi=figure(2);
set(Hi,'position',[50 50 Wi Le/2],'paperpositionmode','auto','color','w');
hi = axes('parent',Hi);
hold(hi,'on')

Ps = 1e4;
y0 = Seq(:,3)/Ps;
y = Sneq(:,3)/Ps;
plot(hi,x,y0,'linewidth',2,'color',myc(1,:));
plot(hi,x,y,'linewidth',2,'color',myc(2,:));
set(hi,'linewidth',2,'fontsize',22,'tickdir','out','ytick',[0,1],'xtick',[0,Tend/Ts])
xlabel('Time [1/k_{-}]')
ylabel('Protein #')
xlim(hi,[T(1),T(end)]/Ts)
ylim(hi,[0,1])

%% Plot accessible regulatory phenotypes (make Fig. 3C)
clc
close all

% Regulatory phenotypes of a 'good' model
UF = SPEC/Smax .* (1-PHI).^2 .* (HILL-min(HILL(:)))./(N-min(HILL(:))) .* (1./(5+TT0));
[~,I] = max(UF(:));
[kl,ka] = ind2sub(size(UF),I);

%%% Hill coefficient vs Residence time, specificity colorcoded
Hi=figure(1);
set(Hi,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
hi = axes('parent',Hi);
hold(hi,'on')

[~,Is] = sort(SPEC(:),'ascend');
scatter(hi,TT0(Is),HILL(Is),3,SPEC(Is)/Smax,'filled')
[~,Is] = sort(TT0(:),'ascend');
plot(hi,TT0(Is),HILL_EQ(Is),'-k','linewidth',3);

xlabel(hi,'Residence time T_{TF} [1/k_{-}]')
ylabel(hi,'Hill coefficient H')
set(hi,'linewidth',2,'fontsize',22,'tickdir','out','xscale','log','ytick',1:3)
cb=colorbar(hi,'LineWidth',2,'xtick',0:1);
set(get(cb,'ylabel'),'String','Specificity S/S_{max}','fontsize',22);
xlim(hi,[min(TT0(:)),1e3])
ylim(hi,[1,N])
caxis(hi,[0,1])
colormap(hi,cmap)

% highlight 'good' model
plot(hi,TT0(kl,ka),HILL(kl,ka),'*k','markersize',12,'linewidth',2)

% Make a grid and interpolation to represent phenotypes as a function of
% specificity and TF residence time
nn = 500;
xx = linspace(log10(min(TT0(:))),3,nn);
yy = linspace(0,1,nn);
[xq,yq] = meshgrid(xx,yy);

In = ~isnan(DS);
vqpp = griddata(log10(TT0(In)),SPEC(In)/Smax,PHI(In),xq,yq);
vqhh = griddata(log10(TT0(In)),SPEC(In)/Smax,HILL(In),xq,yq);

% Remove griddata triangulation artefact
[~,Is] = sort(SPEC_EQ(:),'ascend');
x = log10(TT0(Is));
y = SPEC_EQ(Is)/Smax;
tt = interp1(x(~isnan(y)),y(~isnan(y)),xx);
Ia = false(nn,nn);
for i=1:nn
    Ii = yq(:,1) < tt(i);
    Ia(Ii,i) = true;
end
vqpp(Ia) = nan;
vqhh(Ia) = nan;

%%% Specificity vs Residence time, noise colorcoded (Fig. 3C left)
Hi=figure(2);
set(Hi,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
hi = axes('parent',Hi);
hold(hi,'on')

imagesc(hi,10.^xq(1,:),yq(:,1),vqpp,'AlphaData',~isnan(vqpp))
[~,Is] = sort(TT0(:),'ascend');
plot(hi,TT0(Is),SPEC_EQ(Is)/Smax,'-k','linewidth',3);

xlabel(hi,'Residence time T_{TF} [1/k_{-}]')
ylabel(hi,'Specificity S/S_{max}')
set(hi,'linewidth',2,'fontsize',22,'tickdir','out','xscale','log','ytick',0:1)
cb=colorbar(hi,'LineWidth',2,'xtick',0:1);
set(get(cb,'ylabel'),'String','Propagated noise \Phi','fontsize',22);
xlim(hi,[min(TT0(:)),1e3])
ylim(hi,[0,1])
caxis(hi,[0,1])
colormap(hi,cmap)

% highlight 'good' model
plot(hi,TT0(kl,ka),SPEC(kl,ka)/Smax,'*k','markersize',12,'linewidth',2)

%%% Specificity vs Residence time, Hill coef colorcoded (Fig. 3C right)
Hi=figure(3);
set(Hi,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
hi = axes('parent',Hi);
hold(hi,'on')

imagesc(hi,10.^xq(1,:),yq(:,1),vqhh,'AlphaData',~isnan(vqhh))
[~,Is] = sort(TT0(:),'ascend');
plot(hi,TT0(Is),SPEC_EQ(Is)/Smax,'-k','linewidth',3);

xlabel(hi,'Residence time T_{TF} [1/k_{-}]')
ylabel(hi,'Specificity S/S_{max}')
set(hi,'linewidth',2,'fontsize',22,'tickdir','out','xscale','log','ytick',0:1)
cb=colorbar(hi,'LineWidth',2,'xtick',1:3);
set(get(cb,'ylabel'),'String','Hill coefficient H','fontsize',22);
xlim(hi,[min(TT0(:)),1e3])
ylim(hi,[0,1])
caxis(hi,[min(HILL(:)),N])
colormap(hi,cmap)

% highlight 'good' model
plot(hi,TT0(kl,ka),SPEC(kl,ka)/Smax,'*k','markersize',12,'linewidth',2)

%% Plot pheno (make Fig. 3D)
close all
clc

% Define an utility function UF
% Here the functional form is artificial and only meant as an example
% The given function favors high specificity S, low propageted noise Phi,
% high Hill coefficient and low TF residence time.
UF = SPEC/Smax .* (1-PHI).^2 .* (HILL-min(HILL(:)))./(N-min(HILL(:))) .* (1./(5+TT0));
UF = UF/max(UF(:));

%%% Parameter space and Utility function (Fig. 3D)
Hi=figure(1);
set(Hi,'position',[50 50 Wi Le],'paperpositionmode','auto','color','w');
hi = axes('parent',Hi);
hold(hi,'on')

imagesc(hi,KLINK,ALPHA,UF','AlphaData',~isnan(UF'));
xlabel(hi,'Linking rate k_{link}')
ylabel(hi,'Cooperativity \alpha')
set(hi,'linewidth',2,'fontsize',22,'xscale','log','yscale','log','tickdir','out','xtick',[1e-4,1e1,1e6])
xlim(hi,[1e-4,1e6])
ylim(hi,[1e2,1e6])
colormap(hi,cmap)
caxis(hi,[0,1])
cb=colorbar(hi,'LineWidth',2,'xtick',0:1);
set(get(cb,'ylabel'),'String','Utility F(S,\Phi,H,T_{TF})','fontsize',22);

% Select a 'good' model with respect to UF and higlight its location in
% the parameter space
[~,I] = max(UF(:));
[kl,ka] = ind2sub(size(UF),I);
plot(hi,KLINK(kl),ALPHA(ka),'*k','markersize',12,'linewidth',2)

disp(['TTF=',num2str(TTF(kl,ka)/T0)])
disp(['H=',num2str(HILL(kl,ka))])
disp(['Phi=',num2str(PHI(kl,ka))])
disp(['S=',num2str(SPEC(kl,ka)/Smax)])
disp(['DS=',num2str(DS(kl,ka))])

% Plot contours
% Noise and specificity (left plot)
lPHI = 0.3;
contour(hi,KLINK,ALPHA,PHI',lPHI*[1,1],'color',[0.5,0.5,0.5],'linewidth',3);
lSS = 0.25;
contour(hi,KLINK,ALPHA,SPEC'/Smax,lSS*[1,1],'color',[0.5,0.5,0.5],'linewidth',3);
% Residence time and Hill coefficient (right plot)
lTTF = log10(20);
contour(hi,KLINK,ALPHA,log10(TT0'),lTTF*[1,1],'color',[0.5,0.5,0.5],'linewidth',3);
lHILL = 2;
contour(hi,KLINK,ALPHA,HILL',lHILL*[1,1],'color',[0.5,0.5,0.5],'linewidth',3);

