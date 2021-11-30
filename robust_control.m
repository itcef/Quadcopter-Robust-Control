clc;
close all;
clear all;

m=0.7; % Poids
% m = ureal('m',0.7,'Percentage',20);
Ix=7.5*10^-3; % Moment d'inertie sur l'axe X
% Ix = ureal('Ix',7.5*10^-3,'Percentage',20);
Iy=7.5*10^-3; % Moment d'inertie sur l'axe Y
% Iy = ureal('Iy',7.5*10^-3,'Percentage',20);
Iz=1.3*10^-2; % Moment d'inertie sur l'axe Z
% Iz = ureal('Iz',1.3*10^-2,'Percentage',20);
l=0.23; % Longueur Moteur/CoG
b=3.13*10^-5; % Coefficient de poussée
% b = ureal('b',3.13*10^-5,'Percentage',10);
d=7.5*10^-7; % Coefficient de drag
% d = ureal('d',7.5*10^-7,'Percentage',10);
g=9.81; % Gravité
om=200^2;
om2=100^2;
tint=0.1*7.71; % Intensité de la rafale de vent à 15 knots

A= [ 0 1 0 0 0 0 0 0 ;
     0 0 0 0 0 0 0 0 ;
     0 0 0 0 0 1 0 0 ;
     0 0 0 0 0 0 1 0 ;
     0 0 0 0 0 0 0 1 ;
     0 0 0 0 0 0 0 0 ;
     0 0 0 0 0 0 0 0 ;
     0 0 0 0 0 0 0 0
     ];

B= [0 0 0 0    ;
    1/m 0 0 0  ;
    0 0 0 0    ;
    0 0 0 0    ;
    0 0 0 0    ;
    0 1/Ix 0 0 ;
    0 0 1/Iy 0 ;
    0 0 0 1/Iz
    ];


C=[ 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 0 0
    ];

D=[0 0 0 0;
   0 0 0 0;
   0 0 0 0;
   0 0 0 0
   ];

%% Matrice d'interconnexion pour le modèle de Dryden
W=[ 0 1 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 1 0;
    0 0 1 0 0 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];
%%% u w q v p r

BW = [0 0 0 0 1 0 0  ;
      1/m 0 0 0 0 0 0  ;
      0 0 0 0 0 0 1    ;
      0 0 0 0 0 1 0    ;
      0 0 0 0 0 0 0    ;
      0 1/Ix 0 0 0 0 0 ;
      0 0 1/Iy 0 0 0 0 ;
      0 0 0 1/Iz 0 0 0 ];

DW=[0 0 0 0 0 0 0 ;
   0 0 0 0 0 0 0 ;
   0 0 0 0 0 0 0 ;
   0 0 0 0 0 0 0 ];

G=ss(A,B,C,D); % R.E du système
GW=ss(A,BW,C,DW); % R.E augmentée avec le modèle du vent Dryden
s=zpk('s');

% CONTROL=rank(ctrb(G))
% OBSVERV=rank(obsv(G))

td11=0.999/(0.0625*s^2+0.34*s+1); % Wn=8 et Eps=0.71 -   (Z )
td12=0.999/(0.0156*s^2+0.1775*s+1); % Wn=8 et Eps=0.71 -   ( phi )
td13=0.999/(0.0156*s^2+0.1775*s+1);   % Wn=8 et Eps=0.68 -   ( theta )
td14=0.999/(0.04*s^2+0.2692*s+1);   % Wn=5 et Eps=0.6731 - ( psi )

w11=1/(1-td11);
w12=1/(1-td12);
w13=1/(1-td13);
w14=1/(1-td14);

% Matrice de pondération Wp ( performances )
W1=[w11 0 0 0;
    0 w12 0 0;
    0 0 w13 0 ;
    0 0 0 w14 ];

w21=1/1600;
w22=w21;
w23=w21;
w24=w21;

% Matrice de pondération Wu ( saturation )
W2=[w21 0 0 0; 
    0 w22 0 0; 
    0 0 w23 0;
    0 0 0 w24];



wdZ=80;
wdPhi=80;
wdTheta=wdPhi;
wdPsi=70;
eps=10^-2;
Am=0.1;

w31=(s+wdZ*Am)/(eps*s+wdZ);
w32=(s+wdPhi*Am)/(eps*s+wdPhi);
w33=(s+wdTheta*Am)/(eps*s+wdTheta);
w34=(s+wdPsi*Am)/(eps*s+wdPsi);

% Matrice de pondération Wr ( robustesse )
W3=[w31 0 0 0;
    0 w32 0 0 ;
    0 0 w33 0;
    0 0 0 w34];

P=augw(G,W1,W2,W3); % Système augmenté avec les pondérations
[K,CL,GAM]=hinfsyn(P); % Synthèse du contrôleur K
GAM 

L=G*K; % Système en boucle fermé
I=eye(size(L));
S=feedback(I,L); % Fonction de sensibilité
T=I-S; % Fonction de sensibilité complémentaire
figure(1)
hold on
sigma(S,'b',1/W1,'--r')
legend('S','1/W1')
hold off
figure(2)
hold on
sigma(S*K,'b',1/(ss(W2)),'--r')
legend('S*K','1/W2')
hold off
figure(3)
hold on
sigma(T,'b',1/W3,'--r')
legend('T','1/W3')