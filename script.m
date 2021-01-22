%-------------------------------
% Donnees Entree
%-------------------------------
% Numerique
sauv_graphe = false; % parametre 
nptz=100;

% Physique
z0=0; %debut du domaine : m
L=100; %longueur domaine : m
lambda =0.125;
%-------------------------------
% Initialisation des vecteurs
%-------------------------------
z=zeros(1,nptz);
w_m_1=zeros(1,nptz);
b_1=zeros(1,nptz);

%% MAILLAGE
%Pas dz
dz = L/real(nptz-1);

%initialisation des noeuds
znoeuds(1) = 0;

%Noeuds des mailles
for i=2:nptz
    znoeuds(i) = znoeuds(1)+dz*(i-1);
end

%% Jet sans stratification

w_m_1(1)=1;
b_1(1)=1;




%% Figures
%-------------------------------
%Sauvegarde des figures
%-------------------------------
if (sauv_graphe == true)
    path = pwd ;   % mention your path
    myfolder = 'graphe' ;   % new folder name
    folder = mkdir([path,filesep,myfolder]) ;
    path  = [path,filesep,myfolder] ;
    for i = 1:7
        figure(i);
        temp=[path,filesep,'fig',num2str(i),'.png'];
        saveas(gcf,temp);
    end
end
