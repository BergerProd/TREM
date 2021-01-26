%Transport & Melange
%Resolution Numérique
%Quentin BERGE & Vincent BERGNES
%ENSEEIHT - 3A MFEE - SEE
%-------------------------------
% Donnees Entree
%-------------------------------
% Numerique
sauv_graphe = true; % parametre 
cas=3;
%cas 1 = Jet sans stratification
%cas 2 = Panache force avec stratification
%cas 3 = Hauteur terminale

% Physique
z0=0; %debut du domaine : m
L=250; %longueur domaine : m
lambda =0.125;

%% Cas
if cas==1
%% Jet sans stratification
disp("Jet sans stratification")
% Frequence Brunt Vaissala
N=0;

% Conditions initiales
y0 = [pi pi 0];
%Zspan entre min et max
zspan=[z0 L];

%y(1)=Q
%y(2)=M
%y(3)=F


%Appel ODE45
[z,y] = ode45(@(z,y)odefun(z,y,lambda,N),zspan,y0);

%Creation des grandeurs b,w,g
b=(((y(:,1).^2))/(pi*y(:,2)))^(1/2);
w=y(:,2)./(y(:,1));
g=(y(:,3)./y(:,1));

% Figures
%Q
figure(1)
hold on
plot(z,y(:,1))
xlabel('z')
ylabel('Q')
%title('')
grid on
hold off

%M
figure(2)
hold on
plot(z,y(:,2))
xlabel('z')
ylabel('M')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off

%F
figure(3)
hold on
plot(z,y(:,3))
xlabel('z')
ylabel('F')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off

%b
figure(4)
hold on
plot(z,b)
xlabel('z')
ylabel('b')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off


%w
figure(5)
hold on
plot(z,w)
xlabel('z')
ylabel('w')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off


%g
figure(6)
hold on
plot(z,g)
xlabel('z')
ylabel('g')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off

elseif cas==2
%% Panache force avec stratification
disp("Panache force avec stratification")

N=0.05; % Frequence Brunt Vaissala

%y(1)=Q
%y(2)=M
%y(3)=F

% Conditions initiales
y0 = [pi pi 10*pi];

%valeurs par defaut avant passage dans la boucle while
indicateur_1=0;
indicateur_2=0;

%Tant que pour augmenter la taille du domaine si pas d'annulation des 2
%grandeurs
while(indicateur_1 ==0  || indicateur_2==0) %si indicateur 1=0 ou indicateur 2=0
    L=L+10;%augmentation du domaine
    
    %Zspan entre min et max
    zspan=[z0 L];
    %Appel ODE45
    [z,y] = ode45(@(z,y)odefun(z,y,lambda,N),zspan,y0);
    
    %Creation des grandeurs b,w,g
    b=(((y(:,1).^2))/(pi*y(:,2)))^(1/2);
    w=y(:,2)./(y(:,1));
    g=(y(:,3)./y(:,1));
    
    
    
    % Recherche des annulations
    
    
    % Flux de flottabilité
    indice_annulation=recherche_indice_annulation(y(:,3));
    if indice_annulation ~= 0
        abcisse_annulation=interp1(real(y(:,3)),z,0,'linear');%fonction interpolation lineaire plus facile
        fprintf("le flux de flottabilité s annule entre %6.2f m et %6.2f m \n", z(indice_annulation), z(indice_annulation+1))
        fprintf("avec une interpolation lineaire, s'annule précisement en %6.2f m \n",abcisse_annulation)
        indicateur_1=1;
    else
        disp("pas d'annulation de la flottabilité dans l'intervalle, augmentation du domaine de 10 m")
        indicateur_1=0;
    end
    
    % Vitesse d'ascension
    indice_annulation=recherche_indice_annulation(w(:));
    if indice_annulation ~= 0
        abcisse_annulation=interp1(real(w(:)),z,0,'linear');%fonction interpolation lineaire plus facile
        fprintf("la vitesse d'ascension s annule entre %6.2f m et %6.2f m \n", z(indice_annulation), z(indice_annulation+1))
        fprintf("avec une interpolation lineaire, s'annule précisement en %6.2f m \n",abcisse_annulation)
        indicateur_2=1;
    else
        disp("pas d'annulation de la vitesse d'ascension dans l'intervalle, augmentation du domaine de 10 m")
        indicateur_2=0;
    end
end

% Figures
%Q
figure(7)
hold on
plot(z,y(:,1))
xlabel('z')
ylabel('Q')
%title('')
grid on
hold off

%M
figure(8)
hold on
plot(z,y(:,2))
xlabel('z')
ylabel('M')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off

%F
figure(9)
hold on
plot(z,y(:,3))
xlabel('z')
ylabel('F')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off

%b
figure(10)
hold on
plot(z,b)
xlabel('z')
ylabel('b')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off


%w
figure(11)
hold on
plot(z,w)
xlabel('z')
ylabel('w')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off


%g
figure(12)
hold on
plot(z,g)
xlabel('z')
ylabel('g')
%title('Solutions transitoire de u et v dans le plan u,v')
grid on
hold off

elseif cas==3
%% Influence des conditions initiales sur le cas du panache forcé avec stratification
disp("Hauteur Terminal : Panache force avec stratification")


%y(1)=Q
%y(2)=M
%y(3)=F

% Conditions initiales
y0=zeros(3,3);
y0(:,1) = [pi pi 10*pi];
y0(:,2)=[0.049087 0.01227 0.049087];%w=0,25 b=0,25 g=1
y0(:,3)=[0.3927 0.1963 3.927];%w=0,5,b=0,5,g=10
N=[0.05, 0.25, 0.03]; % Frequence Brunt Vaissala

for i=1:3
    
%valeurs par defaut avant passage dans la boucle while
indicateur_1=0;

    while(indicateur_1 ==0) %si indicateur 1=0 ou indicateur 2=0
        L=L+10;%augmentation du domaine

        %Zspan entre min et max
        zspan=[z0 L];
        %Appel ODE45
        [z,y] = ode45(@(z,y)odefun(z,y,lambda,N(i)),zspan,y0(:,i));
        
        %Creation des grandeurs b,w,g
        b=(((y(:,1).^2))/(pi*y(:,2)))^(1/2);
        w=y(:,2)./(y(:,1));
        g=(y(:,3)./y(:,1));
                       
        % Recherche des annulations
        
        % Flux de flottabilité
        indice_annulation=recherche_indice_annulation(y(:,3));
        if indice_annulation ~= 0
            abcisse_annulation=interp1(real(y(:,3)),z,0,'linear');%fonction interpolation lineaire plus facile
            fprintf("le flux de flottabilité s annule entre %6.2f m et %6.2f m \n", z(indice_annulation), z(indice_annulation+1))
            fprintf("avec une interpolation lineaire, s'annule précisement en %6.2f m \n",abcisse_annulation)
            indicateur_1=1;
        else
            fprintf("pas d'annulation de la flottabilité dans l'intervalle, augmentation du domaine de 10 m, L= %i, i=%i \n",L,i)
            indicateur_1=0;
        end
    end
    ht(i)=abcisse_annulation;
end

% Résolution du systeme

% syms A n1 n2
% eqn1=A==ht(1)/(N(1)^n1*y0(3,1)^n2);
% eqn2=n1==(log(ht(2)/(A*y0(3,2)^n2))/(log(N(2))));
% eqn3=n2==(log(ht(3)/(A*N(3)^n1))/(log(y0(3,3))));
% 
% S = solve(A,n1,n2);
% A= double(S.A) 
% n1 = double(S.n1)
% n2=double(S.n2)

    

else
    disp("erreur saisie cas cas =1,2 ou 3")
end
%% Figures
%-------------------------------
%Sauvegarde des figures
%-------------------------------
if (sauv_graphe == true)
    path = pwd ;   % mention your path
    myfolder = 'graphe' ;   % new folder name
    folder = mkdir([path,filesep,myfolder]) ;
    path  = [path,filesep,myfolder] ;
    if cas==1
        for i = 1:6
            figure(i);
            temp=[path,filesep,'fig',num2str(i),'.png'];
            saveas(gcf,temp);
        end
     elseif cas==2
        for i = 7:12
            figure(i);
            temp=[path,filesep,'fig',num2str(i),'.png'];
            saveas(gcf,temp);
        end
    end
end

%% Fonctions

%Systeme derive pour ODE45
function dydz=odefun(z,y,lambda,N)
dydz=zeros(3,1);
dydz(1)=2*((y(1)^2)/(pi*y(2)))^(1/2)*lambda*(y(2)/y(1));
dydz(2)=(y(3)/y(1))*(y(1)^2/(pi*y(2)));
dydz(3)=-N^2*y(1)/(pi);
end

% renvoit l'indice ou cela passe par 0
function indice_annulation=recherche_indice_annulation(vecteur)
for i=1:length(vecteur)-1
    prod=vecteur(i)*vecteur(i+1);
    if prod < 0
        indice_annulation=i;
        break
    else
        indice_annulation=0;
    end
end
end


