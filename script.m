%-------------------------------
% Donnees Entree
%-------------------------------
% Numerique
sauv_graphe = true; % parametre 
cas=2;
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
     fprintf("le flux de flottabilité s annule entre %f et %f \n", z(indice_annulation), z(indice_annulation+1))
else
     disp("pas d'annulation de la flottabilité dans l'intervalle, augmentation du domaine de 10 m")
     indicateur_1=0;
end     

% Vitesse d'ascension
indice_annulation=recherche_indice_annulation(w(:));
if indice_annulation ~= 0
     fprintf("la vitesse d'ascension s annule entre %f et %f \n", z(indice_annulation), z(indice_annulation+1) )
else
     disp("pas d'annulation de la vitesse d'ascension dans l'intervalle, augmentation du domaine de 10 m")
     indicateur_2=0;
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


