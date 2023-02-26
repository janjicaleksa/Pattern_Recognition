clear
close all
clc

broj_iteracija = [];

for iter=1:1
%% Generisanje odbiraka nelinearno separabilnih klasa
N = 500;

%Klasa 1
C = [6;8];
R = 2*rand(1,N);
Teta = 2*pi*rand(1,N);
K1 = [R.*cos(Teta) ; R.*sin(Teta)] + C*ones(1,N);

%Klasa 2
C = [6;8];
R =2*rand(1,N)+4;
Teta = 2*pi*rand(1,N);
K2 = [R.*cos(Teta);R.*sin(Teta)]+ C*ones(1,N);

figure(1);
scatter(K1(1,:),K1(2,:),'ro');
hold on
scatter(K2(1,:),K2(2,:),'bx');
grid on;
hold off
legend('Klasa 1','Klasa 2','Interpreter','latex');

%% Pocetna klasterizacija 

K = [K1(:,101:500) K2(:,101:500)];
%K1_novo = K1(:,1:100);
%K2_novo = K2(:,1:100);
K = [K1 K2];
K1_novo = []; K2_novo = [];

pom = rand(1,length(K));

for i = 1:length(pom)
    if pom(i) <=0.5
        K1_novo = [K1_novo K(:,i)];
    else
        K2_novo = [K2_novo K(:,i)];
    end
end

figure(2);
scatter(K1_novo(1,:),K1_novo(2,:),'ro');
hold on
scatter(K2_novo(1,:),K2_novo(2,:),'bx');
grid on;
hold off
legend('Klasa 1','Klasa 2','Interpreter','latex');
title('Iteracija 0 - proizvoljna raspodela','Interpreter','latex');

%% Normalna dekompozicija

l = 0; lmax = 100; reklas = 1;

while (l<lmax) && (reklas == 1)
    reklas = 0;
    l = l + 1;
    
    N1 = max(size(K1_novo)); N2 = max(size(K2_novo));
    M1 = mean(K1_novo,2);    M2 = mean(K2_novo,2);
    P1 = N1/(2*N);           P2 = N2/(2*N);
    
    S1 = cov(K1_novo');      S2 = cov(K2_novo');
    
    K1_pom = [];             K2_pom = [];
    
    for i = 1:N1
        J1 = 1/2 * (K1_novo(:,i)-M1)'*inv(S1)*(K1_novo(:,i)-M1)+1/2*log(det(S1))-log(P1);
        J2 = 1/2 * (K1_novo(:,i)-M2)'*inv(S2)*(K1_novo(:,i)-M2)+1/2*log(det(S2))-log(P2);
        if J1 < J2
            K1_pom = [K1_pom K1_novo(:,i)];
        else
            K2_pom = [K2_pom K1_novo(:,i)];
            reklas = 1;
        end
    end
    
    for i = 1:N2
        J1 = 1/2 * (K2_novo(:,i)-M1)'*inv(S1)*(K2_novo(:,i)-M1)+1/2*log(det(S1))-log(P1);
        J2 = 1/2 * (K2_novo(:,i)-M2)'*inv(S2)*(K2_novo(:,i)-M2)+1/2*log(det(S2))-log(P2);
        if J1 < J2
            K1_pom = [K1_pom K2_novo(:,i)];
            reklas = 1;
        else
            K2_pom = [K2_pom K2_novo(:,i)];
        end
    end
    
    clear K1_novo K2_novo
    K1_novo = K1_pom;
    K2_novo = K2_pom;
    
    
end

broj_iteracija = [broj_iteracija l];

figure(3)
scatter(K1_novo(1,:),K1_novo(2,:),'ro');
hold on 
scatter(K2_novo(1,:),K2_novo(2,:),'bx');
grid on;
hold off
legend('Klasa 1','Klasa 2','Interpreter','latex');
title(['Prikaz nakon reklasifikacije - iteracija ' num2str(broj_iteracija(end))],'Interpreter','latex');
end
srednja_vrednost = mean(broj_iteracija);
disp(['Srednja vrendnost potrebnog broja iteracije je: ',num2str(srednja_vrednost)]);
