clear
close all
clc

N = 500;

M1=[0;0];
S1=[0.5,0;0,0.5];  

M2=[0;7];
S2=[0.8,0.2;0.2,0.5];

M3=[7;7];
S3=[0.6,-0.3;-0.3,0.8];

M4=[7;0];
S4=[0.4,0.3;0.3,0.8];


K1 = mvnrnd(M1,S1,N)';
K2 = mvnrnd(M2,S2,N)';
K3 = mvnrnd(M3,S3,N)';
K4 = mvnrnd(M4,S4,N)';

figure(1)
hold all
scatter(K1(1,:),K1(2,:),'ro');
scatter(K2(1,:),K2(2,:),'bx');
scatter(K3(1,:),K3(2,:),'mv');
scatter(K4(1,:),K4(2,:),'yd');
grid on;
hold off
legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4','Interpreter','latex');

% Mesanje klasa

pom = rand(1,4*N);
K1_novo = []; K2_novo = [];
K3_novo = []; K4_novo = []; K5_novo = [];

for i = 1:N %Prolazak kroz klasu K1
    if pom(i)<0.2
        K1_novo = [K1_novo K1(:,i)];
    elseif pom(i)<0.4
        K2_novo = [K2_novo K1(:,i)];
    elseif pom(i)<0.6
        K3_novo = [K3_novo K1(:,i)];
    elseif pom(i)<0.8
        K4_novo = [K4_novo K1(:,i)];
    else
        K5_novo = [K5_novo K1(:,i)];
    end
end

for i = N+1:2*N
    if pom(i)<0.2
        K1_novo = [K1_novo K2(:,i-N)];
    elseif pom(i)<0.4
        K2_novo = [K2_novo K2(:,i-N)];
    elseif pom(i)<0.6
        K3_novo = [K3_novo K2(:,i-N)];
    elseif pom(i)<0.8
        K4_novo = [K4_novo K2(:,i-N)];
    else
        K5_novo = [K5_novo K2(:,i-N)];
    end
end
 
for i = 2*N+1:3*N
    if pom(i)<0.2
        K1_novo = [K1_novo K3(:,i-2*N)];
    elseif pom(i)<0.4
        K2_novo = [K2_novo K3(:,i-2*N)];
    elseif pom(i)<0.6
        K3_novo = [K3_novo K3(:,i-2*N)];
    elseif pom(i)<0.8
        K4_novo = [K4_novo K3(:,i-2*N)];
    else
        K5_novo = [K5_novo K3(:,i-2*N)];
    end
end
 
for i = 3*N+1:4*N
    if pom(i)<0.2
        K1_novo = [K1_novo K4(:,i-3*N)];
    elseif pom(i)<0.4
        K2_novo = [K2_novo K4(:,i-3*N)];
    elseif pom(i)<0.6
        K3_novo = [K3_novo K4(:,i-3*N)];
    elseif pom(i)<0.8
        K4_novo = [K4_novo K4(:,i-3*N)];
    else
        K5_novo = [K5_novo K4(:,i-3*N)];
    end
end
    
figure(2)
hold all
scatter(K1_novo(1,:),K1_novo(2,:),'ro');
scatter(K2_novo(1,:),K2_novo(2,:),'bx');
scatter(K3_novo(1,:),K3_novo(2,:),'mv');
scatter(K4_novo(1,:),K4_novo(2,:),'yd');
scatter(K5_novo(1,:),K5_novo(2,:),'k.');
grid on;
hold off
legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4','Klasa 5','Interpreter','latex');
title('Iteracija 0 - proizvoljna raspodela','Interpreter','latex');
    
%% C-mean metod

N1 = max(size(K1_novo));
N2 = max(size(K2_novo));
N3 = max(size(K3_novo));
N4 = max(size(K4_novo));
N5 = max(size(K5_novo));

M1 = mean(K1_novo,2);    
M2 = mean(K2_novo,2);    
M3 = mean(K3_novo,2);    
M4 = mean(K4_novo,2);
M5 = mean(K5_novo,2);

l = 1; lmax = 100; reklas = 1;

while (l<lmax) && (reklas==1)
    K1pom = []; K2pom = [];
    K3pom = []; K4pom = []; K5pom = [];
    reklas = 0;
    
    for i = 1:N1  %prolazimo kroz K1_novo
        d1 = sum((K1_novo(:,i)-M1).^2);
        d2 = sum((K1_novo(:,i)-M2).^2);
        d3 = sum((K1_novo(:,i)-M3).^2);
        d4 = sum((K1_novo(:,i)-M4).^2);
        d5 = sum((K1_novo(:,i)-M5).^2);
        D = [d1,d2,d3,d4,d5];
        D = sort(D);
        if D(1) == d1
            K1pom = [K1pom K1_novo(:,i)];
        elseif D(1)==d2
            K2pom = [K2pom K1_novo(:,i)];
            reklas = 1;
        elseif D(1)==d3
            K3pom = [K3pom K1_novo(:,i)];
            reklas = 1;
        elseif D(1)==d4
            K4pom = [K4pom K1_novo(:,i)];
            reklas = 1;
        else 
            K5pom = [K5pom K1_novo(:,i)];
            reklas = 1;
        end
    end

    for i = 1:N2  %prolazimo kroz K2_novo
        d1 = sum((K2_novo(:,i)-M1).^2);
        d2 = sum((K2_novo(:,i)-M2).^2);
        d3 = sum((K2_novo(:,i)-M3).^2);
        d4 = sum((K2_novo(:,i)-M4).^2);
        d5 = sum((K2_novo(:,i)-M5).^2);
        D = [d1,d2,d3,d4,d5];
        D = sort(D);
        if D(1) == d1
            K1pom = [K1pom K2_novo(:,i)];
            reklas = 1;
        elseif D(1)==d2
            K2pom = [K2pom K2_novo(:,i)];
        elseif D(1)==d3
            K3pom = [K3pom K2_novo(:,i)];
            reklas = 1;
        elseif D(1)==d4
            K4pom = [K4pom K2_novo(:,i)];
            reklas = 1;
        else
            K5pom = [K5pom K2_novo(:,i)];
            reklas = 1;
        end
     end
    
     for i = 1:N3  %prolazimo kroz K3_novo
        d1 = sum((K3_novo(:,i)-M1).^2);
        d2 = sum((K3_novo(:,i)-M2).^2);
        d3 = sum((K3_novo(:,i)-M3).^2);
        d4 = sum((K3_novo(:,i)-M4).^2);
        d5 = sum((K3_novo(:,i)-M5).^2);
        D = [d1,d2,d3,d4,d5];
        D = sort(D);
        if D(1) == d1
            K1pom = [K1pom K3_novo(:,i)];
            reklas = 1;
        elseif D(1)==d2
            K2pom = [K2pom K3_novo(:,i)];
            reklas = 1;
        elseif D(1)==d3
            K3pom = [K3pom K3_novo(:,i)];
        elseif D(1)==d4
            K4pom = [K4pom K3_novo(:,i)];
            reklas = 1;
        else
            K5pom = [K5pom K3_novo(:,i)];
            reklas = 1;
        end
     end
    
     for i = 1:N4  %prolazimo kroz K4_novo
        d1 = sum((K4_novo(:,i)-M1).^2);
        d2 = sum((K4_novo(:,i)-M2).^2);
        d3 = sum((K4_novo(:,i)-M3).^2);
        d4 = sum((K4_novo(:,i)-M4).^2);
        d5 = sum((K4_novo(:,i)-M5).^2);
        D = [d1,d2,d3,d4,d5];
        D = sort(D);
        if D(1) == d1
            K1pom = [K1pom K4_novo(:,i)];
            reklas = 1;
        elseif D(1)==d2
            K2pom = [K2pom K4_novo(:,i)];
            reklas = 1;
        elseif D(1)==d3
            K3pom = [K3pom K4_novo(:,i)];
            reklas = 1;
        elseif D(1)==d4
            K4pom = [K4pom K4_novo(:,i)];
        else
            K5pom =[K5pom K4_novo(:,i)];
            reklas =1;
        end
     end
     
     
     for i = 1:N5  %prolazimo kroz K4_novo
        d1 = sum((K5_novo(:,i)-M1).^2);
        d2 = sum((K5_novo(:,i)-M2).^2);
        d3 = sum((K5_novo(:,i)-M3).^2);
        d4 = sum((K5_novo(:,i)-M4).^2);
        d5 = sum((K5_novo(:,i)-M5).^2);
        D = [d1,d2,d3,d4,d5];
        D = sort(D);
        if D(1) == d1
            K1pom = [K1pom K5_novo(:,i)];
            reklas = 1;
        elseif D(1)==d2
            K2pom = [K2pom K5_novo(:,i)];
            reklas = 1;
        elseif D(1)==d3
            K3pom = [K3pom K5_novo(:,i)];
            reklas = 1;
        elseif D(1)==d4
            K4pom = [K4pom K5_novo(:,i)];
            reklas = 1;
        else
            K5pom =[K5pom K5_novo(:,i)];
           
        end
     end
    
    clear K1_novo K2_novo K3_novo K4_novo K5_novo
    K1_novo = K1pom; K2_novo = K2pom;
    K3_novo = K3pom; K4_novo = K4pom; K5_novo = K5pom;
    figure()
    hold all
    scatter(K1_novo(1,:),K1_novo(2,:),'ro');
    scatter(K2_novo(1,:),K2_novo(2,:),'bx');
    scatter(K3_novo(1,:),K3_novo(2,:),'mv');
    scatter(K4_novo(1,:),K4_novo(2,:),'yd');
    scatter(K5_novo(1,:),K5_novo(2,:),'k.');
    grid on;
    hold off
    legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4','Klasa 5','Interpreter','latex');
    title(['Iteracija broj ' num2str(l)],'Interpreter','latex');
    pause;
    
    N1 = max(size(K1_novo));    N2 = max(size(K2_novo));
    N3 = max(size(K3_novo));    N4 = max(size(K4_novo));
    N5 = max(size(K5_novo));

    M1 = mean(K1_novo,2);       M2 = mean(K2_novo,2);    
    M3 = mean(K3_novo,2);       M4 = mean(K4_novo,2);
    M5 = mean(K5_novo,2);
    
    l = l+1;
end