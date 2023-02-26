clear
clc
close all

broj_iteracija = [];
N = 500;
for iter = 1:1
    
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

%% Pocetna klasterizacija

pom = rand(1,4*N);
% K1_novo = K1(:,1:125);
% K2_novo = K2(:,1:125);
% K3_novo = K3(:,1:125);
% K4_novo = K4(:,1:125);
K1_novo = []; K2_novo = [];
K3_novo = []; K4_novo = [];

for i = 1:N %Prolazak kroz klasu K1
    if pom(i)<0.25
        K1_novo = [K1_novo K1(:,i)];
    elseif pom(i)<0.5
        K2_novo = [K2_novo K1(:,i)];
    elseif pom(i)<0.75
        K3_novo = [K3_novo K1(:,i)];
    else
        K4_novo = [K4_novo K1(:,i)];
    end
end

for i = N+1:2*N
    if pom(i)<0.25
        K1_novo = [K1_novo K2(:,i-N)];
    elseif pom(i)<0.5
        K2_novo = [K2_novo K2(:,i-N)];
    elseif pom(i)<0.75
        K3_novo = [K3_novo K2(:,i-N)];
    else
        K4_novo = [K4_novo K2(:,i-N)];
    end
end
 
for i = 2*N+1:3*N
    if pom(i)<0.25
        K1_novo = [K1_novo K3(:,i-2*N)];
    elseif pom(i)<0.5
        K2_novo = [K2_novo K3(:,i-2*N)];
    elseif pom(i)<0.75
        K3_novo = [K3_novo K3(:,i-2*N)];
    else
        K4_novo = [K4_novo K3(:,i-2*N)];
    end
end
 
for i = 3*N+1:4*N
    if pom(i)<0.25
        K1_novo = [K1_novo K4(:,i-3*N)];
    elseif pom(i)<0.5
        K2_novo = [K2_novo K4(:,i-3*N)];
    elseif pom(i)<0.75
        K3_novo = [K3_novo K4(:,i-3*N)];
    else
        K4_novo = [K4_novo K4(:,i-3*N)];
    end
end

figure(2)
hold all
scatter(K1_novo(1,:),K1_novo(2,:),'ro');
scatter(K2_novo(1,:),K2_novo(2,:),'bx');
scatter(K3_novo(1,:),K3_novo(2,:),'mv');
scatter(K4_novo(1,:),K4_novo(2,:),'yd');
grid on;
hold off
legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4','Interpreter','latex');
title('Iteracija 0 - proizvoljna raspodela','Interpreter','latex');

%% Odredjivanje pocetnih parametara

N1 = max(size(K1_novo));    N2 = max(size(K2_novo));
N3 = max(size(K3_novo));    N4 = max(size(K4_novo));


P1 = N1/4/N; P2 = N2/4/N;
P3 = N3/4/N; P4 = N4/4/N;
% Nu = N1+N2+N3+N4;
% P1 = N1/Nu; P2 = N2/Nu;
% P3 = N3/Nu; P4 = N4/Nu;

M1 = mean(K1_novo,2);    M2 = mean(K2_novo,2);    
M3 = mean(K3_novo,2);    M4 = mean(K4_novo,2);

S1 = cov(K1');           S2 = cov(K2');
S3 = cov(K3');           S4 = cov(K4');

lmax=100;
l=0; % pocetna iteracija 
L=4; % broj klasa
f1=zeros(1,4*N);         f2=zeros(1,4*N);
f3=zeros(1,4*N);         f4=zeros(1,4*N);
f = zeros(1,4*N);

q1 = zeros(1,4*N);       q2 = zeros(1,4*N);
q3 = zeros(1,4*N);       q4 = zeros(1,4*N);

% f1=zeros(1,Nu);         f2=zeros(1,Nu);
% f3=zeros(1,Nu);         f4=zeros(1,Nu);
% f = zeros(1,Nu);
% 
% q1 = zeros(1,Nu);       q2 = zeros(1,Nu);
% q3 = zeros(1,Nu);       q4 = zeros(1,Nu);

K = [K1_novo K2_novo K3_novo K4_novo];

const1 = 1/(2*pi*det(S1)^0.5);  const2 = 1/(2*pi*det(S2)^0.5);
const3 = 1/(2*pi*det(S3)^0.5);  const4 = 1/(2*pi*det(S4)^0.5);

for i = 1:4*N %kad nije slucajna raspodela na pocetku onda do Nu
    f1(i) = const1*exp(-0.5*(K(:,i)-M1)'*inv(S1)*(K(:,i)-M1));
    f2(i) = const2*exp(-0.5*(K(:,i)-M2)'*inv(S2)*(K(:,i)-M2));
    f3(i) = const3*exp(-0.5*(K(:,i)-M3)'*inv(S3)*(K(:,i)-M3));
    f4(i) = const4*exp(-0.5*(K(:,i)-M4)'*inv(S4)*(K(:,i)-M4));
    f(i) = P1*f1(i)+P2*f2(i)+P3*f3(i)+P4*f4(i);
    
    q1(i) = P1*f1(i)/f(i);          q2(i) = P2*f2(i)/f(i);
    q3(i) = P3*f3(i)/f(i);          q4(i) = P4*f4(i)/f(i);
end

%% Maximum Likelihood 

zavrsi = 0;
T = 0.01;
n = 0;

while (zavrsi == 0)
    n = n+1;
    
    P1 = sum(q1)/length(q1);    N1 = P1*(4*N); 
    P2 = sum(q2)/length(q2);    N2 = P2*(4*N);
    P3 = sum(q3)/length(q3);    N3 = P3*(4*N);
    P4 = sum(q4)/length(q4);    N4 = P4*(4*N);
    
%     N1 = P1*Nu; N2 = P2*Nu; N3 = P3*Nu; N4 = P4*Nu;
    
    M1 = [0;0]; M2 = M1; M3 = M1; M4 = M1;
    
    for i = 1:4*N %inace Nu
        M1 = M1 + q1(i)*K(:,i);
        M2 = M2 + q2(i)*K(:,i);
        M3 = M3 + q3(i)*K(:,i);
        M4 = M4 + q4(i)*K(:,i);
    end
    
    M1 = M1/N1;     M2 = M2/N2;
    M3 = M3/N3;     M4 = M4/N4;
    
    S1 = zeros(2, 2);       S2 = zeros(2, 2);
    S3 = zeros(2, 2);       S4 = zeros(2, 2);
    for i = 1:4*N %inace Nu
       S1 = S1 + q1(i)*(K(:, i) - M1)*(K(:, i) - M1)'; 
       S2 = S2 + q2(i)*(K(:, i) - M2)*(K(:, i) - M2)'; 
       S3 = S3 + q3(i)*(K(:, i) - M3)*(K(:, i) - M3)'; 
       S4 = S4 + q4(i)*(K(:, i) - M4)*(K(:, i) - M4)'; 
    end
    
    S1 = S1/N1;
    S2 = S2/N2;
    S3 = S3/N3;
    S4 = S4/N4;
    
    const1 = 1/(2*pi*det(S1)^0.5);  const2 = 1/(2*pi*det(S2)^0.5);
    const3 = 1/(2*pi*det(S3)^0.5);  const4 = 1/(2*pi*det(S4)^0.5);

for i = 1:4*N %inace Nu
    f1(i) = const1*exp(-0.5*(K(:,i)-M1)'*inv(S1)*(K(:,i)-M1));
    f2(i) = const2*exp(-0.5*(K(:,i)-M2)'*inv(S2)*(K(:,i)-M2));
    f3(i) = const3*exp(-0.5*(K(:,i)-M3)'*inv(S3)*(K(:,i)-M3));
    f4(i) = const4*exp(-0.5*(K(:,i)-M4)'*inv(S4)*(K(:,i)-M4));
    f(i) = P1*f1(i)+P2*f2(i)+P3*f3(i)+P4*f4(i);
    
    q1_novo(i) = P1*f1(i)/f(i);          q2_novo(i) = P2*f2(i)/f(i);
    q3_novo(i) = P3*f3(i)/f(i);          q4_novo(i) = P4*f4(i)/f(i);
end
    
    m1 = max(abs(q1-q1_novo));          m2 = max(abs(q2-q2_novo));
    m3 = max(abs(q3-q3_novo));          m4 = max(abs(q4-q4_novo));
    
    m = max([m1 m2 m3 m4]);
    
    if m<T
        zavrsi = 1;
    else
        clear q1 q2 q3 q4
        q1 = q1_novo;       q2 = q2_novo;
        q3 = q3_novo;       q4 = q4_novo;
    end
    
end
broj_iteracija = [broj_iteracija n];
end

%% Rasporedjivanje odbiraka u nove klase

K1_novo = []; K2_novo = [];
K3_novo = []; K4_novo = [];


for i = 1:4*N %inace Nu
    q = sort([q1_novo(i),q2_novo(i),q3_novo(i) q4_novo(i)],'descend');
    if q(1) == q1_novo(i)
        K1_novo = [K1_novo K(:,i)];
    elseif q(1) == q2_novo(i)
        K2_novo = [K2_novo K(:,i)];
    elseif q(1) == q3_novo(i)
        K3_novo = [K3_novo K(:,i)];
    elseif q(1) == q4_novo(i)
        K4_novo = [K4_novo K(:,i)];
    end
end

figure(3)
hold all
scatter(K1_novo(1,:),K1_novo(2,:),'ro');
scatter(K2_novo(1,:),K2_novo(2,:),'bx');
scatter(K3_novo(1,:),K3_novo(2,:),'mv');
scatter(K4_novo(1,:),K4_novo(2,:),'yd');
grid on;
hold off
legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4','Interpreter','latex');
title(['Prikaz nakon reklasifikacije - iteracija ' num2str(broj_iteracija(end))],'Interpreter','latex');

srednja_vrednost = mean(broj_iteracija);
disp(['Srednja vrednost broja iteracija: ',num2str(srednja_vrednost)]);