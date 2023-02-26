clear all
close all 
clc

%% tacka 1a)
%% Generisanje odbiraka tri linearno separabilne klase

M1 = [-8;-8]; S1 = [1 -0.5;-0.5 1];
M2 = [0;0]; S2 = [1 0.2;0.2 1];
M3 = [8;0]; S3 = [1 0.5;0.5 1];
N = 500;

rng(100)
X1 = mvnrnd(M1,S1,N);
X2 = mvnrnd(M2,S2,N);
X3 = mvnrnd(M3,S3,N);

figure(1)
hold all
plot(X1(:,1),X1(:,2),'c*');
plot(X2(:,1),X2(:,2),'md');
plot(X3(:,1),X3(:,2),'ys');
grid on
legend('Klasa 1','Klasa 2','Klasa 3','Interpreter','latex');

%% Treca numericka metoda - hold out
% Podela podataka na podatke za testiranje i obucavanje

X1_obucavanje = X1 (1:round(0.75*N),:)';
X1_testiranje = X1 (round(0.75*N):N,:)';

X2_obucavanje = X2 (1:round(0.75*N),:)';
X2_testiranje = X2 (round(0.75*N):N,:)';

X3_obucavanje = X3 (1:round(0.75*N),:)';
X3_testiranje = X3 (round(0.75*N):N,:)';

M1_est = mean(X1_obucavanje,2); S1_est = cov(X1_obucavanje');
M2_est = mean(X2_obucavanje,2); S2_est = cov(X2_obucavanje');
M3_est = mean(X3_obucavanje,2); S3_est = cov(X3_obucavanje');

N1 = length(X1_testiranje);
N2 = length(X2_testiranje);
N3 = length(X3_testiranje);

s = 0:0.001:1;
Neps_s = zeros(3,length(s));
V0_opt_s = zeros(3,length(s));

for i = 1:length(s)
    V1 = (s(i)*S1_est + (1-s(i))*S2_est)^(-1)*(M2_est - M1_est);
    V2 = (s(i)*S2_est + (1-s(i))*S3_est)^(-1)*(M3_est - M2_est);
    V3 = (s(i)*S1_est + (1-s(i))*S3_est)^(-1)*(M3_est - M1_est);
    
    Y1_1 = V1'*X1_testiranje; Y2_2 = V2'*X2_testiranje; Y3_3 = V3'*X3_testiranje;
    Y2_1 = V1'*X2_testiranje; Y2_3 = V2'*X3_testiranje; Y3_1 = V3'*X1_testiranje;
    
    Y1 = [Y1_1 Y2_1]; Y1 = sort(Y1);
    Y2 = [Y2_2 Y2_3]; Y2 = sort(Y2);        
    Y3 = [Y3_3 Y3_1]; Y3 = sort(Y3);
    V0 = zeros(3,length(Y1)-1);
    Neps = zeros(3,length(Y1)-1);
    for j =1: (length(Y1)-1)
        V0(1,j) = -(Y1(j) + Y1(j+1))/2;
        V0(2,j) = -(Y2(j) + Y2(j+1))/2;
        V0(3,j) = -(Y3(j) + Y3(j+1))/2;
        
        for k = 1:N1
            if(Y1(k)>-V0(1,j))
                Neps(1,j) = Neps(1,j)+1;
            end
            if (Y2(k) < -V0(1,j))
                Neps(1,j)= Neps(1,j)+1;
            end
        end
        
        for k = 1:N2
           if(Y2(k)>-V0(2,j))
                Neps(2,j) = Neps(2,j)+1;
           end
           if (Y3(k) < -V0(2,j))
                Neps(2,j)= Neps(2,j)+1;
           end
        end
        
        for k = 1:N3
            if(Y1(k)>-V0(3,j))
                Neps(3,j) = Neps(3,j)+1;
            end
            if (Y3(k) < -V0(3,j))
                Neps(3,j)= Neps(3,j)+1;
            end
        end
    
    end
    
    [Neps_s(1,i),index] = min(Neps(1,:));
    V0_opt_s(1,i) = V0(1,index);
    
    [Neps_s(2,i),index] = min(Neps(2,:));
    V0_opt_s(2,i) = V0(2,index);
    
    [Neps_s(3,i),index] = min(Neps(3,:));
    V0_opt_s(3,i) = V0(3,index);
    
end

Neps_opt = zeros(3,1);
index = zeros(3,1);

[Neps_opt(1,1), index(1,1)] = min(Neps_s(1,:));
V0_opt_1 = V0_opt_s(1,index(1,1));
s_opt_1 = s(index(1,1));

[Neps_opt(2,1), index(2,1)] = min(Neps_s(2,:));
V0_opt_2 = V0_opt_s(2,index(2,1));
s_opt_2 = s(index(2,1));

[Neps_opt(3,1), index(3,1)] = min(Neps_s(3,:));
V0_opt_3 = V0_opt_s(3,index(3,1));
s_opt_3 = s(index(3,1));

%% Iscrtavanje klasifikacionih linija

V1 = (s_opt_1*S1_est + (1-s_opt_1)*S2_est)^(-1)*(M2_est - M1_est);
V2 = (s_opt_2*S2_est + (1-s_opt_2)*S3_est)^(-1)*(M3_est - M2_est);
V3 = (s_opt_3*S1_est + (1-s_opt_3)*S3_est)^(-1)*(M3_est - M1_est); 

x = -12:0.01:12;
x1 = -(V0_opt_1 + V1(1)*x)/V1(2);
x2 = -(V0_opt_2 + V2(1)*x)/V2(2);
x3 = -(V0_opt_3 + V3(1)*x)/V3(2);

figure(1)
plot(x,x1,'g--','LineWidth',3);
plot(x,x2,'k--','LineWidth',3);
plot(x,x3,'r--','LineWidth',3);
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
legend('Klasa 1','Klasa 2','Klasa 3','Klasifikator klasa 1 i 2','Klasifikator klasa 2 i 3','Klasifikator klasa 1 i 3','Interpreter','latex');
axis([-12 12 -15 20])
hold off

%% Izracunavanje broja pogresno klasifikovanih

a1 = (x1(20)-x1(10))/(x(20)-x(10));
b1 = -a1*x(10) + x1(10);

br_gresaka_12 = 0;
br_gresaka_21 = 0;

%Racunanje broja odbiraka K1 iznad klasifikacione linije izmedju K1 i K2
for i =1:length(X1_testiranje)
   y1 = X1_testiranje(2,i);
   y0 = a1*X1_testiranje(1,i) + b1;
   if y1>y0
       br_gresaka_12 = br_gresaka_12+1;
   end
end

%Racunanje broja odbiraka K1 ispod klasifikacione linije izmedju K1 i K2
for i =1:length(X2_testiranje)
   y1 = X2_testiranje(2,i);
   y0 = a1*X2_testiranje(1,i) + b1;
   if y1<y0
       br_gresaka_21 = br_gresaka_21+1;
   end
end

Konfuziona_matrica_12 = zeros(2,2);
Konfuziona_matrica_12(1,1) = length(X1_testiranje) - br_gresaka_12;
Konfuziona_matrica_12(2,2) = length(X2_testiranje) - br_gresaka_21;
Konfuziona_matrica_12(1,2) = br_gresaka_12;
Konfuziona_matrica_12(2,1) = br_gresaka_21;
disp('Konfuziona matrica - klase 1 i 2');
disp(Konfuziona_matrica_12)

% Klase 2 i 3
a2 = (x2(20)-x2(10))/(x(20)-x(10));
b2 = -a2*x(10) + x2(10);

br_gresaka_23 = 0;
br_gresaka_32 = 0;

%Racunanje broja odbiraka K2 ispod klasifikacione linije izmedju K2 i K3
for i =1:length(X2_testiranje)
   y1 = X2_testiranje(2,i);
   y0 = a2*X2_testiranje(1,i) + b2;
   if y1<y0
       br_gresaka_23 = br_gresaka_23 + 1;
   end
end

%Racunanje broja odbiraka K3 iznad klasifikacione linije izmedju K2 i K3
for i =1:length(X3_testiranje)
   y1 = X3_testiranje(2,i);
   y0 = a2*X3_testiranje(1,i) + b2;
   if y1>y0
       br_gresaka_32 = br_gresaka_32 + 1;
   end
end

Konfuziona_matrica_23 = zeros(2,2);
Konfuziona_matrica_23(1,1) = 126 - br_gresaka_23;
Konfuziona_matrica_23(2,2) = 126 - br_gresaka_32;
Konfuziona_matrica_23(1,2) = br_gresaka_23;
Konfuziona_matrica_23(2,1) = br_gresaka_32;
disp('Konfuziona matrica - klase 2 i 3');
disp(Konfuziona_matrica_23)

% Klase 1 i 3
a3 = (x3(20)-x3(10))/(x(20)-x(10));
b3 = -a3*x(10) + x3(10);

br_gresaka_13 = 0;
br_gresaka_31 = 0;

%Racunanje broja odbiraka K1 ispod klasifikacione linije izmedju K1 i K3
for i =1:length(X1_testiranje)
   y1 = X1_testiranje(2,i);
   y0 = a3*X1_testiranje(1,i) + b3;
   if y1<y0
       br_gresaka_13 = br_gresaka_13 + 1;
   end
end

%Racunanje broja odbiraka K3 iznad klasifikacione linije izmedju K1 i K3
for i =1:length(X3_testiranje)
   y1 = X3_testiranje(2,i);
   y0 = a3*X3_testiranje(1,i) + b3;
   if y1>y0
       br_gresaka_31 = br_gresaka_31 + 1;
   end
end

Konfuziona_matrica_13 = zeros(2,2);
Konfuziona_matrica_13(1,1) = length(X1_testiranje) - br_gresaka_13;
Konfuziona_matrica_13(2,2) = length(X3_testiranje) - br_gresaka_31;
Konfuziona_matrica_13(1,2) = br_gresaka_13;
Konfuziona_matrica_13(2,1) = br_gresaka_31;
disp('Konfuziona matrica - klase 1 i 3');
disp(Konfuziona_matrica_13)

%% tacka 1b)
%% Iscrtavanje generisanih odbiraka tri linearno separabilne klase iz tacke 1a)

figure(2)
hold all
plot(X1(:,1),X1(:,2),'c*');
plot(X2(:,1),X2(:,2),'md');
plot(X3(:,1),X3(:,2),'ys');
grid on
legend('Klasa 1','Klasa 2','Klasa 3','Interpreter','latex');

%% Linearni klasifikator metodom zeljenog izlaza

U1 = [-1*ones(1,N),ones(1,N);-1*X1',X2'];
U2 = [-1*ones(1,N),ones(1,N);-1*X2',X3'];
U3 = [-1*ones(1,N),ones(1,N);-1*X1',X3'];

Gama1 = ones(2*N,1); %Ako hocemo neku od klasa da bude prioritetna samo pomnozimo Gama 
Gama2 = ones(2*N,1);
Gama3 = ones(2*N,1);

W1 = inv(U1*U1')*U1*Gama1;
W2 = inv(U2*U2')*U2*Gama2;
W3 = inv(U3*U3')*U3*Gama3;

V0a = W1(1);
V0b = W2(1);
V0c = W3(1);
V1a = W1(2);
V1b = W2(2);
V1c = W3(2);
V2a = W1(3);
V2b = W2(3);
V2c = W3(3);

x = -12:0.1:12;
x1 = -(V0a + V1a*x)/V2a;
x2 = -(V0b + V1b*x)/V2b;
x3 = -(V0c + V1c*x)/V2c;

figure(2)
plot(x,x1,'g--','LineWidth',3);
plot(x,x2,'k--','LineWidth',3);
plot(x,x3,'r--','LineWidth',3);
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
legend('Klasa 1','Klasa 2','Klasa 3','Klasifikator klasa 1 i 2','Klasifikator klasa 2 i 3','Klasifikator klasa 1 i 3','Interpreter','latex');
axis([-12 12 -15 20])
hold off

%% Izracunavanje broja pogresno klasifikovanih

X1 = X1'; X2 = X2'; X3 = X3';

a1 = (x1(20)-x1(10))/(x(20)-x(10));
b1 = -a1*x(10) + x1(10);

br_gresaka_12 = 0;
br_gresaka_21 = 0;

%Racunanje broja odbiraka K1 iznad klasifikacione linije izmedju K1 i K2
for i =1:length(X1)
   y1 = X1(2,i);
   y0 = a1*X1(2,i) + b1;
   if y1>y0
       br_gresaka_12 = br_gresaka_12 + 1;
   end
end

%Racunanje broja odbiraka K1 ispod klasifikacione linije izmedju K1 i K2
for i =1:length(X2)
   y1 = X2(2,i);
   y0 = a1*X2(1,i) + b1;
   if y1<y0
       br_gresaka_21 = br_gresaka_21 + 1;
   end
end

Konfuziona_matrica_12 = zeros(2,2);
Konfuziona_matrica_12(1,1) = length(X1) - br_gresaka_12;
Konfuziona_matrica_12(2,2) = length(X2) - br_gresaka_21;
Konfuziona_matrica_12(1,2) = br_gresaka_12;
Konfuziona_matrica_12(2,1) = br_gresaka_21;
disp('Konfuziona matrica - metod zeljenog izlaza - klase 1 i 2');
disp(Konfuziona_matrica_12)

% Klase 2 i 3
a2 = (x2(20)-x2(10))/(x(20)-x(10));
b2 = -a2*x(10) + x2(10);

br_gresaka_23 = 0;
br_gresaka_32 = 0;

%Racunanje broja odbiraka K2 ispod klasifikacione linije izmedju K2 i K3
for i =1:length(X2)
   y1 = X2(2,i);
   y0 = a2*X2(1,i) + b2;
   if y1<y0
       br_gresaka_23 = br_gresaka_23 + 1;
   end
end

%Racunanje broja odbiraka K3 iznad klasifikacione linije izmedju K2 i K3
for i =1:length(X3)
   y1 = X3(2,i);
   y0 = a2*X3(1,i) + b2;
   if y1>y0
       br_gresaka_32 = br_gresaka_32 + 1;
   end
end

Konfuziona_matrica_23 = zeros(2,2);
Konfuziona_matrica_23(1,1) = length(X2) - br_gresaka_23;
Konfuziona_matrica_23(2,2) = length(X3) - br_gresaka_32;
Konfuziona_matrica_23(1,2) = br_gresaka_23;
Konfuziona_matrica_23(2,1) = br_gresaka_32;
disp('Konfuziona matrica - metod zeljenog izlaza - klase 2 i 3');
disp(Konfuziona_matrica_23)

% Klase 1 i 3
a3 = (x3(20)-x3(10))/(x(20)-x(10));
b3 = -a3*x(10) + x3(10);

br_gresaka_13 = 0;
br_gresaka_31 = 0;

%Racunanje broja odbiraka K1 iznad klasifikacione linije izmedju K1 i K3
for i =1:length(X1)
   y1 = X1(2,i);
   y0 = a3*X1(1,i) + b3;
   if y1>y0
       br_gresaka_13 = br_gresaka_13 + 1;
   end
end

%Racunanje broja odbiraka K3 ispod klasifikacione linije izmedju K1 i K3
for i =1:length(X3)
   y1 = X3(2,i);
   y0 = a3*X3(1,i) + b3;
   if y1<y0
       br_gresaka_31 = br_gresaka_31+1;
   end
end

Konfuziona_matrica_13 = zeros(2,2);
Konfuziona_matrica_13(1,1) = length(X1) - br_gresaka_13;
Konfuziona_matrica_13(2,2) = length(X3) - br_gresaka_31;
Konfuziona_matrica_13(1,2) = br_gresaka_13;
Konfuziona_matrica_13(2,1) = br_gresaka_31;
disp('Konfuziona matrica - metod zeljenog izlaza - klase 1 i 3');
disp(Konfuziona_matrica_13)

%% Tacka 2)
%% Generisanje dve nelinearno separabilne klase

N = 500;
phi1 = -1*rand(1,N)*2*pi;
rho1 = rand(1,N);

X1 = zeros(2,N);
X1(1,:) = rho1.*cos(phi1);
X1(2,:) = rho1.*sin(phi1);

phi2 = rand(1,N)*2*pi + 2.5;
rho2 = rand(1,N) + 2.5;
X2 = zeros(2,N);

X2(1,:) = rho2.*cos(phi2);
X2(2,:) = rho2.*sin(phi2);

figure(3)
hold all
scatter(X1(1,:),X1(2,:),'ro');
scatter(X2(1,:),X2(2,:),'bo');
axis([-4 4 -4 4]);
grid on;

%% Kvadratni klasifikator metodom zeljenog izlaza

Gama = [ones(N,1); 1.7*ones(N,1)]; % veci prioritet dajemo drugoj klasi, zato su 500 drugih odbiraka sa vrednoscu 1.7

U = [-1*ones(1,N),ones(1,N);...
    -1*X1,X2 ; -1*(X1(1,:)).^2,(X2(1,:)).^2;...
    -1*(X1(2,:)).^2,(X2(2,:)).^2;...
    -2*X1(1,:).*X1(2,:), 2*X2(1,:).*X2(2,:)];

W = inv(U*U')*U*Gama;
V0 = W(1);
V1 = W(2);
V2 = W(3);
Q11 = W(4);
Q22 = W(5);
Q12 = W(6);

%diskriminaciona kriva
x1 = -3.5:0.1:3.5;
x2 = -3.5:0.1:3.5;
h = zeros(length(x1),length(x2));

for i = 1:length(x1)
    for j = 1:length(x2)
        h(i,j) = V0+V1*x1(i)+V2*x2(j)+Q11*(x1(i))^2 + Q22*(x2(j))^2 + Q12*x1(i)*x2(j);
    end
end
                
figure(3)
contour(x1,x2,h,[0 0],'g','LineWidth',3);
%legend('Klasa 1','Klasa 2','Kvadratni klasifikator - metoda željenog izlaza','Interpreter','latex');
axis([-4 4 -4 4]);