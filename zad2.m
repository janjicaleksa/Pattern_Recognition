clear all
close all
clc

%% Biranje parametara za generisanje raspodela

rng(100)
M11 = [6; 8]; S11 = [0.7 -0.5; -0.5 0.7];
M12 = [5; 5]; S12 = [0.8 0.5; 0.5 0.8];
M21 = [1; 3]; S21 = [1.2 -0.8; -0.8 1.2];
M22 = [0.5; 1.5]; S22 = [0.8 0.5; 0.5 0.8];

N = 500;

% Generisanje dve polimodalne Gausove raspodele 
P11 = 0.6; P12 = 0.4;
P21 = 0.3; P22 = 0.7;

K11 = mvnrnd(M12,S12,N); K12 = mvnrnd(M12,S12,N);
K21 = mvnrnd(M21,S21,N); K22 = mvnrnd(M22,S22,N);
pom1 = rand(N,1);
pom2 = rand(N,1);

K1 = (pom1 < P11).*K11 + (pom1 >=P11).*K12;
K2 = (pom2 < P21).*K21 + (pom2 >=P21).*K22;

% Prikazivanje odbiraka na dijagramu
figure(1)
hold all
scatter(K1(:,1),K1(:,2),'ro');
scatter(K2(:,1),K2(:,2),'bo');
legend('Klasa 1','Klasa 2','Location','NorthWest','Interpreter','latex');
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
grid on

%% Funkcije gustine verovatnoce

x = -1:0.1:9;
y = -1:0.1:12;
f1 = zeros(length(x),length(y));
f2 = zeros(length(x),length(y));

const1 = 1/(2*pi*det(S11)^(0.5));
const2 = 1/(2*pi*det(S12)^(0.5));
const3 = 1/(2*pi*det(S21)^(0.5));
const4 = 1/(2*pi*det(S22)^(0.5));

for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i) y(j)]';
        f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
        f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
        f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
        f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
        f1(i,j) = P11*f11 + P12*f12;
        f2(i,j) = P21*f21 + P22*f22;
   end
end
fgv1 = f1;
fgv2 = f2;

% Prikaz funkcija gustina verovatnoca

figure(2)
hold on
mesh(x, y, fgv1');
mesh(x, y, fgv2');
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
colorbar;
axis([-1 9 -1 10])
hold off

% Prikaz histograma

figure(3)
histogram2(K1,K2,40,'Normalization','pdf')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
set(gca,'YTickLabel',[0 2 4 6 8 10]);
set(gca,'XTickLabel',[0 2 4 6 8]);

%% Bajesov klasifikator minimalne verovatnoce greske

h = -log(f1./f2);
figure(1)
hold on
contour(x, y, h', [0,0], 'k')
legend('Klasa 1','Klasa 2','BKMVG','Location','NorthWest','Interpreter','latex')
axis([-4 8 -2 8])

%Spajamo sve odbirke u jedan vektor, njihove stvarne klase znamo a
%klasifikator nam daje predikciju

K1 = K1';
K2 = K2';
Xs = [K1,K2];
X_true = [ones(1,N), ones(1,N)*2]; %stvarne klase
X_pred = zeros(size(X_true));

for i = 1:length(Xs)
    X = Xs(:,i);
    f11 = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12 = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21 = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22 = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f2 = P21*f21 + P22*f22;
    f1 = P11*f11 + P12*f12;
    if (i<=500)
        h1(i) = -log(f1/f2);
    elseif (i>500)
        h2(i-N) = -log(f1/f2);
    end
    if (f1>f2)
        X_pred(i) = 1;
    else
        X_pred(i) = 2;
    end
end

Konf_matrica = confusionmat(X_true,X_pred)
 
err1 = Konf_matrica(2,1)/N; %procenat pogresno klasifikovanih iz prve klase
err2 = Konf_matrica(1,2)/N;
err = 1/2*(err1+err2);

disp(['Verovatnoca greske 1. tipa: ' num2str(err1*100),'%'])
disp(['Verovatnoca greske 2. tipa: ' num2str(err2*100),'%'])

% Teorijski pristup

f1t = zeros(length(x),length(y));
f2t = zeros(length(x),length(y));

eps1 = 0;
eps2 = 0;

for i = 1:length(x)
    for j = 1:length(y)
        X = [x(i) ; y(j)];
        f11t = const1*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
        f12t = const2*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
        f21t = const3*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
        f22t = const4*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
        f1t(i,j) = P11*f11t + P12*f12t;
        f2t(i,j) = P21*f21t + P22*f22t;
        h = -log(f1t(i,j)/f2t(i,j));
        if h<0 % u oblasti L1
            eps2 = eps2 + 0.1*0.1*f2t(i,j);
        else % u oblasti L2
            eps1 = eps1 + 0.1*0.1*f1t(i,j);
        end

    end
end
 
disp(['Teorijska verovatnoca greske 1. tipa: ' num2str(eps1*100),'%'])
disp(['Teorijska verovatnoca greske 2. tipa: ' num2str(eps2*100),'%'])

%% Bajesov klasifikator minimalne cene

f1 = fgv1;
f2 = fgv2;
c11 = 0; c22 = 0; % Ne penalizujemo tacnu klasifikaciju
c12 = 1; c21 = 5; % Vise penalizujemo pogresnu klasifickaicju odbiraka iz K1

t = (c12-c22)/(c21-c11);
k = f1./f2;
figure(1);
contour(x,y,k',[t t],'g','LineWidth',1.5); %gledamo gde je k jednako pragu
legend('Klasa 1','Klasa 2','BKMVG','BKMC','Location','NorthWest','Interpreter','latex')
axis([-4 8 -2 8])
hold off

%% Neyman-Pearson-ov test

h = -log(f1./f2);
br = 0;

for mi= 0.1:0.1:1
    br = br + 1;
    Eps2(br) = 0;
    for i = 1:length(x)-1
        for j = 1:length(y)-1
            if (h(i,j) < -log(mi))
                Eps2(br) = Eps2(br) + 0.01*f2(i,j);
            end
        end
    end
end
figure(6)
hold all
plot(0.1:0.1:1,Eps2);
plot(0.1:0.1:1, err2*ones(length(0.1:0.1:1)), 'r')
grid on;
xlabel('$\mu$','Interpreter','latex');
ylabel('$\epsilon_{2}$','Interpreter','latex');
legend('$\epsilon_{2}$','$\epsilon_{0}$','Interpreter','latex');
mi = 0.28; % dobijeno mi je priblizno = 0.28

%% Wald-ov test

E1 = 0.000000000001;
E2 = 0.000000000001;
A = (1-E1)/E2; a = -log(A);
B = E1/(1-E2); b = -log(B);

figure(7)
hold on
plot(1:8, zeros(8), 'k')
plot(1:8, a*ones(8), 'k--')
plot(1:8, b*ones(8), 'k--')
grid on
xlabel('n','Interpreter','latex')
ylabel('$S_{m}$','Interpreter','latex')

for j = 1:100 %Crtamo npr. oko 100 odbiraka
    n1 = round(rand*500); n2 = round(rand*500);
    Sm1(1) = h1(n1);
    Sm2(1) = h2(n2);
    i = 2;
    while ((i <= 500) && (Sm1(i-1) >= a) && (Sm1(i-1) <= b))
        Sm1(i) = Sm1(i-1) + h1(round(rand*500));
        i = i + 1;
    end

    figure(7)
    hold on
    plot(1:i-1, Sm1, 'b')
    plot(i-1, Sm1(i-1), 'bd') 
    
    i = 2;
    while ((i <= 500) && (Sm2(i-1) >= a) && (Sm2(i-1) <= b))
        Sm2(i) = Sm2(i-1) + h2(round(rand*500));
        i = i + 1;
    end
    figure(7)
    hold on
    plot(1:i-1, Sm2, 'r')
    plot(i-1, Sm2(i-1), 'rd')
    
    clear Sm1
    clear Sm2   
end

%% Odredjivanje srednjeg broja odbiraka prve klase za donosenje odluke

E1 = [10^-12, 10^-10, 10^-8, 10^-6, 10^-4,10^-3, 10^-2, 10^-1, 0.16, 0.198,0.4,0.98];
E2 = [10^-12, 10^-10, 10^-8, 10^-6, 10^-4,10^-3, 10^-2, 10^-1, 0.16, 0.198,0.4,0.98];  
x_axis= [10^-12, 10^-10, 10^-8, 10^-6, 10^-4,10^-3, 10^-2, 10^-1, 0.16, 0.198,0.4,0.98];

h1_mean = mean(h1);
m1 = (-log((1-E1)/E2(1)).*(1-E1) - log(E1/(1-E2(1))).*E1)/h1_mean;   
figure(8)
semilogx(x_axis, m1)
    
m1 = (-log((1-E1(1))./E2)*(1-E1(1)) - log(E1(1)./(1-E2))*E1(1))/h1_mean;    
figure(8)
hold on
semilogx(x_axis(1:12), m1)
ylabel('$m_{1}$','Interpreter','latex');
xlabel('$\eta_{1_/2}$','Interpreter','latex');
grid on 
legend('$\eta_{2} = const$', '$\eta_{1} = const$','Location','SouthWest','Interpreter','latex')