clear 
close all 
clc

%% Ucitavanje slova
PO_slova = load('PO_slova.mat');

% Izbor prvih 10 slova po 100 odbiraka

A = PO_slova.a(:,1:100);
B = PO_slova.b(:,1:100);
C = PO_slova.c(:,1:100);
D = PO_slova.d(:,1:100);
E = PO_slova.e(:,1:100);
G = PO_slova.g(:,1:100);
H = PO_slova.h(:,1:100);
I = PO_slova.i(:,1:100);
M = PO_slova.m(:,1:100);
N = PO_slova.n(:,1:100);

% Biranje jednog primerka od svakog od slova

A1 = A{45};
B1 = B{45};
C1 = C{45};
D1 = D{45};
E1 = E{45};
G1 = G{45};
H1 = H{45};
I1 = I{45};
M1 = M{45};
N1 = N{45};

[A1_vx,A1_vy,A1_p,A1_x_poz,A1_y_poz] = ekstrakcija_podataka(A1);
[B1_vx,B1_vy,B1_p,B1_x_poz,B1_y_poz] = ekstrakcija_podataka(B1);
[C1_vx,C1_vy,C1_p,C1_x_poz,C1_y_poz] = ekstrakcija_podataka(C1);
[D1_vx,D1_vy,D1_p,D1_x_poz,D1_y_poz] = ekstrakcija_podataka(D1);
[E1_vx,E1_vy,E1_p,E1_x_poz,E1_y_poz] = ekstrakcija_podataka(E1);
[G1_vx,G1_vy,G1_p,G1_x_poz,G1_y_poz] = ekstrakcija_podataka(G1);
[H1_vx,H1_vy,H1_p,H1_x_poz,H1_y_poz] = ekstrakcija_podataka(H1);
[I1_vx,I1_vy,I1_p,I1_x_poz,I1_y_poz] = ekstrakcija_podataka(I1);
[M1_vx,M1_vy,M1_p,M1_x_poz,M1_y_poz] = ekstrakcija_podataka(M1);
[N1_vx,N1_vy,N1_p,N1_x_poz,N1_y_poz] = ekstrakcija_podataka(N1);

%% Iscrtavanje slova 

figure(1)                   
subplot(2,2,[1 3]);
scatter(A1_x_poz,A1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(A1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(A1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo A','Interpreter','latex'); 

figure(2)
subplot(2,2,[1 3]);
scatter(B1_x_poz,B1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(B1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(B1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo B','Interpreter','latex'); 

figure(3)
subplot(2,2,[1 3]);
scatter(C1_x_poz,C1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(C1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(C1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo C','Interpreter','latex'); 

figure(4)
subplot(2,2,[1 3]);
scatter(D1_x_poz,D1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(D1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(D1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo D','Interpreter','latex'); 

figure(5)
subplot(2,2,[1 3]);
scatter(E1_x_poz,E1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(E1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(E1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo E','Interpreter','latex'); 

figure(6)
subplot(2,2,[1 3]);
scatter(G1_x_poz,G1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
title('Brzina po $x$ osi','Interpreter','latex');
plot(G1_vx);
grid on;
subplot(2,2,4)
plot(G1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo G','Interpreter','latex'); 

figure(7)
subplot(2,2,[1 3]);
scatter(H1_x_poz,H1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(H1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(H1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo H','Interpreter','latex'); 

figure(8)
subplot(2,2,[1 3]);
scatter(I1_x_poz,I1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(I1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(I1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo I','Interpreter','latex'); 

figure(9)
subplot(2,2,[1 3]);
scatter(M1_x_poz,M1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(M1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(M1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo M','Interpreter','latex'); 

figure(10)
subplot(2,2,[1 3]);
scatter(N1_x_poz,N1_y_poz)
grid on;
title('Pozicija po $x$ i $y$ osi','Interpreter','latex');
subplot(2,2,2)
plot(N1_vx);
grid on;
title('Brzina po $x$ osi','Interpreter','latex');
subplot(2,2,4)
plot(N1_vy);
grid on;
title('Brzina po $y$ osi','Interpreter','latex');
sgtitle('Slovo N','Interpreter','latex'); 

%% Ekstrakcija obelezja
% Izabrao sam 100 primeraka svakog slova

broj_obelezja = 7;
A_obelezja = zeros(broj_obelezja,100);
B_obelezja = zeros(broj_obelezja,100);
C_obelezja = zeros(broj_obelezja,100);
D_obelezja = zeros(broj_obelezja,100);
E_obelezja = zeros(broj_obelezja,100);
G_obelezja = zeros(broj_obelezja,100);
H_obelezja = zeros(broj_obelezja,100);
I_obelezja = zeros(broj_obelezja,100);
M_obelezja = zeros(broj_obelezja,100);
N_obelezja = zeros(broj_obelezja,100);

for i = 1:100
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(A{i});
    A_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(B{i});
    B_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(C{i});
    C_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(D{i});
    D_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(E{i});
    E_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(G{i});
    G_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(H{i});
    H_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(I{i});
    I_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(M{i});
    M_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
    [o1,o2,o3,o4,o5,o6,o7] = ekstrakcija_obelezja(N{i});
    N_obelezja(1:broj_obelezja,i) = [o1,o2,o3,o4,o5,o6,o7];
end

%% Projektovanje klasifikatora
%Usrednjavanje po vrstama i trazenje kovarijacionih matrica

M_A = mean(A_obelezja,2); KovA = cov(A_obelezja');
M_B = mean(B_obelezja,2); KovB = cov(B_obelezja');
M_C = mean(C_obelezja,2); KovC = cov(C_obelezja');
M_D = mean(D_obelezja,2); KovD = cov(D_obelezja');
M_E = mean(E_obelezja,2); KovE = cov(E_obelezja');
M_G = mean(G_obelezja,2); KovG = cov(G_obelezja');
M_H = mean(H_obelezja,2); KovH = cov(H_obelezja');
M_I = mean(I_obelezja,2); KovI = cov(I_obelezja');
M_M = mean(M_obelezja,2); KovM = cov(M_obelezja');
M_N = mean(N_obelezja,2); KovN = cov(N_obelezja');

% Redukcija dimenzija - LDA metod

Sw = 1/10*(KovA + KovB + KovC + KovD + KovE + KovG + KovH + KovI + KovM + KovN);
M_0 = 1/10*(M_A + M_B + M_C + M_D + M_E + M_G + M_H + M_I + M_M + M_N);

Sb = 1/10*(M_A-M_0)*(M_A-M_0)' + 1/10*(M_B-M_0)*(M_B-M_0)' + 1/10*(M_C-M_0)*(M_C-M_0)'...
    + 1/10*(M_D-M_0)*(M_D-M_0)' + 1/10*(M_E-M_0)*(M_E-M_0)' + 1/10*(M_G-M_0)*(M_G-M_0)'...
    + 1/10*(M_H-M_0)*(M_H-M_0)' + 1/10*(M_I-M_0)*(M_I-M_0)' + 1/10*(M_M-M_0)*(M_M-M_0)'...
    + 1/10*(M_N-M_0)*(M_N-M_0)';

J = Sw^(-1)*Sb;
[F,L] = eig(J);

transformaciona_matrica = F(:,1:3);

mat_Y1 = A_obelezja'*transformaciona_matrica;
mat_Y2 = B_obelezja'*transformaciona_matrica;
mat_Y3 = C_obelezja'*transformaciona_matrica;
mat_Y4 = D_obelezja'*transformaciona_matrica;
mat_Y5 = E_obelezja'*transformaciona_matrica;
mat_Y6 = G_obelezja'*transformaciona_matrica;
mat_Y7 = H_obelezja'*transformaciona_matrica;
mat_Y8 = I_obelezja'*transformaciona_matrica;
mat_Y9 = M_obelezja'*transformaciona_matrica;
mat_Y10 = N_obelezja'*transformaciona_matrica;

figure(11)
hold on
scatter3(mat_Y1(:,1),mat_Y1(:,2),mat_Y1(:,3),'*');
scatter3(mat_Y2(:,1),mat_Y2(:,2),mat_Y2(:,3),'o');
scatter3(mat_Y3(:,1),mat_Y3(:,2),mat_Y3(:,3),'x');
scatter3(mat_Y4(:,1),mat_Y4(:,2),mat_Y4(:,3),'s');
scatter3(mat_Y5(:,1),mat_Y5(:,2),mat_Y5(:,3),'d'); 
scatter3(mat_Y6(:,1),mat_Y6(:,2),mat_Y6(:,3),'.'); 
scatter3(mat_Y7(:,1),mat_Y7(:,2),mat_Y7(:,3),'v');
scatter3(mat_Y8(:,1),mat_Y8(:,2),mat_Y8(:,3),'p');
scatter3(mat_Y9(:,1),mat_Y9(:,2),mat_Y9(:,3),'h');
scatter3(mat_Y10(:,1),mat_Y10(:,2),mat_Y10(:,3),'k<');
hold off
grid on;
legend('A','B','C','D','E','G','H','I','M','N');

%% Izracunavanja srednjih vrednosti i kovarijacionih matrica redukovanih obelezja

MA = mean(mat_Y1); KovA = cov(mat_Y1);
MB = mean(mat_Y2); KovB = cov(mat_Y2);
MC = mean(mat_Y3); KovC = cov(mat_Y3);
MD = mean(mat_Y4); KovD = cov(mat_Y4);
ME = mean(mat_Y5); KovE = cov(mat_Y5);
MG = mean(mat_Y6); KovG = cov(mat_Y6);
MH = mean(mat_Y7); KovH = cov(mat_Y7);
MI = mean(mat_Y8); KovI = cov(mat_Y8);
MM = mean(mat_Y9); KovM = cov(mat_Y9);
MN = mean(mat_Y10); KovN = cov(mat_Y10);

%% Testiranje hipoteze
klase = ['A','B','C','D','E','G','H','I','M','N'];

Konfuziona_matrica = zeros(10,10);

for k = 1:10
   if klase(k) == 'A' 
       mat = mat_Y1;
   elseif klase(k) =='B'
       mat = mat_Y2;
   elseif klase(k) =='C'
       mat = mat_Y3;
   elseif klase(k) =='D'
       mat = mat_Y4;
   elseif klase(k) =='E'
       mat = mat_Y5;
   elseif klase(k) =='G'
       mat = mat_Y6;
   elseif klase(k) =='H'
       mat = mat_Y7;
   elseif klase(k) =='I'
       mat = mat_Y8;
   elseif klase(k) =='M'
       mat = mat_Y9;
   elseif klase(k) =='N'
       mat = mat_Y10;
   end
    
    for i = 1:100
        x = mat(i,:);
        fa = 1/((2*pi)^(3/2)*det(KovA)^0.5)*exp(-1/2 * (x-MA)*inv(KovA)*(x-MA)');
        fb = 1/((2*pi)^(3/2)*det(KovB)^0.5)*exp(-1/2 * (x-MB)*inv(KovB)*(x-MB)');
        fc = 1/((2*pi)^(3/2)*det(KovC)^0.5)*exp(-1/2 * (x-MC)*inv(KovC)*(x-MC)');
        fd = 1/((2*pi)^(3/2)*det(KovD)^0.5)*exp(-1/2 * (x-MD)*inv(KovD)*(x-MD)');
        fe = 1/((2*pi)^(3/2)*det(KovE)^0.5)*exp(-1/2 * (x-ME)*inv(KovE)*(x-ME)');
        fg = 1/((2*pi)^(3/2)*det(KovG)^0.5)*exp(-1/2 * (x-MG)*inv(KovG)*(x-MG)');
        fh = 1/((2*pi)^(3/2)*det(KovH)^0.5)*exp(-1/2 * (x-MH)*inv(KovH)*(x-MH)');
        fi = 1/((2*pi)^(3/2)*det(KovI)^0.5)*exp(-1/2 * (x-MI)*inv(KovI)*(x-MI)');
        fm = 1/((2*pi)^(3/2)*det(KovM)^0.5)*exp(-1/2 * (x-MM)*inv(KovM)*(x-MM)');
        fn = 1/((2*pi)^(3/2)*det(KovN)^0.5)*exp(-1/2 * (x-MN)*inv(KovN)*(x-MN)');
        
        najveci_f = max([fa,fb,fc,fd,fe,fg,fh,fi,fm,fn]);
        
        if najveci_f == fa
            Konfuziona_matrica(k,1)= Konfuziona_matrica(k,1)+1;
        elseif najveci_f == fb
            Konfuziona_matrica(k,2)= Konfuziona_matrica(k,2)+1;
        elseif najveci_f == fc
            Konfuziona_matrica(k,3)= Konfuziona_matrica(k,3)+1; 
        elseif najveci_f == fd
            Konfuziona_matrica(k,4)= Konfuziona_matrica(k,4)+1;
        elseif najveci_f == fe
            Konfuziona_matrica(k,5)= Konfuziona_matrica(k,5)+1;
        elseif najveci_f == fg
            Konfuziona_matrica(k,6)= Konfuziona_matrica(k,6)+1;
        elseif najveci_f == fh
            Konfuziona_matrica(k,7)= Konfuziona_matrica(k,7)+1;
        elseif najveci_f == fi
            Konfuziona_matrica(k,8)= Konfuziona_matrica(k,8)+1;
        elseif najveci_f == fm
            Konfuziona_matrica(k,9)= Konfuziona_matrica(k,9)+1;
        elseif najveci_f == fn
            Konfuziona_matrica(k,10)= Konfuziona_matrica(k,10)+1;
        end
    end
end

disp(Konfuziona_matrica)
ukupna_greska = (sum(sum(Konfuziona_matrica))-trace(Konfuziona_matrica))/sum(sum(Konfuziona_matrica));
disp(ukupna_greska)

%% Primer tacno i netacno klasifikovanog slova

indeks_pogresnog = zeros(1,100);
k = 1;
d_pogresno = 0;
e_pogresno = 0;
h_pogresno = 0;

 for i = 1:100
        x = mat_Y10(i,:); % slovo N
        fa = 1/((2*pi)^(3/2)*det(KovA)^0.5)*exp(-1/2 * (x-MA)*inv(KovA)*(x-MA)');
        fb = 1/((2*pi)^(3/2)*det(KovB)^0.5)*exp(-1/2 * (x-MB)*inv(KovB)*(x-MB)');
        fc = 1/((2*pi)^(3/2)*det(KovC)^0.5)*exp(-1/2 * (x-MC)*inv(KovC)*(x-MC)');
        fd = 1/((2*pi)^(3/2)*det(KovD)^0.5)*exp(-1/2 * (x-MD)*inv(KovD)*(x-MD)');
        fe = 1/((2*pi)^(3/2)*det(KovE)^0.5)*exp(-1/2 * (x-ME)*inv(KovE)*(x-ME)');
        fg = 1/((2*pi)^(3/2)*det(KovG)^0.5)*exp(-1/2 * (x-MG)*inv(KovG)*(x-MG)');
        fh = 1/((2*pi)^(3/2)*det(KovH)^0.5)*exp(-1/2 * (x-MH)*inv(KovH)*(x-MH)');
        fi = 1/((2*pi)^(3/2)*det(KovI)^0.5)*exp(-1/2 * (x-MI)*inv(KovI)*(x-MI)');
        fm = 1/((2*pi)^(3/2)*det(KovM)^0.5)*exp(-1/2 * (x-MM)*inv(KovM)*(x-MM)');
        fn = 1/((2*pi)^(3/2)*det(KovN)^0.5)*exp(-1/2 * (x-MN)*inv(KovN)*(x-MN)');
        
        najveci_f = max([fa,fb,fc,fd,fe,fg,fh,fi,fm,fn]);
        
        if najveci_f == fd
            indeks_pogresnog(k) = i;
            d_pogresno = i;
            k = k+1;
        elseif najveci_f == fe
            indeks_pogresnog(k) = i;
            e_pogresno = i;
            k = k+1;
        elseif najveci_f == fh
            indeks_pogresnog(k) = i;
            h_pogresno = i;
            k = k+1;
        end
end

N_kao_D = N{d_pogresno};
[N_kao_D_vx,N_kao_D_vy,N_kao_D_p,N_kao_D_x_poz,N_kao_D_y_poz] = ekstrakcija_podataka(N_kao_D);
N_kao_E = N{e_pogresno};
[N_kao_E_vx,N_kao_E_vy,N_kao_E_p,N_kao_E_x_poz,N_kao_E_y_poz] = ekstrakcija_podataka(N_kao_E);
N_kao_H = N{h_pogresno};
[N_kao_H_vx,N_kao_H_vy,N_kao_H_p,N_kao_H_x_poz,N_kao_H_y_poz] = ekstrakcija_podataka(N_kao_H);

figure(12)
subplot(2,3,[1 4]);
scatter(N_kao_D_x_poz,N_kao_D_y_poz)
title('Klasifikovano kao D','Interpreter','latex');
grid on;
subplot(2,3,[2 5]);
scatter(N_kao_E_x_poz,N_kao_E_y_poz)
title('Klasifikovano kao E','Interpreter','latex');
grid on;
subplot(2,3,[3 6]);
scatter(N_kao_H_x_poz,N_kao_H_y_poz)
title('Klasifikovano kao H','Interpreter','latex');
grid on;
sgtitle('Slovo N','Interpreter','latex'); 

%% Odabir dva slova i dva obelezja za klasifikaciju

% Slova B i G

figure(13)
hold on
scatter(mat_Y2(:,1),mat_Y2(:,2),'ro');
scatter(mat_Y6(:,1),mat_Y6(:,2),'bo');
legend('Slovo B','Slovo G','Interpreter','latex');
grid on;

%% Linearni klasifikator - druga numericka metoda

N1 = length(mat_Y2(:,1));
N2 = length(mat_Y6(:,1));
mat_Y2 = mat_Y2(:,1:2)';
mat_Y6 = mat_Y6(:,1:2)';
M1_est = mean(mat_Y2,2);
M2_est = mean(mat_Y6,2);
S1_est = cov(mat_Y2');
S2_est = cov(mat_Y6');

s = 0:1e-3:1;
V0_opt_s = []; %optimalno vo za fiksirano s
Neps_s = []; %odgovarajuca greska klasifikacije

for i = 1:length(s)
    %pronalazimo odgovarajuce V
    V = ((s(i)*S1_est + (1-s(i))*S2_est)^(-1))*(M2_est-M1_est);
    %projektujemo odbirke klase na pravac V
    Y1 = V'*mat_Y2;
    Y2 = V'*mat_Y6;
    Y = [Y1 Y2];
    Y = sort(Y);
    %petlja po V0 + brojanje pogresno klasifikovanih
    V0 = [];
    Neps = []; %za jedno s i svako V0 pamtim gresku
    
    for j = 1:(length(Y)-1)
        V0(j) = -(Y(j)+Y(j+1))/2; % tacno izmedju dve susedne projekcije
        Neps(j) = 0;
        for k = 1:N1
            if Y1(k) >-V0(j)
                Neps(j) = Neps(j)+1;
            end
            if Y2(k) < - V0(j)
                Neps(j) = Neps(j)+1;
            end
        end
    end
    
    [Neps_s(i),index] = min(Neps);
    V0_opt_s(i) = V0(index);

end

[Neps_opt, index] = min(Neps_s);
V0_opt = V0_opt_s(index);
s_opt = s(index);   
    
%% Iscrtavanje klasifikacione linije
    
V = (s_opt*S1_est + (1-s_opt)*S2_est)^(-1)*(M2_est - M1_est);

x1 =  0.1:0.01:0.9;
x2 = -(V0_opt +V(1)*x1)/V(2);
    
figure(13)
plot(x1,x2,'k');
hold off
legend('Slovo B','Slovo G','Klasifikaciona linija','Location','NorthWest', 'Interpreter', 'latex');
