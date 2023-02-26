clc 
clear all 
close all 

%% Generisanje 4 linearno separabilne klase 
N=500; 

M1=[0;0];
S1=[0.5,0;0,0.5];  

M2=[0;7];
S2=[0.8,0.2;0.2,0.5];

M3=[7;7];
S3=[0.6,-0.3;-0.3,0.8];

M4=[7;0];
S4=[0.4,0.3;0.3,0.8];

X1 = mvnrnd(M1,S1,N)'; 
X2 = mvnrnd(M2,S2,N)'; 
X3 = mvnrnd(M3,S3,N)'; 
X4 = mvnrnd(M4,S4,N)';

figure(1)
hold on
plot(X1(1,:),X1(2,:),'ro')
plot(X2(1,:),X2(2,:),'bx')
plot(X3(1,:),X3(2,:),'mv')
plot(X4(1,:),X4(2,:),'yd')
grid on;
hold off
legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4','Location','best','Interpreter','latex');
%%
iter_array=[];
for iter=1:1 % promeniti granicu iter za pokretanje vise klasterizacija    
    % Pocetna klasterizacija
    pom= rand(1,4*N);
    X = [X1,X2,X3,X4]; % smestaju se svi odbirci iz sve cetiri klase u jednu klasu 
    K1=[];
    K2=[];
    K3=[];
    K4=[];
    for i=1:4*N
        if(pom(i)<0.25)
            K1=[K1 X(:,i)];
        elseif (pom(i)>= 0.25 && pom(i)<0.5)
            K2=[K2 X(:,i)];
        elseif (pom(i)>=0.5 && pom(i)<0.75)
            K3=[K3 X(:,i)];
        else
            K4=[K4 X(:,i)];
        end
    end
    
    % Prikaz pocetne klasterizacije
    figure(2)
    hold on
    plot(K1(1,:),K1(2,:),'ro')
    plot(K2(1,:),K2(2,:),'bx')
    plot(K3(1,:),K3(2,:),'mv')
    plot(K4(1,:),K4(2,:),'yd')
    grid on;
    hold off
    legend('Klasa 1','Klasa 2','Klasa 3','Klasa 4', 'Location', 'best', 'Interpreter', 'latex');
    title('Iteracija 0 - proizvoljna raspodela','Interpreter','latex')
    
    %% Procena parametara
    N1= length(K1(1,:));
    N2= length(K2(1,:));
    N3= length(K3(1,:));
    N4= length(K4(1,:));

    M1=mean(K1,2);
    M2=mean(K2,2);
    M3=mean(K3,2);
    M4=mean(K4,2);

    % Reklasifikacija

    lmax=100;
    l=1;
    reklas=1;
    while((l<lmax) && (reklas==1))
        K1_pom=[];
        K2_pom=[];
        K3_pom=[];
        K4_pom=[];
        reklas=0;

        for i = 1:N1
            %racunamo rastojanje i trazimo minimum
            d1 = sum((K1(:, i) - M1).^2); 
            d2 = sum((K1(:, i) - M2).^2); 
            d3 = sum((K1(:, i) - M3).^2); 
            d4 = sum((K1(:, i) - M4).^2);
            d_min=min([d1,d2,d3,d4]);
            if d_min==d1
                K1_pom = [K1_pom K1(:, i)];
            elseif d_min == d2
                K2_pom = [K2_pom K1(:, i)];
                reklas=1;
            elseif d_min== d3
                K3_pom = [K3_pom K1(:, i)];
                reklas=1;
            elseif d_min== d4
                K4_pom = [K4_pom K1(:, i)];
                reklas=1;
            end
        end

        for i = 1:N2
            %racunamo rastojanje i trazimo minimum
            d1 = sum((K2(:, i) - M1).^2); 
            d2 = sum((K2(:, i) - M2).^2); 
            d3 = sum((K2(:, i) - M3).^2); 
            d4 = sum((K2(:, i) - M4).^2);
            d_min=min([d1,d2,d3,d4]);
            if d_min==d1
                K1_pom = [K1_pom K2(:, i)];
                reklas=1;
            elseif d_min == d2
                K2_pom = [K2_pom K2(:, i)];
            elseif d_min== d3
                K3_pom = [K3_pom K2(:, i)];
                reklas=1;
            elseif d_min== d4
                K4_pom = [K4_pom K2(:, i)];
                reklas=1;
            end
        end

        for i = 1:N3
            %racunamo rastojanje i trazimo minimum
            d1 = sum((K3(:, i) - M1).^2); 
            d2 = sum((K3(:, i) - M2).^2); 
            d3 = sum((K3(:, i) - M3).^2); 
            d4 = sum((K3(:, i) - M4).^2);
            d_min=min([d1,d2,d3,d4]);
            if d_min==d1
                K1_pom = [K1_pom K3(:, i)];
                reklas=1;
            elseif d_min == d2
                K2_pom = [K2_pom K3(:, i)];
                reklas=1;
            elseif d_min== d3
                K3_pom = [K3_pom K3(:, i)];
            elseif d_min== d4
                K4_pom = [K4_pom K3(:, i)];
                reklas=1;
            end
        end

        for i = 1:N4
            %racunamo rastojanje i trazimo minimum
            d1 = sum((K4(:, i) - M1).^2); 
            d2 = sum((K4(:, i) - M2).^2); 
            d3 = sum((K4(:, i) - M3).^2); 
            d4 = sum((K4(:, i) - M4).^2);
            d_min=min([d1,d2,d3,d4]);
            if d_min==d1
                K1_pom = [K1_pom K4(:, i)];
                reklas=1;
            elseif d_min == d2
                K2_pom = [K2_pom K4(:, i)];
                reklas=1;
            elseif d_min== d3
                K3_pom = [K3_pom K4(:, i)];
                reklas=1;
            elseif d_min== d4
                K4_pom = [K4_pom K4(:, i)];
            end
        end
        
        maxOdb=max([length(K1_pom),length(K2_pom),length(K3_pom),length(K4_pom)]);
        if(isempty(K1_pom))
            if(maxOdb==length(K2_pom))
                K1_pom=K2_pom(:,1:round(maxOdb/2));
                K2_pom=K2_pom(:,round(maxOdb/2)+1:maxOdb);
            else if (maxOdb==length(K3_pom))
                K1_pom=K3_pom(:,1:round(maxOdb/2));
                K3_pom=K3_pom(:,round(maxOdb/2)+1:maxOdb);
            else 
                K1_pom=K4_pom(:,1:round(maxOdb/2));
                K4_pom=K4_pom(:,round(maxOdb/2)+1:maxOdb);
            end
            end
        end
        if(isempty(K2_pom))
            if(maxOdb==length(K1_pom))
                K2_pom=K1_pom(:,1:round(maxOdb/2));
                K1_pom=K1_pom(:,round(maxOdb/2)+1:maxOdb);
            else if (maxOdb==length(K3_pom))
                K2_pom=K3_pom(:,1:round(maxOdb/2));
                K3_pom=K3_pom(:,round(maxOdb/2)+1:maxOdb);
            else 
                K2_pom=K4_pom(:,1:round(maxOdb/2));
                K4_pom=K4_pom(:,round(maxOdb/2)+1:maxOdb);
            end
            end
        end
        
        if(isempty(K3_pom))
            if(maxOdb==length(K1_pom))
                K3_pom=K1_pom(:,1:round(maxOdb/2));
                K1_pom=K1_pom(:,round(maxOdb/2)+1:maxOdb);
            else if (maxOdb==length(K2_pom))
                K3_pom=K2_pom(:,1:round(maxOdb/2));
                K2_pom=K2_pom(:,round(maxOdb/2)+1:maxOdb);
            else 
                K3_pom=K4_pom(:,1:round(maxOdb/2));
                K4_pom=K4_pom(:,round(maxOdb/2)+1:maxOdb);
            end
            end
        end
        
        if(isempty(K4_pom))
            if(maxOdb==length(K1_pom))
                K4_pom=K1_pom(:,1:round(maxOdb/2));
                K1_pom=K1_pom(:,round(maxOdb/2)+1:maxOdb);
            else if (maxOdb==length(K2_pom))
                K4_pom=K2_pom(:,1:round(maxOdb/2));
                K2_pom=K2_pom(:,round(maxOdb/2)+1:maxOdb);
            else 
                K4_pom=K3_pom(:,1:round(maxOdb/2));
                K3_pom=K3_pom(:,round(maxOdb/2)+1:maxOdb);
            end
            end
        end
        
        clear K1 K2 K3 K4
        K1=K1_pom;
        K2=K2_pom;
        K3=K3_pom;
        K4=K4_pom;
        
%         figure()
%         hold on
%         plot(K1(1,:),K1(2,:),'ro');
%         plot(K2(1,:),K2(2,:),'bx');
%         plot(K3(1,:),K3(2,:),'mv');
%         plot(K4(1,:),K4(2,:),'yd');
%         grid on;
%         hold off
%         legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4', 'Location', 'Best', 'Interpreter', 'latex')
%         title(['Iteracija broj ' num2str(l)], 'Interpreter', 'latex');
%         pause;
        
        N1= length(K1(1,:));
        N2= length(K2(1,:));
        N3= length(K3(1,:));
        N4= length(K4(1,:));

        M1=mean(K1,2);
        M2=mean(K2,2);
        M3=mean(K3,2);
        M4=mean(K4,2);
        
        if reklas==0
            iter_array=[iter_array l-1];
        end
        
        l=l+1;
    end
    
    figure()
    hold on
    plot(K1(1,:),K1(2,:),'ro');
    plot(K2(1,:),K2(2,:),'bx');
    plot(K3(1,:),K3(2,:),'mv');
    plot(K4(1,:),K4(2,:),'yd');
    grid on;
    hold off
    legend('Klasa 1', 'Klasa 2', 'Klasa 3', 'Klasa 4', 'Interpreter', 'latex')
    title(['Iteracija broj ' num2str(l-2)],'Interpreter', 'latex');
    iter_array=[iter_array l];
    %{
    if(flag)
        figure()
        hold on
        plot(K1(1,:),K1(2,:),'ro');
        plot(K2(1,:),K2(2,:),'bx');
        plot(K3(1,:),K3(2,:),'mv');
        plot(K4(1,:),K4(2,:),'yd');
        hold off
        legend('K_1', 'K_2', 'K_3', 'K_4')
        title(['Iteracija broj' num2str(l)]);
        iter_array=[iter_array l];
    end
    %}
    %}
end

mean_iter = sum(iter_array)/length(iter_array);
disp(['Prosecan broj iteracija je: ', num2str(mean_iter)]);