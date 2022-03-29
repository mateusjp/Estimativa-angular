%% Implementacao do Algoritmo descrito no artigo - Seel (2014)
clc 
close all
clear all

%% Adicionando dados

%nova_coleta = load('data_samuel.mat');

movimento = load('dadosBTSEXOIMU.mat');

% Amostragem
dt = 1/50;
%% Criando variáveis

% Selecao de método

metodo = 1; % Médodo 1 = Gauss, Método = 2 Kalman
fk = 1;
% Encoder

enc = movimento.angExo(:,2);

% Definindo o numero de segmentos
seg = 2;

% Definindo quais segmentos

% No caso: Coxa e Canela
primeiro_seg = 2;
ultimo_seg   = primeiro_seg + (seg -1);

% Acelerometro
acel = cell(seg, 1);

% Giro
giro = cell(seg, 1);

%% Estimativa do vetor
j1Hatv = zeros(1,3);
j2Hatv = zeros(1,3);

% Valores para analise e plot das superfícies

step_batch_v = [10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 500];
ordem_v = [1, 2, 3, 4, 6, 8, 10];

% Valor para plot da articulação

step_batch_v = [200];
ordem_v = [6];

rmse_j = zeros(length(step_batch_v),length(ordem_v));
rmse_o = zeros(length(step_batch_v),length(ordem_v));
cor_v = zeros(length(step_batch_v),length(ordem_v));
for b=1:length(step_batch_v) 
    for o =1:length(ordem_v)
        step_batch = step_batch_v(b);
        ordem = ordem_v(o);
        phi1   =  0.5; %unifrnd(-pi/2, pi/2, 1); % Phi 1
        theta1 =  0.5; %unifrnd(0, 2*pi, 1); % Theta 1
        phi2   =  0.5; %unifrnd(-pi/2, pi/2, 1); % Phi 2
        theta2 =  0.5; %unifrnd(0, 2*pi, 1); % Theta 2
        x_j    =  [phi1 theta1 phi2 theta2]'; % Vetor de variáveis


        o1Hat = zeros(3,1); %unifrnd(-1,1,[3,1]);
        o2Hat = zeros(3,1); %unifrnd(-1,1,[3,1]);

        x_o = [o1Hat; o2Hat];

        alpha_acc_giro = zeros(1,length(enc));
        alpha_acc_giro_corr = zeros(1,length(enc));
        alpha_giro = zeros(1,length(enc));
        alpha_acc = zeros(1,length(enc));

        % Normal

        cont_seg = 1;
        for i = primeiro_seg : ultimo_seg % 1 - Body, 2 - Thigh, 3 - Shank, 4 - Foot
            acel{cont_seg} = movimento.dataIMU_P{i}.outAcce;
            giro{cont_seg} = movimento.dataIMU_P{i}.outGyro;
            cont_seg = cont_seg+1;
        end

        %%%%%% Disturbio %%%%%%%
        for k = 9*50:length(enc)
            % rotacionando somente o sensor da coxa
            rot_ruido = angle2dcm(0,50*pi/180,0); 
            giro{1}(k,:) = (rot_ruido*giro{1}(k,:)')'; 
            acel{1}(k,:) = (rot_ruido*acel{1}(k,:)')';    
        end
        cont_seg = 1;

        int_aux = 0;
        lambda = 0.01;

        for i = 1:seg-1
            cont_batch = 1;
            for k=3:length(enc(:,1))-2
                    if cont_batch == step_batch
                        cont_batch = 1;
                        % Selecione os valores referentes ao segmento
                        g1 = giro{i}(k-(step_batch-1)-2:k+2,:);
                        g2 = giro{i+1}(k-(step_batch-1)-2:k+2,:);
                        a1 = acel{i}(k-(step_batch-1)-2:k+2,:);
                        a2 = acel{i+1}(k-(step_batch-1)-2:k+2,:);

                        % Calc da derivada de g
                        g1Dot =  (g1(1:step_batch,:)-8.*g1(2:step_batch+1,:)+8.*g1(4:step_batch+3,:)-g1(5:step_batch+4,:))/(12*dt);
                        g2Dot =  (g2(1:step_batch,:)-8.*g2(2:step_batch+1,:)+8.*g2(4:step_batch+3,:)-g2(5:step_batch+4,:))/(12*dt);
                        g1 = g1(3:length(g1)-2,:);
                        g2 = g2(3:length(g2)-2,:);

                        % Calc J1 e J2
                        [j1Hat,j2Hat,x_j,rmse_j(b,o)] = estimative_j_kalman(x_j, g1, g2, step_batch, ordem);
                        % Cal O1 e O2
                        [o1,o2,~,rmse_o(b,o)] = estimative_o_kalman(x_o, g1, g2, g1Dot, g2Dot, j1Hat, j2Hat, a1, a2, step_batch, ordem);     

                        % c é qualquer vetor que não seja paralelo a j1 e j2 - Seel 2014 - Eq. 13
                        c = [1,0,0]';

                        % Seel 2014 - Eq. 13
                        x1 = cross(j1Hat,c);
                        y1 = cross(j1Hat,x1);
                        x2 = cross(j2Hat,c);
                        y2 = cross(j2Hat,x2);
                    end
                    if k >= step_batch+2
                        %% Estimativa pelo giroscópio
                        if i == 2 && primeiro_seg == 1
                            alpha_giro(k)= int_aux + (giro{i}(k,:)*j1Hat-giro{i+1}(k,:)*j2Hat)*dt;
                        else
                            alpha_giro(k)= int_aux + (giro{i}(k,:)*j1Hat-giro{i+1}(k,:)*j2Hat)*dt;
                        end
                        int_aux = alpha_giro(k);

                        %% Estimativa pelo acelerômetro
                        g1_dot = (giro{i}(k-2,:)-8*giro{i}(k-1,:)+8*giro{i}(k+1,:)-giro{i}(k+2,:))/(12*dt);
                        g2_dot = (giro{i+1}(k-2,:)-8*giro{i+1}(k-1,:)+8*giro{i+1}(k+1,:)-giro{i+1}(k+2,:))/(12*dt);
                        gamma1 = cross(giro{i}(k,:),cross(giro{i}(k,:),o1))+cross(g1_dot,o1);
                        gamma2 = cross(giro{i+1}(k,:),cross(giro{i+1}(k,:),o2))+cross(g2_dot,o2);

                        % Seel 2014 - Eq. 12
                        a1Til = acel{i}(k,:)-gamma1;
                        a2Til = acel{i+1}(k,:)-gamma2;

                        % Seel 2014 - Eq. 14
                        % u é a matriz do primeiro termo
                        % v é a matriz do segundo termo
                        u = [a1Til*x1, a1Til*y1,0]';
                        v = [a2Til*x2, a2Til*y2,0]';

                        % Cálculo do angulo articular usando o acelerômetro
                        % Acc usando arcotangente
                        alpha_acc(k) = atan2(norm(cross(u,v)),dot(u,v));

                        if k == step_batch+2
                            alpha_acc_giro(k) = alpha_acc(k);
                            alpha_acc_giro(k) = lambda*alpha_acc(k)+(1 - lambda)*(alpha_giro(k));
                        else
                            alpha_acc_giro(k) = lambda*alpha_acc(k)+(1-lambda)*(alpha_acc_giro(k-1)+alpha_giro(k)-alpha_giro(k-1));
                        end
                        
                        %%%%% Filtro de Kalman %%%%%
                        if k == step_batch+2
                            xp_joint = 0;
                            Pp_joint = 1e-0;
                            A_joint = 1;
                            Q_joint = 0; % Foi considerado o sistema livre de ruídos
                            H_joint = 1;
                            R_joint = 2;                    
                        end
                        if fk == 1
                            % Predição
                            xp_joint = A_joint*xp_joint;
                            Pp_joint = A_joint*Pp_joint*A_joint' + Q_joint;

                            % Ganho de Kalman
                            K_joint = Pp_joint*H_joint*(H_joint*Pp_joint*H_joint' + R_joint)^(-1);

                            % Medição (erro)
                            z_joint = -alpha_acc(k) + alpha_giro(k);

                            % Atualização
                            xu_joint =   xp_joint + K_joint*(z_joint - H_joint*xp_joint);
                            Pu_joint =  Pp_joint -  K_joint*H_joint*Pp_joint;

                            xp_joint =  xu_joint;
                            Pp_joint =  Pu_joint;

                            % Correção Final do joint angle
                            alpha_acc_giro(k) = alpha_giro(k) - xp_joint;
                        end
                    end
                    cont_batch = cont_batch + 1;                      
            end

            alpha_acc_v(i,:) = alpha_acc;
            alpha_giro_v(i,:) = alpha_giro;
            alpha_acc_giro_v(i,:) = alpha_acc_giro;  
        end
    
    cor = corrcoef(enc(2*step_batch+2:end-2),rad2deg(alpha_acc_giro(2*step_batch+2:end-2)));
    cor_v(b,o) = cor(1,2);
  
    end
end

% Identificando os valores para a maior correlação

[maxValue, linearIndexesOfMaxes] = max(cor_v(:));
[i j] = find(cor_v == maxValue,1,'first');

%[valor, ind] = max(cor_v(:));
%[i,j] = ind2sub(size(cor_v),ind);

disp('    _____________________________________')
disp('    Valores para a maior correlação')
disp(['    Batch:  ', int2str(step_batch_v(i))])
disp(['    Ordem:  ', int2str(ordem_v(j))])
disp(['    Correlação: ', num2str(cor_v(i,j)*100)])
disp('    _____________________________________')
%% Gráficos
if length(step_batch_v) > 1            
    figure
    surf(ordem_v,step_batch_v,cor_v)
    grid minor
    xlabel('Ordem')
    ylabel('Tamanho Batch')
    zlabel('Correlação')
    figure
    surf(ordem_v,step_batch_v,rmse_o)
    grid minor
    xlabel('Ordem')
    ylabel('Tamanho Batch')
    zlabel('Erro do offset')
    figure
    surf(ordem_v,step_batch_v,rmse_j)
    grid minor
    xlabel('Ordem')
    ylabel('Tamanho Batch')
    zlabel('Erro do versor')
else
tempo = dt*(1:length(enc(:,1))-(3 + step_batch));
for i = 1:seg-1
    figure;hold on
    plot(tempo, rad2deg(alpha_acc_v(i,step_batch+2:end-2)),'color',[0,0,0]+0.8,'LineWidth',2) % Acelerometro
    plot(tempo, rad2deg(alpha_giro_v(i,step_batch+2:end-2)),'color',[0.4940 0.1840 0.5560],'LineWidth',2) % Giroscópio
    plot(tempo, rad2deg(alpha_acc_giro_v(i,step_batch+2:end-2)),'color',[39, 143, 108]/255,'LineWidth',2) % Resultado
    plot(tempo, enc(step_batch+2:end-2,i),'k--','LineWidth',2) % Referência
    grid minor
    title('Estimativa ângulo articular')
    legend('Acel.','Giro.','Filtro')
    xlim([0,31])
    ylim([-30,60])
    ylabel(' ângulo articular do joelho [graus]')
    xlabel('tempo [s]')
    % divisão dos batches
    j = step_batch+2;

    while j < length(enc)-step_batch
        xline(dt*j,'--b','HandleVisibility','off');
        j = j + step_batch;
    end
end

end
rms_angulo = rms((rad2deg(alpha_acc_giro_v(i,2*step_batch+2:end-2))-enc(2*step_batch+2:end-2,i)'));
% figure(10)
% subplot(2,1,1)
% hold on
% xline(500,'--b')
% plot(acel{1})
% subplot(2,1,2)
% plot(giro{1})
% hold on
% xline(500,'--b')
% 
% 
%cor = corrcoef(enc(step_batch+2:end-2),rad2deg(alpha_acc_giro(step_batch+2:end-2))) % Cálculo da correlação

% Resultados para as estimativas
% disp('     Phi 1     Phi 2    Theta 1   Theta 2')
% disp('    _____________________________________')
% disp([x_j(1) x_j(3) x_j(2) x_j(4)])
% disp('    _____________________________________')
% fprintf('    Correlação: %.2f \n', cor(1,2)*100)