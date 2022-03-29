function [j1Hat,j2Hat,x_j,rmse] = estimative_j_kalman(x_j, g1, g2, step_batch, ordem)

    F = eye(4);
    Q = eye(4).*0.001;
    R = eye(step_batch).*0.001;
    P = eye(4);
    
    persistent firstRun_j

    H = zeros(step_batch,4);
    r = zeros(step_batch,1);
    
    j1Hat_ant = [cos(x_j(1))*cos(x_j(2)),cos(x_j(1))*sin(x_j(2)),sin(x_j(1))]';
    j2Hat_ant = [cos(x_j(3))*cos(x_j(4)),cos(x_j(3))*sin(x_j(4)),sin(x_j(3))]';
    
    for l = 1:ordem        
        % Predição
        x_p = F*x_j;
        Pp = F*P*F'+ Q;     
        
        % J1 e J2 - Seel 2014 - Eq. 4 e 5
        j1Hat = [cos(x_p(1))*cos(x_p(2)),cos(x_p(1))*sin(x_p(2)),sin(x_p(1))]';
        j2Hat = [cos(x_p(3))*cos(x_p(4)),cos(x_p(3))*sin(x_p(4)),sin(x_p(3))]';   
        
        for j = 1:step_batch
            % Auxiliar para calcular de/dj - Seel 2012 - Eq. 2
            normJ1 = (norm( cross(g1(j,:),j1Hat)));
            normJ2 = (norm(cross(g2(j,:),j2Hat)));
            
            if normJ1 == 0
                normJ1 = 1e-5;
            end
            
            if normJ2 == 0
                normJ2 = 1e-5;
            end
            
            % Derivadas do Erro em relação J1 e J2 - Seel 2012 - Eq. 2
            dj1 = (cross((cross(g1(j,:),j1Hat)),g1(j,:)))/normJ1;
            dj2 = -(cross((cross(g2(j,:),j2Hat)),g2(j,:)))/normJ2;

            % Derivadas de J1 e J2 em relação a X
            % X é o vator com os valores de phi1,phi2,theta1 e theta2
            dj1dx = [-sin(x_p(1))*cos(x_p(2)) -cos(x_p(1))*sin(x_p(2)) 0 0
                     -sin(x_p(1))*sin(x_p(2))  cos(x_p(1))*cos(x_p(2)) 0 0
                      cos(x_p(1))              0                       0 0];
            dj2dx = [0 0 -sin(x_p(3))*cos(x_p(4)) -cos(x_p(3))*sin(x_p(4))
                     0 0 -sin(x_p(3))*sin(x_p(4))  cos(x_p(3))*cos(x_p(4))
                     0 0  cos(x_p(3))              0];

           % Regra da cadeia derivada do Erro em relação a X
           H(j,:) = (dj1*dj1dx + dj2*dj2dx);   
           % Calculo do erro

           r(j) = -(norm(cross(g1(j,:),j1Hat))-norm(cross(g2(j,:),j2Hat)));
        end 
        % Correção
        K = Pp*H'*pinv(H*Pp*H'+R);
        x_j = x_p + K*r;
        P = (eye(4) - K*H)*Pp;
    end
    
    j1Hat = [cos(x_j(1))*cos(x_j(2)),cos(x_j(1))*sin(x_j(2)),sin(x_j(1))]';
    j2Hat = [cos(x_j(3))*cos(x_j(4)),cos(x_j(3))*sin(x_j(4)),sin(x_j(3))]';


    
    if isempty(firstRun_j)
        firstRun_j = 1;
        rmse = sqrt(sum(r.^2)/length(r));
    else
        for j = 1:step_batch
            r_ant(j) = (norm(cross(g1(j,:),j1Hat_ant))-norm(cross(g2(j,:),j2Hat_ant)));
        end
        rmse = sqrt(sum(r.^2)/length(r));
        rmse_ant = sqrt(sum(r_ant.^2)/length(r_ant));
        if rmse_ant <= rmse
            j1Hat = j1Hat_ant;
            j2Hat = j2Hat_ant;
            rmse = rmse_ant;
        end
        
    end
    x_j = inclinacao_orientacao(x_j);
    j1Hat = [cos(x_j(1))*cos(x_j(2)),cos(x_j(1))*sin(x_j(2)),sin(x_j(1))]';
    j2Hat = [cos(x_j(3))*cos(x_j(4)),cos(x_j(3))*sin(x_j(4)),sin(x_j(3))]';
end