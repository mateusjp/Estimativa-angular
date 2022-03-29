function [o1,o2,x_o,rmse] = estimative_o_kalman(x_o, g1, g2, g1Dot, g2Dot, j1Hat, j2Hat, a1, a2, step_batch, ordem)

    F = eye(6);
    Q = eye(6).*0.001;
    R = eye(step_batch).*0.001;
    P = eye(6);
    
    persistent firstRun_o
    
    H = zeros(step_batch,6);
    r = zeros(step_batch,1);
   
    gamma1 = zeros(step_batch,3);
    gamma2 = zeros(step_batch,3);
    
    o1_ant = x_o(1:3,1);
    o2_ant = x_o(4:6,1);
    
    for l = 1:ordem

        % Predição
        x_p = F*x_o;
        Pp = F*P*F'+ Q; 
        
        o1 = x_p(1:3,1);
        o2 = x_p(4:6,1);
        
                
        for j = 1:step_batch
            
            gamma1(j,:) = cross(g1(j,:),cross(g1(j,:),o1))+cross(g1Dot(j,:),o1);
            gamma2(j,:) = cross(g2(j,:),cross(g2(j,:),o2))+cross(g2Dot(j,:),o2);

            a1Til = a1(j,:)-gamma1(j,:);
            a2Til = a2(j,:)-gamma2(j,:);

            norm_o1 = norm(a1(j,:)-gamma1(j,:));
            norm_o2 = norm(a2(j,:)-gamma2(j,:));

            Mg1 = [ 0       -g1(j,3)  g1(j,2)
                    g1(j,3)  0       -g1(j,1)
                   -g1(j,2)  g1(j,1)  0      ];
            Mg2 = [ 0       -g2(j,3)  g2(j,2)
                    g2(j,3)  0       -g2(j,1)
                   -g2(j,2)  g2(j,1)  0      ];
            Mg1Dot = [ 0          -g1Dot(j,3)  g1Dot(j,2)
                       g1Dot(j,3)  0          -g1Dot(j,1)
                      -g1Dot(j,2)  g1Dot(j,1)  0      ];
            Mg2Dot = [ 0          -g2Dot(j,3)  g2Dot(j,2)
                       g2Dot(j,3)  0          -g2Dot(j,1)
                      -g2Dot(j,2)  g2Dot(j,1)  0      ];

            % Calculo derivada parcial em relacao a o1

            dedo1 = ((Mg1^2-Mg1Dot)*a1Til'/norm_o1)';
            dedo2 = ((Mg2^2-Mg2Dot)*a2Til'/norm_o2)';

           % Regra da cadeia derivada do Erro em relação a X
           H(j,:) = [-dedo1 dedo2];   

           % Calculo do erro
           r(j) = -(norm(a1(j,:)-gamma1(j,:))-norm(a2(j,:)-gamma2(j,:)));
        end
        % Correção
        K = Pp*H'*pinv(H*Pp*H'+R);
        x_o = x_p + K*r;
        P = (eye(6) - K*H)*Pp;
    end
    
    o1 = x_o(1:3,1);
    o2 = x_o(4:6,1);
        
    o1 = o1 - j1Hat*(dot(o1,j1Hat)+dot(o2,j2Hat))/2;
    o2 = o2 - j2Hat*(dot(o1,j1Hat)+dot(o2,j2Hat))/2;   
    
    if isempty(firstRun_o)
        firstRun_o = 1;
        rmse = sqrt(sum(r.^2)/length(r));
    else
        for j = 1:step_batch
            
            gamma1(j,:) = cross(g1(j,:),cross(g1(j,:),o1_ant))+cross(g1Dot(j,:),o1_ant);
            gamma2(j,:) = cross(g2(j,:),cross(g2(j,:),o2_ant))+cross(g2Dot(j,:),o2_ant);
            
            r_ant(j) = (norm(a1(j,:)-gamma1(j,:))-norm(a2(j,:)-gamma2(j,:)));
        end
        rmse = sqrt(sum(r.^2)/length(r));
        rmse_ant = sqrt(sum(r_ant.^2)/length(r_ant));
        if rmse_ant <= rmse
            o1 = o1_ant;
            o2 = o2_ant;
            rmse = rmse_ant;
        end
        
    end

end