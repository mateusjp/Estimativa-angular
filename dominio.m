function Ksi = dominio(Ksi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% "Garantindo" o domínio da função
while Ksi(1) > pi/2 %%%%% p1 %%%%
    Ksi(1) = Ksi(1) - pi;
    %cont1 = cont1 + 1
end

while Ksi(1) < -pi/2 %%%% p1 %%%%
    Ksi(1) = Ksi(1) + pi;
    %cont2 = cont2 + 1
end

while Ksi(3) > pi/2 %%%%% p2 %%%%
    Ksi(3) = Ksi(3) - pi;
    %cont3 = cont3 + 1
end

while Ksi(3) < -pi/2 %%%% p2 %%%%
    Ksi(3) = Ksi(3) + pi;
    %cont4 = cont4 + 1
end

while Ksi(2) > 2*pi %%%%% t1 %%%%%
    Ksi(2) = Ksi(2) - 2*pi;
    %cont5 = cont5 + 1
end

while Ksi(2) < 0 %%%%%% t1 %%%%%%%
    Ksi(2) = Ksi(2) + 2*pi;
    %cont6 = cont6 + 1
end

while Ksi(4) > 2*pi %%%%%% t2 %%%%%%%
    Ksi(4) = Ksi(4) - 2*pi;
    %cont7 = cont7 + 1
end

while Ksi(4) < 0 %%%%%% t2  %%%%%%%
    Ksi(4) = Ksi(4) + 2*pi;
    %cont8 = cont8 + 1
end
end

