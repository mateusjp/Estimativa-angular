function Ksi = sentidoj1j2(Ksi)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% "Garante" os versores no mesmo sentido

if Ksi(1) < 0
    Ksi(2) = Ksi(2) + pi;
    Ksi(1) = - Ksi(1);
end
if Ksi(3) < 0
    Ksi(3) = - Ksi(3);
    Ksi(4) = Ksi(4) + pi;
end

end


% if Ksi(1)/abs(Ksi(1)) ~= Ksi(3)/abs(Ksi(3))
%     sel = randi(2); % roleta 50/50
%     if sel == 1
%         Ksi(1) = - Ksi(1);
%         if Ksi(1) < 0
%             Ksi(2) = Ksi(2) + pi;
%         else
%             Ksi(2) = Ksi(2) - pi;
%         end
%     else
%         Ksi(3) = - Ksi(3);
%         if Ksi(3) < 0
%             Ksi(4) = Ksi(4) + pi;
%         else
%             Ksi(4) = Ksi(4) - pi;
%         end
%     end
% end
% end

