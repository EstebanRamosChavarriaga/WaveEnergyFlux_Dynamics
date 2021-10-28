function [k,iter,e_a,omega] = jfpa_dispersionNewtonRaphson(g,h,T)
% % % T=8;
% % % g=9.81;
% % % h=16;

omega = 2*pi/T;

mu = (omega^2)/g; % Eq. 11

% Newton-Raphson to quantify k_0 and k_n
% [Chapra & Canale 2010 Section 6.2]
k_01 = 1;
iter = 0;
e_a = 100;

while abs(e_a) > 0.1
        f_k = k_01*tanh(k_01*h) - mu;
        fp_k = k_01*h*(1-(tanh(k_01*h))^2) + tanh(k_01*h);
        k_02 = k_01 - f_k/fp_k; % Eq. 6.6
        e_a = (k_02 - k_01)/k_02*100; % Percent relative error, Eq. 3.5
        
        iter = iter+1;
        k_01 = k_02;
end

k = k_02;

% Check...
% % % [k,~,~] = jfpa_dispersionGuo(g,h,T)

end