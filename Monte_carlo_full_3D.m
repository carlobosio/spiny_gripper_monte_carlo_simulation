% monte carlo method for expected F_pull
% we assume to have 2 fingers

clc, clear, close all

N = 12;         % spines on each side
alpha = pi/6;   % [rad] angle of the tangent to the surface
R = 20;         % [cm] rock radius of curvature
x = 2;          % [cm] distance between wrist and rock
mu = 0.39;      % friction coefficient
Fint = 20;      % [N] gripper preload
Fmax = 12;      % [N] maximum transmittable force by a single spine (cap due to strength)
kn = 15;        % [N/m] normal stiffness (spine compliance)
kt = 30;        % [N/m] tangential stiffness (spine compliance)
kc = 1;         % [N/m] tentative compliance
n_mc = 200;     % amount of monte carlo iterations
n_samples = 10000;

% function handles for linear equations from statics

A = @(alpha,kn,kt,kc) [sin(alpha), cos(alpha), 0, -0.5*sin(alpha), -0.5*cos(alpha), -sqrt(3)/2, -0.5*sin(alpha), -0.5*cos(alpha), sqrt(3)/2;
                       0, 0, 1, sqrt(3)/2*sin(alpha), sqrt(3)/2*cos(alpha), -0.5, -sqrt(3)/2*sin(alpha), -sqrt(3)/2*cos(alpha), -0.5;
                       cos(alpha), -sin(alpha), 0, cos(alpha), -sin(alpha), 0, cos(alpha), -sin(alpha), 0;
                       0, 0, 0, sin(alpha)*cos(alpha)*sqrt(3)/2, -(sin(alpha))^2*sqrt(3)/2, 0, -sin(alpha)*cos(alpha)*sqrt(3)/2, (sin(alpha))^2*sqrt(3)/2, 0;
                       -sin(alpha)*cos(alpha), (sin(alpha))^2, 0, 0.5*sin(alpha)*cos(alpha), -0.5*(sin(alpha))^2, 0, 0.5*sin(alpha)*cos(alpha), -0.5*(sin(alpha))^2, 0;
                       0, 0, 1, 0, 0, 1, 0, 0, 1;
                       sin(alpha)/kn, cos(alpha)/kt, 0, sin(alpha)/kn, cos(alpha)/kt, 0, sin(alpha)/kn, cos(alpha)/kt, 0;
                       sin(alpha)/kn, cos(alpha)/kt, 0, 0, 0, 1/(kc*sqrt(3)), 0, 0, -1/(kc*sqrt(3));
                       0.5*sin(alpha)/kn, 0.5*cos(alpha)/kt, 0.5*sqrt(3)/kc, 0, 0, 0, sin(alpha)/kn, cos(alpha)/kt, 0];

b = @(Fboom,beta,phi,alpha,R,x) [-Fboom*sin(beta)*cos(phi); 
                                 -Fboom*sin(beta)*sin(phi);
                                 -Fboom*cos(beta); 
                                 Fboom*sin(beta)*sin(phi)*(1-cos(alpha) + x/R); 
                                 -Fboom*sin(beta)*cos(phi)*(1-cos(alpha) + x/R); 
                                 0; 
                                 0; 
                                 0; 
                                 0];

Forces = @(alpha,kn,kt,kc,Fboom,beta,phi,R,x) A(alpha,kn,kt,kc)\b(Fboom,beta,phi,alpha,R,x);    % solves contact forces at the finger tips
                                                                    % in the order
                                                                    % [Fln; Flt; Frn; Frt]
compute_mu = @(Fn,Ft,Fp,alpha,psi) (Ft*cos(psi) - Fn*sin(psi) + Fp*cos(alpha+psi))/(Ft*sin(psi) + Fn*cos(psi) + Fp*sin(alpha+psi));
detachment = @(Fn,Ft,Fp,alpha) Fn*sin(alpha) + Ft*cos(alpha) + Fp < 0;

%%

F = zeros(2,n_samples);
beta = pi/2*rand(1,n_samples);
phi = 2*pi*rand(1,n_samples);

parfor j = 1:n_samples
    Fpull = zeros(1,n_mc);
    for i = 1:n_mc

        psi = (pi/2 - (atan(1/mu) - alpha))*rand(N,3) + (atan(1/mu) - alpha)*ones(N,3);     % generate the local asperity angle for each spine from a uniform distribution
    %     coeff = 2*sin(alpha)./(cos(psi) - mu*sin(psi)).*(sin(alpha)*(sin(psi) + mu*cos(psi)) + cos(alpha)*(mu*sin(psi) - cos(psi)));

        Fp = Fint/N*ones(1,3);            % [N] initial preload on single spines
        ok = 1;
        engaged = ones(N,3);    % we assume almost all the spines to be engaged at the beginning
    %     engaged(coeff<0) = 0;

        while ok

            n = sum(engaged,1);     % computes how many spines are engaged on each side
            Ff = Forces(alpha,kn,kt,kc,Fpull(i),beta(j),phi(j),R,x);       % computes the average load on spines (different on each side)
            Ffs = [sqrt(Ff(1)^2 + Ff(2)^2 + Ff(3)^2), sqrt(Ff(4)^2 + Ff(5)^2 + Ff(6)^2),sqrt(Ff(7)^2 + Ff(8)^2 + Ff(9)^2)];
            Fe = Ffs./n;
            Fp = Fint./n;

            % computes the required friction coefficient there should be on
            % each spine contact to be able to transmit the contact force
            % resulting from equilibrium
            mus_1 = arrayfun(compute_mu,ones(N,1)*Ff(1),ones(N,1)*Ff(2),ones(N,1)*Fint,ones(N,1)*alpha,psi(:,1));
            mus_2 = arrayfun(compute_mu,ones(N,1)*Ff(4),ones(N,1)*Ff(5),ones(N,1)*Fint,ones(N,1)*alpha,psi(:,2));
            mus_3 = arrayfun(compute_mu,ones(N,1)*Ff(7),ones(N,1)*Ff(8),ones(N,1)*Fint,ones(N,1)*alpha,psi(:,3));
            mus = [mus_1, mus_2, mus_3].*engaged;
            mus(mus == 0) = NaN;
            [max_mus, index] = max(mus);
            disengagement = max_mus > mu;

            detach_1 = detachment(Ff(1),Ff(2),Fint,alpha);
            detach_2 = detachment(Ff(4),Ff(5),Fint,alpha);
            detach_3 = detachment(Ff(7),Ff(8),Fint,alpha);

            if detach_1 || detach_2 || detach_3 > 0
                ok = 0;
            end

    %         ref_coeff = Fe./Fp; 
    %         actual_coeff = coeff.*engaged;
    %         actual_coeff(actual_coeff == 0) = NaN;
    %         [risk_coeff,index] = min(actual_coeff);  % takes the minimum nonzero ratio between Fe and Fp
    %         disengagement = risk_coeff < ref_coeff;

            if sum(disengagement > 0)
                k = find(disengagement);
                engaged(index(k),k) = 0;
            end

            if min(sum(engaged)) == 0 | sum(Fe > Fmax) > 0
                ok = 0;
            end

            Fpull(i) = Fpull(i) + 0.1;
        end
        Fpull(i) = Fpull(i) - 0.1;
    end
    F(:,j) = [mean(Fpull), std(Fpull)];
end

%save('output',"F",'alpha','beta','phi','Fint','Fmax','kn','kt','kc','N','R','x')