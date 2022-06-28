% monte carlo method for expected F_pull
% we assume to have 2 fingers

clc, clear, close all

N = 10;         % spines on each side
alpha = pi/6;   % [rad] angle of the tangent to the surface
mu = 0.39;      % friction coefficient
Fint = 10;       % [N] gripper preload
Fmax = 10;      % [N] maximum transmittable force by a single spine (cap due to strength)
kn = 15;        % [N/m] normal stiffness (spine compliance)
kt = 30;        % [N/m] tangential stiffness (spine compliance)
n_mc = 400;     % amount of monte carlo iterations
n_beta = 90/5 + 1;

% function handles for linear equations from statics

A = @(alpha,kn,kt) [sin(alpha), cos(alpha), -sin(alpha), -cos(alpha);
-cos(alpha), sin(alpha), -cos(alpha), sin(alpha);
cos(alpha), -sin(alpha),  -cos(alpha), sin(alpha);
-sin(alpha)/kn, -cos(alpha)/kt, -sin(alpha)/kn, -cos(alpha)/kt];

b = @(Fboom,beta) [-Fboom*sin(beta); Fboom*cos(beta); Fboom*tan(alpha)*sin(beta); 0];

Forces = @(alpha,kn,kt,Fboom,beta) A(alpha,kn,kt)\b(Fboom,beta);    % solves contact forces at the finger tips
                                                                    % in the order
                                                                    % [Fln; Flt; Frn; Frt]
% compute_mu calculates the friction coefficient that would be needed for the spine not to slip (Eq. 4.5)                                                  
compute_mu = @(Fn,Ft,Fp,alpha,psi) (Ft*cos(psi) - Fn*sin(psi) + Fp*cos(alpha+psi))/(Ft*sin(psi) + Fn*cos(psi) + Fp*sin(alpha+psi));
% detachment check whole finger detachment 
detachment = @(Fn,Ft,Fp,alpha) Fn*sin(alpha) + Ft*cos(alpha) + Fp < 0;

F = zeros(2,n_beta);
beta = linspace(0,pi/2,n_beta);

parfor j = 1:n_beta
    Fpull = zeros(1,n_mc); % at the end of each beta iteration will store mean and std dev of this vector in vector F
    for i = 1:n_mc

        psi = (pi/2 - (atan(1/mu) - alpha))*rand(N,2) + (atan(1/mu) - alpha);     % generates the local asperity angle for each spine from a uniform distribution
                                                                                  % the range is set in order for the spines to not slip when only the preload is applied
    %     coeff = 2*sin(alpha)./(cos(psi) - mu*sin(psi)).*(sin(alpha)*(sin(psi) + mu*cos(psi)) + cos(alpha)*(mu*sin(psi) - cos(psi)));

        Fp = Fint/N*ones(1,2);            % [N] initial preload on single spines
        ok = 1;     % variable used in the while loop below ok=1 when grasp is still active, ok=0 when failure happens
        engaged = ones(N,2);    % we assume all the spines to be engaged at the beginning
    %     engaged(coeff<0) = 0;

        while ok

            n = sum(engaged,1);     % computes how many spines are currently engaged on each side
            Ff = Forces(alpha,kn,kt,Fpull(i),beta(j));       % computes the average load on spines (different on each side)
            Ffs = [sqrt(Ff(1)^2+Ff(2)^2), sqrt(Ff(3)^2+Ff(4)^2)];   % vector sum of normal and tangent contact force on each side
            Fe = Ffs./n;    % force applied on single spines
            Fp = Fint./n;   % preload applied on single spines (internal force divided by number of active spines)

            % computes the required friction coefficient there should be on
            % each spine contact to be able to transmit the contact force
            % resulting from equilibrium without slipping
            left_mus = arrayfun(compute_mu,ones(N,1)*Ff(1),ones(N,1)*Ff(2),ones(N,1)*Fint,ones(N,1)*alpha,psi(:,1));
            right_mus = arrayfun(compute_mu,ones(N,1)*Ff(3),ones(N,1)*Ff(4),ones(N,1)*Fint,ones(N,1)*alpha,psi(:,2));
            mus = [left_mus, right_mus].*engaged;
            mus(mus == 0) = NaN; % weird thing to get around the spines that already failed
            [max_mus, index] = max(mus); % computes the maximum (across spines) required mu on both sides
            disengagement = max_mus > mu;

            left_detach = detachment(Ff(1),Ff(2),Fint,alpha);
            right_detach = detachment(Ff(3),Ff(4),Fint,alpha);

            if (left_detach + right_detach) > 0
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

%%
close all
% beta = linspace(0,pi/2,n_beta);

polarplot(pi/2-beta,F(1,:),'-*b')
hold on
polarplot(pi/2-beta,F(1,:)+F(2,:),'-*r')
polarplot(pi/2-beta,F(1,:)-F(2,:),'-*r')
polarplot(pi/2+beta,F(1,:),'-*b')
polarplot(pi/2+beta,F(1,:)+F(2,:),'-*r')
polarplot(pi/2+beta,F(1,:)-F(2,:),'-*r')
% xlabel('Fx [N]')
% ylabel('Fy [N]')
title('Fint=15')