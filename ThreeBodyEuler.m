function ThreeBodyEuler(m1,S1,m2,S2,P,tspan)
% 
% Author: Erik Lamoureux
% Date: Novemeber 17, 2015
% 
% Purpose: Approximates a solution to Euler's Three Body Problem
% Inputs:
%   m1 is the mass of body 1 (in solar mass)
%   S1 is a vector of length 2:
%       S1(1) is the x-position of Star 1
%       S1(2) is the y-position of Star 1
%   m2 is the mass of body 2 (in solar mass)
%   S2 is a vector of length 2:
%       S2(1) is the x-position of Star 2
%       S2(2) is the y-position of Star 2
%   P is a vector of length 4:
%       P(1) is the initial x-position of the planet at time tspan(1)
%       P(2) is the initial x-velocity of the planet at time tspan(1)
%       P(3) is the initial y-position of the planet at time tspan(1)
%       P(4) is the initial y-velocity of the planet at time tspan(1)
%   tspan is the interval of integration [t0,tf]
% 

% Define the function defining the system of equations
    function dudt = odefun(t,u)
        dudt = zeros(4,1);
        % Calculate square root of distance from S1
        d1 = sqrt( (S1(1) - u(1))^2 + (S1(2) - u(3))^2 );
        % calculate square root of distance from S2
        d2 = sqrt( (S2(1) - u(1))^2 + (S2(2) - u(3))^2 );
        
        % Calculate gravitational constant
        G = 4*pi^2;
        
        % Calculate system of equations
        dudt(1) = u(2);
        dudt(2) = ( G * m1 * (S1(1) - u(1)) ) / (d1^3) + ( G * m2 * (S2(1) - u(1)) ) / (d2^3);
        dudt(3) = u(4);
        dudt(4) = ( G * m1 * (S1(2) - u(3)) ) / (d1^3) + ( G * m2 * (S2(2) - u(3)) ) / (d2^3);
        
    end

% Define an array of evenly spaced time values
t = linspace( tspan(1), tspan(2), (tspan(2)-tspan(1))*100);

% Solve using ode23
% Output T is the same as input t
% Output U is an aray of values at the time values in T
[T,U] = ode23(@odefun,t,P);

dmax = max(max([abs(U(:,1)),abs(U(:,3))]));

for i = 1:length(T)
    % Plot planet and its trajectory
    plot(U(:,1),U(:,3),'b','LineWidth',2);
    hold on;
    plot(U(i,1),U(i,3),'ro','MarkerSize',7,'MarkerEdgeColor','r','MarkerFaceColor',[1,0,0]);
    
    % Plot Star 1 and 2
    plot(S1(1),S1(2),'go','MarkerSize',m1*10,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,0]);
    plot(S2(1),S2(2),'go','MarkerSize',m2*10,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,0]);
    hold off
    axis equal;
    grid on;
    xlim(1.2*[-dmax,dmax]);
    ylim(1.2*[-dmax,dmax]);
    % Animate
    drawnow;
end

end