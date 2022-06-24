% Discrete-Time Fast Terminal Sliding Mode Control for Permanent Magnet
% Linear Motor without Disturbance compensation

clc; clear all; close all;
m = 5.4; R = 16.8; kf = 130; ke = 123; h = 0.005;
a = kf*ke/(R*m); b = kf/(R*m); c1 = 1.5; c2 = 1.5; alpha = 1/2;
u = zeros(2001,1); e1 = u; e2 = u; xr_dot = u; xr_ddot = u; t = u; F = u;
xr = u; x1 = u; x2 = u; Ff = u; Fr = u;
fc = 10; fs = 20; fv = 10; xsd = 0.1; A1 = 8.5; A2 = 4.25; A3 = 2; w = 3.14;

for k = 1:length(u)-1
    % current time
    t(k+1) = t(k) + h;

    % desired value and its derivatives
    xr(k+1) = 5*sin(t(k+1));
    xr_dot(k+1) = 5*cos(t(k+1));
    xr_ddot(k+1) = -5*sin(t(k+1));

    % control signal
    u(k) = 1/(h*b)*((1+c1*h-h*a)*e2(k)+c1*e1(k)+h*(a*xr_dot(k)+xr_ddot(k))+...
        c2*sig(e1(k)+h*e2(k),alpha));

    % state error caclulation
    e1(k+1) = e1(k) + h*e2(k);
    e2(k+1) = e2(k) - h*b*u(k) - h*a*e2(k) + h*(a*xr_dot(k)+xr_ddot(k)) + h*F(k);

    % states update
    x1(k+1) = xr(k+1) - e1(k+1);
    x2(k+1) = xr_dot(k+1) - e2(k+1);

    % disturbance
    Ff(k+1) = (fc + (fs - fc)*exp(-(x2(k+1)/xsd)^2) + fv*x2(k+1)*sign(x2(k+1)));
    Fr(k+1) = A1*sin(w*x1(k+1)) + A2*sin(3*w*x1(k+1)) + A3*sin(5*w*x1(k+1));
    F(k+1) = (Ff(k+1) + Fr(k+1))/m;

end

% plot results
figure(1);
plot(t,x1,t,xr,':r');
grid on;
axis([0 10 -5 5]);
figure(2);
plot(t,e1);
grid on;
axis([min(t) max(t) min(e1) max(e1)]);

function y=sig(x,a)
    y = sign(x)*abs(x)^a;
end