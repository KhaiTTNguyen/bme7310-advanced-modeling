format long
center = pi/6;
center_val = sin(pi/6);
analytical_deriv = 0.866025403784;
h1 = pi/20;
h2 = pi/40;
h3 = pi/80;

center_err = zeros(1,3);
backward_err = zeros(1,3);
forward_err = zeros(1,3);
[center_err(1), backward_err(1), forward_err(1)] = get_u_from_h(center, h1, analytical_deriv);
[center_err(2), backward_err(2), forward_err(2)] = get_u_from_h(center, h2, analytical_deriv);
[center_err(3), backward_err(3), forward_err(3)] = get_u_from_h(center, h3, analytical_deriv);

h = 1./[h1,h2,h3];
figure(1)
loglog(h, center_err,...
       h, backward_err,...
       h, forward_err)
legend('Cntr', 'Back', 'Fwd')
grid on
ylabel('Error')
xlabel('1/h')
title('Loglog plot of Error w.r.t to reciprocal of respective step size.')

function [center_error, backward_error, forward_error] =  get_u_from_h(center, h, analytical_deriv)

u = sin(center-2*h : h : center+2*h);

center_diff = (u(1) - 8*u(2) + 8*u(4) - u(5)) / (12*h);
backward_diff = (u(1) - 4*u(2) + 3*u(3)) / (2*h);
forward_diff = (u(4) - u(3))/h;

center_error = abs(center_diff - analytical_deriv);
backward_error = abs(backward_diff - analytical_deriv);
forward_error = abs(forward_diff - analytical_deriv);
fprintf('---------Values for derivatives at %.12f with %.12f--------\n', sin(center), h);
fprintf('Center diff   : %.12f, Error %.12f \n', center_diff, center_error);
fprintf('Backward diff : %.12f, Error %.12f \n', backward_diff, backward_error);
fprintf('Forward diff  : %.12f, Error %.12f \n', forward_diff, forward_error);


end




