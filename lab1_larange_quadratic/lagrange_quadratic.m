function [xr,L20,L21,L22]=lagrange_quadratic(x,yvals,xi)
% Function assumes x0 < x1 < x2
% Inputs:
%         x      :    the three x values [x0 x1 x2]
%         yvals  :    the specific values of your function at each x0,x1,x2
%         xi     :    the specific value for x that you want to interpolate
%                     your function to
x0=x(1);
x1=x(2);
x2=x(3);

dx=(x2-x0)/100;
xr=[x0:dx:x2];
h01=x1-x0;
h02=x2-x0;
h12=x2-x1;

nn=length(xr);

% In this loop, I just want you to evaluate the 3 Lagrange quadratic
% polynomials, L20, L21, L22 for the range of xvalues between x0 and and
% x2.  The values have been stored in 'xr' above.
for i=1:nn
    L20(i)=(xr(i)-x1)*(xr(i)-x2)/((-h01)*(-h02));
    L21(i)=(xr(i)-x0)*(xr(i)-x2)/((h01)*(-h12));
    L22(i)=(xr(i)-x0)*(xr(i)-x1)/((h02)*(h12));
end
% ---------------------- Part a ------------------------
figure(1)
plot(xr, L20)
hold on
plot(xr,L21)
hold on
plot(xr,L22)
hold off
legend('L20','L21','L22')

% ---------------------- Part b ------------------------
% Now let's evaluate a sample interpolated value using the quadratic
% Lagrange basis.  Now you will evaluate the 3 basis functions at
% a specific value x, your 'xi' input above.  You will calculate the 
% interpolated value of your function at that specific value.
L20v=(xi-x1)*(xi-x2)/((-h01)*(-h02));
L21v=(xi-x0)*(xi-x2)/((h01)*(-h12));
L22v=(xi-x0)*(xi-x1)/((h02)*(h12));

% Here is where you will calculate the interpolated value of your function
% at the specific value of x, 'xi'
yv = yvals(1)*L20v + yvals(2)*L21v + yvals(3)*L22v;

fprintf('Your Values for the Basis functions at %f\n', xi);
fprintf('L20, L21, L22 = %f, %f, %f, repectively\n',L20v, L21v, L22v);
fprintf('f(Xi) %f\n',yv);

% ---------------------- Part c ------------------------

for i=1:nn
    y_interpolated(i) = yvals(1)*L20(i) + yvals(2)*L21(i) + yvals(3)*L22(i);
end
figure(2)
plot(xr, y_interpolated)
legend('f(x)')

figure(3)
plot(xr, y_interpolated)
hold on
plot(xr,L20)
hold on
plot(xr,L21)
hold on
plot(xr,L22)
hold off
legend('f(x)','L20','L21','L22')
end