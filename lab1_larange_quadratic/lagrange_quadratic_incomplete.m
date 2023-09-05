function [xr,L20,L21,L22]=lagrange_quadratic_incomplete(x,yvals,xi)
% Function assumes x0 < x1 < x2
% Inputs:
%         x      :    the three x values [x0 x1 x2]
%         yvals  :    the specific values of your function at each x0,x1,x2
%         xi     :    the specific value for x that you want to interpolate
%                     your function too
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
    L21(i)=???;
    L22(i)=???;
end

% Now let's evaluate a sample interpolated value using the quadratic
% Lagrange basis.  Now you will not only evaluate the 3 basis functions at
% a specific value x, your 'xi' input above.  You will calculate the 
% interpolated value of your function at that specific value.
L20v=???
L21v=???
L22v=???

% Here is where you will calculate the interpolated value of your function
% at the specific value of x, 'xi'
yv=???

fprintf('Your Values for the Basis functions at the Xi\n');
fprintf('L20, L21, L22 = %f, %f, %f, repectively\n',L20v, L21v, L22v);
fprintf('f(Xi) %f\n',yv);

