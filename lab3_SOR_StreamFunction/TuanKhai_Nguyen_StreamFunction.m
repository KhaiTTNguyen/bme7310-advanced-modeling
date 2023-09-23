figure(1); subplot(1,2,1);subplot(1,2,2);
figure(2); 

% Gauss-Seidel
clear;
h=0.05;
x=[0:h:2];
y=[0:h:1];
A=zeros(length(y), length(x));
%% ---- BC -------
A(2:end-1,1) = (1/4) * ( A(1:end-2, 1) + A(3:end, 1) + 2*A(2:end-1, 2));
A(2:end-1,end) = 1/(4+2*h) * ( A(1:end-2, end) + A(3:end, end) + 2*A(2:end-1, end-1) - 2*h); % type3;
A(1,:)=1;
A(end,:)=3;
% ---- obstruction -----
A(6:14,16:24) = 0; 
       
error=1;
itr=0;
pitr=0;
while (error > 1e-5 && itr < 10000)
    itr=itr+1;
    Aold=A;
    for i=2:(length(y)-1)
        for j=1:length(x)
            if j==1
                A(i,j) = 1/4 * (A(i-1, j) + Aold(i+1, j) + 2*Aold(i, j+1));
            elseif j==length(x)
                A(i,j) = 1/(4+2*h) * ( A(i-1, j) + Aold(i+1, j) + 2*A(i, j-1) - 2*h);
            elseif 6<=i && i<=14 && 16<=j && j<=24
                A(i,j) = 0;
            elseif j~=length(x) && j~=1   
                A(i,j)=1/4*(A(i-1,j)+Aold(i+1,j)+A(i,j-1)+Aold(i,j+1));
            end
        end
    end
    errorold=error;
    error=max(max(abs(A-Aold)));
    errornew=error;
    spectral(itr)=errornew/errorold;
    pitr=pitr+1;
    if pitr==5
       figure(1);subplot(1,2,1);
       contour(A);
       pitr=0;
    end   
end
 
fprintf('Gauss Seidel Iterations %d\n',itr);
figure(1); subplot(1,2,1);
contour(A);
xlabel('x')
ylabel('y')
title(['GaussSeidel: # of Iterations ' num2str(itr)]);
figure(1); subplot(1,2,2);
[X,Y]=meshgrid(x, y);
surf(X,Y,A)
title(['GaussSeidel: # of Iterations ' num2str(itr)]);

figure(2), clf
plot(spectral);
gauss_spectral = mean(spectral(itr-5:itr));
fprintf('Spectral Radius %f\n',gauss_spectral);

fprintf('Value at point  x=0.5, y=0.5 %.12f \n', A(10,10))
%%
Vx=zeros(length(y), length(x));
Vy=zeros(length(y), length(x));
 
Vx(2:end-1, 2:end-1) = ( A(3:end, 2:end-1) - A(1:end-2, 2:end-1) ) ./ (2*h);
Vy(:, 2:end-1) = - ( A(:,3:end) - A(:,1:end-2) ) ./ (2*h);
figure(3), clf
quiver(x, y, Vx, Vy)
title('Velocity vector')

% % SOR
% clearvars -except gauss_spectral 
% 
% h=0.05;
% y=[0:.05:1]';
% A=zeros(21,21);
% for i=1:21
%   A(i,21)=1;
%   A(1,i)=1;
%   A(21,i)=1;
%   A(i,1)=cos(2*pi*(i*h-0.05));
% end
%  
% omega = 2 / (1 + sqrt(1 - gauss_spectral));
% error=1;
% itr=0;
% pitr=0;
% while (error > 1e-5 & itr < 10000)
%    itr=itr+1;
%    Aold=A;
%    for i=2:20
%        for j=2:20
%            A(i,j)=omega * 1/4 * ( A(i-1,j) + Aold(i+1,j) + A(i,j-1) + Aold(i,j+1)) ...
%                + (1-omega)*Aold(i,j);
%        end
%    end
% 
%    errorold=error;
%    error=max(max(abs(A-Aold)));
%    errornew=error;
%    spectral(itr)=errornew/errorold;
%    pitr=pitr+1;
%    if pitr==5
%        figure(1); subplot(1,2,2);
%        contour(A);
%        pitr=0;
%    end   
% end
%  
%  
% fprintf('SOR Iterations %d\n',itr);
% figure(1); subplot(1,2,2);
% contour(A);
% title(['SOR: # of Iterations ' num2str(itr)]);
% figure(2); subplot(1,2,2);
% plot(spectral);
% fprintf('Spectral Radius %f\n',mean(spectral(itr-5:itr)));
% fprintf('Analytical Spectral Radius %f\n',omega-1);

