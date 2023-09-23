figure(1); subplot(1,3,1);subplot(1,3,2);subplot(1,3,3);

figure(2); subplot(1,3,1);subplot(1,3,2);subplot(1,3,3);

% JACOBI METHOD
clear;
h=0.05;
y=[0:.05:1]';
 
Aold=zeros(21,21);
for i=1:21
  Aold(i,21)=1;
  Aold(1,i)=1;
  Aold(21,i)=1;
  Aold(i,1)=cos(2*pi*(i*h-0.05));
end
A=zeros(21,21);
A=Aold;
error=1;
errorold=1;
itr=0;
pitr=0;
while (error > 1e-5 & itr < 10000)
   itr=itr+1;
   A(2:20,2:20)=1/4*(Aold(1:19,2:20)+Aold(3:21,2:20)+Aold(2:20,1:19)+Aold(2:20,3:21));
   errorold=error;
   error=max(max(abs(A-Aold)));
   spectral(itr)=error/errorold;
   Aold=A;
   pitr=pitr+1;
   if pitr==5
       figure(1); subplot(1,3,3);
       contour(A);
       pitr=0;
   end
end
 
fprintf('Jacobi Iterations %d\n',itr);
figure(1); subplot(1,3,1);
contour(A);
title(['Jacobi: # of Iterations ' num2str(itr)]);
figure(2); subplot(1,3,1);
plot(spectral);
fprintf('Spectral Radius %f\n',mean(spectral(itr-5:itr)));

theorretical_spectral_r_Jacobi = 1 - (pi^2 * h^2 /(2 * y(end)^2));
fprintf('Theroretical_spectral_r_Jacobi %.12f \n', theorretical_spectral_r_Jacobi)
index_x_V = find(y==0.7);
fprintf('Jacobi V(x==0.7, x==0.7): %.12f\n', A(index_x_V, index_x_V))
fprintf('-------------------------------\n')

% Gauss-Seidel
clear;
h=0.05;
y=[0:.05:1]';
A=zeros(21,21);
for i=1:21
  A(i,21)=1;
  A(1,i)=1;
  A(21,i)=1;
  A(i,1)=cos(2*pi*(i*h-0.05));
end
 
error=1;
itr=0;
pitr=0;
while (error > 1e-5 & itr < 10000)
   itr=itr+1;
   Aold=A;
   for i=2:20
        for j=2:20
         A(i,j)=1/4*(A(i-1,j)+Aold(i+1,j)+A(i,j-1)+Aold(i,j+1));
     end
    end
   errorold=error;
   error=max(max(abs(A-Aold)));
   errornew=error;
   spectral(itr)=errornew/errorold;
   pitr=pitr+1;
   if pitr==5
       figure(1);subplot(1,3,2);
       contour(A);
       pitr=0;
   end   
end
 
fprintf('Gauss Seidel Iterations %d\n',itr);
figure(1); subplot(1,3,2);
contour(A);
title(['GaussSeidel: # of Iterations ' num2str(itr)]);
figure(2); subplot(1,3,2);
plot(spectral);
gauss_spectral = mean(spectral(itr-5:itr));
fprintf('Spectral Radius %f\n',gauss_spectral);
theorretical_spectral_r_Gauss = 1 - (pi^2 * h^2 /(y(end)^2)); % OR (theroretical_spectral_r_Jacobi)^2;
fprintf('Theroretical_spectral_r_Gauss %.12f \n', theorretical_spectral_r_Gauss)
index_x_V = find(y==0.7);
fprintf('Gauss Seidel V(x==0.7, x==0.7): %.12f\n', A(index_x_V, index_x_V))
fprintf('-------------------------------\n')

% SOR
clearvars -except gauss_spectral 

h=0.05;
y=[0:.05:1]';
A=zeros(21,21);
for i=1:21
  A(i,21)=1;
  A(1,i)=1;
  A(21,i)=1;
  A(i,1)=cos(2*pi*(i*h-0.05));
end
 
omega = 2 / (1 + sqrt(1 - gauss_spectral));
error=1;
itr=0;
pitr=0;
while (error > 1e-5 & itr < 10000)
   itr=itr+1;
   Aold=A;
   for i=2:20
       for j=2:20
           A(i,j)=omega * 1/4 * ( A(i-1,j) + Aold(i+1,j) + A(i,j-1) + Aold(i,j+1)) ...
               + (1-omega)*Aold(i,j);
       end
   end

   errorold=error;
   error=max(max(abs(A-Aold)));
   errornew=error;
   spectral(itr)=errornew/errorold;
   pitr=pitr+1;
   if pitr==5
       figure(1); subplot(1,3,3);
       contour(A);
       pitr=0;
   end   
end
 
%% 
fprintf('SOR Iterations %d\n',itr);
figure(1); subplot(1,3,3);
contour(A);
title(['SOR: # of Iterations ' num2str(itr)]);
figure(2); subplot(1,3,3);
plot(spectral);
fprintf('Spectral Radius %f\n',mean(spectral(itr-5:itr)));
fprintf('Analytical Spectral Radius %.12f\n',omega-1);
index_x_V = find(y==0.7);
fprintf('SOR V(x==0.7, x==0.7): %.12f\n', A(index_x_V, index_x_V))
