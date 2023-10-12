
% JACOBI METHOD
clear;
h=0.05;
y=[0:.05:1]';
 
Aold=zeros(21,21);
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
       figure(1);
       contour(A);
       pitr=0;
   end
end
 
fprintf('Jacobi Iterations %d\n',itr);
figure(1);
contour(A);
figure(2);
plot(spectral);
fprintf('Spectral Radius %f\n',mean(spectral(itr-5:itr)));
