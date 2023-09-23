
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
       figure(3);
       contour(A);
       pitr=0;
   end   
end
 
 
fprintf('Gauss Seidel Iterations %d\n',itr);
figure(3);
contour(A);
figure(4);
plot(spectral);
fprintf('Spectral Radius %f\n',mean(spectral(itr-5:itr)));


