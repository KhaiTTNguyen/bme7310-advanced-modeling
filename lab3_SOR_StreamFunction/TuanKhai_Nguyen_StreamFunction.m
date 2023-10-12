figure(1); subplot(1,2,1);subplot(1,2,2); 

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
    % for i=2:(length(y)-1)
    %     for j=1:length(x)-1
    %         if j==1
    %             A(i,j) = 1/4 * (A(i-1, j) + Aold(i+1, j) + 2*Aold(i, j+1));
    %         elseif j==length(x)
    %             % A(i,j) = 1/(4+2*h) * ( A(i-1, j) + Aold(i+1, j) + 2*A(i, j-1) - 2*h);
    %             A(i,j)= w/(4-2*h)*(2*A(i,j-1)+A(i-1,j)+Aold(i+1,j)+2*h);
    %         elseif 7<=i && i<=15 && 17<=j && j<=25
    %             A(i,j) = 0;
    %         else  
    %             A(i,j)=1/4*(A(i-1,j)+Aold(i+1,j)+A(i,j-1)+Aold(i,j+1));
    %         end
    %     end
    % end
    for i=2:20 
     j=1;
     A(i,j)=1/4*(2*Aold(i,2)+Aold(i+1,j)+A(i-1,j));
     for j=2:40
         if i>=7 & i<=15 & j>=17 & j <=25
         A(i,j)=0;
         else 
         A(i,j)=1/4*(Aold(i+1,j)+Aold(i,j+1)+A(i-1,j)+A(i,j-1)); 
        end
     end
     j=41;
     A(i,j)=1/(4-2*h)*(2*A(i,j-1)+A(i-1,j)+Aold(i+1,j)+2*h);
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

% figure(2), clf
% plot(spectral);
% gauss_spectral = mean(spectral(itr-5:itr));
% fprintf('Spectral Radius %f\n',gauss_spectral);

fprintf('Value at point  x=0.5, y=0.5 %.12f \n', A(10,10))
%%
Vx=zeros(length(y), length(x));
Vy=zeros(length(y), length(x));

Vy(:, 2:end-1)= -( A(:, 3:end) - A(:, 1:end-2) ) ./ (2*h);
Vx(2:end-1, 2:end-1) = ( A(3:end, 2:end-1) - A( 1:end-2, 2:end-1) ) ./ (2*h);

% Vy=-(A(2:20,3:41)-A(2:20,1:39))/2/h;
% Vx=(A(3:21,2:40)-A(1:19,2:40))/2/h;
% 
% c1=0; 
% for i=1:19
%  for j=1:39
%  c1=c1+1; 
%  xpos(c1)=j*.05;
%  ypos(c1)=i*.05;
%  vx(c1)=Vx(i,j);
%  vy(c1)=Vy(i,j); 
%  end
%  end
figure(3);
% quiver(xpos,ypos,vx,vy)
quiver(x, y, Vx, Vy)
title('Velocity vector')

