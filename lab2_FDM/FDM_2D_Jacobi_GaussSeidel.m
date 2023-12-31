% ---------------- 2-D FDM with Jacobi ----------------
figure(1); subplot(1,2,1);subplot(1,2,2);
Tol=1e-5;
x=0:.05:1;
n = length(x);
h=0.05;
V=zeros(n,n);
V(1,:)=1;
V(n,:)=1;
V(:,n)=1;
i=1:1:21;
V(:,1)=cos(2*pi*(i*h-0.05));
Vnew=V; V_prev=V;

Err=max(abs(Vnew));
itr=0;
Pitr=0;
spectral_Jacobi = 0;
while Err>Tol
  itr=itr+1;
  Pitr=Pitr+1;
  Vnew(2:n-1, 2:n-1)= 0.25 * ( ...
      V(1:n-2, 2:n-1) + V(3:n, 2:n-1) + ...
      V(2:n-1, 3:n) + V(2:n-1, 1:n-2));
    if itr > 1
       Err = norm(V-Vnew, Inf);
       Err_prev = norm(V_prev-V, Inf);
       spectral_Jacobi = max(spectral_Jacobi, Err/Err_prev  );
    end
    V_prev=V;
    V=Vnew; 
%   V_prev(2:n-1, 2:n-1) = V(2:n-1, 2:n-1);
%   V(2:n-1, 2:n-1)=Vnew(2:n-1, 2:n-1);
  if Pitr==20
      Pitr=0;
      figure(1);subplot(1,2,1);
        [X,Y] = meshgrid(x,x);
        contourf(X,Y,V);xlabel('Position X (m)');ylabel('Position Y (m)');drawnow;
  end
end
title(['Jacobi: # of Iterations ' num2str(itr)]);
%%
index_x_V = find(x==0.7);
fprintf('Jacobi V(x==0.7, x==0.7): %.12f\n', V(index_x_V, index_x_V))
theorretical_spectral_r_Jacobi = 1 - (pi^2 * h^2 /(2 * x(end)^2));
fprintf('theroretical_spectral_r_Jacobi %.12f \n', theorretical_spectral_r_Jacobi)
fprintf('numerical_spectral_Jacobi %.12f \n', spectral_Jacobi)
fprintf('Theoretical number of iterations for Jacobi at Tol=%f: %.12f \n', Tol, 1/abs(log10(theorretical_spectral_r_Jacobi)))
fprintf('Estimated number of iterations for Jacobi at Tol=%f: %.12f \n', Tol, 1/abs(log10(spectral_Jacobi)))

% --------------- Gauss Seidel --------------- 
%%
x=0:.05:1;
n = length(x);
h=0.05;
V=zeros(n,n);

V(1,:)=1;
V(n,:)=1;
V(:,n)=1;
V(:,1)=cos(2*pi*x);

Vnew=V; V_prev=V;

Err=max(abs(Vnew));
itr_g=0;
Pitr=0;
spectral_Gauss = 0;

while Err>Tol
  itr_g=itr_g+1;
  Pitr=Pitr+1;
  Vnew = V;
  for i=2:n-1
      for j=2:n-1
        Vnew(i,j)= 0.25 *( Vnew(i-1, j) + Vnew(i+1, j) + Vnew(i, j+1)+ Vnew(i, j-1) );
      end
  end
    if itr_g > 1
       Err = norm(V-Vnew, Inf);
       Err_prev = norm(V_prev-V, Inf);
       spectral_Gauss = max(spectral_Gauss, Err/Err_prev  );
    end 
  V_prev(2:n-1, 2:n-1) = V(2:n-1, 2:n-1);
  V(2:n-1, 2:n-1)=Vnew(2:n-1, 2:n-1);
  if Pitr==20
      Pitr=0;
      figure(1);subplot(1,2,2);
        [X,Y] = meshgrid(x,x);
        contourf(X,Y,V);xlabel('Position X (m)');ylabel('Position Y (m)');drawnow;
  end
end
title(['Gauss Seidel: # of Iterations ' num2str(itr_g)]);
%%
index_x_V = find(x==0.7);
fprintf('Gauss Seidel V(x==0.7, x==0.7): %.12f\n', V(index_x_V, index_x_V))
theorretical_spectral_r_Gauss = 1 - (pi^2 * h^2 /(x(end)^2)); % OR (theroretical_spectral_r_Jacobi)^2;
fprintf('theroretical_spectral_r_Gauss %.12f \n', theorretical_spectral_r_Gauss)
fprintf('numerical_spectral_r_Gauss %.12f \n', spectral_Gauss)
%%
fprintf('Theoretical number of iterations for GaussSeidel at Tol=%f: %.12f \n', Tol, 1/abs(log10(theorretical_spectral_r_Gauss)))
fprintf('Estimated number of iterations for GaussSeidel at Tol=%f: %.12f \n', Tol, 1/abs(log10(spectral_Gauss)))

