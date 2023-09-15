% ---------------- 2-D FDM with Jacobi ----------------
figure(1); subplot(1,2,1);subplot(1,2,2);
Tol=1e-5;
x=0:.05:1;
n = length(x);
h=0.05;
V=zeros(n,n);
Vnew=V; V_prev=V;

V(1,:)=1;
V(n,:)=1;
V(:,n)=1;
V(:,1)=cos(2*pi*x);

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
       Err = norm(abs(V-Vnew), Inf);
       Err_prev = norm(abs(V_prev-V), Inf);
       spectral_Jacobi = max(spectral_Jacobi, Err/Err_prev  );
    end 
  V_prev(2:n-1, 2:n-1) = V(2:n-1, 2:n-1);
  V(2:n-1, 2:n-1)=Vnew(2:n-1, 2:n-1);
  if Pitr==20
      Pitr=0;
      figure(1);subplot(1,2,1);
        [X,Y] = meshgrid(x,x);
        contour(X,Y,V);xlabel('Position X (m)');ylabel('Position Y (m)');drawnow;
  end
end
title(['# of Iterations ' num2str(itr)]);
%%
index_x_V = find(x==0.7);
fprintf(' V(x==0.7, x==0.7): %.12f\n', V(index_x_V, index_x_V))
theroretical_spectral_r_Jacobi = 1 - (pi* 0.05/(2 * x(end)^2));
fprintf('theroretical_spectral_r_Jacobi %.12f \n', theroretical_spectral_r_Jacobi)
fprintf('numerical_spectral_Jacobi %.12f \n', spectral_Jacobi)
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
        Vnew(i,j)= 0.25 * ( Vnew(i-1, j) + Vnew(i+1, j) + Vnew(i, j+1)+ Vnew(i, j-1) );
      end
  end
    if itr_g > 1
       Err = norm(abs(V-Vnew), Inf);
       Err_prev = norm(abs(V_prev-V), Inf);
       spectral_Gauss = max(spectral_Gauss, Err/Err_prev  );
    end 
  V_prev(2:n-1, 2:n-1) = V(2:n-1, 2:n-1);
  V(2:n-1, 2:n-1)=Vnew(2:n-1, 2:n-1);
  if Pitr==20
      Pitr=0;
      figure(1);subplot(1,2,2);
        [X,Y] = meshgrid(x,x);
        contour(X,Y,V);xlabel('Position X (m)');ylabel('Position Y (m)');drawnow;
  end
end
title(['# of Iterations ' num2str(itr_g)]);
%%
index_x_V = find(x==0.7);
fprintf(' V(x==0.7, x==0.7): %.12f\n', V(index_x_V, index_x_V))
theroretical_spectral_r_Gauss = (theroretical_spectral_r_Jacobi)^2;
fprintf('theroretical_spectral_r_Gauss %.12f \n', theroretical_spectral_r_Gauss)
fprintf('numerical_spectral_r_Gauss %.12f \n', spectral_Gauss)
