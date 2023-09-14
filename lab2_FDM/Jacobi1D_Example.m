% ---------------- 2-D FDM with Jacobi ----------------
figure(1); subplot(1,2,1);subplot(1,2,2);
Tol=1e-6;
x=0:.05:1;
n = length(x);
h=0.05;
V=zeros(n,n);
Vnew=V;

V(1,:)=1;
V(n,:)=1;
V(:,n)=1;
V(:,1)=cos(2*pi*x);

Vnew(1,:)=1;
Vnew(n,:)=1;
Vnew(:,n)=1;
Vnew(:,1)=cos(2*pi*x);

Err=max(abs(Vnew));
itr=0;
Pitr=0;

while Err>Tol
  itr=itr+1;
  Pitr=Pitr+1;
  Vnew(2:n-1, 2:n-1)= 0.25 * ( ...
      V(1:n-2, 2:n-1) + V(3:n, 2:n-1) + ...
      V(2:n-1, 3:n) + V(2:n-1, 1:n-2));
  Err=max(sum(abs(V-Vnew)));
  V(2:n-1, 2:n-1)=Vnew(2:n-1, 2:n-1);  
  if Pitr==5
      Pitr=0;
      figure(1);subplot(1,2,1);
        [X,Y] = meshgrid(x,x);
        contour(X,Y,V);xlabel('Position X (m)');ylabel('Position Y (m)');drawnow;
  end
end
title(['# of Iterations ' num2str(itr)]);


% --------------- Gauss Seidel --------------- 

x=0:.05:1;
n = length(x);
h=0.05;
V=zeros(n,n);
Vnew=V;

V(1,:)=1;
V(n,:)=1;
V(:,n)=1;
V(:,1)=cos(2*pi*x);

Vnew(1,:)=1;
Vnew(n,:)=1;
Vnew(:,n)=1;
Vnew(:,1)=cos(2*pi*x);

Err=max(abs(Vnew));
itr=0;
Pitr=0;

while Err>Tol
  itr=itr+1;
  Pitr=Pitr+1;
  for i=2:n-1
      for j=2:n-1
        Vnew(i,j)= 0.25 * ( V(i-1, j) + V(i+1, j) + V(i, j+1)+ V(i, j-1) );
      end
  end
  Err=max(sum(abs(V-Vnew)));
  V(2:n-1, 2:n-1)=Vnew(2:n-1, 2:n-1);  
  if Pitr==20
      Pitr=0;
      figure(1);subplot(1,2,2);
        [X,Y] = meshgrid(x,x);
        contour(X,Y,V);xlabel('Position X (m)');ylabel('Position Y (m)');drawnow;
  end
end
title(['# of Iterations ' num2str(itr)]);
