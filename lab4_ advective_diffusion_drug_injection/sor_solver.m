function [x]=sor_solver(A,b,xo,omega,itrmax,TOL)

k=1;
n=length(b);
x=zeros(n,1);
while k <= itrmax
    for i=1:n
        S1=0;
        S2=0;
        if i>1
          S1=sum(A(i,1:i-1)'.*x(1:i-1));
        end
        if i<n
          S2=sum(A(i,i+1:n)'.*xo(i+1:n));
        end
        x(i)=(1.0-omega)*xo(i)+ (omega*(-S1-S2+b(i)))/A(i,i);
    end
    error=max(abs(x-xo));
    if error < TOL
       fprintf('Solution convergence in %d iterations\n',k);        
       return;           
    end    
    k=k+1;
    xo=x;
end
fprintf('Maximum number of iterations exceeded');
