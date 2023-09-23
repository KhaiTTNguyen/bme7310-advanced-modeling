clear
% ------------------------------------
%perform a series of simulations that enable you to plot iteration count as a function of omega to confirm whether 
% the theoretical optimum is the same as that found in practice for this problem 
interations_v = zeros(1);
itr=0;

for omega=0.8:0.02:1.95
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
           pitr=0;
       end   
    end
    if omega==0.8
        interations_v(end)=itr;
    else
        interations_v(end+1)=itr;
    end
end

figure(5)
plot([0.8:0.02:1.95], interations_v)
xlabel('omega')
ylabel('SOR num iter')
title('SOR number of iterations w.r.t omega')