h=0.05;
y=[0:.05:1]';
Aold=zeros(21,41);
for i=1:41
 Aold(21,i)=3;
 Aold(1,i)=1;
end
A=Aold;
w=1.785;
error=1;
itr=0;
while (error > 1e-5 & itr < 10000)
 itr=itr+1;
 for i=2:20 
     j=1;
     A(i,j)=w/4*(2*Aold(i,2)+Aold(i+1,j)+A(i-1,j))+(1-w)*Aold(i,j);
     for j=2:40
         if i>=7 & i<=15 & j>=17 & j <=25
         A(i,j)=0;
         else 
         A(i,j)=w/4*(Aold(i+1,j)+Aold(i,j+1)+A(i-1,j)+A(i,j-1))+(1-w)*Aold(i,j); 
        end
     end
     j=41;
     A(i,j)=w/(4-2*h)*(2*A(i,j-1)+A(i-1,j)+Aold(i+1,j)+2*h)+ (1-w)*Aold(i,j);
 end 
 errorold=error;
 error=max(max(abs(A-Aold)))
 errornew=error;
 spectral(itr)=errornew/errorold;
 Aold=A;
end
xi=[0:.05:2];
yi=[0:0.05:1];
figure(1);
contour(xi,yi,A,30)



Vy=-(A(2:20,3:41)-A(2:20,1:39))/2/h; % column x
Vx=(A(3:21,2:40)-A(1:19,2:40))/2/h;  % row y
 
c1=0; 
for i=1:19
 for j=1:39
 c1=c1+1; 
 xpos(c1)=j*.05;
 ypos(c1)=i*.05;
 vx(c1)=Vx(i,j);
 vy(c1)=Vy(i,j); 
 end
 end
figure(2);
quiver(xpos,ypos,vx,vy)