%Project ADI%
clc
clear all


% initializing the values


ax =0;bx=2*pi;
ay=0; by=2*pi;
nx =30;ny=30;
D=1;
nt=140;
dx = (bx-ax)/nx;
dy = (by-ay)/ny;
dt=0.06;


x=ax:dx:bx;
y=ay:dy:by;
nx=nx+1;
ny=ny+1;



U=zeros(ny,nx,nt);%3d matrix to save time data for y direction
U_x=zeros(ny,nx);
lamdax= D*dt/(dx*dx);%calcualting Sx
lamday= D*dt/(dy*dy);%calcualting Sy


%%%setting boundary conditions
U_x(1,:) = (y-ay).^(2).*sin(pi.*(y-ay)./(2*(by-ay)));
U_x(nx,:) = ((cos(pi.*(y-ay))-1).*cosh(by-y));


%making x Tridiagonal
a=lamdax/2;
b=-(1+lamdax);
c=lamdax/2;

TriDiagx=b*diag(ones(nx,1))+a*diag(ones(nx-1,1),-1)+c*diag(ones(nx-1,1),+1);

TriDiagx(1,1:2)=[1 0];

TriDiagx(nx,nx-1:nx)=[ 0 1];

%making y Tridiagonal
a1=lamday/2;
b1=-(1+lamday);
c1=lamday/2;

TryDiagy=b1*diag(ones(ny,1))+a1*diag(ones(ny-1,1),-1)+c1*diag(ones(ny-1,1),+1);

TryDiagy(1,1:2)=[1 0];

TryDiagy(ny,ny-1:ny)=[ 0 1];

d1=-lamdax/2;
e1=-(1-lamdax);
f1=-lamdax/2;



%Time loop
for n=1:nt
    
    %Visualisation
    h = surf(x,y,U_x','EdgeColor','none');
    colorbar
    shading interp
    %axis ([0 2*pi 0 2*pi -180 70])
    title({['time (\itt) = ',num2str(n*dt)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('Solution to Diffution Equation \rightarrow')
    drawnow;
    refreshdata(h)
    
    fx=zeros(nx,1);
    fx(1,1)=0;
    fx(nx,1)=0;
    
    % x sweep<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    for j=2:ny-1
        %coefficient of right hand side vector
        d=-lamday/2;
        e=-(1-lamday);
        f=-lamday/2;
        for i=2:nx-1
            fx(i)=d*U_x(j-1,i)+e*U_x(j,i)+f*U_x(j+1,i);   
        end
        % Solve the diagonal matrix.
        ut = TriDiagx\fx;
        U_x(j,:)=ut;
    end
 
    %y sweep<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   
    %declaring the vector
    fy=zeros(ny,1);
    fy(1,1)= 0;
    fy(ny,1)=0;
    U_y=zeros(ny,nx);
    %setting boundary condition
    U_y(1,:) = (y-ay).^(2).*sin(pi.*(y-ay)./(2*(by-ay)));
    U_y(nx,:) = ((cos(pi.*(y-ay))-1).*cosh(by-y));
    for i=2:nx-1
        for j=2:ny-1
            fy(j)=d*U_x(j,i-1)+e*U_x(j,i)+f*U_x(j,i+1);
        end
        % Solve the diagonal matrix.
        ut = TryDiagy\fy;
        U_y(:,i)=ut;
    end
    U(:,:,n)=U_y;
end


