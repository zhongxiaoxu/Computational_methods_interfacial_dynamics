% ------------------- Project background --------------------
% A droplet (rho1) settles down in an other liquid (rho2)
% Small particles distributed in both liquid with same diffusion coefficient
% Insoluble surfactants distributed at the droplet surface
% The surface tension is described by Langmuir-relation


clear all; clc; close all;
tic

% base_path = '/Users/zhongxiaoxu/OneDrive - purdue.edu/study/course/Computational_methods_for_interface_dynamics/project/code';
% folder_name = 'case11';
% mkdir case11

% --------------------- dimensionless parameters -----------------
gamma0 = 10; % dimensionless surface tension coefficient = gamma0 / L / rho2 / U^2
Pe_Gamma = 1000; % Peclet number for surfactant concentration
c0 = 0.3; % initial particle concentration
Gamma0 = 0.01; % initial surfactant concentration

Lx = 1.0; Ly = 1.0;                                                  % domain size 
gx = 0.0; gy = -100.0; rho1 = 2; rho2 = 1; 
Re = 100; 
Pe_c = 100; % Peclet number for particle concentration
SIGMA = 1; % surfactant activity = Gamma_m * Rg * T / gamma0
unorth = 0; usouth = 0; veast = 0; vwest = 0;                        % boundary conditions
rad = 0.15; xc = 0.5; yc = 0.7*Ly;   % initial drop size and location

% ---------------------- initialization --------------------
time = 0.0; 

nx = 80; ny = 80; dx = Lx/nx; dy = Ly/ny; dt = 1e-4;
nstep = round(0.4/dt); maxit = 200; maxError = 0.001; omg = 1.5; Nf = 100;
export_fre = round(0.01 / dt); plot_freq = round(0.01/dt);

dropCenter = zeros(3, nstep+1); % [t; xc; yc]
dropCenter(2, 1) = xc; dropCenter(3,1) = yc;

u=zeros(nx+1,ny+2); ut = u   ; uplot = zeros(nx+1,ny+1);
v=zeros(nx+2,ny+1); vt = u   ; vplot = zeros(nx+1,ny+1);
p=zeros(nx+2,ny+2); tmp1 = p ; tmp2  = p; r = p; chi = p; 

xf=zeros(1,Nf+2); yf=zeros(1,Nf+2); 
uf=zeros(1,Nf+2); vf=zeros(1,Nf+2);

xh = linspace(0,Lx,nx+1)         ; yh = linspace(0,Ly,ny+1);                % velocity points
x  = linspace(-dx/2,Lx+dx/2,nx+2); y  = linspace(-dy/2,Ly+dy/2,ny+2);       % pressure points

% --------------------- initialize c -------------------
c = p;
XC = zeros(nx, ny); YC = zeros(nx, ny); CC = zeros(nx, ny);              % concentration points
% c(2:nx+1, 2:ny+1) = c(2:nx+1, 2:ny+1) + c0;
for i = 2:nx+1
    for j = 2 : ny+1
%         c(i,j) = c0;
        c(i,j) = 0.1 + 0.7 * (j-0.5) * dy;
    end
end
c(:, 1) = c(:,2); c(1,:) = c(2, :); c(:, ny+2) = c(:, ny+1); c(nx+2, :) = c(nx+1, :);

for i = 1 : nx
    for j = 1 : ny
        XC(i,j) = (i-0.5) * dx;  YC(i,j) = (j-0.5) * dy;
    end
end

% ---------------------- initialize surfactant concentration Gamma ------------
Gamma = zeros(1, Nf+2) + Gamma0; % Gamma locates at the middle point of the front points
delta_s = zeros(1, Nf+2); % arc length occupied by Gamma


% ---------------------------------------------------------------------------------
r = zeros(nx+2,ny+2) + rho2;                                                % initial density 
for i = 2:nx+1; for j = 2:ny+1                   
  if((x(i)-xc)^2+(y(j)-yc)^2 < rad^2); r(i,j) = rho1; chi(i,j)=1.0; end; 
end; end
                                          
for l=1:Nf+2; xf(l) = xc - rad*sin(2*pi*(l-1)/Nf); yf(l) = yc + rad*cos(2*pi*(l-1)/Nf); end   % initialize front    

 % ---------------- initialize the arc length between front point ----------
 delta_s(1:Nf+1) = sqrt((xf(2:Nf+2) - xf(1:Nf+1)).^2 + (yf(2:Nf+2) - yf(1:Nf+1)).^2);
 delta_s(Nf+2) = delta_s(2);
 
 % ------------ calculate total surfactant on the interface ------
%  disp(['total surfactant is: ', num2str(sum(Gamma(1:Nf) .* delta_s(1:Nf)))]);
            
for is=1:nstep
    
    delta_sold = delta_s;
    Gammaold = Gamma;
    cold = c;
    
    % ----------------- update c ----------------   
    for i = 2: nx+1
        for j = 2 : ny+1
            c(i,j) = cold(i,j)-(0.5*dt/dx)*(u(i,j)*(cold(i+1,j)...
                   + cold(i,j))-u(i-1,j)*(cold(i-1,j)+cold(i,j)))...
                   - (0.5* dt/dy)*(v(i,j)*(cold(i,j+1)...
                   + cold(i,j))-v(i,j-1)*(cold(i,j-1)+cold(i,j))  )...
                   + (1/Pe_c*dt/dx/dx)*(cold(i+1,j)-2.0*cold(i,j)+cold(i-1,j))...
                   + (1/Pe_c*dt/dy/dy)*(cold(i,j+1)-2.0*cold(i,j)+cold(i,j-1));
        end
    end
    c(:, 1) = c(:,2); c(1,:) = c(2, :); c(:, ny+2) = c(:, ny+1); c(nx+2, :) = c(nx+1, :);
    % ---------------------------------------------

 for l=2:Nf+1                       % interpolate front velocities
   ip = floor(xf(l)/dx)+1; jp = floor((yf(l)+0.5*dy)/dy)+1;
   ax = xf(l)/dx-ip+1; ay = (yf(l)+0.5*dy)/dy-jp+1;	   
   uf(l) = (1.0-ax)*(1.0-ay)*u(ip,jp) + ax*(1.0-ay)*u(ip+1,jp) + (1.0-ax)*ay*u(ip,jp+1) + ax*ay*u(ip+1,jp+1);
   ip = floor((xf(l)+0.5*dx)/dx)+1; jp = floor(yf(l)/dy)+1;
   ax = (xf(l)+0.5*dx)/dx-ip+1; ay = yf(l)/dy-jp+1;
   vf(l) = (1.0-ax)*(1.0-ay)*v(ip,jp) + ax*(1.0-ay)*v(ip+1,jp) + (1.0-ax)*ay*v(ip,jp+1) + ax*ay*v(ip+1,jp+1);
 end     

 xf(2:Nf+1) = xf(2:Nf+1) + dt*uf(2:Nf+1); yf(2:Nf+1) = yf(2:Nf+1) + dt*vf(2:Nf+1);  % move front
 xf(1) = xf(Nf+1); yf(1) = yf(Nf+1); xf(Nf+2) = xf(2); yf(Nf+2) = yf(2); 
 
 % ---------------- calculate the arc length between front point ----------
 
 delta_s(1:Nf+1) = sqrt((xf(2:Nf+2) - xf(1:Nf+1)).^2 + (yf(2:Nf+2) - yf(1:Nf+1)).^2); % update delta_s
 delta_s(Nf+2) = delta_s(2);
 
 % ------------------ update surfactant concentration Gamma ------------------
 
 for l = 2 : Nf + 1
     Gamma(l) = 1 / delta_s(l) * (Gammaold(l) * delta_sold(l) + 2 * dt / Pe_Gamma * (...
         (Gammaold(l+1) - Gammaold(l)) / (delta_sold(l) + delta_sold(l+1)) - (Gammaold(l) - Gammaold(l-1)) / (delta_sold(l) + delta_sold(l-1))));
 end
 Gamma(1) = Gamma(Nf+1); Gamma(Nf+2) = Gamma(2);
 
 
 % --------------------------------------------------------------------------------
  	
 d(2:nx+1,2:ny+1)=4*Lx;                                        % placeholder distance to front
 for l=2:Nf+1
   nfx = -(yf(l+1)-yf(l)); nfy =  (xf(l+1)-xf(l));             % outer normal vector
   ds  = sqrt(nfx*nfx+nfy*nfy); nfx = nfx/ds; nfy = nfy/ds;    % unit outer normal
   xfront = 0.5*(xf(l)+xf(l+1)); yfront = 0.5*(yf(l)+yf(l+1)); % find midpoint in the segment
   ip=floor((xfront+0.5*dx)/dx)+1; jp=floor((yfront+0.5*dy)/dy)+1;

   d1=sqrt(((xfront-  x(ip)))^2+((yfront-  y(jp)))^2);
   d2=sqrt(((xfront-x(ip+1)))^2+((yfront-  y(jp)))^2);
   d3=sqrt(((xfront-x(ip+1)))^2+((yfront-y(jp+1)))^2);
   d4=sqrt(((xfront-  x(ip)))^2+((yfront-y(jp+1)))^2);

   if d1<d(ip,jp)
    d(ip,jp) = d1; dn1 = (x(ip) - xfront)*nfx/dx + (y(jp) - yfront)*nfy/dy;
    chi(ip,jp) = 0.5*(1 + sign(dn1)); 
    if abs(dn1) < 0.5; chi(ip,jp) = 0.5 + dn1; end;
   end
   if d2<d(ip+1,jp)
    d(ip+1,jp) = d2; dn2 = (x(ip+1) - xfront)*nfx/dx + (y(jp) - yfront)*nfy/dy;
    chi(ip+1,jp) = 0.5*(1 + sign(dn2));
    if abs(dn2) < 0.5; chi(ip+1,jp) = 0.5 + dn2; end;
   end
   if d3<d(ip+1,jp+1)
    d(ip+1,jp+1) = d3; dn3 = (x(ip+1)-xfront)*nfx/dx+(y(jp+1)-yfront)*nfy/dy;
    chi(ip+1,jp+1)=0.5*(1 + sign(dn3)); 
    if abs(dn3) < 0.5; chi(ip+1,jp+1) = 0.5 + dn3; end;
   end
   if d4<d(ip,jp+1)
    d(ip,jp+1) = d4; dn4 = (x(ip) - xfront)*nfx/dx + (y(jp+1) - yfront)*nfy/dy;
    chi(ip,jp+1) = 0.5*(1 + sign(dn4)); 
    if abs(dn4) < 0.5; chi(ip,jp+1) = 0.5 + dn4; end;
   end
  end
              
  ro = r;
  r  = rho1*chi + rho2*(1-chi);  % obtain density from charact func

  for l=1:Nf+1
    ds=sqrt((xf(l+1)-xf(l))^2+(yf(l+1)-yf(l))^2);
    tx(l)=(xf(l+1)-xf(l))/ds; ty(l)=(yf(l+1)-yf(l))/ds; % unit tangent vectors
  end
  tx(Nf+2)=tx(2); ty(Nf+2)=ty(2);
  

  fgx=zeros(nx+2,ny+2); fgy=zeros(nx+2,ny+2);
  for l=2:Nf+1                                               % distribute to the fixed grid
    fglx = gamma0 * (1 + SIGMA * log(1- Gamma(l))) * tx(l) - gamma0* (1 + SIGMA * log(1- Gamma(l-1)) ) * tx(l-1);
    fgly = gamma0 * (1 + SIGMA * log(1- Gamma(l))) * ty(l) - gamma0* (1 + SIGMA * log(1- Gamma(l-1)) )  * ty(l-1);
    ip   = floor(xf(l)/dx)+1;     jp   = floor((yf(l)+0.5*dy)/dy)+1;
    ax   = xf(l)/dx-ip+1;         ay   = (yf(l)+0.5*dy)/dy-jp+1;
    fgx(ip,jp)     = fgx(ip,jp)     + (1.0-ax)*(1.0-ay)*fglx/dx/dy;
    fgx(ip+1,jp)   = fgx(ip+1,jp)   + ax*(1.0-ay)*fglx/dx/dy;
    fgx(ip,jp+1)   = fgx(ip,jp+1)   + (1.0-ax)*ay*fglx/dx/dy;
    fgx(ip+1,jp+1) = fgx(ip+1,jp+1) + ax*ay*fglx/dx/dy;

    ip = floor((xf(l)+0.5*dx)/dx)+1; jp = floor(yf(l)/dy)+1;
    ax = (xf(l)+0.5*dx)/dx-ip+1;     ay = yf(l)/dy-jp+1;	  
    fgy(ip,jp)     = fgy(ip,jp)     + (1.0-ax)*(1.0-ay)*fgly/dx/dy;
    fgy(ip+1,jp)   = fgy(ip+1,jp)   + ax*(1.0-ay)*fgly/dx/dy;
    fgy(ip,jp+1)   = fgy(ip,jp+1)   + (1.0-ax)*ay*fgly/dx/dy;
    fgy(ip+1,jp+1) = fgy(ip+1,jp+1) + ax*ay*fgly/dx/dy;  
  end

  fgx(1:nx+2,2) = fgx(1:nx+2,2) + fgx(1:nx+2,1); fgx(1:nx+2,ny+1) = fgx(1:nx+2,ny+1) + fgx(1:nx+2,ny+2);  % bring all forces to interior
  fgy(2,1:ny+2) = fgy(2,1:ny+2) + fgy(1,1:ny+2); fgy(nx+1,1:ny+2) = fgy(nx+1,1:ny+2) + fgy(nx+2,1:ny+2);  
  
  u(1:nx+1,1) = 2*usouth-u(1:nx+1,2); u(1:nx+1,ny+2) = 2*unorth-u(1:nx+1,ny+1); % tangential vel BC
  v(1,1:ny+1) = 2*vwest -v(2,1:ny+1); v(nx+2,1:ny+1) = 2*veast -v(nx+1,1:ny+1); % tangential vel BC

  for i=2:nx; for j=2:ny+1     % temporary u-velocity (boundary values are not touched)
    ut(i,j) = (2.0/(r(i+1,j)+r(i,j)))*(0.5*(ro(i+1,j)+ro(i,j))*u(i,j)+ dt* (...
            - (0.25/dx)*(ro(i+1,j)*(u(i+1,j)+u(i,j))^2-ro(i,j)*(u(i,j)+u(i-1,j))^2)...
            - (0.0625/dy)*( (ro(i,j)+ro(i+1,j)+ro(i,j+1)+ro(i+1,j+1))*(u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) ...
            - (ro(i,j)+ro(i+1,j)+ro(i+1,j-1)+ro(i,j-1))*(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))...
            + 1/Re / (1- 0.5 * (c(i,j) + c(i+1,j)))^2 *((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+ (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2)...
            + 0.5*(ro(i+1,j)+ro(i,j))*gx + fgx(i,j) ) );
  end; end

  for i=2:nx+1; for j=2:ny       % temporary v-velocity (boundary values are not touched)
    vt(i,j) = (2.0/(r(i,j+1)+r(i,j)))*(0.5*(ro(i,j+1)+ro(i,j))*v(i,j)+ dt* (...     
            - (0.0625/dx)*( (ro(i,j)+ro(i+1,j)+ro(i+1,j+1)+ro(i,j+1))*(u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) ...
            - (ro(i,j)+ro(i,j+1)+ro(i-1,j+1)+ro(i-1,j))*(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )...                                 
            - (0.25/dy)*(ro(i,j+1)*(v(i,j+1)+v(i,j))^2-ro(i,j)*(v(i,j)+v(i,j-1))^2 )...
            + 1/Re / (1 - 0.5 * (c(i,j) + c(i,j+1)))^2 *((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2)...
            + 0.5*(ro(i,j+1)+ro(i,j))*gy + fgy(i,j) ) );    
  end; end     

  for i = 2:nx+1; for j = 2:ny+1
    tmp1(i,j) = (0.5/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy );
    tmp2(i,j) =1/( (1/dx)*(1/(dx*(r(i+1,j)+r(i,j)))+ 1/(dx*(r(i-1,j)+r(i,j))) )+ ...
                   (1/dy)*(1/(dy*(r(i,j+1)+r(i,j)))+ 1/(dy*(r(i,j-1)+r(i,j))) )   );
  end; end

  for it = 1:maxit	               % solve for pressure by SOR
    pold   = p;
    p(1,:) = p(2,:); p(nx+2,:) = p(nx+1,:); p(:,1) = p(:,2); p(:,ny+2) = p(:,ny+1); % set gosht values
    for i=2:nx+1; for j=2:ny+1
      p(i,j) = (1.0-omg)*p(i,j) + omg*tmp2(i,j)*(        ...
      (1/dx)*( p(i+1,j)/(dx*(r(i+1,j)+r(i,j)))+ p(i-1,j)/(dx*(r(i-1,j)+r(i,j))) )+    ...
      (1/dy)*( p(i,j+1)/(dy*(r(i,j+1)+r(i,j)))+ p(i,j-1)/(dy*(r(i,j-1)+r(i,j))) ) - tmp1(i,j));
    end; end
    if max(max(abs(pold-p))) < maxError; break;  end
  end
                                      
  for i=2:nx; for j=2:ny+1   % correct the u-velocity 
    u(i,j)=ut(i,j)-dt*(2.0/dx)*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j));
  end; end
      
  for i=2:nx+1; for j=2:ny   % correct the v-velocity 
    v(i,j)=vt(i,j)-dt*(2.0/dy)*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j));
  end; end

  xfold=xf; yfold=yf; j=1;   % add or/and delete points in the front
  Gammaold = Gamma; delta_sold = delta_s;
  sum_Gamma = 0; sum_delta_s = 0;
  for l=2:Nf+1
    ds=sqrt( ((xfold(l)-xf(j))/dx)^2 + ((yfold(l)-yf(j))/dy)^2);
    if (ds > 0.5)                                    % add point
      j=j+1; xf(j) = 0.5*(xfold(l)+xf(j-1)); yf(j) = 0.5*(yfold(l)+yf(j-1)); 
      Gamma(j-1) = (sum_Gamma + Gammaold(l-1) * delta_sold(l-1)) / (sum_delta_s + delta_sold(l-1));
      j=j+1; xf(j) = xfold(l); yf(j) = yfold(l); 
      Gamma(j-1) = Gamma(j-2);
      sum_delta_s = 0; sum_Gamma = 0;    
    elseif (ds < 0.25)                               % skip point
        sum_delta_s = sum_delta_s + delta_sold(l-1);
        sum_Gamma = sum_Gamma + Gammaold(l-1) * delta_sold(l-1);
    else
      j = j+1; xf(j) = xfold(l); yf(j) = yfold(l); 
      Gamma(j-1) = (sum_Gamma + Gammaold(l-1) * delta_sold(l-1)) / (sum_delta_s + delta_sold(l-1));  % copy point  
      sum_delta_s = 0; sum_Gamma = 0;
    end    
  end
  Nf = j-1; xf(1) = xf(Nf+1); yf(1) = yf(Nf+1); xf(Nf+2) = xf(2); yf(Nf+2) = yf(2);
  Gamma(Nf+1) = Gamma(1); Gamma(Nf+2) = Gamma(2);
  delta_s(1:Nf+1) = sqrt((xf(2:Nf+2) - xf(1:Nf+1)).^2 + (yf(2:Nf+2) - yf(1:Nf+1)).^2); % update delta_s
 delta_s(Nf+2) = delta_s(2);

  time = time+dt     
  
  % get center of drop
  
      
dropCenter(1, is +1) = time;
temp1 = sum(sum(chi(2:nx+1, 2:ny+1)));
dropCenter(2, is+1) = sum(sum(chi(2:nx+1, 2:ny+1) .* XC(1:nx, 1:ny)))  / temp1; % x coordinate of center
dropCenter(3, is+1) = sum(sum(chi(2:nx+1, 2:ny+1) .* YC(1:nx, 1:ny)))  / temp1; % y coordinate of center
    
    % get surfactant concentration and its position
    sc = zeros(2, Nf);
    for i = 1 : Nf
        temp1 = (yf(i) + yf(i+1)) * 0.5 - dropCenter(3, is+1);
        temp2 = (xf(i) + xf(i+1)) * 0.5 - dropCenter(2, is+1);
        THETA = acos(temp1 / sqrt(temp1^2 + temp2^2));
        if (xf(i) + xf(i+1)) * 0.5 - dropCenter(2, is+1) > 0
            sc(1, i) = THETA;
        else
            sc(1, i) = 2 * pi - THETA;
        end
        
        sc(2, i) = Gamma(i);
    end
  
  
  if (mod(is,plot_freq)==0) | (is==1);                         % plot solution
      
%----------------- save parameters --------------------
%         if is == 1
%             para = [c0, Gamma0, Pe_Gamma, gamma0, dt, nx, ny];
%             para_name = 'para.mat';
%             save(fullfile(base_path, folder_name, para_name), 'para');
%         end
        
%         data_file_name = sprintf('Gamma_%d.mat', round(is / plot_freq));
%         save(fullfile(base_path, folder_name, data_file_name), 'sc');

%         data_file_name = sprintf('xf_yf_%d.mat', round(is / plot_freq));
%         frontshape = zeros(2, Nf); frontshape(1,:) = xf(1:Nf); frontshape(2,:) = yf(1:Nf);
%         save(fullfile(base_path, folder_name, data_file_name), 'frontshape');
      
      % -------------- get velocity field ----------------
      uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
      vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));

      for i = 1 : nx
          for j = 1 : ny
              CC(i, j) = c(i+1, j+1);
          end
      end
      
%       figure(1);
%       contour(x,y,r'); axis equal; axis([0 Lx 0 Ly]); hold on;
%       title([sprintf('t=%0.2f', time)], 'FontSize', 20);
%       quiver(xh,yh,uu',vv','r'); hold on; plot(xf(1:Nf),yf(1:Nf),'k','linewidth',3);
%       xlabel('x', 'FontSize', 16); ylabel('y', 'FontSize', 16); set(gca, 'FontSize', 20);
%       fig= gcf;  fig.Color = 'white';
%       
%           
%       drawnow;
%       hold off
      
      % ---------------- plot surfactant concentration along interface ----
%       figure(1);
%       plot(sc(1,:), sc(2,:));
%       xlim([0 2*pi]); ylim([Gamma0 - 0.1, Gamma0 + 0.1]);
%       title([sprintf('t=%0.2f', time)], 'FontSize', 20);
%       xlabel('theta', 'FontSize', 16); ylabel('Gamma', 'FontSize', 16); set(gca, 'FontSize', 20);
%       drawnow;
      
      
      
   % ------------- contour plot of particle concentration -------
      contourf(XC, YC, CC); 
      xlabel('x', 'FontSize', 16); ylabel('y', 'FontSize', 16); set(gca, 'FontSize', 20);
      title(['Concentration at ', sprintf('t=%0.2f', time)], 'FontSize', 16);
      colormap jet; colorbar; axis tight;
      caxis([0.1, 0.8]);
      fig = gcf;  fig.Color = 'white';
      drawnow;
      
%       fig_name = sprintf('shape_%d.png', round(is / plot_freq));
%       saveas(gcf, fullfile(base_path, folder_name, fig_name));
      
  end

end

% data_file_name = sprintf('centroid%d.mat');
% dropCenter1 = dropCenter(:, 1 : round(0.001 / dt) : end);
% save(fullfile(base_path, folder_name, data_file_name), 'dropCenter1');



toc
