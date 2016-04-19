function PWP_Byron_Earle(met_input_file, profile_input_file, pwp_output_file)
%--------------------------------------------------------------------------
% MATLAB implementation of the Price Weller Pinkel mixed layer model
% 19 May 2011
% Byron Kilbourne
% University of Washington
% School of Oceanography
%--------------------------------------------------------------------------
% input (.mat format)
% met_input_file -> path and file name for MET forcing
% profile_input_file -> path and file name for intial density profile
% output (.mat format)
% pwp_output -> path and file name for output
%--------------------------------------------------------------------------
% The 2001 version of PWP was having some tricky bugs related to the timing
% of respnse to the wind stress.  I will rewrite the matlab version for
% MATLAB 2001a
% I will: remove obsolete function calls (eg sw_den0)
%         change structure to remove global variables
%         preallocate array size memory
% Effort has been made to preserve the original format for input and output
% for compatability with the previous version.  The variable names, though
% obtuse, have been preserved as well.
%--------------------------------------------------------------------------
% this code has sections
% 1) set parameters
% 2) read in forcing data
% 3) preallocate variables
% 4) main model loop
%   4.1) heat and salt fluxes
%   4.2) rotate, adjust to wind, rotate
%   4.3) bulk Richardson number mixing
%   4.4) gradient Richardson number mixing
% 5) save results to output file

%--------------------------------------------------------------------------
% diagnostic plots
% 0 -> no plots
% 1 -> shows depth integrated KE and momentum and current,T, and S profiles
% 2 -> plot momentum after model run
% 3 -> plot mixed layer depth after model run
diagnostics = 1;
%--------------------------------------------------------------------------
tic % log runtime
% set parameters

dt			= 3600*3;%900;          %time-step increment (seconds)
dz			= 1;%5;            %depth increment (meters)
% the days parameter has been deprecated
%days 		= 17;           %the number of days to run (max time grid)
depth		= 100;          %the depth to run (max depth grid)
dt_save     = 1;            %time-step increment for saving to file (multiples of dt)
lat 		= 74;%-58.2;        %latitude (degrees)
g			= 9.8;          %gravity (9.8 m/s^2)
cpw			= 4183.3;       %specific heat of water (4183.3 J/kgC)
rb			= 0.65;         %critical bulk richardson number (0.65)
rg			= 0.25;         %critical gradient richardson number (0.25)
rkz			= 0;            %background vertical diffusion (0) m^2/s
beta1   	= 0.6;          %longwave extinction coefficient (0.6 m)
beta2   	= 20;           %shortwave extinction coefficient (20 m)

f = sw_f(lat);              %coriolis term (rad/s)
ucon = (.1*abs(f));         %coefficient of inertial-internal wave dissipation (0) s^-1
%--------------------------------------------------------------------------
% load forcing data

load(met_input_file)
load(profile_input_file)
dtd = dt/86400; 
% days parameter deprecated
time = met.time(1):dtd:met.time(end);
nmet = length(time);
clear dtd
qi = interp1(met.time,met.sw,time); 
qo = interp1(met.time,(met.lw + met.qlat + met.qsens),time); 
tx = interp1(met.time,met.tx,time);
ty = interp1(met.time,met.ty,time);
precip = interp1(met.time,met.precip,time);
% make depth grid
zmax = max(profile.z);
if zmax < depth
    depth = zmax;
    disp(['Profile input shorter than depth selected, truncating to ' num2str(zmax)])
end
clear zmax
z = 0:dz:depth;
nz = length(z);
% Check the depth-resolution of the profile file
profile_increment = (profile.z(end)-profile.z(1))/(length(profile.z)-1);
if dz < profile_increment/5
    yorn = input('Depth increment, dz, is much smaller than profile resolution. Is this okay? (y/n)','s');
    if yorn == 'n'
        error(['Please restart PWP.m with a new dz >= ' num2str(profile_increment/5)])
    end
end
t = zeros(nz,nmet);
s = zeros(nz,nmet);
t(:,1)	= interp1(profile.z,profile.t,z.','linear','extrap');
s(:,1)	= interp1(profile.z,profile.s,z.','linear','extrap');
d	= sw_dens0(s(:,1),t(:,1));
% Interpolate evaporation minus precipitaion at dt resolution
evap = (0.03456/(86400*1000))*interp1(met.time,met.qlat,floor(time),'nearest');
emp	= evap - precip;
emp(isnan(emp)) = 0;
%--------------------------------------------------------------------------
% preallocate space
u = zeros(nz,nmet); % east velocity m/s
v = zeros(nz,nmet); % north velocity m/s
% set initial conditions in PWP 18 Mar 2013
%u(1:13,1) = -5.514e-2;v(1:13,1) = -1.602e-1; 
mld = zeros(1,nmet); % mized layer depth m
absrb = absorb(beta1,beta2,nz,dz); % absorbtion of incoming rad. (units?)
dstab = dt*rkz/dz^2; % Courant number
if dstab > 0.5
    disp('!unstable CFL condition for diffusion!')
    pause
end
% output space
pwp_output.dt = dt;
pwp_output.dz = dz;
pwp_output.lat = lat;
pwp_output.z = z;
pwp_output.time	= zeros(nz,floor(nmet/dt_save));
pwp_output.t = zeros(nz,floor(nmet/dt_save));
pwp_output.s = zeros(nz,floor(nmet/dt_save));
pwp_output.d = zeros(nz,floor(nmet/dt_save));
pwp_output.u = zeros(nz,floor(nmet/dt_save));
pwp_output.v = zeros(nz,floor(nmet/dt_save));
pwp_output.mld	= zeros(nz,floor(nmet/dt_save));

%--------------------------------------------------------------------------
% model loop
for n = 2:nmet
     disp(['loop iter. ' sprintf('%2d',n)])
    % pwpgo function does the "math" for fluxes and vel
    [s(:,n), t(:,n), u(:,n), v(:,n), mld(n)] = pwpgo(qi(n-1),qo(n-1),emp(n-1),tx(n-1),ty(n-1), ...
        dt,dz,g,cpw,rb,rg,nz,z,t(:,n-1),s(:,n-1),d,u(:,n-1),v(:,n-1),absrb,f,ucon,n);
    
    % vertical (diapycnal) diffusion
    if rkz > 0
        diffus(dstab,t);
        diffus(dstab,s);
        d = sw_dens0(s,t);
        diffus(dstab,u);
        diffus(dstab,v);
    end % diffusion

    % Diagnostic plots
    switch diagnostics
        case 1
            figure(1)
            subplot(211)
            plot(time(n)-time(1),trapz(z,.5.*d.*(u(:,n).^2+v(:,n).^2)),'b.')
            if n == 2
                set(gcf,'position',[1 400 700 400])
                hold on
                grid on
                title('depth int. KE')
            end
            subplot(212)
            plot(time(n)-time(1),trapz(z,d.*sqrt(u(:,n).^2+v(:,n).^2)),'b.')
            if n == 2
                hold on
                grid on
                title('depth int. momentum')
            end
            figure(2)
            if n == 2
                set(gcf,'position',[700 400 700 400])
            end
            
            subplot(121)
            plot(u(:,n),z,'b',v(:,n),z,'r')
            axis ij
            grid on
            xlabel('u (b) v (r)')
            subplot(143)
            plot(t(:,n),z,'b')
            axis ij
            grid on
            xlabel('temp')
            subplot(144)
            plot(s(:,n),z,'b')
            axis ij
            grid on
            xlabel('salinity')
            
            pause(.2)
        case 2
            mom = zeros(1,nmet);
            if n == nmet
                for k = 1:nmet
                    mom(k) = trapz(z,d.*sqrt(u(:,k).^2+v(:,k).^2));
                    mli = find(sqrt(u(:,k).^2+v(:,k).^2) < 10^-3,1,'first');
                    mld(k) = z(mli);
                end
                w = u(1,:)+1i*v(1,:);
                tau = tx+1i*ty;
                T = tau./mld;
                dTdt = gradient(T,dt);
                pi_w = d(1).*mld.*real((w./1i*f).*dTdt);
                
                figure(1)
                plot(time-time(1),mom,'b.-')

            end
        case 3
            if n == nmet
                plot(time-time(1),mld,'b.')
                axis ij
                pause
            end
    end
    
end % model loop
% save ouput
pwp_output.s = s(:,1:dt_save:end);
pwp_output.t = t(:,1:dt_save:end);
pwp_output.u = u(:,1:dt_save:end);
pwp_output.v = v(:,1:dt_save:end);
pwp_output.d = sw_dens0(pwp_output.s,pwp_output.t);
time = repmat(time,nz,1);
pwp_output.time = time(:,1:dt_save:end);
mld(1)=mld(2);
pwp_output.mld = mld(:,1:dt_save:end);

save(pwp_output_file,'pwp_output')
toc
end % PWP driver routine
%--------------------------------------------------------------------------

function [s t u v mld] = pwpgo(qi,qo,emp,tx,ty,dt,dz,g,cpw,rb,rg,nz,z,t,s, ...
    d,u,v,absrb,f,ucon,n)
    % pwpgo is the part of the model where all the dynamics "happen"
  
t_old = t(1); s_old = s(1); 
 t(1) = t(1) + (qi*absrb(1)-qo)*dt./(dz*d(1)*cpw); 
s(1) = s(1)/(1-emp*dt/dz); 

Tf = sw_fp(s_old,1);
if t(1) < Tf;
    t(1) = Tf;
end

%  Absorb solar radiation at depth.
if size(dz*d(2:nz)*cpw,2) > 1
    t(2:nz) = t(2:nz)+(qi*absrb(2:nz)*dt)./(dz*d(2:nz)*cpw)'; %THIS IS THE ORIG ONE!
else
    t(2:nz) = t(2:nz)+(qi*absrb(2:nz)*dt)./(dz*d(2:nz)*cpw);
end

d = sw_dens0(s,t);
[t s d u v] = remove_si(t,s,d,u,v); %relieve static instability

% original ml_index criteria
ml_index = find(diff(d)>1E-4,1,'first'); %1E
%ml_index = find(diff(d)>1E-3,1,'first');
%ml_index = find( (d-d(1)) > 1e-4 ,1,'first');

% debug MLD index
%{
figure(86)
clf
plot(d,-z,'b.-')
hold on
plot(d(ml_index),-z(ml_index),'r*')
pause
%}
if isempty(ml_index)
	error_text = 'Mixed layer depth is undefined';
	error(error_text)
end

%  Get the depth of the surfacd mixed-layer.
ml_depth = z(ml_index+1);
mld = ml_depth; % added 2013 03 07 

% rotate u,v do wind input, rotate again, apply mixing
ang = -f*dt/2; % should be moved to main driver
[u v] = rot(u,v,ang);
du = (tx/(ml_depth*d(1)))*dt;
dv = (ty/(ml_depth*d(1)))*dt;
u(1:ml_index) = u(1:ml_index)+du;
v(1:ml_index) = v(1:ml_index)+dv;

% Apply drag to the current (this is a horrible parameterization of
% inertial-internal wave dispersion).
if ucon > 1E-10
	u = u*(1-dt*ucon);
	v = v*(1-dt*ucon);
end
[u v] = rot(u,v,ang);

% Bulk Richardson number instability form of mixing (as in PWP).
if rb > 1E-5
	[t s d u v] = bulk_mix(t,s,d,u,v,g,rb,nz,z,ml_index);
end

% Do the gradient Richardson number instability form of mixing.
if rg > 0
	[t,s,~,u,v] = grad_mix(t,s,d,u,v,dz,g,rg,nz);
end

% Debugging plots
%{
figure(1)
subplot(211)
plot(1,1)
hold on
plot(u,z,'b')
plot(v,z,'r')
axis ij
grid on
xlim([-.2 .2])

subplot(212)
plot(1,1)
hold on
plot(s-30,z,'b')
plot(t,z,'r')
grid on
axis ij
xlim([2 6])
pause(.05)
clf
%}
end % pwpgo

%--------------------------------------------------------------------------
function absrb = absorb(beta1,beta2,nz,dz)
%  Compute solar radiation absorption profile. This
%  subroutine assumes two wavelengths, and a double
%  exponential depth dependence for absorption.
%
%  Subscript 1 is for red, non-penetrating light, and
%  2 is for blue, penetrating light. rs1 is the fraction
%  assumed to be red.

rs1 = 0.6;
rs2 = 1.0-rs1;
%absrb = zeros(nz,1);
z1 = (0:nz-1)*dz;
z2 = z1 + dz;
z1b1 = z1/beta1;
z2b1 = z2/beta1;
z1b2 = z1/beta2;
z2b2 = z2/beta2;
absrb = (rs1*(exp(-z1b1)-exp(-z2b1))+rs2*(exp(-z1b2)-exp(-z2b2)))';
end % absorb

%--------------------------------------------------------------------------

function [t s d u v] = remove_si(t,s,d,u,v)
% Find and relieve static instability that may occur in the
% density array d. This simulates free convection.
% ml_index is the index of the depth of the surface mixed layer after adjustment,

while 1
	ml_index = find(diff(d)<0,1,'first');
	if isempty(ml_index)
		break
    end
    %%{
    figure(86)
    clf
    plot(d,'b-')
    hold on
	[t s d u v] = mix5(t,s,d,u,v,ml_index+1);
    plot(d,'r-')
    %pause
    %}
	%[t s d u v] = mix5(t,s,d,u,v,ml_index+1);
    
end

end % remove_si
%--------------------------------------------------------------------------

function [t s d u v] = mix5(t,s,d,u,v,j)
%  This subroutine mixes the arrays t, s, u, v down to level j.
t(1:j) = mean(t(1:j));
s(1:j) = mean(s(1:j));
d(1:j) = sw_dens0(s(1:j),t(1:j));
u(1:j) = mean(u(1:j));
v(1:j) = mean(v(1:j));
end % mix5

%--------------------------------------------------------------------------

function [u v] = rot(u,v,ang)
%  This subroutine rotates the vector (u,v) through an angle, ang
r = (u+1i*v)*exp(1i*ang);
u = real(r);
v = imag(r);

end %rot
%--------------------------------------------------------------------------

function [t s d u v] = bulk_mix(t,s,d,u,v,g,rb,nz,z,ml_index)

rvc = rb;
for j = ml_index+1:nz
	h 	= z(j);
	dd 	= (d(j)-d(1))/d(1);
	dv 	= (u(j)-u(1))^2+(v(j)-v(1))^2;
	if dv == 0
		rv = Inf;
	else
		rv = g*h*dd/dv;
	end
	if rv > rvc
		break
	else
		[t s d u v] = mix5(t,s,d,u,v,j);
	end
end

end % bulk_mix
%--------------------------------------------------------------------------

function [t s d u v] = grad_mix(t,s,d,u,v,dz,g,rg,nz)

%  This function performs the gradeint Richardson Number relaxation
%  by mixing adjacent cells just enough to bring them to a new
%  Richardson Number.

rc 	= rg;

%  Compute the gradeint Richardson Number, taking care to avoid dividing by
%  zero in the mixed layer.  The numerical values of the minimum allowable
%  density and velocity differences are entirely arbitrary, and should not
%  effect the calculations (except that on some occasions they evidnetly have!)

j1 = 1;
j2 = nz-1;

while 1
    r = zeros(size(j1:j2));
	for j = j1:j2
		if j <= 0
			keyboard
		end
		dd = (d(j+1)-d(j))/d(j);
		dv = (u(j+1)-u(j))^2+(v(j+1)-v(j))^2;
		if dv == 0
			r(j) = Inf;
		else
			r(j) = g*dz*dd/dv;
		end
	end

	%  Find the smallest value of r in profile
	[rs js] = min(r);
	%js = find(r==rs,1,'first');

	%  Check to see whether the smallest r is critical or not.
	if rs > rc
		return
	end

	%  Mix the cells js and js+1 that had the smallest Richardson Number
	[t s d u v] = stir(t,s,d,u,v,rc,rs,js);

	%  Recompute the Richardson Number over the part of the profile that has changed
	j1 = js-2;
	if j1 < 1
		 j1 = 1;
	end
	j2 = js+2;
	if j2 > nz-1
		 j2 = nz-1;
	end
end

end % grad_mix

%--------------------------------------------------------------------------
function [t s d u v] = stir(t,s,d,u,v,rc,r,j)

%  This subroutine mixes cells j and j+1 just enough so that
%  the Richardson number after the mixing is brought up to
%  the value rnew. In order to have this mixing process
%  converge, rnew must exceed the critical value of the
%  richardson number where mixing is presumed to start. If
%  r critical = rc = 0.25 (the nominal value), and r = 0.20, then
%  rnew = 0.3 would be reasonable. If r were smaller, then a
%  larger value of rnew - rc is used to hasten convergence.

%  This subroutine was modified by JFP in Sep 93 to allow for an
%  aribtrary rc and to achieve faster convergence.

rcon 			= 0.02+(rc-r)/2;
rnew 			= rc+rcon/5;
f 				= 1-r/rnew;
dt				= (t(j+1)-t(j))*f/2;
t(j+1)		= t(j+1)-dt;
t(j) 			= t(j)+dt;
ds				= (s(j+1)-s(j))*f/2;
s(j+1)		= s(j+1)-ds;
s(j) 			= s(j)+ds;
d(j:j+1)	= sw_dens0(s(j:j+1),t(j:j+1));
du				= (u(j+1)-u(j))*f/2;
u(j+1)		= u(j+1)-du;
u(j) 			= u(j)+du;
dv				= (v(j+1)-v(j))*f/2;
v(j+1)		= v(j+1)-dv;
v(j) 			= v(j)+dv;

end % stir
