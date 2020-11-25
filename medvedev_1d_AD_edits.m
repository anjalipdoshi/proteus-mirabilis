%% trying to store the values of u & v @ every time pt
% (i.e. arrays that are x by t)
% 11-14-20

%% Implementing Medvedev Model of Mirabilis Patterning %%
% https://epubs.siam.org/doi/abs/10.1137/S0036139998344635
% This model does not use age component like Esipov Shapiro but instead
% only uses reaction-diffusion eqns.
% Authors say they used "explicit Euler scheme" for discretization-- using
% a 1st order forward diff method for du/dt and 2nd order centered diff for
% laplacian operators
% u =  swarm cell surface density (# swarm cells per unit area)
% v = swimmer cell surface density (# swimmer cells per unit area)
clear all;
cla reset;

%% Define Mesh
N = 50; %100; %N must be even!
h = 10/N; %10/N; % step size in x and y
x = h*(0:N); % x coordinates of grid

% for time-stepping loop below
dt = 0.03*h*h; % time step - usually small
t = 0; 
tmax = 10; %1; %10; 
nsteps = round(tmax/dt); % number of time steps
time_points = 1*(0:nsteps);

%space-time mesh:
[xx,tt] = meshgrid(x,time_points);

%% Define Parameters

% Start with values from fig 1 caption
k = 5; % in diffusion equation
mu = 0.6; % cell division
nu = 0.6; % differentiation to sw cells-> only occurs when concentration of swimmers v is between vflat and vtilda
alpha = 0.7; % cell growth
sigma = 1.5; %1.5;% very small, not sure what value?
D0 = 1; %5; %.05; %1; % diffusion coefficient
vflat = 12; %12;
vtilda = 16; %16; % saturating concentration above which colony development stops

%% Define ICs

% initially, there should be no swarm & swimmer cells outside inoculum
u = zeros(size(xx));
v = zeros(size(xx));

% add a small perturbation Inecessary for pattern to form?)
%u = u + .01; 
%v = u + 1;

% initial inoculum (should just be swimmers, but must add small amount of swarmers)
radius = 5; %(looks like it was 10 in paper fig 1?), %5; %initial colony radius
v(1,[N/2-radius+1:N/2+radius+1]) = 1; %5; %1; 
u(1,[N/2-radius+1:N/2+radius+1]) = 0; %.05;

%load('matlab.mat')
%load('matlab2.mat')
%pause; 

D = zeros(size(xx)); % D is a function of u (initially 0 everywhere)--> change to D=u perturbation?
%D = u;
%D = D0*(u./(u+k*v));
%D = D0.*ones(size(xx));

%nu(v), mu(v), and alpha(v) are functions of v (swimmer cell density)
%nu = differentiation to swarmer cells
%mu = cellular division
%alpha = cellular growth
nuv = zeros(size(xx)); % nuv is initially 0 wherever v is not in range (vflat,vtilda)
muv = zeros(size(xx)); % when v>vtilda, muv = 0
alphav= zeros(size(xx)); % when v>vtilda, alpha = 0

%% Main time stepping loop

t_vec = t:dt:tmax; % a vector from t=0 to tmax with steps of dt

for n = 1:nsteps 
%      if mod(n,10)==0
%          n
%      end
%     
    t = t+dt;
    
    % Save previous loop's values
    uo = u(n,:); % for nth time-point, initial u is u(nth row, all cols)
    vo = v(n,:);
    
  %  pause; 
   
    % Will calculate D & nuv based on previous loop's values before updating
    
    %Calculate Laplacian
    uE = u(n,[2:N+1 N]);%uE= u(:,[2:N+1 1]); %I don't really get how we write uE, uW,N,S - the commented formulas make sense for me
    uW = u(n,[2 1:N]);%uW = u(:,[N+1 1:N]); 
    Lu = uE+uW-2*u; %approximation
    
    % Calculate diffusion coefficient from previous step's values; Eq. 1.2
    D(uo==0 & vo==0)=0;
    Dtemp=D0*(uo./(uo+k*vo));
    D(~(uo==0 & vo==0))=Dtemp(~(uo==0 & vo==0));
    
    %D(n,:) = D0*(uo./(uo+k*vo));
    %    D(n,:) = D0*(uo./(uo+k*vo));
   
    
    D_E = D(n,[2:N+1 N]);
    D_W = D(n,[2 1:N]);
  
    % Determine where differentiation is 0 or nu
     locs1 = ((vflat+sigma) < v(n,:)) & (v(n,:) < (vtilda-sigma)); 
     locs2_a = (v(n,:) > vtilda);
     locs2_b = (v(n,:) < vflat);
     locs3_a = ((v(n,:) > vflat) & (v(n,:) < (vflat+sigma)));
     locs3_b = ((v(n,:) > (vtilda-sigma)) & (v(n,:) < vtilda));
      
     nuv(n,locs1) = nu; % differentiation to swarmers where swimmers are in correct concentration range
     nuv(n,locs2_a) = 0; % 0 differentiation where swim cells concentration are not in range
     nuv(n,locs2_b) = 0; % 0 differentiation where swim cells concentration are not in range
     nuv(n,locs3_a) = (v(n,locs3_a)-vflat).*(nu/sigma); %"continous function of v" 
     %nuv(locs3_a) = exp(nu); %TRYING OTHER CONTINUOUS FUNCTIONS
     %nuv(locs3_a) = v(locs3_a);
     nuv(n,locs3_b) = nu - (v(n,locs3_b)-(vtilda - sigma)).*(nu/sigma); %"continous function of v" 
     %nuv(locs3_b) = exp(-nu); %TRYING OTHER CONTINUOUS FUNCTIONS
     %nuv(locs3_b) = v(locs3_b);
     
     locs4 = (v(n,:) > vtilda);
     locs5 = (v(n,:) < (vtilda - sigma));
     locs6 = (v(n,:) > (vtilda-sigma) & (v(n,:) < vtilda));
     % don't know how alpha and mu look like when v is between vtilda-sigma
     % and vtilda; assume linear relationship
     alphav(n,locs4) = 0; % no growth when v>vtilda
     alphav(n,locs5) = alpha;
     alphav(n,locs6) = alpha - (v(n,locs6)-(vtilda-sigma)).*(alpha/sigma);  %"continous function of v" btwn locs4 & locs5
     %alphav(locs6) = exp(alpha); %TRYING OTHER CONTINUOUS FUNCTIONS
     %alphav(locs6) = v(locs6);
     
     muv(n,locs4) = 0; % no division when v>vtilda
     muv(n,locs5) = mu;
     muv(locs6) = mu - (v(n,locs6) -(vtilda-sigma)).*(mu/sigma); %"continous function of v" btwn locs4 & locs5
     %muv(locs6) = exp(mu); %TRYING OTHER CONTINUOUS FUNCTIONS
     %muv(n,locs6) = v(n,locs6);

    % calculate "centered" D's
    D_c_E=(D_E+D(n,:))/2;
    D_c_W=(D_W+D(n,:))/2;
    
    D_term = (D_c_E.*(uE-uo)-D_c_W.*(uo-uW))/h^2;
    
    u(n+1,:) = uo + dt*(vo.*nuv(n,:) + uo.*(alphav(n,:)-muv(n,:))+ D_term);
    
    v(n+1,:) = vo + dt*(vo.*(alphav(n,:)-nuv(n,:)) + uo.*muv(n,:));
   

    if (mod(n,200)==0)
        colormap('jet')
        caxis([0 20]);
        figure(1);
        subplot(1, 2, 1);
        imagesc(u);% displays the data in array u as an image that uses the full range of colors in the colormap. 
           %Each element of u specifies the color for one pixel of the image
        title(['swarm cell density (u)'],'fontsize',16);
        xlabel('position, x');
        ylabel('time step, dt');
        caxis([0 20]);
        colorbar;

        subplot(1, 2, 2);
        imagesc(v);% displays the data in array u as an image that uses the full range of colors in the colormap. 
           %Each element of u specifies the color for one pixel of the image
        title(['swim cell density (v)'],'fontsize',16);
        xlabel('position, x');
        ylabel('time step, dt');
        caxis([0 20]);
        colorbar;
        
        sgtitle(sprintf('time step %d', n));

    end
    
end
%%
%determining colony radius over time
%currently these thresholds only work the parameters set in this file
for i = 1:(nsteps+1)
    uzeros = find((u(i,[1:(N/2+1)]))<=0.012); % get the locations where u is 0.01 (~0) on just the left half of the plate up to the inoculum
    urad = N/2-max(uzeros); % radius of swarmer cells
    vzeros = find((v(i,[1:(N/2+1)]))<.0091); % get the locations where v is 0 on just the left half of the plate up to the inoculum
    vrad = N/2-max(vzeros); % radius of swimmer cells
    colrad(i) = max(urad,vrad); % colony radius is max radius of swarm vs swim cells
    %colrad(i)=urad; %if we're only considered swarmer cells
end
    
    colormap('jet')
    figure(6);
    imagesc(u);% displays the data in array u as an image that uses the full range of colors in the colormap. 
       %Each element of u specifies the color for one pixel of the image
    title(['swarm cell density (u) over space and time'],'fontsize',16);
    xlabel('position, x');
    ylabel('time step, dt');
    colorbar;

figure(2);
imagesc(v);% displays the data in array u as an image that uses the full range of colors in the colormap. 
       %Each element of u specifies the color for one pixel of the image
title(['swimmer cell density (v) over space and time'],'fontsize',16);
xlabel('position, x');
ylabel('time step, dt');
colorbar;
      
figure(3);
plot(x,u(:,:));
xlabel('position, x');
ylabel('Swarm cell surface density, u');
title('Spatial swarm cell development');

figure(4);
plot(x,v(:,:));
xlabel('position, x');
ylabel('Swimmer cell surface density, v');
title('Spatial swimmer cell development');

figure(5);
plot(x,u(length(time_points),:));
hold on;
plot(x,v(length(time_points),:));
lgd = legend('u = swarm cells', 'v = swimmer cells');
xlabel('position, x');
ylabel('cell surface density');
title('Final cell surface densities');

figure(6);
plot(t_vec,colrad);
xlabel('time, t');
ylabel('Colony radius');
title('Colony radius vs. time');

%fourier = fft(u);
%https://www.mathworks.com/help/matlab/ref/fft.html

