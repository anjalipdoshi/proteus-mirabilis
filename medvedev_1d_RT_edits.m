%% trying to store the values of u & v @ every time pt
% (i.e. arrays that are x by t)


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
h = 0.2; %10/N; % step size in x and y
x = h*(0:N); % x coordinates of grid

% for time-stepping loop below
dt = 0.03*h*h; % time step - usually small
t = 0; 
tmax = 15; %1; %10; 
nsteps = round(tmax/dt); % number of time stepss

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
u = zeros(nsteps,length(x));
v = zeros(nsteps,length(x));


% initial inoculum (should just be swimmers, but must add small amount of swarmers)
radius = 3; %initial colony radius
v(1,[N/2-radius+1:N/2+radius+1]) = 1; %5; %1; 
u(1,[N/2-radius+1:N/2+radius+1]) = 0; %.05;

D = zeros(nsteps,length(x)); % D is a function of u (initially 0 everywhere)


%nu(v), mu(v), and alpha(v) are functions of v (swimmer cell density)
%nu = differentiation to swarmer cells
%mu = cellular division
%alpha = cellular growth
nuv = zeros(nsteps,length(x)); % nuv is initially 0 wherever v is not in range (vflat,vtilda)
muv = zeros(nsteps,length(x)); % when v>vtilda, muv = 0
alphav= zeros(nsteps,length(x)); % when v>vtilda, alpha = 0

%create video for u and v in space and time
video_1 = VideoWriter('v_u_movie.avi');
video_2 = VideoWriter('D_movie.avi');
open(video_1);
open(video_2);

%% Main time stepping loop

t_vec = t:dt:tmax; % a vector from t=0 to tmax with steps of dt

for n = 1:nsteps    
    t = t+dt;
    
    % Save previous loop's values
    uo = u(n,:); % for nth time-point, initial u is u(nth row, all cols)
    vo = v(n,:);
    
  %  pause; 
   
    % Will calculate D & nuv based on previous loop's values before updating
    
    %Calculate Laplacian
    uE = u(n,[2:N+1 N]);
    uW = u(n,[2 1:N]);
    Lu = uE+uW-2*u; 
    
    % Calculate diffusion coefficient from previous step's values; Eq. 1.2
    %D(n,:)=D0+zeros(1,length(x));%linear D
    
    D(n,uo==0 & vo==0)=0;
    Dtemp=D0*(uo./(uo+k*vo));
    D(n,~(uo==0 & vo==0))=Dtemp(~(uo==0 & vo==0));
    
   % D(n,:) = D0*(uo./(uo+k*vo));

    D_E = D(n,[2:N+1 N]);
    D_W = D(n,[2 1:N]);
  
    %Sigmoid functions for alpha, mu and nu
    alphav(n,:) = (alpha*(1+zeros(size(x))))./(1+exp(7*(v(n,:)-15.25)));
      
    muv(n,:) = (mu*(1+zeros(size(x))))./(1+exp(7*(v(n,:)-15.25)));
      
    nuv(n,:) = (nu*(1+zeros(size(x))))./((1+exp((-7)*(v(n,:)-12.75))).*(1+exp(7*(v(n,:)-15.25))));
      


    % calculate "centered" D's
    D_c_E=(D_E+D(n,:))/2;
    D_c_W=(D_W+D(n,:))/2;
    
    D_term = (D_c_E.*(uE-uo)-D_c_W.*(uo-uW))/h^2;
    
    u(n+1,:) = uo + dt*(vo.*nuv(n,:) + uo.*(alphav(n,:)-muv(n,:))+ D_term);
    
    v(n+1,:) = vo + dt*(vo.*(alphav(n,:)-nuv(n,:)) + uo.*muv(n,:));
   
    %Make u and v on the edges equal to 0
    u(n+1,1)=0;  u(n+1,length(x))=0;
    v(n+1,1)=0;  v(n+1,length(x))=0;

  
    
    %plot u v and D 
    if (mod(n,200)==0)
        colormap('jet')
        caxis([0 20]);
        figure(1);
        subplot(1, 2, 1);
        imagesc(u);% displays the data in array u as an image that uses the full range of colors in the colormap. 
           %Each element of u specifies the color for one pixel of the image
        title(['swarm cell density (u)'],'fontsize',10);
        xlabel('position, x');
        ylabel('time step, dt');
        caxis([0 20]);
        colorbar;

        subplot(1, 2, 2);
        imagesc(v);% displays the data in array u as an image that uses the full range of colors in the colormap. 
           %Each element of u specifies the color for one pixel of the image
        title(['swim cell density (v)'],'fontsize',10);
        xlabel('position, x');
        ylabel('time step, dt');
        caxis([0 20]);
        colorbar;
        
        sgtitle(sprintf('time step %d', n));
        
        frame_1=getframe(gcf);
        writeVideo(video_1,frame_1);
        
        figure(2);
        imagesc(D);
        title(['Diffusivity (D) over space and time'],'fontsize',16);
        xlabel('position, x');
        ylabel('time step, dt');
        colorbar;
        
        frame_2=getframe(gcf);
        writeVideo(video_2,frame_2);
        
        
    end
    
end

close(video_1);
close(video_2);

%%
%determining colony radius over time
%currently these thresholds only work the parameters set in this file
% colrad = zeros(1,nsteps+1);
% for i = 1:(nsteps+1)
%     uzeros = find((u(i,1:(N/2+1)))<=0.012); % get the locations where u is 0.01 (~0) on just the left half of the plate up to the inoculum
%     urad = N/2+1-max(uzeros); % radius of swarmer cells
%     vzeros = find((v(i,1:(N/2+1)))<=.0091); % get the locations where v is 0 on just the left half of the plate up to the inoculum
%     vrad = N/2+1-max(vzeros); % radius of swimmer cells
%     colrad(i) = max(urad,vrad); % colony radius is max radius of swarm vs swim cells
%     %colrad(i)=urad; %if we're only considered swarmer cells
% end
    


figure(3);
imagesc(v);% displays the data in array u as an image that uses the full range of colors in the colormap. 
       %Each element of u specifies the color for one pixel of the image
title(['swimmer cell density (v) over space and time'],'fontsize',16);
xlabel('position, x');
ylabel('time step, dt');
colorbar;

figure(4);
imagesc(u);
title(['swarmer cell density (u) over space and time'],'fontsize',16);
xlabel('position, x');
ylabel('time step, dt');
colorbar;



figure(5);
plot(x,u(:,:));
xlabel('position, x');
ylabel('Swarm cell surface density, u');
title('Spatial swarm cell development');

figure(6);
plot(x,v(:,:));
xlabel('position, x');
ylabel('Swimmer cell surface density, v');
title('Spatial swimmer cell development');

figure(7);
plot(x,u(nsteps,:));
hold on;
plot(x,v(nsteps,:));
lgd = legend('u = swarm cells', 'v = swimmer cells');
xlabel('position, x');
ylabel('cell surface density');
title('Final cell surface densities');


% figure(8);
% plot(t_vec,colrad);
% xlabel('time, t');
% ylabel('Colony radius');
% title('Colony radius vs. time');

%fourier = fft(u);
%https://www.mathworks.com/help/matlab/ref/fft.html

