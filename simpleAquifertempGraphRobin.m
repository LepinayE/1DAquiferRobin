% simpleAquifertempGraphRobin - This script analyses the temperature profiles of 
% the aquifer from the solveAquifertemp.m file
%
%--------------------------------------------------------------------------
% Syntax of solveAquifertemp.m :
% [results] = solveAquifertempRobin(K_d,ncyc,years,K_r, U0, b)

%
% Inputs:
%   K_d - Effective thermal dispersivity of fracture ( 10^-5 to 10^-7)
%   ncyc - Total number of cycles 
%   years - Number of years system is run
%   K_r - Molecular heat diffusivity of rock (10^-7)
%   U0 - Fluid speed in fracture (10^-5) 
%   Tinj - Temperature of injected fluid
%   Taq - Initial aquifer temperature
%   b - Height of fracture (model is symmetrical and takes b/2)
%    
%
% Outputs:
%    results - structure with NON-DIM fields:
%                       .frac - T_f(X,t), fracture temperature
%                       .rock = T_r(X,Y,t), rock temperature
%                       .velocity = U, flow velocity
%                       .dt - time step
%                       .t_vec - vectorised time s.t t_vec(k) = k*dt
%                       .X_vec - vectorised horizontal length space
%                       .dX - horizontal step
%                       .Y_vec - vectorised vertical length space
%                       .dY - vertical step
%
% Other m-files required: solveAquifertemp.m 
% Subfunctions: flowVelocity.m
% See also: analyseAquifertemp.m
%
%--------------------------------------------------------------------------
% Analysis includes:
%                   - Temperature profile along fracture for injections and
%                   extractions: figure (1) & (2)
%                   - Contour of rock temperature just above fracture
%                   through out the cycles: figure(3)
%                   - Temperature profile in rock and fracture at 2 pts
%                   down fracture through the cycles: figure(4) & (5)
%                   - Temperature contour of rock after last injection:
%                   figure (6)
%                   - Extraction Temperature at Well during each cycle:
%                   figure (7)
%                       
%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 09/09/2022; Last revision: 20/09/2022
% Version: R2022a

addpath '/Users/lepinay/Desktop/Aquifer Matlab'/FluxbcStepAquifer
%------------- BEGIN CODE -------------------------------------------------
% Aquifer Conditions

K_d = [10^(-5), 10^(-6), 10^(-7)];

K_r = 10^(-7); % Molecular

ncyc = [1,10,50,100];

years = [1, 10]; 

U0 =  10^(-5); % m/s^2

Q =  10^(-5); % Make Qhat 1 with Tinj = 1, Taq = 0

Tinj = 1; 

Taq = 0; 

b = 2; % meters

%-------------
%Parameter wanted

iyears = 1;

incyc = 2 ;

iKd = 1;

% Call function
results = solveAquifertempRobin(K_d(iKd),ncyc(incyc),years(iyears),K_r, U0, b);

%% Temperature along Fracture and Rock during Injection and Extraction

% Find Index of Injections and Extractions
idinjEnd = islocalmax(results.frac(1,:)); % End of injections
idextEnd = islocalmin(results.frac(1,:));% End of extractions
idextEnd(end) = 1;

idfinishes = idextEnd|idinjEnd; % Logical of ends of injections and ...
%                                   extractions

injEnd = find(idinjEnd); % Time index of end of injections
extEnd = find(idextEnd); % Time index of end of extraction
 
finishes = find(idfinishes);% Time index of end of injections and...
%                               extractions

nlines = 20; % Number of lines wanted on graph

indBigtimept = 1 : floor((length(results.t_vec)) /(2*nlines)):...
           (length(results.t_vec)); % Time step index

indBigtimeptyears = indBigtimept * results.dt ;
idinjVelocity = results.velocity(indBigtimept) >=0 ; %Index of time step 
%                                                       during injections


%-------------
% Fracture Temperature during Injection
figure(1)

TfInj = results.frac(:,indBigtimept(idinjVelocity));

plot(results.X_vec, TfInj ,'linewidth',2)
drawnow
legend(sprintfc('Time = %.2f yr', indBigtimeptyears(idinjVelocity), false))
xlabel('X')
xlim([0 1])
ylabel('$\hat{T}_f$' ,'Interpreter','latex')

%-------------
% Fracture Temperature during Extraction
figure(2)

TfExt = results.frac(:,indBigtimept(~idinjVelocity));

plot(results.X_vec(1:end), TfExt ,'linewidth',2)
drawnow
legend(sprintfc('Time = %.2f yr', indBigtimeptyears(~idinjVelocity), false))
xlabel('X')
ylabel('$\hat{T}_f$' ,'Interpreter','latex')

ylim([0 1])

%% Contour of rock temperature just above fracture through out cycles

figure(3)


fAty = squeeze(results.rock(:,1,:));

contour(results.t_vec,results.X_vec,fAty)
xlabel('Time (yr)')
ylabel('X-direction')
colorbar
xlim([0 0.5])
%{

% Zoom in on graph
y1 = 0; % Start of y axis label 
y2 = 0.6 ; % End of y axis label
ylim([y1+results.dX ,y2]) 

% Changing ticks location
ticksYaxis = y1:y2 /6:y2 ;
ticksYaxis(1) = results.dX;
yticks(ticksYaxis)

%------------- 
% Run after running this section
% Changing labels of ticks
yl = yticklabels;
yl(1) = {'0'};
yticklabels(yl);
%}

%% Temperature profile of rock and fracture through the cycles
% at 2 different pts down the fracture (vary K_eff)

% X coordinates
indx1 = 1;
indx2 = 15;

% Legend
xdim = ((b/2)^(2))* (U0)/(K_r);

x1 = xdim*results.X_vec(indx1+1);
x1 = round (x1,1);
x1 = sprintf("%.1f",x1);

x2 = xdim*results.X_vec(indx2+1);
x2 = round (x2,1);
x2 = sprintf("%.1f",x2);

txt = ['x = ',num2str(x1) , 'm'];
txt2 = ['x = ',num2str(x2) , 'm'];

%-------------
% Changes in Fracture temp through the cycles

figure(4)

plot(results.t_vec,results.frac(indx1,:),'LineWidth', 2, 'DisplayName', txt)
hold on 
plot(results.t_vec,results.frac(indx2,:),'LineWidth', 2, 'DisplayName', txt2)
ylabel('$\hat{T}_f$' ,'Interpreter','latex')
xlabel('Time (yr)')
legend ('Location', 'northeastoutside')
xlim([0, results.t_vec(end)])

%-------------
% Changes in Rock temp through the cycles

TrAty = results.rock(:,10,:);
TrAty = squeeze(TrAty);

Y = results.Y_vec(10)

% Plot
figure(5)

plot(results.t_vec,TrAty(indx1,:),'LineWidth', 2, 'DisplayName', txt)
hold on 
plot(results.t_vec,TrAty(indx2,:),'LineWidth', 2, 'DisplayName', txt2)
ylabel('$\hat{T}_r$' ,'Interpreter','latex')
xlabel('Time (yr)')
legend ('Location', 'northeastoutside')
xlim([0, results.t_vec(end)])


%% Contour plot of Rock temp at end of Last Injection

% Find Index of Injections
idinjEnd = islocalmax(results.frac(1,:)); % Logical, End of injections

injEnd = find(idinjEnd); % Time index of end of injections
lastinjEnd = injEnd(end); % Time index of end of last injection

% Rock Temperature at end of Last Injection
Trcontour = results.rock (:,:,lastinjEnd);
Trcontour = squeeze(Trcontour);
Trcontour = Trcontour .' ;
 
%-------------
figure(6)

contour(results.X_vec,results.Y_vec, Trcontour)
xlabel('X')
ylabel('Y')

ylim([0 2.1])
xlim ([0, 0.5])
colorbar


%% Extraction Temperature at Well during each cycle 

%Index of time step during extraction
indxTime = 1:length(results.t_vec); 
idextrVelocity = results.velocity(indxTime) <0 ; 


% Find index Start and End of Extractions
idextStart = islocalmax(results.frac(1,:)); % Logical, Start of extractions
idextEnd = islocalmin(results.frac(1,:));% Logical, End of extractions
idextEnd(end) = 1;

extStart = find(idextStart); % Time index of end of injections
extEnd = find(idextEnd); % Time index of end of extraction

%-------------
figure(7)

for n = 1:ncyc(incyc)
    txt3 = ['Cycle number = ',num2str(n)];

    extTf = [];
    
    extTf = results.frac(1, extStart(n): extEnd(n));
    
    plot(extTf, 'o', 'DisplayName', txt3);
    hold on

   
end

xlabel('Time during Extraction (days)')
ylabel('Extraction Temperarture at Well')
legend('Location','northeastoutside')


%-------------
%Reshaping Time axis

xt = xticks;
xt = xt * results.dt; % Time in years 
xt = xt* 12; %Time in months
xt = floor (xt *31); % Rough time in days 
xticklabels(xt); % Change ticks labels 

%% Change in Extraction Temperature at Well over each cycle 

%Index of time step during extraction
indxTime = 1:length(results.t_vec); 
idextrVelocity = results.velocity(indxTime) <0 ; 


% Find index Start and End of Extractions
idextStart = islocalmax(results.frac(1,:)); % Logical, Start of extractions
idextEnd = islocalmin(results.frac(1,:));% Logical, End of extractions
idextEnd(end) = 1;

extStart = find(idextStart); % Time index of end of injections
extEnd = find(idextEnd); % Time index of end of extraction

% Difference of temp vector
n = 1:ncyc(incyc);
DiffextTf(n) = results.frac(1, extStart(n)) - results.frac(1, extEnd(n));

%-------------
figure(8)
plot(DiffextTf, 'o');


xlabel('Cycle Number')
ylabel('Change in Temperarture')




