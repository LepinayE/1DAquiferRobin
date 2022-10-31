% HeatconsAquifertempRobin.m - Heat Conservation check & Heat ratio 
%
% During injection flux in due to advection and diffusion must equal total
% energy in Aquifer at each time step. During extraction, internal energy
% of Aquifer must equal advective heat flux out.
% Fracture is taken as the Aquifer.
% i = 1 grid point is not part of Aquifer.
%
%--------------------------------------------------------------------------
% Syntax of solveAquifertempRobin.m :
% [results] = solveAquifertempRobin(K_d,ncyc,years,K_r, U0, b)
%
% Inputs:
%   K_d - Effective thermal dispersivity of fracture ( 3*10^-5 to 3*10^-7)
%   ncyc - Total number of cycles 
%   years - Number of years system is run
%   K_r - Molecular heat diffusivity of rock (10^-7)
%   U0 - Fluid speed in fracture (10^-5) 
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
%
%--------------------------------------------------------------------------
% Analysis includes:
%                   - Figure(1): Thermal Flux in and out of Aquifer
%                   - Figure(2): Total Thermal Flux in and out of Aquifer
%                                   (J)
%                   - Figure(3): Internal Energy of Aquifer (J)
%                   - Figure(4): Energy changes overview of System (J)
%                   - Figure(5): Thermal Flux in and out of Aquifer
%                   - Figure(6): Heat ratio of Flux in to Internal Aquifer
%                               during each Injection cycle                               
%                   - Figure(7): Heat ratio of Flux in to Internal Aquifer
%                               at each time step 
%               
%--------------------------------------------------------------------------
%
% Other m-files required: solveAquifertempRobin.m
% Subfunctions: flowVelocity.m
% MAT-files required: none
%
% See also: analyseAquifertemp.m
%
%--------------------------------------------------------------------------
% Author: Emma Lepinay
% Email: el547@cam.ac.uk
% Date: 18/10/2022; Last revision: 31/10/2021
% Version: R2022a

addpath '/Users/lepinay/Desktop/Aquifer Matlab'/FluxbcStepAquifer

%------------- BEGIN CODE -------------------------------------------------

% Aquifer Conditions

K_d = [10^(-5), 10^(-6), 10^(-7)];

K_r = 10^(-7); % Molecular

ncyc = [1, 10,50,100];

years = [1, 10]; 

U0 =  10^(-5); % m/s^2

b = 2; % meters

Tinj = 10; % Injection temp

Taq = 2; % Aquifer temp (celsius)

%-------------
% Run solveAquifertemp.m

incyc = 1;
iyears = 1;
iKd = 3;

results = solveAquifertempRobin(K_d(iKd),ncyc(incyc),years(iyears),K_r,...
                                  U0, b);
%-------------
% Useful parameters
time = results.t_vec; % Time vector
timeStep = 1: length(time); % Time steps vector

%%
%-------------
% Dimensionalise variables
%Tf = (Tinj - Taq).*(results.frac) +Taq;
%Tr = (Tinj - Taq).*(results.rock) +Taq;
Tf = results.frac;
Tr = results.rock;
%-------------
% Thermal Flux into Aquifer (Flow of energy per unit area per unit time)

cp = 4000 ; % Heat capacity of water J/kgdegrees
row = 1000 ; % Density of water kg/m^3

bigQinj = zeros (1, timeStep(end));


% Build Thermal Flux vector
dX = results.dX;
idInjvelocity = results.velocity(timeStep) >= 0; % Logical 1 for injection
    
%bigQinj(idInjvelocity) = cp*row*Tinj*U0 .*results.velocity(idInjvelocity); % J/sm^2
bigQinj(idInjvelocity) = cp*row*U0 .*results.velocity(idInjvelocity);

dTfdx(idInjvelocity) = Tf(2,idInjvelocity)- Tf(1,idInjvelocity);
bigQinjtwo(idInjvelocity) = cp*row*U0.* results.velocity(idInjvelocity)...
                  .*(Tf(1,idInjvelocity)) + ...
                    cp*row* ((4 *K_d(iKd) * K_r)/((b^2) *U0 * dX))...
                    .*(dTfdx(idInjvelocity));

% Total flux into Aquifer
totalFluxin = cumsum(bigQinj);
totalFluxin = ((b/2)^2)* (1/K_r)*results.dt * totalFluxin; % J/m^2
totalFluxin = (b/2)* totalFluxin; %J

totalFluxin2 = cumsum(bigQinjtwo);
totalFluxin2 = ((b/2)^2)* (1/K_r)*results.dt * totalFluxin2; % J/m^2
totalFluxin2 = (b/2)* totalFluxin2;

%-------------
% Thermal Flux out of Aquifer (Flow of energy per unit area per unit time)

bigQout = zeros (1, timeStep(end));

bigQout(~idInjvelocity) = row* cp* U0 .* results.velocity(~idInjvelocity) ...
                            .* (Tf(1,~idInjvelocity) ) ;
bigQout = - bigQout; % J/sm^2 

% Total flux out of Aquifer
totalFluxout = cumsum(bigQout);
totalFluxout = ((b/2)^2)* (1/K_r)*results.dt * totalFluxout;
totalFluxout = (b/2)* totalFluxout;

%-------------
% Internal Energy of Aquifer

% Thermal Energy in fracture
fracEnergy = Tf;
fracEnergy(1,:) = 0;
fracEnergy =  sum (fracEnergy); % Vector of sum of temp at each time step
fracEnergy = row * cp * ((b/2)^3)* U0* (1/K_r) * results.dX * fracEnergy;


% Thermal Energy in rock
rockEnergy = sum(Tr, [1 2]);% Sum over the X and Y variables
rockEnergy = squeeze (rockEnergy);
rockEnergy = rockEnergy.' ;

rockEnergy = row * cp * ((b/2)^3)* U0* (1/K_r)*results.dX * results.dY* ...
            rockEnergy; %Vector of sum of temp at each time step

% Thermal Energy at t = 0
%TempAq = Taq * ones(length(results.X_vec) , length(results.Y_vec)); % T(x,y)
%originalIntEnergy = sum(TempAq, 'all');
%originalIntEnergy = results.dX * results.dY * originalIntEnergy;

% Internal Energy at each time k
intEnergy =  fracEnergy + rockEnergy; % + originalIntEnergy
 % J

%-------------
%%
% Heat Ratio, each time step

heatDiffratio = diff(totalFluxin)./ diff(intEnergy); % Energy injected to...
%                                       energy in aquifer at each time step

heatDiffratio (heatDiffratio <0) = NaN ; % Zeroes values below zero due to
%                                       discontinuous Internal Energy
heatDiffratio(end +1) = heatDiffratio(end); % Allows plot against time

%-------------
% Heat Ratio at time t

%Index for start and end of Injections
idInjEnd = islocalmax(intEnergy); % Logical 1 at index of end injection
indexInjEnd = find(idInjEnd);% Index of end of each injection 
indexInjfirstEnd = indexInjEnd(1);

idInjStart = islocalmin(intEnergy); % Logical 1 at start of injection
idInjStart(1) = 1 ;
indexInjStart = find(idInjStart); % Index of start of each injection 

% Thermal Flux in per cycle
Fluxinpercycle = cumsum(bigQinj);
for n = 1: ncyc(incyc)
   
    totalFluxinpercycle(indexInjStart(n):indexInjEnd(n)) = ...
                        Fluxinpercycle(indexInjStart(n):indexInjEnd(n)) ...
                        - Fluxinpercycle(indexInjStart(n));

end
totalFluxinpercycle = ((b/2)^2)* (1/K_r)*results.dt .* totalFluxinpercycle; 
totalFluxinpercycle = (b/2)* totalFluxinpercycle; %J

% Internal Energy  per cycle
for n = 1: ncyc(incyc)
intEnergypercyc(indexInjStart(n):indexInjEnd(n)) = ...
                        intEnergy((indexInjStart(n):indexInjEnd(n))) ...
                        - intEnergy(indexInjStart(n));
end

% Energy in to energy in aquifer at each time t
heatRatio = totalFluxinpercycle./ intEnergypercyc ; 


%% Energy Plots 
%
%-------------
% Thermal Flux in and out Aquifer 
figure(1)

plot(time, bigQinj)
hold on
plot(time, bigQout)
ylabel('Thermal Flux (J/sm^2)')
xlabel('Time (yrs)')
legend ('During Injection', 'During Extraction', 'Location',...
            'northeastoutside')
grid on

% Total Total Thermal Flux In and Out Aquifer
figure(2)

plot(time, totalFluxin)
hold on 
plot(time, totalFluxout)
ylabel('Total Thermal Flux (J)')
xlabel('Time (yrs)')
legend ('During Injection', 'During Extraction', 'Location',...
            'northeastoutside')
grid on

%-------------
% Internal Energy of Aquifer
figure(3)

plot(time,intEnergy)
ylabel('Internal Energy of Aquifer (J)')
xlabel('Time (yrs)')
grid on

%-------------
% Heat Conservation
figure(4)
plot(time, totalFluxin)
hold on 
plot(time, totalFluxout)
hold on 
plot(time,intEnergy)
hold on 
plot(time, fracEnergy)
hold on 
plot(time, rockEnergy)

ylabel('Energy of Aquifer (J)')
xlabel('Time (yrs)')
legend ('Total Flux in', 'Total Flux out', 'Internal Energy','Int Fracture Energy','Int Rock Energy', 'Location',...
            'northeastoutside')
grid on

%-------------
% Heat Ratio
%{
figure(5)
plot(time, heatRatio)
xlabel('Time (yrs)')
ylabel(['Ratio of Total Flux in to Internal Energy of Aquifer during' ...
    'each injection cycle'])
%}

figure(6)

plot(time, heatDiffratio , '*')
xlabel('Time (yrs)')
ylabel('Ratio of Total Flux in to Internal Energy of Aquifer per time step')
grid on 
ylim([0.9 1])

%------------- END OF CODE ------------------------------------------------
%

%