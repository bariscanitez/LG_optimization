% EE553-Optimization, Final Project Landing Gear Response Optimization  Fall 2023-2024 METU
% Objective Function Fitting
% Baris Canitez - 2033207

% Note-1: The r-square value is displayed in the command window and the
% optimum point is saved to the workspace with the handle "x_iter"
% Note-2: The code takes approximately 200 to 300 seconds to run due to large
% size and depending on the other software running. The variance of the timer is
%due to the randomized combinations of experiment matrix.

tic %Initiating timer

clear variables

global D_low D_flapper F_preload %Defining desicion variables as global variables

%% Forming the experiment matrix
exp_vec = linspace(0,7,8); %Expressing the normalized sampling intervals of the desicion variables with 8 points
CodedValue = combvec(exp_vec,exp_vec,exp_vec)'; %Taking all combinations of the normalized sampling points
runorder = randperm(size(CodedValue,1)); %Randomizing the normalized sample points
exp_D_low = [1 4]; %Sampling interval for the D_low parameter [mm]
exp_D_flapper = [6 10]; %Sampling interval for the D_flapper parameter [mm]
exp_F_preload = [100 400]; %Sampling interval for the F_preload parameter [mm]
bounds = [exp_D_low;exp_D_flapper;exp_F_preload]; %Expressing desicion variable intervals in a matrix

RealValue = zeros(size(CodedValue)); %Pre-allocation for the experiment matrix of actual values

%For-loop for converting normalized experiment matrix to the experiment matrix with actual values
for i = 1:size(CodedValue,2)
    zmax = max(CodedValue(:,i));
    zmin = min(CodedValue(:,i));
    RealValue(:,i) = interp1([zmin zmax],bounds(i,:),CodedValue(:,i));
end
%% Solution of the differential equation system

%Initial conditions for the differential equation system
sink_speed = 8.2; %Impact velocity of the landing gear [ft/s]
sink_speed = sink_speed*0.3048; %Unit conversion for the impact velocity [m/s]
ic1 = 0; %IC for cylinder (drop mass) displacement [m]
ic2 = sink_speed; %IC for cylinder (drop mass) velocity [m/s]
ic3 = 0; %IC for flapper displacement [m]
ic4 = sink_speed; %IC for flapper velocity [m/s]
ic5 = 0; %IC for piston displacement [m]
ic6 = sink_speed; %IC for piston velocity [m/s]

ics = [ic1 ic2 ic3 ic4 ic5 ic6]; %Initial condition vector

%Pre-allocation for shock absorber efficieny recorded for each model run
SA_efficiency = zeros(size(RealValue,1),1); 

%For-loop for solving differential equation system for all sample points
%and recording shock absorber efficiency values for each run
for i = 1:size(RealValue,1)
    D_low = RealValue(i,1); %Selecting current value of D_low
    D_flapper = RealValue(i,2); %Selecting current value of D_flapper
    F_preload = RealValue(i,3); %Selecting current value of F_preload
    
    %Solving the system
    [t,y] = ode23s(@(t,y) odefun(t,y),[0 0.18], ics);
    stroke = y(:,1)-y(:,5); %Post-processing for obtaining stroke [m]
    k_tire = 700000; %Defining tire stiffness for calculating ground force [N/m]
    F_tire = y(:,5)*k_tire; %Post-processing for ground force calculation [N]
    SA_energy = cumtrapz(stroke,F_tire); %Integrating ground force w.r.t stroke to obtain total work done by shock absorber [N*m]
    SA_efficiency(i) = SA_energy(end)/(max(F_tire)*max(stroke)); %Calculating shock absorber efficiency as in Eqn. 13 in report [-]
end
%% Fitting of the Objective Function
syms q11 q12 q13 q22 q23 q33 b1 b2 b3 %Defining elements of Q matrix & b vector symbolically
Q = [q11 q12 q13;q12 q22 q23;q13 q23 q33]; %Forming the symmetric Q matrix
b = [b1;b2;b3]; %Forming the b vector
% For loop for generating r vector
for i = 1:size(SA_efficiency,1)
    x_exp = RealValue(i,:).';
    r_store = SA_efficiency(i)-(0.5*x_exp.'*Q*x_exp-b.'*x_exp);
    r(i,1) =r_store;
end

%Initial point where Q is 3x3 identity matrix & b is vector of ones with size of 3
x_initial = [1 0 0 1 0 1 1 1 1].';
vars = [q11 q12 q13 q22 q23 q33 b1 b2 b3].'; %Defining vector of variables for ease of substitution
x_iter = x_initial; %Setting x value as the initial point
eps = 1e-7; %Tolearance for stopping condition
x_old = zeros(size(vars,2),1); %Preallocation for old x value


J = zeros(size(SA_efficiency,1),size(vars,1)); %Preallocation for Jacobian matrix
% For-loop for calculating the Jacobian matrix which is constant
for i=1:size(SA_efficiency,1)
    for j = 1:size(vars,1)
        J(i,j) = diff(r(i),vars(j));
    end
end

% While loop for Newton's Method
while norm(x_iter-x_old)/max(1,norm(x_old))>=eps %Stopping condition considering norm of x(k+1)- x(k)
    x_old = x_iter; %Storing x(k) value for checking stopping condition
    x_iter = x_iter-inv(transpose(J)*J)*transpose(J)*(eval(subs(r,vars,x_iter))); %Calculating the new x value by Gauss-Newton method
    r_val = eval(subs(r,vars,x_iter)); %Calculating the function value at the new x value
    r_val.'*r_val %r-square value calculated for display
end

toc %Stopping Timer
%% Differential Equation system as a MATLAB function
function dydt = odefun(t,y)
global D_low D_flapper F_preload %Globally defined desicion variables
% Distributing states
x1 = y(1); %drop mass displacement [m]
x1d = y(2); %drop mass velocity [m/s]
x2 = y(3);% flapper displacement [m]
x2d = y(4);  %flapper velocity [m/s]
x3 = y(5);% unsprung mass displacement [m]
x3d = y(6); % unsprung mass velocity [m/s]

% System Parameters
P_inf = 1.75; %Inflation pressure [MPa]
V_inf = 644000; %Inflation volume [mm^3]
POD = 60.27; %Piston outer diameter [mm]
D_nozzle = 5.8; %Nozzle diameter [mm]
% D_low = 3.2; %Low-speed hole diameter [mm]
N_low = 4; %Number of low-speed holes [-]
% D_flapper = 5.9; %Flapper diameter [mm]
Cd_low = 0.7; %Discharge coefficient of low-speed holes [-]
Cd_nozzle = 0.7; %Discharge coefficient of nozzle [-]
k_spring = 130; %Flapper spring stiffness [N/mm]
% F_preload = 300; %Spring preload [N]
gamma = 1.4; %Polytropic constent [-]
rho = 842.35;%Hydraulic fluid Density [kg/m^3]
k_tire = 700000; %Tire Stiffness
k_contact = 50000000000; %contact stiffness [N/m]
c_contact = 10000000; %Contact damping [Ns/m]
g = 9.80665; %Gravitational acceleration [m/s^2]

%Masses
m1 = 1780; %Drop mass [kg]
m2 = 0.05; %Flapper [kg]
m3 = 5; %Unsprung mass [kg]

% Relevant pneumatic & hydraulic area calculations
A_pod = (pi*POD^2)/4; %Piston outer diameter area [mm^2]

x_lim=(pi*(D_nozzle^2)/4)/(pi*D_flapper); %Limit flapper lift due to nozzle area
flp_lift = min([(x2-x3)*1000 x_lim]); %Instantaneou flapper lift [mm]

A_nozzle=pi*D_flapper*(flp_lift+abs(flp_lift))/2; %Nozzle area considering flapper movement [mm^2]


A_low = N_low*(pi*D_low^2)/4; %Total low-speed (taxi) hole area [mm^2]

% Force/Pressure definitions
P2 = P_inf*(V_inf/(V_inf-1000*(x1-x3)*A_pod))^gamma; %Pneumatic chamber pressure [MPa]

deltaP_overall = (rho*A_pod^2)/(2*((Cd_nozzle^2)*A_nozzle^2+(Cd_low^2)*A_low^2))*(x1d-x3d)*abs(x1d-x3d); %Overall pressure difference along the orifice [Pa]

deltaP_tot =deltaP_overall*(1e-006); %Unit conversion from Pa to MPa [MPa]

P1 = P2+deltaP_tot; % Pressure in the main hydraulic chamber [MPa]
F_spring = k_spring*1000*(x2-x3)+F_preload; %Spring force including preload [N]

% Contact force between flapper and flapper seat 
if x3>=x2 %Condition when the flapper penetrates the flapper seating
    F_contact = k_contact*(x3-x2)+c_contact*(x3d-x2d);
else %No contact force when the flapper is lifted
    F_contact = 0;
end

F_tire = x3*k_tire; %Calculating tire (ground) force [N]

dydt = zeros(size(y)); %Preallocation for the vector of derivative of states
dydt(1) = x1d; %1st differential equation for displacement of the drop mass displacement
dydt(2) = (m1*g-P1*A_pod)/m1; %2nd differential equation for drop mass velocity
dydt(3) = x2d; %3rd differential equation for flapper displacement
dydt(4) = (F_contact + (P1-P2)*((pi*D_nozzle^2)/4)-F_spring)/m2; %4th differential equation for flapper velocity
dydt(5) = x3d; %5th differential equation for unsprung mass displacement
dydt(6) = (P1*(A_pod-(pi*D_nozzle^2)/4)-F_tire-F_contact+F_spring)/m3; %6th differential equation for unsprung mass velocity
end