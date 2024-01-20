% EE553-Optimization, Final Project Landing Gear Response Optimization  Fall 2023-2024 METU
% Minimization by Fitted Functions
% Baris Canitez - 2033207

% Note-1: The "Efficiency_fit.mat" & "max_stroke_fit.mat" files should be present
% in the current folder to run the code. If they are not available, the
% values presented in the report can also be used directly
% Note-2: The optimum value and the optimum point is saved to the 
%workspace with the handles "min_val" and "x_iter" respectively
clear variables

syms x1 x2 x3 %Defining Desicion variables symbolically

x = [x1;x2;x3]; %Defining variables as a vector for ease of use

% Loading .mat files obtained as the output of
% The necessary files for loading into the workspace has been provided
%Loading elements of Q matrix and b vector representing the quadratic model fitted for objective function
efficiency_fit = load("Efficiency_fit.mat");
efficiency_fit = efficiency_fit.x_iter;

%Loading elements of the b vector representing the quadratic model fitted for objective function
max_stroke_fit = load("max_stroke_fit.mat");
max_stroke_fit = max_stroke_fit.x_iter;

% Substituting the loaded values into the Q matrix (objective function)
q11 = efficiency_fit(1);
q12 = efficiency_fit(2);
q13 = efficiency_fit(3);
q22 = efficiency_fit(4);
q23 = efficiency_fit(5);
q33 = efficiency_fit(6);
Q=[q11 q12 q13;q12 q22 q23;q13 q23 q33];

% Substituting the loaded values into the b vector (objective function)
b1 = efficiency_fit(7);
b2 = efficiency_fit(8);
b3 = efficiency_fit(9);
b = [b1;b2;b3];

% Substituting the loaded values into the a vector (constraint function)
a1 = max_stroke_fit(1);
a2 = max_stroke_fit(2);
a3 = max_stroke_fit(3);
a = [a1;a2;a3];
% Forming the A matrix including the additional constraints imposed by the domain of x
A = [a1 a2 a3;1 0 0;0 1 0;0 0 1;-1 0 0;0 -1 0;0 0 -1];

% Forming the c matrix consisting of allowable max. stroke and the limits of the domain of x
c = [0.19;4;10;400;-1;-6;-100];

f_o = -(0.5*transpose(x)*Q*x-transpose(b)*x); %Expressing objective function symbolically


phi = 0; %Preallocation for the barrier function
for i = 1:size(A,1)
phi = phi -log(-(A(i,:)*x-c(i))); %Obtaining barrier function symbolically
end

nabla_f_o = -Q*x+b; % Gradient of the objective function expressed symbolically
hess_f = -Q; %Hessian matrix of the objective function (quadratic !)

nabla_phi=zeros(size(x)); %Pre-allocation for obtaining gradient of the log. barrier function
% For-loop for obtaining the gradient of the log. barrier function symbolically
for i = 1:size(A,1)
nabla_phi = nabla_phi +(1/(c(i)-A(i,:)*x))*transpose(A(i,:)); %Gradient expression in Eqn. 22 in report
end 

hessian_phi=zeros(size(x,1),size(x,1));%Pre-allocation for obtaining Hessian of the log. barrier function
% For-loop for obtaining the Hessian of the log. barrier function symbolically
for i = 1:size(A,1)
hessian_phi = hessian_phi + 1/((c(i)-A(i,:)*x)^2)*transpose(A(i,:))*A(i,:); %Hessian expression in Eqn. 23 in report
end 

mu =10; %Selected mu value for magnifying t
t = 2; %Selected t value which defines accuracy of the barrier function
m = size(c); %Obtaining m for calculating duality gap
epsilon = 1e-7; %Tolerance for the duality gap value
x0 = [1.01;6.01;101]; %Initial value selected as the lower values of the domain of x

x_iter = x0; %Setting current x value as the initial point

while m/t>epsilon %Outer loop for the barrier method
x_old= [0;0;0]; %Resetting the previous x value befor starting each Newton's loop
while norm(x_old-x_iter)>epsilon %Inner loop (Newton's loop)
x_old=x_iter; %Setting previous x value to the current x value for checking stopping condition
grad_phi = eval(subs(nabla_phi,x,x_iter)); %Calculating gradient of log. the barrier function at current x value
hess_phi = eval(subs(hessian_phi,x,x_iter)); %Calculating Hessian of log. the barrier function at current x value
grad_f_o = eval(subs(nabla_f_o,x,x_iter)); %Calculating gradient of the objective function at current x value
delta_x_nt = -inv(t*hess_f+hess_phi)*(t*grad_f_o+grad_phi); %Calculating Newton's step at the current x value
x_iter = x_iter+delta_x_nt; %Updating current x value by the Newton's step calculated
min_val = eval(subs(-f_o,x,x_iter)) %Displaying value of the objective function
s_max = A(1,:)*x_iter  %Displaying value of the constraint function
end
t = mu*t; %Updating t value by the specified mu value
end
min_val = eval(subs(-f_o,x,x_iter)); %Optimum value of the shock absorber efficency