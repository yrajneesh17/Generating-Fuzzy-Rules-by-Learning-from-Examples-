clc;
close all;
clear all ;
 
%Loading the Data sequences generated from 14 initial states
load tb1.dat 
load tb2.dat 
load tb3.dat 
load tb4.dat 
load tb5.dat 
load tb6.dat 
load tb7.dat 
load tb8.dat 
load tb9.dat 
load tb10.dat 
load tb11.dat 
load tb12.dat 
load tb13.dat 
load tb14.dat
 
table=[tb1; tb2; tb3; tb4; tb5; tb6; tb7; tb8; tb9; tb10; tb11; tb12; tb13; tb14];
 
%b=inputdlg('Enter the length of the Truck','Value of b',1,{'4'});
b=4;      %Length of the truck (str2double(b))
iteration = 50;     % No. of iterations for trajectory plotting
 
%%Step 1:- Divide the input output space of the given numerical data into
%%fuzzy regions 
% Divide mx,mphi mtheta into 2N+1 regions
 
x=0:0.1:20 ; %Universe of discourse of input x
mx(:,1)=trapmf(x,[0 0 1.5 7]);
mx(:,2)=trimf(x,[4 7 10]);
mx(:,3)=trimf(x,[9 10 11]);
mx(:,4)=trimf(x,[10 13 16]);
mx(:,5)=trapmf(x,[13 18.5 20 20]);
%plot the membership functions of x
figure(1)
plot(x,mx,'Linewidth',2.0);
xlabel('x');
ylabel('m(x)');
title('Input Memebrship Function - 1');
legend('S2','S1','CE','B1','B2');
legend boxoff;
ylim([0 1]);
set(gca,'Fontsize',13,'FontName','Times New Roman');
 
phi = -115:0.1:295; % Universe of discourse of input phi.
mphi(:,1) = trimf(phi,[-115 -65 -15]);
mphi(:,2) = trimf(phi,[-45 0 45]);
mphi(:,3) = trimf(phi,[15 52.5 90]);
mphi(:,4) = trimf(phi,[80 90 100]);
mphi(:,5) = trimf(phi,[90 127.5 165]);
mphi(:,6) = trimf(phi,[135 180 225]);
mphi(:,7) = trimf(phi,[195 245 295]);
 
 % plot the MFs of phi.
figure(2)
plot(phi,mphi,'Linewidth',2.0);
xlabel('\phi');
ylabel('m(\phi)');
title('Input Memebrship Function - 2');
legend('S3','S2','S1','CE','B1','B2','B3');
legend boxoff;
ylim([0 1]);
set(gca,'Fontsize',13,'FontName','Times New Roman');
 
theta=-40:0.1:40 ;%universe of discourse of output theta
mtheta(:,1)=trimf(theta,[-40 -40 -20]);
mtheta(:,2)=trimf(theta,[-33 -20 -7]);
mtheta(:,3)=trimf(theta,[-14 -7 0]);
mtheta(:,4)=trimf(theta,[-4 0 4]);
mtheta(:,5)=trimf(theta,[0 7 14]);
mtheta(:,6)=trimf(theta,[7 20 33]);
mtheta(:,7)=trimf(theta,[20 40 40]);
% plot the MFs of theta.
figure(3)
plot(theta,mtheta,'Linewidth',2.0);
xlabel('\theta');
ylabel('m(\theta)');
title('Output Memebrship Function - 1');
legend('S3','S2','S1','CE','B1','B2','B3');
legend boxoff;
ylim([0 1]);
set(gca,'Fontsize',13,'FontName','Times New Roman');
 
%% Step 2: Generating Fuzzy Rules from Given 14 Data Pairs.
% Calculating the Degree for given datasets
% Input x
mu_xS2 = trapmf(table(:,1),[0 0 1.5 7]);
mu_xS1 = trimf(table(:,1),[4 7 10]);
mu_xCE = trimf(table(:,1),[9 10 11]);
mu_xB1 = trimf(table(:,1),[10 13 16]);
mu_xB2 = trapmf(table(:,1),[13 18.5 20 20]);
mu_x = [mu_xS2 mu_xS1 mu_xCE mu_xB1 mu_xB2];
 
% Input phi
mu_phiS3 = trimf(table(:,2),[-115 -65 -15]);
mu_phiS2 = trimf(table(:,2),[-45 0 45]);
mu_phiS1 = trimf(table(:,2),[15 52.5 90]);
mu_phiCE = trimf(table(:,2),[80 90 100]);
mu_phiB1 = trimf(table(:,2),[90 127.5 165]);
mu_phiB2 = trimf(table(:,2),[135 180 225]);
mu_phiB3 = trimf(table(:,2),[195 245 295]);
mu_phi = [mu_phiS3 mu_phiS2 mu_phiS1 mu_phiCE mu_phiB1 mu_phiB2 mu_phiB3];
 
%Output theta
mu_thetaS3 = trimf(table(:,3),[-40 -40 -20]);
mu_thetaS2 = trimf(table(:,3),[-33 -20 -7]);
mu_thetaS1 = trimf(table(:,3),[-14 -7 0]);
mu_thetaCE = trimf(table(:,3),[-4 0 4]);
mu_thetaB1 = trimf(table(:,3),[0 7 14]);
mu_thetaB2 = trimf(table(:,3),[7 20 33]);
mu_thetaB3 = trimf(table(:,3),[20 40 40]);
mu_theta = [mu_thetaS3 mu_thetaS2 mu_thetaS1 mu_thetaCE mu_thetaB1 mu_thetaB2 mu_thetaB3];
 
%degree and index
[avg_degree_x,avg_index_x] = max(mu_x,[],2);  
%Returning column vector ,max value in each row and their corresponding row
%location
%max(mu_x_test,[],1) then max value ,row vector,max value in each column
%and their corresponding column location 
[avg_degree_phi,avg_index_phi] = max(mu_phi,[],2);
[avg_degree_theta,avg_index_theta] = max(mu_theta,[],2);
%Defining Rules
rules = [avg_index_x avg_index_phi avg_index_theta];
rules_degree = avg_degree_x.*avg_degree_phi.*avg_degree_theta;
 
%% Step 3: Assign a Degree to Each Rule.
% For conflicting rule in rule base, assign the maximum degree for the
% conflict group.
R_temp=[rules rules_degree];
n=length(R_temp);
for i=1:n
    for j=1:n
        if((R_temp(i,1)==R_temp(j,1))&&(R_temp(i,2)==R_temp(j,2)))
            if(R_temp(i,4)>=R_temp(j,4))
            R_temp(j,3)=R_temp(i,3);
            R_temp(j,4)=R_temp(i,4);
            else
               R_temp(i,3)=R_temp(j,3);
               R_temp(i,4)=R_temp(j,4); 
            end
        end
    end
end
final_rules=unique(R_temp,'rows','sorted');%It removes redundent rows 
                                           % and sort the Data
fprintf('The rule base is ready. \n');
 
%Step 4 :-Create a Combined Fuzzy Rule base
%Not requires in this Truck Backer-upper because this step is used only where or
%connectors are used
 
%% Step 5: Determine a Mapping Based on the Combined Fuzzy Rule Base.
%Defuzzification is also done in this step ,Centroid defuzzification is
%used to determine the output .
figure;
axis([0 20 -10 100]); % Defining x-axis and y-axis limits
iter = 65; % number of iterations for trajectory plotting.
y_t = [-40 -20 -7 0 7 20 40]; % the points where centroid is lying .
 
%% Test Input 1
x_input = 10; % Input x (by user)
phi_input = 220; % Input phi (by user)
y_sample = 2; % initial y position
 
y_final = zeros(1,1);
x_final = zeros(1,1);
x_final(1,1) = x_input;
y_final(1,1) = y_sample;
Input1 = strcat('(',num2str(x_input),',',num2str(phi_input),')');
 
for i = 1:iter % no. of iterations for trajectory tracking.
 
    % value of x.
    mu_x1_test = trapmf(x_input,[0 0 1.5 7]);
    mu_x2_test = trimf(x_input,[4 7 10]);
    mu_x3_test = trimf(x_input,[9 10 11]);
    mu_x4_test = trimf(x_input,[10 13 16]);
    mu_x5_test = trapmf(x_input,[13 18.5 20 20]);
 
    % value of phi.
    mu_phi1_test = trimf(phi_input,[-115 -65 -15]);
    mu_phi2_test = trimf(phi_input,[-45 0 45]);
    mu_phi3_test = trimf(phi_input,[15 52.5 90]);
    mu_phi4_test = trimf(phi_input,[80 90 100]);
    mu_phi5_test = trimf(phi_input,[90 127.5 165]);
    mu_phi6_test = trimf(phi_input,[135 180 225]);
    mu_phi7_test = trimf(phi_input,[195 245 295]);
 
    mu_x_test = [mu_x1_test mu_x2_test mu_x3_test mu_x4_test mu_x5_test];
    mu_p_test = [mu_phi1_test mu_phi2_test mu_phi3_test mu_phi4_test mu_phi5_test mu_phi6_test mu_phi7_test];
    mo = mu_x_test(final_rules(:,1)).*mu_p_test(final_rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(final_rules(:,3));
    theta_output = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
 
    % Approximate kinematics of truck backer upper control for calculating the next states.
    x_input = x_input + cosd(phi_input + theta_output) + sind(theta_output)*sind(phi_input);
    x_final(i+1,1) = x_input;
    y_sample = y_sample + sind(phi_input + theta_output) - sind(theta_output)*cosd(phi_input); % displacement on y axis.
    y_final(i+1,1) = y_sample;
    phi_input = phi_input - asind(2*sind(theta_output)/b);
 
end
 
fprintf('Plotting the trajectory for %s....... \n',Input1);
 
% Plot the trajectory.
plot(x_final,y_final,'.-','MarkerSize',12.0);
text(x_final(1,1),y_final(1,1)+1,Input1);
hold on;
fprintf('Trajectory for %s has been plotted. \n',Input1);
 
%% Test Input 2
x_input = 3; % Input x (by user)
phi_input = -30; % Input phi (by user)
y_sample = 2; % initial y position
 
y_final = zeros(1,1);
x_final = zeros(1,1);
x_final(1,1) = x_input;
y_final(1,1) = y_sample;
Input2 = strcat('(',num2str(x_input),',',num2str(phi_input),')');
 
for i = 1:iter % no. of iterations for trajectory tracking.
 
    % value of x.
    mu_x1_test = trapmf(x_input,[0 0 1.5 7]);
    mu_x2_test = trimf(x_input,[4 7 10]);
    mu_x3_test = trimf(x_input,[9 10 11]);
    mu_x4_test = trimf(x_input,[10 13 16]);
    mu_x5_test = trapmf(x_input,[13 18.5 20 20]);
 
    % value of phi.
    mu_phi1_test = trimf(phi_input,[-115 -65 -15]);
    mu_phi2_test = trimf(phi_input,[-45 0 45]);
    mu_phi3_test = trimf(phi_input,[15 52.5 90]);
    mu_phi4_test = trimf(phi_input,[80 90 100]);
    mu_phi5_test = trimf(phi_input,[90 127.5 165]);
    mu_phi6_test = trimf(phi_input,[135 180 225]);
    mu_phi7_test = trimf(phi_input,[195 245 295]);
 
    mu_x_test = [mu_x1_test mu_x2_test mu_x3_test mu_x4_test mu_x5_test];
    mu_p_test = [mu_phi1_test mu_phi2_test mu_phi3_test mu_phi4_test mu_phi5_test mu_phi6_test mu_phi7_test];
    mo = mu_x_test(final_rules(:,1)).*mu_p_test(final_rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(final_rules(:,3));
    theta_output = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
 
    % Approximate kinematics of truck backer upper control for calculating the next states.
    x_input = x_input + cosd(phi_input + theta_output) + sind(theta_output)*sind(phi_input);
    x_final(i+1,1) = x_input;
    y_sample = y_sample + sind(phi_input + theta_output) - sind(theta_output)*cosd(phi_input); % displacement on y axis.
    y_final(i+1,1) = y_sample;
    phi_input = phi_input - asind(2*sind(theta_output)/b);
 
end
 
fprintf('Plotting the trajectory for %s....... \n',Input2);
 
% Plot the trajectory.
plot(x_final,y_final,'.-','MarkerSize',12.0);
text(x_final(1,1),y_final(1,1)+1,Input2);
hold on;
fprintf('Trajectory for %s has been plotted. \n',Input2);
 
%% Test Input 3
x_input = 13; % Input x (by user)
phi_input = 30; % Input phi (by user)
y_sample = 2; % initial y position
 
y_final = zeros(1,1);
x_final = zeros(1,1);
x_final(1,1) = x_input;
y_final(1,1) = y_sample;
Input3 = strcat('(',num2str(x_input),',',num2str(phi_input),')');
 
for i = 1:iter % no. of iterations for trajectory tracking.
 
    % value of x.
    mu_x1_test = trapmf(x_input,[0 0 1.5 7]);
    mu_x2_test = trimf(x_input,[4 7 10]);
    mu_x3_test = trimf(x_input,[9 10 11]);
    mu_x4_test = trimf(x_input,[10 13 16]);
    mu_x5_test = trapmf(x_input,[13 18.5 20 20]);
 
    % value of phi.
    mu_phi1_test = trimf(phi_input,[-115 -65 -15]);
    mu_phi2_test = trimf(phi_input,[-45 0 45]);
    mu_phi3_test = trimf(phi_input,[15 52.5 90]);
    mu_phi4_test = trimf(phi_input,[80 90 100]);
    mu_phi5_test = trimf(phi_input,[90 127.5 165]);
    mu_phi6_test = trimf(phi_input,[135 180 225]);
    mu_phi7_test = trimf(phi_input,[195 245 295]);
 
    mu_x_test = [mu_x1_test mu_x2_test mu_x3_test mu_x4_test mu_x5_test];
    mu_p_test = [mu_phi1_test mu_phi2_test mu_phi3_test mu_phi4_test mu_phi5_test mu_phi6_test mu_phi7_test];
    mo = mu_x_test(final_rules(:,1)).*mu_p_test(final_rules(:,2)); % product operation to determine the degree of output control.
    y_bar = y_t(final_rules(:,3));
    theta_output = sum(mo.*y_bar)/sum(mo); % Value of theta from the rule base.
 
    % Approximate kinematics of truck backer upper control for calculating the next states.
    x_input = x_input + cosd(phi_input + theta_output) + sind(theta_output)*sind(phi_input);
    x_final(i+1,1) = x_input;
    y_sample = y_sample + sind(phi_input + theta_output) - sind(theta_output)*cosd(phi_input); % displacement on y axis.
    y_final(i+1,1) = y_sample;
    phi_input = phi_input - asind(2*sind(theta_output)/b);
 
end
 
fprintf('Plotting the trajectory for %s....... \n',Input3);
 
% Plot the trajectory.
plot(x_final,y_final,'.-','MarkerSize',12.0);
text(x_final(1,1),y_final(1,1)+1,Input3);
hold on;
fprintf('Trajectory for %s has been plotted. \n',Input3);
 
xlabel('x (in meters)');
ylabel('y (in meters)');
ylim([0 65]);
title('Truck Trajectories using the Truncated Data Pairs');
legend('Trajectory for (10,220)','Trajectory for (3,-30)','Trajectory for (13,30)');
legend boxoff;
set(gca,'Fontsize',13,'FontName','Times New Roman');            


