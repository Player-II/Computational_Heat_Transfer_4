%% CA4 Assignment
% MEGN 471 - Spring 2023
%Michael Allen, Cullen Hirstius, Hayeden Payne

format compact;
clc; clear; close all;

%% fluid properties 
%temperatures
T_s=27+273; %surface temperature [K]
T_inf=4+273; %fluid temperature [K]
T_f=0.5*(T_s+T_inf); %film temperature [K]
%properties from Table A4 (air)
alpha=20.982e-6; %thermal diffusivity [m^2/s]
rho=1.215; %density [kg/m^3]
mu=178.85e-7; %dynamic viscocity [Ns/m^2]
v=14.867e-6; %kinematic viscocity [m^2/s]
k=25.38e-3; %thermal conductivity [W/mK]
Pr=0.71; %Prandtl number [-]
Beta=1/T_f; %volumetric thermal expansion coeff. [1/K]
g=9.81; %acc. due to gravity [m/s^2]

%% can properties
H=(20:10:300)./1000; %heights of the cylinders [m]
D=(20:5:100)./1000; %diameters of the cylinders [m]


%% vertical orientation
for i=1:length(H) %loop through the heights for the vertical orientation
    %%%%% Calculate the Grashoff number for the vertical
    Gr_v(i)=( (g*Beta*(T_s-T_inf)*H(i)^3) / (v^2) ); %Grashoff number [-]
    %check bounds of validity for vertical cylinder
    Ra_v(i)=Gr_v(i)*Pr;
    Check_v(:,i)=D./H(i)-35./Gr_v(i)^.25; %checks the condition for each height/diameter if D/H>35/Gr^.25 (value should be positive if this is true)
    %Calculate Nusslt Number and corresponding h
    Nu_v(i) = (4/3)*(Gr_v(i)/4)^(1/4) * (.75*Pr^(1/2) / (.609+1.221*Pr^(1/2)+1.238*Pr)^(1/4)); %assuming laminar flow
    h_v(i)=Nu_v(i) *  k  / H(i); %convective heat transfer coeff. [W/m^2K]
end


%% horizontal orientation

for j=1:length(D) %loop through the diameters for the horzontal orientation
    %%%%% Calculate the Grashoff number for the horizontal
    Gr_h(j)=g*Beta*(T_s-T_inf)*D(j)^3 / (v^2); %Grashoff number [-]
    %check bounds
    Ra_h(j)=Gr_h(j)*Pr;
    %Check_v(:,i)=D./H(i)-35./Gr_v(i)^.25; %checks the condition for each height/diameter if D/H>35/Gr^.25 (value should be positive if this is true)
    %Calculate Nusslt Number and corresponding h
    Nu_h(j)=( .6 + .387*(Ra_h(j)^(1/6)) / ((1+(.559/Pr)^(9/16))^(8/27)) )^2; %Nusselt number [-]
    h_h(j)=Nu_h(j) * k /  D(j); %convective heat transfer coeff. [W/m^2K]
end

%% look up for specific dimensions
D_in=65/1000; %diameter for lookup [m]
H_in=120/1000; %height for lookup [m]
ind_D=find(D==D_in); %finds index of desired diameter
ind_H=find(H==H_in); %finds index of desired height
fprintf('At the diameter D=%0.3f m and H=%0.3f m, the vertical orientation has an h of %0.3f W/m^2K \n and the horizontal orientation had an h of %0.3f W/m^2K', D_in, H_in, h_v(ind_H), h_h(ind_D));

%% plot results
figure('color', 'w','DefaultAxesFontSize',16)
set(gcf, 'Position', [100, 100, 600, 800]) %opens figure to larger size
subplot(2,1,1) %vertical orientation subplot
hold on; grid on; box on;
plot(H, h_v, 'LineWidth', 2);
xlabel('Can Height [m]')
ylabel('$\overline{h}\ [W/m^2K]$', 'interpreter', 'LaTeX')
title('Vertical Orientation')
xlim([min(H), max(H)])
subplot(2,1,2) %horizontal orientation subplot
hold on; grid on; box on;
plot(D, h_h, 'LineWidth', 2);
xlabel('Can Diameter [m]')
ylabel('$\overline{h}\ [W/m^2K]$', 'interpreter', 'LaTeX')
title('Horizontal Orientation')
xlim([min(D), max(D)])

%% Vertical Again
Gr_v=[]
Check_v=[]
Nu_v = []
h_v = []
for i=1:length(D) %loop through the heights for the vertical orientation
    %%%%% Calculate the Grashoff number for the vertical
    Gr_v(i)=( (g*Beta*(T_s-T_inf)*H_in.^3) / (v^2) ); %Grashoff number [-]
    %check bounds of validity for vertical cylinder
    Ra_v(i)=Gr_v(i)*Pr;
    Check_v(:,i)=D(i)./H-35./Gr_v(i)^.25; %checks the condition for each height/diameter if D/H>35/Gr^.25 (value should be positive if this is true)
    %Calculate Nusslt Number and corresponding h
    Nu_v(i) = (4/3)*(Gr_v(i)/4)^(1/4) * (.75*Pr^(1/2) / (.609+1.221*Pr^(1/2)+1.238*Pr)^(1/4)); %assuming laminar flow
    h_v(i)=Nu_v(i) *  k  / H_in; %convective heat transfer coeff. [W/m^2K]
end
%% Horizontal Again

for j=1:length(H) %loop through the diameters for the horzontal orientation
    %%%%% Calculate the Grashoff number for the horizontal
    Gr_h(j)=g*Beta*(T_s-T_inf)*D_in^3 / (v^2); %Grashoff number [-]
    %check bounds
    Ra_h(j)=Gr_h(j)*Pr;
    %Check_v(:,i)=D./H(i)-35./Gr_v(i)^.25; %checks the condition for each height/diameter if D/H>35/Gr^.25 (value should be positive if this is true)
    %Calculate Nusslt Number and corresponding h
    Nu_h(j)=( .6 + .387*(Ra_h(j)^(1/6)) / ((1+(.559/Pr)^(9/16))^(8/27)) )^2; %Nusselt number [-]
    h_h(j)=Nu_h(j) * k /  D_in; %convective heat transfer coeff. [W/m^2K]
end

%% More Plots
figure('color', 'w','DefaultAxesFontSize',16)
set(gcf, 'Position', [100, 100, 600, 800]) %opens figure to larger size
subplot(2,1,1) %vertical orientation subplot
hold on; grid on; box on;
plot(D, h_v, 'LineWidth', 2);
xlabel('Can Width [m]')
ylabel('$\overline{h}\ [W/m^2K]$', 'interpreter', 'LaTeX')
title('Vertical Orientation')
xlim([min(D), max(D)])
subplot(2,1,2) %horizontal orientation subplot
hold on; grid on; box on;
plot(H, h_h, 'LineWidth', 2);
xlabel('Can Height [m]')
ylabel('$\overline{h}\ [W/m^2K]$', 'interpreter', 'LaTeX')
title('Horizontal Orientation')
xlim([min(H), max(H)])
