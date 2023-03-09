%%CalculatingCohesion from We
%%SERA
theta_nb1 = 79.8;
theta_na2 = 82.1;
theta_nanb = 27.2;
H= -2444.1+2483.57; %Measured dif (Big Block: -2441.97+2518.64)
phi_i = atand(76.67/99.58); %rise over run or phi_i = 36
theta_24 = 59.3;
theta_35 = 61.6;
theta_45 = 42.2; 
theta_13 = 83.2;
phi_a = 35.6;
phi_b = 39.2;
zeta = 27.3;
beta = 49.55;

%%FIRSOFF
% theta_nb1 = 86.6;
% theta_na2 = 87.4;
% theta_nanb = 73.4;
% H= -2655.16+2679.10; %Measured dif (Big Block: -2441.97+2518.64)
% phi_i = atand((-2655.16+2679.10)/63.72); %rise over run or phi_i = 35.6
% theta_24 = 47.3;
% theta_35 = 45.4;
% theta_45 = 44.1; 
% theta_13 = 49.4;
% phi_a = 50.4;
% phi_b = 48.3;
% zeta = 73.6;
% beta = 85.1;


B = (cosd(phi_b)-cosd(phi_a)*cosd(theta_nanb))/(sind(phi_i)*sind(theta_nanb)^2);
A = (cosd(phi_a)-cosd(phi_b)*cosd(theta_nanb))/(sind(phi_i)*sind(theta_nanb)^2);
X = sind(theta_24)/(sind(theta_45)*cosd(theta_na2));
Y = sind(theta_13)/(sind(theta_35)*cosd(theta_nb1));

density = 1200:100:2000;
int_frict = 10:5:60;
for i = 1:length(density)
    for j = 1:length(int_frict)
        Cohesion(i,j) = (1-A.*B.*tand(int_frict(j)))./(X+Y).*density(i).*3.711.*76.67./3000;
    end
end


FRICTION = atand(tand(phi_i)*sind(zeta/2)/sind(beta));
for i = 1:length(density)
    Cohesion_F(i) = (1-A.*B.*tand(FRICTION))./(X+Y).*density(i).*3.711.*76.67./3000;
end
