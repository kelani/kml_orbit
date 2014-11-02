clc
close all
clear all

GOCE = 1;
GPS = 2;
MOLNIYA = 3;
GEO = 4;
MICHIBIKI = 5;

deg2rad = @(x) (x/180)*pi;

orbit(GOCE).a = 6629E3;
orbit(GOCE).e = 0.004;
orbit(GOCE).i = deg2rad (96.6);
orbit(GOCE).bo = deg2rad (257.7);
orbit(GOCE).so = deg2rad (144.2);

orbit(GPS).a = 26560E3;
orbit(GPS).e = 0.01;
orbit(GPS).i = deg2rad (55);
orbit(GPS).bo = deg2rad (60);
orbit(GPS).so = deg2rad (0);

orbit(MOLNIYA).a = 26554E3;
orbit(MOLNIYA).e = 0.7;
orbit(MOLNIYA).i = deg2rad (63);
orbit(MOLNIYA).bo = deg2rad (245);
orbit(MOLNIYA).so = deg2rad (270);

orbit(GEO).a = 42164E3;
orbit(GEO).e = 0;
orbit(GEO).i = deg2rad (0);
orbit(GEO).bo = deg2rad (0);
orbit(GEO).so = deg2rad (0);

orbit(MICHIBIKI).a = 42164E3;
orbit(MICHIBIKI).e = 0.075;
orbit(MICHIBIKI).i = deg2rad (41);
orbit(MICHIBIKI).bo = deg2rad (195);
orbit(MICHIBIKI).so = deg2rad (270);

GM = 398.6005*(10^12);
t = 1:1:86400;
T0 = 0;
error = 1E-6;


do_kep_2_orb_standalone = true;
plot_3D_pos = true;
%%

if (do_kep_2_orb_standalone == true)
    
    [n, M, E, r, v, x, y] = kep2orb(orbit, GM, t, T0, error);
    
    
    figure (20000);
    hold off;
    grid on;
    plot(x(1,:),y(1,:), x(2,:),y(2,:), x(3,:),y(3,:), x(4,:),y(4,:), x(5,:),y(5,:));
    legend('GOCE', 'GPS', 'MOLNIYA', 'GEO', 'MICHIBIKI');
    title('Orbit of 5 Satellites in 2-D plane');
    xlabel('x(m)');
    ylabel('y(m)');
    hold off;
    
    figure (30000);
    hold off;
    subplot(2,1,1);
    plot_end = 43000;
    d = (v-M);
    plot (t(1:plot_end), M (2,1:plot_end), t (1:plot_end), E (2,1:plot_end), t (1:plot_end), v (2,1:plot_end), t (1:plot_end), d (2,1:plot_end));
    legend('M', 'E', 'v', 'v-M');
    title('M,E,v,v-M of the GPS satellite');
    xlabel('time(s)');
    ylabel('angle(rad)');
    
    subplot(2,1,2);
    plot (t(1:plot_end), M (3,1:plot_end), t (1:plot_end), E (3,1:plot_end), t (1:plot_end), v (3,1:plot_end), t (1:plot_end), d (3,1:plot_end));
    legend('M', 'E', 'v', 'v-M');
    title('M,E,v,v-M of the MOLNIYA satellite');
    xlabel('time(s)');
    ylabel('angle(rad)');
    
end
%% calc 3D position and velocity
[r_3D, r_3D_dt] = kep2cart (orbit, GM, t, T0, error);

%% plot position in 3D
if(plot_3D_pos == true)
    figure(123123);
    hold off;
    axis equal;
    
    plot3(r_3D(GOCE).r_vec(1,:), r_3D(GOCE).r_vec(2,:), r_3D(GOCE).r_vec(3,:),...
        r_3D(GPS).r_vec(1,:), r_3D(GPS).r_vec(2,:), r_3D(GPS).r_vec(3,:),...
        r_3D(MOLNIYA).r_vec(1,:), r_3D(MOLNIYA).r_vec(2,:), r_3D(MOLNIYA).r_vec(3,:),...
        r_3D(GEO).r_vec(1,:), r_3D(GEO).r_vec(2,:), r_3D(GEO).r_vec(3,:),...
        r_3D(MICHIBIKI).r_vec(1,:), r_3D(MICHIBIKI).r_vec(2,:), r_3D(MICHIBIKI).r_vec(3,:));
    
    legend('GOCE', 'GPS', 'MOLNIYA', 'GEO','MICHIBIKI');
    title('Orbit of 5 Satellites in 3-D space');
    xlabel('x(m)');
    ylabel('y(m)');
    zlabel('z(m)');
    
end
%% plot Velocity in 2D
figure(3243564);
hold off;
axis equal;

plot3(r_3D_dt(GOCE).r_vec(1,:), r_3D_dt(GOCE).r_vec(2,:), r_3D_dt(GOCE).r_vec(3,:),...
    r_3D_dt(GPS).r_vec(1,:), r_3D_dt(GPS).r_vec(2,:), r_3D_dt(GPS).r_vec(3,:),...
    r_3D_dt(MOLNIYA).r_vec(1,:), r_3D_dt(MOLNIYA).r_vec(2,:), r_3D_dt(MOLNIYA).r_vec(3,:),...
    r_3D_dt(GEO).r_vec(1,:), r_3D_dt(GEO).r_vec(2,:), r_3D_dt(GEO).r_vec(3,:),...
    r_3D_dt(MICHIBIKI).r_vec(1,:), r_3D_dt(MICHIBIKI).r_vec(2,:), r_3D_dt(MICHIBIKI).r_vec(3,:));

legend('GOCE', 'GPS', 'MOLNIYA', 'GEO','MICHIBIKI');
title('Velocity of 5 Satellites in 3-D space');
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)');