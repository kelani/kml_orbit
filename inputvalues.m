clc
close all
clear all

GOCE = 1;
GPS = 2;
MOLNIYA = 3;
GEO = 4;
MICHIBIKI = 5;

deg2rad = @(x) (x/180)*pi;
rad2deg = @(x) (x*180)/pi;
GM = 398.6005*(10^12);

orbit(GOCE).a = 6629E3;
orbit(GOCE).e = 0.004;
orbit(GOCE).i = deg2rad (96.6);
orbit(GOCE).bo = deg2rad (257.7);
orbit(GOCE).so = deg2rad (144.2);
orbit(GOCE).T = (2*pi)*sqrt((orbit(GOCE).a^3)/GM);

orbit(GPS).a = 26560E3;
orbit(GPS).e = 0.01;
orbit(GPS).i = deg2rad (55);
orbit(GPS).bo = deg2rad (60);
orbit(GPS).so = deg2rad (0);
orbit(GPS).T = (2*pi)*sqrt((orbit(GPS).a^3)/GM);


orbit(MOLNIYA).a = 26554E3;
orbit(MOLNIYA).e = 0.7;
orbit(MOLNIYA).i = deg2rad (63);
orbit(MOLNIYA).bo = deg2rad (245);
orbit(MOLNIYA).so = deg2rad (270);
orbit(MOLNIYA).T = (2*pi)*sqrt((orbit(MOLNIYA).a^3)/GM);


orbit(GEO).a = 42164E3;
orbit(GEO).e = 0;
orbit(GEO).i = deg2rad (0);
orbit(GEO).bo = deg2rad (0);
orbit(GEO).so = deg2rad (0);
orbit(GEO).T = (2*pi)*sqrt((orbit(GEO).a^3)/GM);


orbit(MICHIBIKI).a = 42164E3;
orbit(MICHIBIKI).e = 0.075;
orbit(MICHIBIKI).i = deg2rad (41);
orbit(MICHIBIKI).bo = deg2rad (195);
orbit(MICHIBIKI).so = deg2rad (270);
orbit(MICHIBIKI).T = (2*pi)*sqrt((orbit(MICHIBIKI).a^3)/GM);


t = 1:1:(86400*2);
T0 = 0;
error = 1E-6;
RotR = ((2*pi) / 86164);
th0 = (RotR * t);


do_kep_2_orb_standalone = false;
plot_3D_pos = false;
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
%% Calculate 3D position and velocity
[r_3D, r_3D_dt] = kep2cart (orbit, GM, t, T0, error);

%% Plot position in 3D
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
    title('Orbit of 5 Satellites in 3-D space (Space-fixed System');
    xlabel('x(m)');
    ylabel('y(m)');
    zlabel('z(m)');
    
end
%% Plot Velocity in 2D
figure(3243564);
hold off;
grid on;
plot (t (1,:), sqrt ((r_3D_dt(GOCE).r_vec(1,:)).^2 + (r_3D_dt(GOCE).r_vec(2,:)).^2 + (r_3D_dt(GOCE).r_vec(3,:)).^2),...
    t (1,:), sqrt ((r_3D_dt(GPS).r_vec(1,:)).^2 + (r_3D_dt(GPS).r_vec(2,:)).^2 + (r_3D_dt(GPS).r_vec(3,:)).^2),...
    t (1,:), sqrt ((r_3D_dt(MOLNIYA).r_vec(1,:)).^2 + (r_3D_dt(MOLNIYA).r_vec(2,:)).^2 + (r_3D_dt(MOLNIYA).r_vec(3,:)).^2),...
    t (1,:), sqrt ((r_3D_dt(GEO).r_vec(1,:)).^2 + (r_3D_dt(GEO).r_vec(2,:)).^2 + (r_3D_dt(GEO).r_vec(3,:)).^2),...
    t (1,:), sqrt ((r_3D_dt(MICHIBIKI).r_vec(1,:)).^2 + (r_3D_dt(MICHIBIKI).r_vec(2,:)).^2 + (r_3D_dt(MICHIBIKI).r_vec(3,:)).^2));
axis([0 86400 1000 10000])
legend('GOCE', 'GPS', 'MOLNIYA', 'GEO','MICHIBIKI');
title('Velocity of 5 Satellites');
xlabel('t(s)');
ylabel('Velocity Magnitude(m/s)');
hold off;

%% From space-fixed system to earth-fixed system

[ r_3D_ef, lat, long ] = cart2efix ( r_3D, th0); % Position
[ r_3D_ef_dt ] = cart2efix ( r_3D_dt, th0); % Velocity


GOCE_2_rev_time = round(orbit(GOCE).T)*2;
GPS_2_rev_time = round(orbit(GPS).T)*2;
MOLNIYA_2_rev_time = round(orbit(MOLNIYA).T)*2;
GEO_2_rev_time = round(orbit(GEO).T)*2;
MICHIBIKI_2_rev_time = round(orbit(MICHIBIKI).T)*2;

%% Draw Orbit of 5 Satellites in 3-D space (Earth-fixed System)
figure(124123);
hold off;
axis equal;
plot3(r_3D_ef(GOCE).r_vec(1,1:GOCE_2_rev_time), r_3D_ef(GOCE).r_vec(2,1:GOCE_2_rev_time), r_3D_ef(GOCE).r_vec(3,1:GOCE_2_rev_time),...
    r_3D_ef(GPS).r_vec(1,1:GPS_2_rev_time), r_3D_ef(GPS).r_vec(2,1:GPS_2_rev_time), r_3D_ef(GPS).r_vec(3,1:GPS_2_rev_time),...
    r_3D_ef(MOLNIYA).r_vec(1,1:MOLNIYA_2_rev_time), r_3D_ef(MOLNIYA).r_vec(2,1:MOLNIYA_2_rev_time), r_3D_ef(MOLNIYA).r_vec(3,1:MOLNIYA_2_rev_time),...
    r_3D_ef(GEO).r_vec(1,1:GEO_2_rev_time), r_3D_ef(GEO).r_vec(2,1:GEO_2_rev_time), r_3D_ef(GEO).r_vec(3,1:GEO_2_rev_time),'*',...
    r_3D_ef(MICHIBIKI).r_vec(1,1:MICHIBIKI_2_rev_time), r_3D_ef(MICHIBIKI).r_vec(2,1:MICHIBIKI_2_rev_time), r_3D_ef(MICHIBIKI).r_vec(3,1:MICHIBIKI_2_rev_time));
axis([-4e7 6e7 -4e7 4e7 -4e7 6e7])
legend('GOCE', 'GPS', 'MOLNIYA', 'GEO','MICHIBIKI');
title('Orbit of 5 Satellites in 3-D space (Earth-fixed System)');
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)');
hold off;
%% Draw satellite ground-tracks on the Earth surface
figure (244000)
Earth_coast (2);
hold all;
ts = 60;
plot (lat(GOCE).r_vec(1:ts:GOCE_2_rev_time), long(GOCE).r_vec(1:ts:GOCE_2_rev_time), '*',...
    lat(GPS).r_vec(1:ts:GPS_2_rev_time), long(GPS).r_vec(1:ts:GPS_2_rev_time), '*',...
    lat(MOLNIYA).r_vec(1:ts:MOLNIYA_2_rev_time), long(MOLNIYA).r_vec(1:ts:MOLNIYA_2_rev_time), '*',...
    lat(GEO).r_vec(1:ts:GEO_2_rev_time), long(GEO).r_vec(1:ts:GEO_2_rev_time), '*',...
    lat(MICHIBIKI).r_vec(1:ts:MICHIBIKI_2_rev_time), long(MICHIBIKI).r_vec(1:ts:MICHIBIKI_2_rev_time), '*');
legend('EARTH','GOCE', 'GPS', 'MOLNIYA', 'GEO','MICHIBIKI', 'location','NorthEastOutside');
%%