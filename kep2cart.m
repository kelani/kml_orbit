function [ r_3D, r_3D_dt ] = kep2cart ( orbit_X, GM, t, T0, error)
R1 = @ (x) [1,0,0;0,cos(x),sin(x);0,-sin(x),cos(x)];
R2 = @ (x) [cos(x),0,-sin(x);0,1,0;sin(x),0,cos(x)];
R3 = @ (x) [cos(x),sin(x),0;-sin(x),cos(x),0;0,0,1];

[n, M, E, r, v, x, y] = kep2orb(orbit_X, GM, t, T0, error);
temp(1).rb = zeros(3,length(x));
temp(2).rb = zeros(3,length(x));
temp(3).rb = zeros(3,length(x));
temp(4).rb = zeros(3,length(x));
temp(length(orbit_X)).rb = zeros(3,length(x));


for it = 1:(length(orbit_X))
    temp(it).rb (1,:) = x (it,:);
    temp(it).rb (2,:) = y (it, :);
    temp(it).rb (3, :) = zeros (1, length (x));
end

r_3D(1).r_vec = zeros(3,length(x) );
r_3D(2).r_vec = zeros(3,length(x) );
r_3D(3).r_vec = zeros(3,length(x) );
r_3D(4).r_vec = zeros(3,length(x) );
r_3D(5).r_vec = zeros(3,length(x) );

for it = 1:length (orbit_X)
    for k = 1: length (x)
        r_3D(it).r_vec(:,k) = (R3 (-orbit_X(it).bo) * (R1 (-orbit_X(it).i) * (R3 (-orbit_X(it).so) * temp(it).rb(:,k)) ));
    end
end

%% Velocity
rb_dt(1).vec = zeros(3,length(v));
rb_dt(2).vec = zeros(3,length(v));
rb_dt(3).vec = zeros(3,length(v));
rb_dt(4).vec = zeros(3,length(v));
rb_dt(5).vec = zeros(3,length(v));

for it = 1:length(orbit_X)
    rb_dt(it).vec (1,:) = sqrt (GM / (orbit_X(it).a*(1-(orbit_X(it).e)^2)))*(-sin(v(it,:)));
    rb_dt(it).vec (2,:) = sqrt (GM / (orbit_X(it).a*(1-(orbit_X(it).e)^2)))*(orbit_X(it).e + cos(v(it,:)));
    rb_dt(it).vec (3,:) = zeros(1,length(v));
end

r_3D_dt(1).r_vec = zeros(3,length(v) );
r_3D_dt(2).r_vec = zeros(3,length(v) );
r_3D_dt(3).r_vec = zeros(3,length(v) );
r_3D_dt(4).r_vec = zeros(3,length(v) );
r_3D_dt(5).r_vec = zeros(3,length(v) );

for it = 1:length (orbit_X)
    for k = 1: length (v)
        r_3D_dt(it).r_vec(:,k) = (R3 (-orbit_X(it).bo) * (R1 (-orbit_X(it).i) * (R3 (-orbit_X(it).so) * rb_dt(it).vec(:,k)) ));
    end
end

end

