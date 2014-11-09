function [ A, E ] = efix2topo2( sat_3d_ef, wlong, wlat )
R3 = @ (x) [cos(x),sin(x),0;-sin(x),cos(x),0;0,0,1];

ref = R3 (wlong)*sat_3d_ef;
R = @ (x,y) 6378 .*[cos(x) * cos(y); cos(x) * sin(y); sin(y)];
sef = ref - R (wlat, wlong);
eE = @ (x) [-sin(x);cos(x);0];
eN = @ (x,y) [-sin(x)*cos(y); -sin(x)*sin(y); cos(x)];
eZ = @ (x,y) [cos(x)*cos(y); cos(x)*sin(y); sin(x)];
E = [eE(wlat), eN(wlat, wlong), eZ(wlat, wlong)]';
s = E*sef;
A = atan2 (s(1),s(2));
E = atan2 (s(3),sqrt(s(1)^2+s(2)^2));
end

