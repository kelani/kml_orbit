function [ r_3D_ef ] = cart2efix( r_3D_sf, th0 )
%CART2EFIX Rotates space-fixed system to earth-fixed system
%

R3 = @ (x) [cos(x),sin(x),0;-sin(x),cos(x),0;0,0,1];

r_3D_ef(1).r_vec = zeros(3,length(th0));
r_3D_ef(2).r_vec = zeros(3,length(th0));
r_3D_ef(3).r_vec = zeros(3,length(th0));
r_3D_ef(4).r_vec = zeros(3,length(th0));
r_3D_ef(5).r_vec = zeros(3,length(th0));

for it = 1:length (r_3D_sf)
    for k = 1: length (th0)
        r_3D_ef(it).r_vec(:,k) = (R3 (th0(k)) * r_3D_sf(it).r_vec(:,k));
    end
end
end
