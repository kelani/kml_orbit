function [ r_3D_ef, lat, long ] = cart2efix( r_3D_sf, th0 )
%CART2EFIX Rotates space-fixed system to earth-fixed system
%
rad2deg = @(x) (x*180)/pi;

R3 = @ (x) [cos(x),sin(x),0;-sin(x),cos(x),0;0,0,1];

r_3D_ef(1).r_vec = zeros(3,length(th0));
r_3D_ef(2).r_vec = zeros(3,length(th0));
r_3D_ef(3).r_vec = zeros(3,length(th0));
r_3D_ef(4).r_vec = zeros(3,length(th0));
r_3D_ef(5).r_vec = zeros(3,length(th0));
lat(1).r_vec = zeros (1, length (r_3D_ef));
lat(2).r_vec = zeros (1, length (r_3D_ef));
lat(3).r_vec = zeros (1, length (r_3D_ef));
lat(4).r_vec = zeros (1, length (r_3D_ef));
lat(5).r_vec = zeros (1, length (r_3D_ef));
long(1).r_vec = zeros (1, length (r_3D_ef));
long(2).r_vec = zeros (1, length (r_3D_ef));
long(3).r_vec = zeros (1, length (r_3D_ef));
long(4).r_vec = zeros (1, length (r_3D_ef));
long(5).r_vec = zeros (1, length (r_3D_ef));


for it = 1:length (r_3D_sf)
    for k = 1: length (th0)
        r_3D_ef(it).r_vec(:,k) = (R3 (th0(k)) * r_3D_sf(it).r_vec(:,k));
        lat(it).r_vec(:,k) = rad2deg (atan2 ((r_3D_ef(it).r_vec(2,k)),(r_3D_ef(it).r_vec(1,k))));
        long(it).r_vec(:,k) = rad2deg (atan2 ((r_3D_ef(it).r_vec(3,k)),(sqrt ((r_3D_ef(it).r_vec(1,k))^2 + (r_3D_ef(it).r_vec(2,k))^2))));
    end
end


end
