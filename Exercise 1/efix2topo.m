function [azim,elev] = efix2topo(r_satellite, rw, RE, REP)
%EFIX2TOPO Transfoms position and velocity in an Earth-fixed system into
%position and velocity in a topocentric system centred at the station Wettzell

deg2rad = @(x) (x/180)*pi;
R2 = @ (x) [cos(x),0,-sin(x);0,1,0;sin(x),0,cos(x)];
R3 = @ (x) [cos(x),sin(x),0;-sin(x),cos(x),0;0,0,1];
Q1 = [-1, 0, 0; 0, 1, 0; 0, 0, 1];

e_topo = sqrt (RE^2-REP^2)/REP;
etopo = sqrt (RE^2-REP^2)/RE;
mtopo = atan2 ((rw(3)*RE),(sqrt(rw(1)^2 + rw(2)^2)*REP));
rwlat = atan2 ((rw(3)+e_topo^2*REP*sin(mtopo)^3), (sqrt(rw(1)^2 + rw(2)^2)-etopo^2*RE*cos(mtopo)^3));
rwlong = acos (rw(1)/ sqrt(rw(1)^2 + rw(2)^2));

azim = zeros (length (r_satellite),1);
elev = zeros (length (r_satellite),1);

for it = 1:length (r_satellite);
    azim (it,1:length (r_satellite(it).r_vec)) = zeros(1,length (r_satellite(it).r_vec));
    elev (it,1:length (r_satellite(it).r_vec)) = zeros(1,length (r_satellite(it).r_vec));
    for k = 1:length (r_satellite(1).r_vec);
        d = r_satellite(it).r_vec(k) - rw;
        dhor = (Q1*(R2((pi/2)-rwlat)*(R3(rwlong)*d)));             
        %azim (it,k) = atan2 (dhor(2), dhor(1));
        %elev (it,k) = (pi/2)- (acos (dhor(3)/sqrt(dhor(1)^2 + dhor(2)^2 + dhor(3)^2)));
        z = acos (dhor(3)/d(3));
        elev (it,k) = (pi/2)-z;
        azim (it,k) = asin ((dhor(2)/d(2))*sin(z));
    end
end

end

