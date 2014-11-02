function [n, M, E, r, v, x, y] = kep2orb( orbit_X, GM_x, t, T0, error )
n = sqrt (GM_x./(([orbit_X.a].').^3));
M = n * (t - T0);

E = zeros(length(orbit_X), length(M));
r = zeros(length(orbit_X), length(M));
v = zeros(length(orbit_X), length(M));
for i = 1:length(orbit_X)
    for k = 1:length(M)
        E(i,k) = kepler_E(orbit_X(i).e, M(i,k), error);
        r(i,k) = ([orbit_X(i).a].')*(1-([orbit_X(i).e])*cos(E(i,k)));
        v(i,k) = 2 * atan2 (sqrt (1 + orbit_X(i).e) * sin (E(i,k)/2), sqrt (1 - orbit_X(i).e) * cos (E(i,k)/2) );
        
    end
end

x = r .* cos(v);
y = r .* sin(v);

end