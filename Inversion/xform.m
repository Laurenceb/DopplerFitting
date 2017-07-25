function p2 = xform( op, p1, c1, c2 )

%XFORM Transform between earth coordinate systems.
% xform( 'init', N ) allocate N NED coord systems
% xform( 'new', P1, C1, Cn ) new Cn NED at P1(C1)
% P2 = xform( 'pos', P1, C1, C2 ) P2(C2) == P1(C1)
% V2 = xform( 'vel', V1, C1, C2 ) V2(C2) == V1(C1)
% K2 = xform( 'cov', K1, C1, C2 ) K2(C2) == K1(C1)
%
% where Pi are point vectors (1x3), Vi are velocity vectors (1x3)
% and Ki are covariance matrices (3x3). Ci are coordinate system
% references with the following numbering scheme:
% Cn = -1 LLA p = [ lat lon alt ] in radians and meters
% = 0 ECEF p = [ x y z ] in meters
% = 1 NED # 1 p = [ x y z ] in meters
% = 2 NED # 2 p = [ x y z ] in meters
% = N NED # N p = [ x y z ] in meters
%
% For example, place the third NED coordinate system origin at
% N30 00 00.0, W90 00 00.0 and 100m elevation and then find
% the ECEF coordinates of a point 200 meters to the North:
% xform( 'new', [pi/6 -pi/2 100], -1, 3 );
% p = xform( 'pos', [200 0 0], 3, 0 );

global XFORM


switch op


% P2 = xform( 'pos', P1, C1, C2 ) P2(C2) <== P1(C1)
% Transform position P1 in coord system C1
% to position P2 in coord system C2
%
case 'pos'

   % transform P1(C1) to ECEF
   %
   if c1 < 0, % p1 is LLA
      p0 = lla2ecef(p1)';
   elseif c1 == 0, % p1 is ECEF
      p0 = p1';
   else % p1 is NED
      p0 = XFORM(c1).Cen * p1' + ... % XFORM(c1).org;
            repmat( XFORM(c1).org, 1, size(p1,1) );
   end
   
   % transform P0 to C2
   %
   if c2 < 0, % p2 is LLA
      p2 = ecef2lla(p0');
   elseif c2 == 0, % p2 is ECEF
      p2 = p0';
   else % p2 is NED
      p2 = XFORM(c2).Cen' * ...
         ( p0 - repmat( XFORM(c2).org, 1, size(p0,2) ) );
      p2 = p2';
   end
   

   
% V2 = xform( 'pos', V1, C1, C2 ) V2(C2) <== V1(C1)
% Transform velocity V1 in coord system C1
% to velocity V2 in coord system C2
% Velocities can be in any units
%
case 'vel'

   % transform V1(C1) to ECEF
   %
   if c1 < 0,
      error('LLA not valid for ''vel''');
   elseif c1 == 0, % v1 is ECEF
      p0 = p1';
   else % v1 is NED
      p0 = XFORM(c1).Cen * p1';
   end
   
   % transform P0 to C2
   %
   if c2 < 0,
      error('LLA not valid for ''vel''');
   elseif c2 == 0, % v2 is ECEF
      p2 = p0';
   else % v2 is NED
      p2 = XFORM(c2).Cen' * p0;
      p2 = p2';
   end



% P2 = xform( 'cov', P1, C1, C2 ) P2(C2) <== P1(C1)
% Transform covariance P1 in coord system C1
% to covariance P2 in coord system C2
%
case 'cov'
   % transform P1(C1) to ECEF
   %
   if c1 < 0,
      error('LLA not valid for ''cov''');
   elseif c1 == 0, % p1 is ECEF
      p0 = p1;
   else % p1 is NED
      p0 = XFORM(c1).Cen * p1 * XFORM(c1).Cen';
   end
   
   % transform P0 to C2
   %
   if c2 < 0,
      error('LLA not valid for ''cov''');
   elseif c2 == 0, % p2 is ECEF
      p2 = p0;
   else % p2 is NED
      p2 = XFORM(c2).Cen' * p0 * XFORM(c2).Cen;
   end

   

% xform( 'new', P1, C1, Ci ) -- define new NED # Ci at P1(C1)
%
case 'new'
   if isempty(XFORM),
      xform init;
   end
   
   % transform P1(C1) to ECEF
   %
   if c1 < 0, % p1 is LLA
      lat = p1(1);
      lon = p1(2);
      XFORM(c2).org = lla2ecef(p1)';
   elseif c1 == 0, % p1 is ECEF
      XFORM(c2).org = p1';
      x = ecef2lla(p1);
      lat = x(1);
      lon = x(2);
   else % p1 is NED
      p0 = XFORM(c1).Cen * p1' + XFORM(c1).org;
      XFORM(c2).org = p0;
      lat = p0(1);
      lon = p0(2);
   end
   
   % Transformation matrix:
   % For column vectors: p.ecef = Cen * p.ned
   % For row vectors: p.ecef = p.ned * Cen'
   %
   XFORM(c2).Cen = [ -sin(lat)*cos(lon), -sin(lon), -cos(lat)*cos(lon);
                  -sin(lat)*sin(lon), cos(lon), -cos(lat)*sin(lon);
                   cos(lat), 0, -sin(lat) ];


% xform( 'init', N ) -- allocate N NED coord systems
%
case 'init'
   if nargin < 2,
      p1 = 10; % default is 10 NED coord systems
   end
   XFORM = struct( 'org', cell(p1,1) ,'Cen', eye(3) );
   [XFORM.org] = deal( zeros(3,1) );
   xform( 'new', [0 0 0], -1, 1 );


otherwise
   error( 'Bad op code' );

end



% LLA to ECEF using WGS-84
% from "Understanding GPS Principles and Applications"
% section 2.2.3
%
function ecef = lla2ecef( lat_lon_alt )

a = 6378.137e3; % earth radius @ equator
f = 1 / 298.257223563; % flattening
b = a * (1-f); % earth radius @ poles
e2 = 1 - (1-f)^2; % eccentricity squared

[n,m] = size( lat_lon_alt );
if m ~= 3,
   error( 'Input must have three columns' );
end

ecef = zeros( n, m );

lat = lat_lon_alt(:,1);
lon = lat_lon_alt(:,2);
alt = lat_lon_alt(:,3);

ecef = [ a*cos(lon) ./ sqrt( 1 + (1-e2)*tan(lat).^2 ) + alt.*cos(lon).*cos(lat), ...
        a*sin(lon) ./ sqrt( 1 + (1-e2)*tan(lat).^2 ) + alt.*sin(lon).*cos(lat), ...
        a*(1-e2)*sin(lat) ./ sqrt( 1 - e2*sin(lat).^2 ) + alt.*sin(lat) ];


% ECEF to LLA using WGS-84
% from "Understanding GPS Principles and Applications"
% section 2.2.3
%
function lla = ecef2lla( ecef )

a = 6378.137e3; % earth radius @ equator
f = 1 / 298.257223563; % flattening
b = a * (1-f); % earth radius @ poles
e2 = 1 - (1-f)^2; % eccentricity squared
ep = sqrt((a/b)^2-1); % second eccentricity

[n,m] = size( ecef );
if m ~= 3,
   error( 'Input must have three columns' );
end

lat_lon_alt = zeros( n, m );

x = ecef(:,1);
y = ecef(:,2);
z = ecef(:,3);

r = sqrt( x.^2 + y.^2 );

E2 = a^2 - b^2;

F = 54 * ( b*z).^2;

G = r.^2 + (1-e2) * z.^2 - e2 * E2;

c = ( e2^2 * F .* r.^2 ) ./ G.^3;

s = ( 1 + c + sqrt(c.^2 + 2*c) ).^(1/3);

P = F ./ ( 3 * (s + 1./s + 1).^2 .* G.^2 );

Q = sqrt( 1 + 2*e2^2 * P );

temp = a^2/2*(1+1./Q ) - (1-e2)*P.*z.^2 ./ (Q.*(1+Q)) - P.*r.^2/2;
r0 = -e2*P.*r ./ (1+Q) + sqrt( max(temp,0) );

U = sqrt( ( r - e2*r0).^2 + z.^2 );
V = sqrt( ( r - e2*r0).^2 + (1-e2)*z.^2 );

z0 = b^2 * z ./ (a*V);

lla = [ atan2( ( z + ep^2*z0 ), r ), ... % lat
      atan2( y, x ), ... % long
      U.*( 1 - b^2./(a*V) ) ]; % alt
return


% Ellipsoid name Semi-major axis Inverse flattenning
% (a) (1/f)
%
% Airy 6377563.3960 299.324964600
% Airy_modified 6377340.1890 299.324964600
% Australian_National 6378160.0000 298.250000000
% Bessel 1841 6377397.1550 299.152812800
% Bessel 1841 in Namibia 6377483.8650 299.152812800
% Clarke 1866 6378206.4000 294.978698200
% Clarke 1880 6378249.1450 293.465000000
% Everest (Sabah & Sarawak) 6377298.5560 300.801700000
% Everest 1830 6377276.3450 300.801700000
% Everest 1948 6377304.0630 300.801700000
% Everest 1956 6377301.2430 300.801700000
% Everest_Modified 6377304.0630 300.801700000
% GRS-80 6378137.0000 298.257222101
% Helmert 1906 6378200.0000 298.300000000
% Hough 6378270.0000 297.000000000
% International 6378388.0000 297.000000000
% Krassovsky 6378245.0000 298.300000000
% Modified Fisher 1960 6378155.0000 298.300000000
% SGS 85 6378136.0000 298.257000000
% South America 1969 6378160.0000 298.250000000
% WGS-72 6378135.0000 298.260000000
% WGS-84 6378137.0000 298.257223563


% e2 = 1 - (b/a)^2; % eccentricity squared
% f = 1 - b/a; % flattening
% ep = sqrt((a/b)^2-1); % second eccentricity

% ep = 0.0820944379496; % second eccentricity
% e2 = ep^2 / (ep^2+1); % eccentricity squared
% b = a / sqrt( ep^2 + 1 ); % earth radius @ poles
% invf = a / (a-b); % inverse flattening 
