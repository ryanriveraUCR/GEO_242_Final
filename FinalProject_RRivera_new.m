%clear all;

%Final Project GEO242 - Ryan Rivera

%Estimating the fault structure of the Maacama Fault using seismicity
%catalog, first using the equation of a plane and then a quadratic surface.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in the double difference seismicity catalog downloaded from https://www.ldeo.columbia.edu/~felixw/NCAeqDD/
doubdiff_cata = load('NCAeqDD.v201112.1_new.txt');
%assign longitudes from column 8 of catalog
data_lon_dd = doubdiff_cata(:,8);
%assign latitudes from column 7 of catalog
data_lat_dd = doubdiff_cata(:,7);
%assign depths from column 9 of catalog
data_depth_dd = doubdiff_cata(:,9);

%now we need to create a polygon around the Maacama Fault in Google Earth
%and use inpolygon to select the earthquakes within the polygon
%copy the polygon from Google Earth to text notepad to see the coordinates

%the y coordinates for the Google Earth polygon - longitudes
yv_m = [-123.0366857199684,-122.7536801368026,-122.6757295355889,-122.9737312594132,-123.0366857199684];
%the x coordinates for the Google Earth polygon - latitudes
xv_m = [38.85986929693475,38.57243282723844,38.61121501445241,38.88008279030988,38.85986929693475];

%using the function inpolygon to return earthquake 
[in,on] = inpolygon(data_lat_dd,data_lon_dd,xv_m,yv_m);

%initiate figure
figure(1)
%plot the earthquakes within the polygon with red crosses
plot(data_lon_dd(in),data_lat_dd(in),'r+')
hold on
%plot the earthquakes outside the polygon with blue circles
plot(data_lon_dd(~in),data_lat_dd(~in),'bo')
%plot the Google Earth polygon
plot(yv_m,xv_m)
hold off

%isolate the earthquake location data within the polygon
%latitudes of earthquakes within polygon
poly_lat = [data_lat_dd(in)];
%longitudes of earthquakes within polygon
poly_lon = [data_lon_dd(in)];
%depths of earthquakes within polygon
poly_depth = [data_depth_dd(in)];

%using deg2utm to convert lat/lon to UTM coordinates
[poly_lat_utm,poly_lon_utm,utmzone] = deg2utm(poly_lat,poly_lon);
%Using the equation of a plane: z = ax + by + c
%where x are the longitudes converted to UTM
%where y are the latitudes converted to UTM
%where z is the depths in km
%set  matrices to preform least squares problem [x]=[A]\[b]
%where [x] and [b] are vectors in order to solve for a, b, and c

%making the [A] matrix with [lon,lat,1s]
A_mac = [poly_lat_utm, poly_depth, ones(length(poly_depth),1)];
%making [b] vector with earthquake depths
b_mac = [poly_lon_utm];
%backslash operator to solve least squares problem to solve for a,b,c
x = A_mac\b_mac;

%creating a mesh grid using poly_lon_utm and poly_lat_utm
[X,Y] = meshgrid(linspace(min(poly_lat_utm),max(poly_lat_utm)), linspace(min(poly_depth),max(poly_depth)));
%solving for the best fitting plane Z
Z = x(1)*X + x(2)*Y + x(3)*ones(size(X));

%initiate figure
figure(2)
%plotting the 3D data
plot3(poly_lat_utm, poly_depth, poly_lon_utm,'.')
hold on
%plotting the meshgrid
meshc(X, Y, Z)
hold off
grid on
%turn on colorbar
colorbar
%axis labels
xlabel('longitude (UTM)'); ylabel('depth (km)'); zlabel('latitude (UTM)');
%title
title('Best Fitting Plane for the Maacama Fault near The Geysers')
text(-20, 50, 450, sprintf('Z= %.3f\\cdotX %+.3f\\cdotY %+3.0f', x))

%solving for the strike and dip of the fault plane
%to solve for strike, we set z=0 in the equation z=ax+by+c, where a,b,c are
%the constants we solved for, and we get the z,x, and y from any of our
%earthquakes
strike = poly_lon_utm(1) - x(1)*poly_lat_utm(1) - x(3);