% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       CUAUV Multipath Simulator
%                           by: Michael Nguyen
% Modelling the pool's geometry, we can calulate the multipath
% contributions to the overall signal.
%
% Essentially models loss of power through water, multipathing
% and models reflections as additional "point sources"
% assumptions: single reflections only contribute to the multipath problem
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

res = 0.1; % e.g. each point in the field refers to 10 cm difference
speed = 1497; %m/s under freshwater at 25 degrees Celsius

%Locations of the transmitter/reciever  (X, Y, Z)
transmit = [10; 10; 10];
receiver = [600; 450; 90];

%http://www.sal2000.com/ds/ds3/Acoustics/Wave%20Reflection.htm
%Using concepts similar to EM wave propagation, we can extend the 
%definition of reflection to sound waves:
%R = (z2 - z1)/(z2 + z1)
z1 = 1480000; %Freshwater
z2 = 8000000; %concerete

R = (z2-z1)/(z2+z1);

%transdec dimensions
majX = uint16(300*0.3048)/(2*res);
minY = uint16(200*0.3048)/res;
depth = uint16(38*0.3048)/res;

%define surfaces of reflection
field = zeros(majX,minY); %array of distances
link = zeros(majX,minY); %array of link budgets
zfield = 0;

%Absorption of sound in water
%http://resource.npl.co.uk/acoustics/techguides/seaabsorption/
alpha = 1.26; % @ 15kHz, 25 degrees C, 38 ft depth

%populate the field with multipath distances
dim = size(field);
for i = 1:dim(1)
    for j = 1:dim(2)
        temp = [j; i; zfield];
        dtrans = sqrt( (temp - transmit)'*(temp - transmit) );
        drecei = sqrt( (temp - receiver)'*(temp - receiver) );
        dist = res*(dtrans + drecei);
        
        %note that spherical rad attenuation and the reflection have
        %negative coefficients because we wish for all the terms within the
        %transmission loss budget to have positive values
        link(i,j) = -1*( -20*log10(1/(4*pi*dist)) + alpha*dist - 20*log10(R));
        link(i,j) = 10^(link(i,j)/10);
        field(i,j) = dist;
    end
end

LOSdist = sqrt( (receiver - transmit)'*(receiver - transmit) )*res;
LOSlink = 10^(-1*( -20*log10(1/(4*pi*LOSdist)) + alpha*LOSdist)/10);

link = link/LOSlink;

timeOfArrival = field/speed; %time of arrival of the signal
invsquare = 1./(field.^2); %Power loss is proportional

%visualization of the multipaths
figure
hold on
imagesc(link);
points = [transmit, receiver];
scatter3(points(1,:)', points(2,:)', points(3,:)',5,'black');
axis([0 minY 0 majX 0 depth]);

c = colorbar;

xlabel('X');
ylabel('Y');
ylabel(c,'LOS to Reflection signal ratio');
zlabel('Z');
hold off

%Signal Calculations
f = 150; %15 kHz
t = 1:100000;
timeres = 1e-6;
t = t*timeres; %seconds
y = sin(2*pi*f*t);

%figure
%plot(t,y);

%max(max(timeOfArrival))
