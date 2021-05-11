nodes = input('Enter no. of nodes');
t = 1.5;
Tinfi = 523;
qdot = 10^8;
l = 0.01;
deltax = l/(nodes-1);
T = linspace(523,523,nodes);
T = T';
alpha = 5*10^(-6);
h = 1100;
k = 30;
deltat = (fix(((k*deltax^2)/(2*alpha*(k + (h*deltax))))*10))/10;
time = fix(t/deltat) + 1;
if (fix(t/deltat) ~= (t/deltat))
    time = time + 1;
end
lastdeltat = t - ((time-2)*deltat);
Tf = [time nodes];
for a = 1:nodes
    Tf(1, a) = T(a);
end
for j = 2:(time)
    if (j ~= time)
        % 1st node
        Tf(j, 1) = ((alpha*deltat/(deltax^2))*(2*Tf((j-1), 2)-2*Tf((j-1), 1))) + Tf((j-1), 1) + qdot*deltat*alpha/k;
        % 2 to nodes - 1
        for i = 2:(nodes-1)
            Tf(j, i) = ((alpha*deltat/(deltax^2))*(Tf((j-1), (i-1))+Tf((j-1), (i+1))-2*Tf((j-1), i))) + Tf((j-1), i) + qdot*deltat*alpha/k;
        end
        % last node
        Tf(j, nodes) = ((2*alpha*deltat/(deltax^2))*(Tf((j-1), (i-1)) + (h*deltax*Tinfi/k) + ((qdot*deltax^2)/(2*k)))) + (1-(2*alpha*deltat/deltax^2) - (2*h*alpha*deltat/(k*deltax)))*Tf((j-1), i);
    else 
        % 1st node
        Tf(j, 1) = ((alpha*lastdeltat/(deltax^2))*(2*Tf((j-1), 2)-2*Tf((j-1), 1))) + Tf((j-1), 1) + qdot*lastdeltat*alpha/k;
        % 2 to nodes - 1
        for i = 2:(nodes-1)
            Tf(j, i) = ((alpha*lastdeltat/(deltax^2))*(Tf((j-1), (i-1))+Tf((j-1), (i+1))-2*Tf((j-1), i))) + Tf((j-1), i) + qdot*lastdeltat*alpha/k;
        end
        % last node
        Tf(j, nodes) = ((2*alpha*lastdeltat/(deltax^2))*(Tf((j-1), (i-1)) + (h*deltax*Tinfi/k) + ((qdot*deltax^2)/(2*k)))) + (1-(2*alpha*lastdeltat/deltax^2) - (2*h*alpha*lastdeltat/(k*deltax)))*Tf((j-1), i);
    end
end

Tf

dist = 0:1:nodes-1;
dist = deltax*dist;
plot(dist, Tf(time, :), 'Color', 'r');
hold on;
plot(-dist, Tf(time, :), 'Color', 'r');
xlabel('Distance in m');
label = ylabel('Temperature in K');
label.Position(1) = -0.1; % change horizontal position of ylabel
title('Variation of Temperature with Length of the element at time 1.5s')
hold off;
figure;

timespan = 0:1:time-1;
timespan = deltat*timespan;
timespan(1, time) = lastdeltat + timespan(1, time-1);
plot(timespan, Tf(:, 1), 'Color', 'g', 'DisplayName', 'Temperature of Mid-plane node');
hold on;
plot(timespan, Tf(:, nodes), 'Color', 'b', 'DisplayName', 'Temperature of Surface node');
xlabel('Time in s');
label_h = ylabel('Temperature in K');
label_h.Position(1) = -0.15; % change horizontal position of ylabel
hold off;
title('Variation of Temperature with Time of Mid-plane node and surface node')
