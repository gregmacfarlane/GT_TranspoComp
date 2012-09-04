% Problem 5.6 from Meyer and Miller
% Gravity Model
clear

%------- Trip Generation Inputs -----------%

% Productions Vector
Prod = [100, 200, 100];

% Attractions Vector;
Att = [200, 50, 150];

% travel time Matrix;
time = [2, 5, 4;...
        5, 2, 3;...
        4, 3, 2];


%------- Trip Distribution ---------------%
    
% b coefficient
b = 1.5;

% tolerance
TOLERANCE = 1;

% Convergence Conditions
difference = Inf;
Trips1 = zeros([3,3]);
Astar = Att; %initialize


% PROGRAM LOOP %
disp('Adjusted Gravity Model')


while difference > TOLERANCE
 
    % Calculate Trips from i to j
    for i = 1:3
        for j = 1:3
            Trips(i,j) = ( Prod(i) * Astar(j)*( time(i,j) )^(-b))/...
                    ( sum(Astar .* ( time (i,:) ).^-b ));
        end
    end

% Adjust Attractions
Astar1 = Astar .* (Att ./ sum( Trips ));
Astar = Astar1;

% Calculate Convergence condition
difference = sum( sum( abs( Trips1 - Trips )));
Trips1 = Trips;

disp(difference);
end

%Output
disp('Final Trip Matrix')
disp(Trips)