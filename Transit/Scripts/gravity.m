% Problem 5.6 from Meyer and Miller
% Adapted for TRB Competition

% Gravity Model
clear

%------- Trip Generation Inputs -----------%

% Station [X1,X2,X3,X4,AM,Y1,Y2,Y3,Y4]
% Productions Vector
Prod = [353,325,308,240,693,407,812,684,327];
% Attractions Vector;
Att = [341,334,289,240,692,423,787,714,329];

% travel time Matrix;
time = [Inf,1,3,4,2,5,4,4,5;...
        1,Inf,2,3,1,4,3,3,4;...
        3,2,Inf,1,1,4,3,3,4;...
        4,3,1,Inf,2,5,4,4,5;...
        2,1,1,2,Inf,2,1,1,2;...
        5,4,4,5,2,Inf,1,3,4;...
        4,3,3,4,1,1,Inf,2,3;...
        4,3,3,4,1,3,2,Inf,1;...
        5,4,4,5,2,4,3,1,Inf];

    % Assumptions:  Everyone travels (t[X1,X1]=Inf)
    %               Adjacent nodes have t=1
    %               transfer penalty    x=1
    

%------- Trip Distribution ---------------%
    
% b coefficient
b = 1.5;

% tolerance
TOLERANCE = 1;

% Convergence Conditions
difference = Inf;
Trips1 = zeros([9,9]);
Astar = Att; %initialize


% PROGRAM LOOP %
disp('Adjusted Gravity Model')


while difference > TOLERANCE
 
    % Calculate Trips from i to j
    for i = 1:9
        for j = 1:9
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