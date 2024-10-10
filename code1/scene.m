function [bs_pos,ris_pos,CV_points,V2V_points,V2V_dist]=scene(r,M)
    %% Communication scenarios

    h=50;
    res=100;
%     parameter_settings;
    %% Mapping of base station locations
    bs_pos=[0,0,25];
    scatter3(bs_pos(1),bs_pos(2),bs_pos(3),'ro');
    hold on;

    %% Mapping of RIS locations
    ris_pos=[60,0,30];
    scatter3(ris_pos(1),ris_pos(2),ris_pos(3),'go');
    hold on;

    %% Randomly generated vehicle locations
    num_cars=50;
    numPoints = poissrnd(num_cars);
    cnt=10+M;

    % Generate random points in a circular area with a radius of 500
    points = zeros(numPoints, 3);
    count = 0;
    while count < numPoints
        randX = (rand * 2 - 1) * 500;
        randY = (rand * 2 - 1) * 500;
        if hypot(randX, randY) <= 500
            count = count + 1;
            points(count, :) = [randX, randY,0];
        end
    end

    points = points(1:cnt, :);
    %% Mapping vehicle locations
    % scatter3(x_car,y_car,z_car,'bo');
    % hold on;

    %% Mapping the space scene
    theta = linspace(0, 2*pi, res); % Construct angular ranges

    x=r * cos(theta);
    y=r * sin(theta);
    z=zeros(size(theta));

    z_end=ones(size(theta))*h;

    surf([x;x],[y;y],[z;z_end],'FaceAlpha',0.3);

    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('空间场景')
    axis equal;

    %% Generate V2I and V2V
    % 从中随机选取五个点
    selectedIndices = randperm(cnt, 5);
    selectedPoints = points(selectedIndices, :);

    % Initialize the list of point pairs
    pairs = [];

    % Find the nearest point to each selected point
    for i = 1:5
        selectedPoint = selectedPoints(i, :);

        % Initialize minimum distance and nearest point indexes
        minDist = inf;
        nearestIndex = 0;

        for j = 1:cnt
            % Skip if it is a selected point or already in another point pair
            if ismember(j, selectedIndices) || ismember(j, pairs)
                continue;
            end

            % Calculate the distance between the selected point and the current point
            dist = norm(selectedPoint - points(j, :));

            % Update minimum distance and nearest point index if current distance is smaller
            if dist < minDist
                minDist = dist;
                nearestIndex = j;
            end
        end

        % Add selected and nearest points to the list of point pairs
        pair = [selectedIndices(i), nearestIndex];
        pairs = [pairs; pair];
    end

    for i=1:5
        startCoord=points(pairs(i,1),1:3);
        endCoord=points(pairs(i,2),1:3);
        plot3([startCoord(1),endCoord(1)],[startCoord(2),endCoord(2)],[startCoord(3),endCoord(3)],'r');
        scatter3(points(pairs(i,1), 1), points(pairs(i,1), 2),points(pairs(i,1),3), 'filled','b');
        scatter3(points(pairs(i,2), 1), points(pairs(i,2), 2),points(pairs(i,2),3), 'filled','b');
    end
   
    remainingIndices = setdiff(1:cnt, pairs);
    scatter3(points(remainingIndices, 1), points(remainingIndices, 2),points(remainingIndices,3), 'filled','g');
    %% Save randomly distributed vehicle coordinates to a file
    CV_points=[points(remainingIndices,1:3)];
    V2V_points=[points(pairs(:,1),1:3),points(pairs(:,2),1:3)];
    
    V2V_dist=zeros(1,5);
    for i=1:5
        V2V_dist(i)=pdist2(points(pairs(i,1),:),points(pairs(i,2),:));
    end
    save('CV_location.mat','CV_points');
    save('V2V_location.mat','V2V_points');
    save('V2V_dist.mat','V2V_dist');
end