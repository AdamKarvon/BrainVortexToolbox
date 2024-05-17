function [relevantArea] = contour_area(xPoint,yPoint,Data,spiral_type)
    % if strcmp(spiral_type,'modulated_pos ') || strcmp(spiral_type,'generated_pos ')       
    %     level = 0.5; 
    % elseif strcmp(spiral_type,'modulated_neg ') || strcmp(spiral_type,'generated_neg ') 
    %     level = -0.5;
    % end
    level = 0;

    contourMatrix = contourc(Data, [level level]);  
    idx = 1;    
    relevantArea = NaN; % Initialize as NaN, will be updated if the point is within a contour

    while idx < size(contourMatrix, 2)
        nPoints = contourMatrix(2, idx);
        contourX = contourMatrix(1, idx+1:idx+nPoints);
        contourY = contourMatrix(2, idx+1:idx+nPoints);

        % Check if the point is inside this contour
        isInside = inpolygon(xPoint, yPoint, contourX, contourY);

        if isInside
            % Calculate area using the shoelace formula
            area = 0;
            for i = 1:nPoints-1
                area = area + contourX(i)*contourY(i+1) - contourY(i)*contourX(i+1);
            end
            area = area + contourX(nPoints)*contourY(1) - contourY(nPoints)*contourX(1);
            relevantArea = abs(area) / 2; % This is the area of the contour surrounding the point

            % Optionally, break the loop if you're only interested in the first relevant contour
            break;
        end

        idx = idx + nPoints + 1; % Move to the next contour in the matrix
    end
end