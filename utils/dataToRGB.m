function RGB = dataToRGB(data, indexValue, typeIdx,suffix)
    % dataToRGB - Generate an RGB matrix transitioning from blue to white to red,
    %             with black color for data values equal to zero.
    % 
    % Syntax: RGB = dataToRGB(data, indexValue, typeIdx)
    % 
    % Inputs:
    %    data - A vector (M x 1) of data values.
    %    indexValue - The value at which the color should be white.
    %    typeIdx - The type of colormap:
    %              1 = Blue to white to red (normal gradient)
    %              2 = Blue to red with white at indexValue (alternative gradient)
    % 
    % Outputs:
    %    RGB - An M x 3 matrix of RGB colors.
    
    % Ensure data is a column vector
    if isrow(data)
        data = data';
    end

    % Define colors
    topColor = [1 0 0];      % Red
    indexColor = [1 1 1];    % White
    bottomColor = [0 0 1];   % Blue
    blackColor = [1 1 1];    % Black

    % Calculate the range of the data
    largest = max(data);                      % Find the largest value in the data
    smallest = min(data(data ~= 0));          % Find the smallest non-zero value in the data
    
    % Initialize the RGB matrix
    M = length(data);
    RGB = zeros(M, 3);

    % Color assignment logic based on typeIdx
    if typeIdx == 1 || typeIdx == 2 && strcmp(suffix ,'std')
        % Type 1: Blue to white to red gradient
        for i = 1:M
            value = data(i);
            if value == 0
                % Set the color to black if the data value is 0
                RGB(i, :) = blackColor;
            elseif value <= indexValue
                % Transition from blue to white
                ratio = (value - smallest) / (indexValue - smallest);
                RGB(i, :) = (1 - ratio) * bottomColor + ratio * indexColor;
            else
                % Transition from white to red
                ratio = (value - indexValue) / (largest - indexValue);
                RGB(i, :) = (1 - ratio) * indexColor + ratio * topColor;
            end
        end
    elseif typeIdx == 2 && strcmp(suffix ,'mean')
        % Type 2: White at indexValue and blue-to-red for values > indexValue
        for i = 1:M
            value = data(i);
            if value == 0
                % Set the color to black if the data value is 0
                RGB(i, :) = blackColor;
            elseif value == indexValue
                % Set the color to white for the indexValue
                RGB(i, :) = indexColor;
            else
                % Transition from blue to red directly for values > indexValue
                ratio = (value - indexValue) / (largest - indexValue);
                RGB(i, :) = (1 - ratio) * bottomColor + ratio * topColor;
            end
        end
    else
        error('Invalid typeIdx. Use 1 for blue-to-white-to-red or 2 for white at indexValue and blue-to-red for higher values.');
    end
end

