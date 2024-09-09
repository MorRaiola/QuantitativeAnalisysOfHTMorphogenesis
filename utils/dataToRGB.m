function RGB = dataToRGB(data, indexValue)
    % dataToRGB - Generate an RGB matrix transitioning from blue to white to red,
    %             with black color for data values equal to zero.
    % 
    % Syntax: RGB = dataToRGB(data, indexValue)
    % 
    % Inputs:
    %    data - A vector (M x 1) of data values.
    %    indexValue - The value at which the color should be white.
    % 
    % Outputs:
    %    RGB - An M x 3 matrix of RGB colors.
    
    % Ensure data is a column vector
    if isrow(data)
        data = data';
    end

    % Define colors
    topColor = [1 0 0]; % Red
    indexColor = [1 1 1]; % White
    bottomColor = [0 0 1]; % Blue
    blackColor = [0 0 0]; % Black

    % Calculate the range of the data
    largest = max(data);
    smallest = min(data(data~=0));
    
    % Initialize the RGB matrix
    M = length(data);
    RGB = zeros(M, 3);

    % Compute the index for the color transition
    index = (indexValue - smallest) / (largest - smallest);

    % Create gradients
    % Gradients from blue to white and from white to red
    for i = 1:M
        value = data(i);
        if value == 0
            % Set the color to black if the data value is 0
            RGB(i, :) = blackColor;
        elseif value <= indexValue
            % Blue to white
            ratio = (value - smallest) / (indexValue - smallest);
            RGB(i, :) = (1 - ratio) * bottomColor + ratio * indexColor;
        else
            % White to red
            ratio = (value - indexValue) / (largest - indexValue);
            RGB(i, :) = (1 - ratio) * indexColor + ratio * topColor;
        end
    end
end
