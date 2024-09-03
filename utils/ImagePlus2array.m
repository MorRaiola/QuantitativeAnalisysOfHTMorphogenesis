function I = ImagePlus2array(imp)
% ImagePlus2array retrieves an arary of image data from ImagePlus object of
% ImageJ1. IJM.getDataset and IJM.getDatasetAs provided by
% net.imagej.matlab.ImageJMATLABCommands class of ImageJ-MATLAB offers the
% similar function. Unlike these methods that rely on ImageJ2 code,
% ImagePlus2array used ImageJ1 API only.
%
% Advantages over getDataset and getDatasetAs:
% ImagePlus2array can work without an image being shown in ImageJ. In other
% words, it can work in headless mode. Although the computation is slower
% than getDataset, getDataset may take longer when you're dealing with a
% large image, whitch takes a long time to show on screen.
%
% In addition, because ImagePlus2array does not rely on IJM
% (net.imagej.matlab.ImageJMATLABCommands Java object), which only exist in
% the base workspace, ImagePlus2array may be easier to be used in a MATLAB
% function.
%
% Disadvantages:
% ImagePlus2array relys on ImageJ1 and can be slower than IJM.getDataset.
% ImagePlus2array only supports up to 5 dimensions.
%
% SYNTAX
% I = ImagePlus2array(imp)
%
% INPUT ARGUMENTS
% imp         ij.ImagePlus Java object
%             ImageJ1 ImagePlus object. It does not have to be shown in screen.
%
%
% OUTPUT ARGUMENTS
% I           uint8 | uint16 
%             Multidimensional (2 to 5) array of pixel data in the
%             appropriate data type. Note that the first and second
%             dimensions are swapped from ImageJ's representation; X and Y
%             axes are represented by the first and second dimensions in
%             ImageJ, and by the second and first in MATLAB.
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 29-Sep-2018 10:20:48
%
% See also
% net.imagej.matlab.ImageJMATLABCommands.getDataset,
% net.imagej.matlab.ImageJMATLABCommands.getDatasetAs
% copytoImagePlus
% 



dim = double(imp.getDimensions);

switch imp.getBitDepth
    case 8
        I = zeros(dim(2),dim(1),dim(3),dim(4),dim(5),'uint8');
        
    case 16
        I = zeros(dim(2),dim(1),dim(3),dim(4),dim(5),'uint16');
        
    case 24 %RGB
        I = zeros(dim(2),dim(1),dim(3),'uint8');
        
        channels = ij.plugin.ChannelSplitter.split(imp);
        
        for c = 1:3
            for z = 1:dim(4)
                channels(c).setZ(z);
                for t = 1:dim(5)
                    channels(c).setT(t);
                    I(:,:,c,z,t) = local_getPixels(channels(c),dim);
                end
            end
        end

        return
    case 32
        error('not supported yet')
    case 64
        error('not supported yet')
    otherwise
        error('not supported yet')
end


if nnz(dim(3:5) == 1) == 2
    
    DIM = dim(find(dim(3:5) ~= 1)+2);
    for i = 1:DIM
        imp.setSlice(i);
        
        switch find(dim(3:5) ~= 1)
            case 1
                c = i;
                z = 1;
                t = 1;
            case 2
                c = 1;
                z = i;
                t = 1;
            case 3
                c = 1;
                z = 1;
                t = i;
        end
        
        P = local_getPixels(imp,dim);
        
        I(:,:,c,z,t) = P;
        
    end
    
    
elseif imp.isHyperStack
    
    for c = 1:dim(3) %TODO parfor? I cannot be used as is
        imp.setC(c);

        for z = 1:dim(4)
            imp.setZ(z);
            for t = 1:dim(5)    
                imp.setT(t);
                
                P = local_getPixels(imp,dim);
               
                I(:,:,c,z,t) = P;

            end
        end
    end
    
elseif imp.isComposite
    
    
    for c= 1:dim(3)
        imp.setSlice(c); % MUCH faster than setC(c) for big images

        for z = 1:dim(4)
            imp.setZ(z);
            for t = 1:dim(5)
                imp.setZ(t);
                
                P = local_getPixels(imp,dim);
                
                I(:,:,c,z,t) = P;
            end
        end
    end
    
else

    P = local_getPixels(imp,dim);
    
    I = P;
        
end



end

%--------------------------------------------------------------------------

function P = local_getPixels(imp,dim)

ip = imp.getChannelProcessor;

% f = ip.getFloatArray; % single, SLOW
% c = ip.getIntArray; % int32, SLOW

p = ip.getPixels; % vector, int16, very quick

switch class(p)
    case 'int32' %RGB or 32-bit
        
    case 'int16'
        ptc = typecast(p,'uint16'); %NOTE need to convert to uint16 from int16
    case 'int8'
        ptc = typecast(p,'uint8'); %NOTE need to convert to uint16 from int16
        
end

P = reshape(ptc,dim(1),dim(2))';

end