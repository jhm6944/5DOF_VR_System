function [LF01D, LF02D, LF03D, LF04D] = LFDConstruction()
%*******************************************************
% Visesh Chari (2020). High accuracy optical flow
% (https://www.mathworks.com/matlabcentral/fileexchange/17500-high-accuracy-optical-flow),
% MATLAB Central File Exchange. Retrieved August 13, 2020.
%
% Depth sampled at 1/10 is used.
%*******************************************************
    global Params;
    
    [TX, TY] = meshgrid(1:Params.Width, 1:Params.Height);
    TX = TX(:);
    TY = TY(:);
    
    LF01D = zeros(Params.N_D * Params.Height * Params.Width, 1);
    LF02D = zeros(Params.N_D * Params.Height * Params.Width, 1);
    LF03D = zeros(Params.N_D * Params.Height * Params.Width, 1);
    LF04D = zeros(Params.N_D * Params.Height * Params.Width, 1);

    for n = 1:Params.N_D
        file = sprintf('LF01_D/%04d.mat', n);
        load(file);
        LF01D((n - 1) * (Params.Height * Params.Width) + (TX - 1) * (Params.Height) + (TY - 1) + 1) = dist(:)/2 + Params.W;
        
        file = sprintf('LF02_D/%04d.mat', n);
        load(file);
        LF02D((n - 1) * (Params.Height * Params.Width) + (TX - 1) * (Params.Height) + (TY - 1) + 1) = dist(:)/2 + Params.W;

        file = sprintf('LF03_D/%04d.mat', n);
        load(file);
        LF03D((n - 1) * (Params.Height * Params.Width) + (TX - 1) * (Params.Height) + (TY - 1) + 1) = dist(:)/2 + Params.W;
        
        file = sprintf('LF04_D/%04d.mat', n);
        load(file);
        LF04D((n - 1) * (Params.Height * Params.Width) + (TX - 1) * (Params.Height) + (TY - 1) + 1) = dist(:)/2 + Params.W;
        
        fprintf('%03d/%03d\n', n, Params.N_D);
    end
end