function [LF01, LF02, LF03, LF04] = LFConstruction()
    global Params;
    
    [TX, TY] = meshgrid(1:Params.Width, 1:Params.Height);
    TX = TX(:);
    TY = TY(:);
    
    LF01 = zeros(Params.N * Params.Height * Params.Width * 3, 1, 'uint8');
    LF02 = zeros(Params.N * Params.Height * Params.Width * 3, 1, 'uint8');
    LF03 = zeros(Params.N * Params.Height * Params.Width * 3, 1, 'uint8');
    LF04 = zeros(Params.N * Params.Height * Params.Width * 3, 1, 'uint8');

    for n=1:Params.N
        file = sprintf('LF01/%04d.jpg', n);
        im = imread(file);
        img_C1 = im(:, :, 1);
        img_C2 = im(:, :, 2);
        img_C3 = im(:, :, 3);
        LF01((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 3) = img_C1(:);
        LF01((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 2) = img_C2(:);
        LF01((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 1) = img_C3(:);

        file = sprintf('LF02/%04d.jpg', n);
        im = imread(file);
        img_C1 = im(:, :, 1);
        img_C2 = im(:, :, 2);
        img_C3 = im(:, :, 3);
        LF02((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 3) = img_C1(:);
        LF02((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 2) = img_C2(:);
        LF02((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 1) = img_C3(:);

        file = sprintf('LF03/%04d.jpg', n);
        im = imread(file);
        img_C1 = im(:, :, 1);
        img_C2 = im(:, :, 2);
        img_C3 = im(:, :, 3);
        LF03((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 3) = img_C1(:);
        LF03((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 2) = img_C2(:);
        LF03((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 1) = img_C3(:);

        file = sprintf('LF04/%04d.jpg', n);
        im = imread(file);
        img_C1 = im(:, :, 1);
        img_C2 = im(:, :, 2);
        img_C3 = im(:, :, 3);
        LF04((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 3) = img_C1(:);
        LF04((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 2) = img_C2(:);
        LF04((n - 1) * (Params.Height * Params.Width * 3) + (TX - 1) * (Params.Height * 3) + (TY - 1) * 3 + 1) = img_C3(:);
        
        fprintf('%03d/%03d\n', n, Params.N);
    end
end