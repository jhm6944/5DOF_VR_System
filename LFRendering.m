clc; clear all;

global Params;

Params.N = 150;
Params.N_D = Params.N/10;
Params.S = Params.N/3;
Params.Height = 1024;
Params.Width  = 2048;

Params.K = 75;
Params.W = 10;
mode = 'Hybrid';
%% Construct 3D LFs
fprintf('Construct LFs\n');
[LF01, LF02, LF03, LF04] = LFConstruction();
%%
fprintf('Load depth data\n');
[LF01D, LF02D, LF03D, LF04D] = LFDConstruction();

%% View generation
curr_posx = 0; %-20 ~ 20
curr_posy = 20; %-20 ~ 20

for curr_posx = -20:10:20
    CurrPos = [curr_posx, curr_posy];
    % Allocated angles for each LF
    [Alloc01, Alloc02, Alloc03, Alloc04, Pos01, Pos02, Pos03, Pos04] = getAngles(CurrPos, mode);

    % View generation
    % View generation of LF01 is divided into the first and second half.
    [Img01A, Over01A] = LFViewGen(LF01, LF01D, Pos01, Alloc01, 1, 0);
    [Img02, Over02] = LFViewGen(LF02, LF02D, Pos02, Alloc02, 2, 2);
    [Img03, Over03] = LFViewGen(LF03, LF03D, Pos03, Alloc03, 3, 2);
    [Img04, Over04] = LFViewGen(LF04, LF04D, Pos04, Alloc04, 4, 2);
    [Img01B, Over01B] = LFViewGen(LF01, LF01D, Pos01, Alloc01, 1, 1);

    % Blending overlapping views
    % Blending between LF01-LF02
    Img01A_out = Img01A(:, 1:size(Img01A, 2) - (Over01A(1, 2) + Over02(1, 1)), :);
    blend_0 = Img01A(:, size(Img01A, 2) - (Over01A(1, 2) + Over02(1, 1)) + 1:size(Img01A, 2), :);
    blend_1 = Img02(:, 1:(Over01A(1, 2) + Over02(1, 1)), :);
    Img01_02 = blend(blend_0, blend_1);

    % Blending between LF02-LF03
    Img02_out = Img02(:, (Over01A(1, 2) + Over02(1, 1))+1:size(Img02, 2) - (Over02(1, 2) + Over03(1, 1)), :);
    blend_0 = Img02(:, size(Img02, 2) - (Over02(1, 2) + Over03(1, 1)) + 1:size(Img02, 2), :);
    blend_1 = Img03(:, 1:(Over02(1, 2) + Over03(1, 1)), :);
    Img02_03 = blend(blend_0, blend_1);

    % Blending between LF03-LF04
    Img03_out = Img03(:, (Over02(1, 2) + Over03(1, 1))+1:size(Img03, 2) - (Over03(1, 2) + Over04(1, 1)), :);
    blend_0 = Img03(:, size(Img03, 2) - (Over03(1, 2) + Over04(1, 1)) + 1:size(Img03, 2), :);
    blend_1 = Img04(:, 1:(Over03(1, 2) + Over04(1, 1)), :);
    Img03_04 = blend(blend_0, blend_1);

    % Blending between LF04-LF01
    Img04_out = Img04(:, (Over03(1, 2) + Over04(1, 1))+1:size(Img04, 2) - (Over04(1, 2) + Over01B(1, 1)), :);
    blend_0 = Img04(:, size(Img04, 2) - (Over04(1, 2) + Over01B(1, 1)) + 1:size(Img04, 2), :);
    blend_1 = Img01B(:, 1:(Over04(1, 2) + Over01B(1, 1)), :);
    Img04_01 = blend(blend_0, blend_1);

    Img01B_out = Img01B(:, (Over04(1, 2) + Over01B(1, 1))+1:size(Img01B, 2), :);

    out = horzcat(Img01A_out, Img01_02, Img02_out, Img02_03, Img03_out, Img03_04, Img04_out, Img04_01, Img01B_out);
    out = imresize(out, [Params.Height, Params.Width]);
    imshow(out);
    fprintf('Pos : %d, %d\n', curr_posx, curr_posy);
    pause;
end

%% Functions
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

function [Alloc01, Alloc02, Alloc03, Alloc04, Pos01, Pos02, Pos03, Pos04] = getAngles(CurrPos, mode)
    global Params;
    
    posx = CurrPos(1, 1);
    posy = CurrPos(1, 2);
    
    Pos01 = [posx, posy];
    Pos02 = [-1 * posy, posx];
    Pos03 = [-1 * posx, -1 * posy];
    Pos04 = [posy, -1 * posx];
    
    if(strcmp(mode, 'RaySpace360'))
        Alloc01(1, 1) = atan2((-1.0 * Params.S / 2 - Pos01(1, 1)), (Params.S / 2 - Pos01(1, 2)));
        Alloc01(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos01(1, 1)), (Params.S / 2 - Pos01(1, 2)));
        
        Alloc02(1, 1) = atan2((-1.0 * Params.S / 2 - Pos02(1, 1)), (Params.S / 2 - Pos02(1, 2)));
        Alloc02(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos02(1, 1)), (Params.S / 2 - Pos02(1, 2)));
        
        Alloc03(1, 1) = atan2((-1.0 * Params.S / 2 - Pos03(1, 1)), (Params.S / 2 - Pos03(1, 2)));
        Alloc03(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos03(1, 1)), (Params.S / 2 - Pos03(1, 2)));
        
        Alloc04(1, 1) = atan2((-1.0 * Params.S / 2 - Pos04(1, 1)), (Params.S / 2 - Pos04(1, 2)));
        Alloc04(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos04(1, 1)), (Params.S / 2 - Pos04(1, 2)));
        
    elseif(strcmp(mode, 'LFU'))
        Alloc01 = [-1 * deg2rad(45), deg2rad(45)];
        Alloc02 = [-1 * deg2rad(45), deg2rad(45)];
        Alloc03 = [-1 * deg2rad(45), deg2rad(45)];
        Alloc04 = [-1 * deg2rad(45), deg2rad(45)];
        
    elseif(strcmp(mode, 'Hybrid'))
        Alloc04(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos04(1, 1)), (Params.S / 2 - Pos04(1, 2)));
        Alloc01(1, 1) = atan2((-1.0 * Params.S / 2 - Pos01(1, 1)), (Params.S / 2 - Pos01(1, 2)));

        if(abs(Alloc01(1, 1)) > deg2rad(Params.K))
            Alloc01(1, 1) = -1 * deg2rad(Params.K);
            Alloc04(1, 2) = deg2rad(90 - Params.K);
        elseif(abs(Alloc04(1, 2)) > deg2rad(Params.K))
            Alloc01(1, 1) = -1 * deg2rad(90 - Params.K);
            Alloc04(1, 2) = deg2rad(Params.K);
        end

        Alloc01(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos01(1, 1)), (Params.S / 2 - Pos01(1, 2)));
        Alloc02(1, 1) = atan2((-1.0 * Params.S / 2 - Pos02(1, 1)), (Params.S / 2 - Pos02(1, 2)));
        
        if(abs(Alloc02(1, 1)) > deg2rad(Params.K))
            Alloc02(1, 1) = -1 * deg2rad(Params.K);
            Alloc01(1, 2) = deg2rad(90 - Params.K);
        elseif(abs(Alloc01(1, 2)) > deg2rad(Params.K))
            Alloc02(1, 1) = -1 * deg2rad(90 - Params.K);
            Alloc01(1, 2) = deg2rad(Params.K);
        end
        
        Alloc02(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos02(1, 1)), (Params.S / 2 - Pos02(1, 2)));
        Alloc03(1, 1) = atan2((-1.0 * Params.S / 2 - Pos03(1, 1)), (Params.S / 2 - Pos03(1, 2)));
        
        if(abs(Alloc03(1, 1)) > deg2rad(Params.K))
            Alloc03(1, 1) = -1 * deg2rad(Params.K);
            Alloc02(1, 2) = deg2rad(90 - Params.K);
        elseif(abs(Alloc02(1, 2)) > deg2rad(Params.K))
            Alloc03(1, 1) = -1 * deg2rad(90 - Params.K);
            Alloc02(1, 2) = deg2rad(Params.K);
        end
        
        Alloc03(1, 2) = atan2(( 1.0 * Params.S / 2 - Pos03(1, 1)), (Params.S / 2 - Pos03(1, 2)));
        Alloc04(1, 1) = atan2((-1.0 * Params.S / 2 - Pos04(1, 1)), (Params.S / 2 - Pos04(1, 2)));

        if(abs(Alloc04(1, 1)) > deg2rad(Params.K))
            Alloc04(1, 1) = -1 * deg2rad(Params.K);
            Alloc03(1, 2) = deg2rad(90 - Params.K);
        elseif(abs(Alloc03(1, 2)) > deg2rad(Params.K))
            Alloc04(1, 1) = -1 * deg2rad(90 - Params.K);
            Alloc03(1, 2) = deg2rad(Params.K);
        end
    end
end

function [OutImg, Over] = LFViewGen(LF, LFD, Pos, Alloc, dir, LFflag)
    global Params;
    
    % Allocate overlapping area
    if(LFflag == 0)
        Alloc_s = 0;
        Alloc_e = Alloc(1, 2) + deg2rad(5);
    elseif(LFflag == 1)
        Alloc_s = Alloc(1, 1) - deg2rad(5);
        Alloc_e = 0;
    else
        Alloc_s = Alloc(1, 1) - deg2rad(5);
        Alloc_e = Alloc(1, 2) + deg2rad(5);
    end
    
    if(abs(Alloc_s) > deg2rad(Params.K)),         Alloc_s = -1 * deg2rad(Params.K);    end
    if(abs(Alloc_e) > deg2rad(Params.K)),         Alloc_e = deg2rad(Params.K);         end
    
    Over(1, 1) = floor((Alloc(1, 1) - Alloc_s) / deg2rad(360/Params.Width));
    Over(1, 2) = floor((Alloc_e - Alloc(1, 2)) / deg2rad(360/Params.Width));
    
    Alloc(1, 1) = Alloc_s;
    Alloc(1, 2) = Alloc_e;
    
    OutWidth = floor((Alloc(1, 2) - Alloc(1, 1)) / deg2rad(360/Params.Width));
    OutImg = zeros(Params.Height, OutWidth, 3, 'uint8');
    PelWidth = (1:OutWidth)';
    PelAngles = Alloc(1, 1) + deg2rad(360/Params.Width) .* (PelWidth - 1);
    
    LF_x = ((Params.S / 2) - Pos(1, 2)) .* tan(PelAngles) + Pos(1, 1);
%     dist = 0.02 * ((Params.S / 2) - Pos(1, 2));
    dist = 0.02 * sqrt((LF_x - Pos(1, 1)) .^ 2 + ((Params.S / 2) - Pos(1, 2)) .^ 2);
    dist = repmat(dist', Params.Height, 1);

    LF_u = PelAngles .* (180.0 / pi) * (1.0 / 180.0) * Params.Width / 2 + Params.Width / 2;
    LF_v = (1:Params.Height)';
    LF_v = repmat(LF_v, 1, OutWidth);
    
    LF_x = repmat(LF_x', size(LF_v, 1), 1);
    LF_x_0 = floor((LF_x + (Params.N / 2))/10);
    LF_x_1 = ceil((LF_x + (Params.N / 2))/10);
    LF_x_r = ((LF_x + (Params.N / 2))/10) - LF_x_0;
    
    LF_x_0(LF_x_0 < 1) = 1;
    LF_x_0(LF_x_0 > Params.N_D) = Params.N_D;
    
    LF_x_1(LF_x_1 < 1) = 1;
    LF_x_1(LF_x_1 > Params.N_D) = Params.N_D;
        
    LF_u = repmat(LF_u', size(LF_v, 1), 1);
    LF_u_0 = floor(LF_u);
    LF_u_1 = ceil(LF_u);
    LF_u_r = LF_u - LF_u_0;
    
    LF_u_0(LF_u_0 < 1) = 1;
    LF_u_0(LF_u_0 > Params.Width) = Params.Width;
    
    LF_u_1(LF_u_1 < 1) = 1;
    LF_u_1(LF_u_1 > Params.Width) = Params.Width;
    
    LF_d = interD_mat(LFD, LF_x_r, LF_x_0-1, LF_x_1-1, LF_u_r, LF_u_0-1, LF_u_1-1, LF_v);
    
    LF_x_0 = floor(LF_x + (Params.N / 2));
    LF_x_1 = ceil(LF_x + (Params.N / 2));
    LF_x_r = (LF_x + (Params.N / 2)) - LF_x_0;
    
    LF_x_0(LF_x_0 < 1) = 1;
    LF_x_0(LF_x_0 > Params.N) = Params.N;
    
    LF_x_1(LF_x_1 < 1) = 1;
    LF_x_1(LF_x_1 > Params.N) = Params.N;
    
    LF_pi = (LF_v-(Params.Height/2))/((Params.Height/2)) * pi/2;
    
    % This function is not exactly the same as that of the paper due to backward.
    LF_pi_corr = atan2((LF_d .* sin(LF_pi)), (LF_d .* cos(LF_pi) - dist));
    LF_v_corr  = (LF_pi_corr / (pi/2)) * (Params.Height/2) + (Params.Height/2);
        
    LF_u(LF_v_corr < 1) = LF_u(LF_v_corr < 1) + Params.Width/2;
    LF_v_corr(LF_v_corr < 1) = -1 * LF_v_corr(LF_v_corr < 1);
    LF_u(LF_v_corr > Params.Height) = LF_u(LF_v_corr > Params.Height) + Params.Width/2;
    LF_v_corr(LF_v_corr > Params.Height) = 2 * Params.Height - LF_v_corr(LF_v_corr > Params.Height);
    LF_v = LF_v_corr;
    
    LF_u(LF_u > Params.Width) = LF_u(LF_u > Params.Width) - Params.Width;
    
    if(dir == 2)
        LF_u = LF_u + Params.Width / 4;
    elseif(dir == 3)
        LF_u = LF_u + Params.Width / 2;
    elseif(dir == 4)
        LF_u = LF_u - Params.Width / 4;
    end
    LF_u(LF_u < 1) = LF_u(LF_u < 1) + Params.Width;
    LF_u(LF_u > Params.Width) = LF_u(LF_u > Params.Width) - Params.Width;

    LF_u_0 = floor(LF_u);
    LF_u_1 = ceil(LF_u);
    LF_u_r = LF_u - LF_u_0;
    
    LF_u_0(LF_u_0 < 1) = 1;
    LF_u_0(LF_u_0 > Params.Width) = Params.Width;
    
    LF_u_1(LF_u_1 < 1) = 1;
    LF_u_1(LF_u_1 > Params.Width) = Params.Width;
    
    LF_v_0 = floor(LF_v);
    LF_v_1 = ceil(LF_v);
    LF_v_r = LF_v - LF_v_0;
    
    LF_v_0(LF_v_0 < 1) = 1;
    LF_v_0(LF_v_0 > Params.Height) = Params.Height;
    
    LF_v_1(LF_v_1 < 1) = 1;
    LF_v_1(LF_v_1 > Params.Height) = Params.Height;
    
    OutImg(:, :, 3) = (inter8_mat(LF, LF_x_r, LF_x_0-1, LF_x_1-1, LF_u_r, LF_u_0-1, LF_u_1-1, LF_v_r, LF_v_0-1, LF_v_1-1, 1));
    OutImg(:, :, 2) = (inter8_mat(LF, LF_x_r, LF_x_0-1, LF_x_1-1, LF_u_r, LF_u_0-1, LF_u_1-1, LF_v_r, LF_v_0-1, LF_v_1-1, 2));
    OutImg(:, :, 1) = (inter8_mat(LF, LF_x_r, LF_x_0-1, LF_x_1-1, LF_u_r, LF_u_0-1, LF_u_1-1, LF_v_r, LF_v_0-1, LF_v_1-1, 3));
end

function IMAGE_MAT = inter8_mat(LF, P_r, P_1, P_2, U_r, U_1, U_2, H_r, H_1, H_2, c)
    global Params;

    P_1(P_r == 1) = 0;    P_2(P_r == 0) = 0;
    U_1(U_r == 1) = 0;    U_2(U_r == 0) = 0;
    H_1(H_r == 1) = 0;    H_2(H_r == 0) = 0;
    
    if(c == 1)
        IMAGE_MAT = ((1.0 - P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_1 * 3 + 1)) + ...
                                                    H_r .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_2 * 3 + 1))) + ...
                         ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_1 * 3 + 1)) + ...
                                                    H_r .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_2 * 3 + 1)))))) + ...
                        ((P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_1 * 3 + 1)) + ...
                                       H_r .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_2 * 3 + 1))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_1 * 3 + 1)) + ...
                                       H_r .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_2 * 3 + 1))))));
    elseif(c == 2)
        IMAGE_MAT = ((1.0 - P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_2 * 3 + 2))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_2 * 3 + 2)))))) + ...
                        ((P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_2 * 3 + 2))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_1 * 3 + 2)) + ...
                                       H_r .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_2 * 3 + 2))))));
    elseif(c == 3)
        IMAGE_MAT = ((1.0 - P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_2 * 3 + 3))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_1) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_2 * 3 + 3)))))) + ...
                        ((P_r) .* ...
                ((1.0 - U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_1 * (Params.Height * 3) + H_2 * 3 + 3))) + ...
                      ((U_r) .* ((1.0 - H_r) .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_1 * 3 + 3)) + ...
                                       H_r .* im2double(LF((P_2) * (Params.Height * Params.Width * 3) + U_2 * (Params.Height * 3) + H_2 * 3 + 3))))));
    end
    IMAGE_MAT = im2uint8(IMAGE_MAT);
end

function Depth_MAT = interD_mat(LFD, P_r, P_1, P_2, U_r, U_1, U_2, H)
    global Params;

    P_1(P_r == 1) = 0;    P_2(P_r == 0) = 0;
    U_1(U_r == 1) = 0;    U_2(U_r == 0) = 0;
    
    Depth_MAT = ((1.0 - P_r) .* ...
                ((1.0 - U_r) .* (im2double(LFD((P_1) * (Params.Height * Params.Width) + U_1 * (Params.Height) + H + 1))) + ...
                      ((U_r) .* (im2double(LFD((P_1) * (Params.Height * Params.Width) + U_2 * (Params.Height) + H + 1)))))) + ...
                      ((P_r) .* ...
                ((1.0 - U_r) .* (im2double(LFD((P_2) * (Params.Height * Params.Width) + U_1 * (Params.Height) + H + 1))) + ...
                      ((U_r) .* (im2double(LFD((P_2) * (Params.Height * Params.Width) + U_2 * (Params.Height) + H + 1))))));
end

function blended_img = blend(im0, im1)
height = size(im0, 1);
width = size(im0, 2);

map0 = repmat((width-1:-1:0)/(width-1), [height, 1]);
map1 = repmat((0:width-1)/(width-1), [height, 1]);

blended_img = im2uint8(map0 .* im2double(im0) + map1 .* im2double(im1));
end
