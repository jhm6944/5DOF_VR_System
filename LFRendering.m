clc; clear all;

global Params;

Params.N = 150;
Params.N_D = Params.N/10; % 1/10 sampled depth is used
Params.S = Params.N/3;
Params.Height = 1024;
Params.Width  = 2048;

Params.K = 75;
Params.W = 10;
mode = 'Hybrid';
%% Construct 3D LFs
fprintf('Construct LFs\n');
[LF01, LF02, LF03, LF04] = LFConstruction();
fprintf('Load depth data\n');
[LF01D, LF02D, LF03D, LF04D] = LFDConstruction();

%% View generation
curr_posx = 0; %-20 ~ 20
curr_posy = -20; %-20 ~ 20

for curr_posx = -20:5:20
    CurrPos = [curr_posx, curr_posy];
    % Allocated angles for each LF
    [Alloc01, Alloc02, Alloc03, Alloc04, Pos01, Pos02, Pos03, Pos04] = getAngles(CurrPos, mode);

    % View generation
    % View generation of LF01 is divided into the first and second half.
    % The middle view of the 360 degrees generated view is from LF03.
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
