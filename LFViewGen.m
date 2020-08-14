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