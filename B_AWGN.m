% Note: Needs Optimization toolbox
function B_AWGN

    % B_AWGN.m
    clc;
    close all;
    clear;

    N_bits = 64;% Defining that 64 bit solution is considered

    x_uniform = rand(100); % Uniformly distributed numbers
    x = x_uniform(1); % one instance considered for code functionality testing
    x_scaled = floor((2^64)*x); % internal scalling to arrive at binary representation without foregoing the uniformity of distribution
    x_bin = dec2bin(x_scaled, N_bits); % binary representation 
    
    x_sign = find_sign(x_bin, N_bits); % obtaining the sign (from the algorithm in the reference from the assignment)
    x_offset = find_offset(x_bin, N_bits); % obtaining the offset (from the algorithm in the reference from the assignment)

    x_lz = lzd(x_bin, N_bits); % obtaining the leading zero position (from the algorithm in the reference from the assignment)
    mul_inpt = find_mul_inpt(x_bin, x_lz, N_bits); % extract the leading zero position (from the algorithm in the reference from the assignment)
    seg_sel = find_seg_sel(x_bin, x_lz, N_bits); % extract the segment for selection (from the algorithm in the reference from the assignment)
    
    ROM_trans = find_ROM_trans(N_bits-x_lz-1); % Create the trans output (from the algorithm in the reference from the assignment)
    
% Could not complete the remaining steps for lack of Optimization toolbox 
    
    % ROM_trans = find_ROM_coef(N_bits-x_lz-1); % Create the trans output (from the algorithm in the reference from the assignment)
    % ICDF_x = sqrt(2)./erf(2*x-1);

    disp('AWGN done!');
end

%%
function x_lz = lzd(x_bin, N_bits)
% Obtaining the leading zero position (from the algorithm in the reference from the assignment)

    x_lz = -1;
    for i = 1:length(x_bin)
        if x_bin(i) == '0'
            x_lz = N_bits-i; % zero base
            break
        end
    end
end

%%
function mul_inpt = find_mul_inpt(x_bin, x_lz, N_bits)
% Extract the leading zero position (from the algorithm in the reference from the assignment)    

    mul_inpt = x_bin(N_bits-17:N_bits-3);
    if x_lz <= 17
        for i = 1:17-x_lz
            mul_inpt(i) = '0'; % Mask block with '0's
        end
    end
end

%%
function seg_sel = find_seg_sel(x_bin, x_lz, N_bits)
% Extract the segment for selection (from the algorithm in the reference from the assignment)    
    seg_sel = x_bin(1 : N_bits-x_lz);
    
end

%%
function x_sign = find_sign(x_bin, N_bits)
% Obtaining the sign (from the algorithm in the reference from the assignment)
    x_sign = 1;
    if x_bin(N_bits) == '1'
        x_sign = -1;
    end
end

%%
function x_offset = find_offset(x_bin, N_bits)
% Obtaining the offset (from the algorithm in the reference from the assignment)
    x_offset = 0;
    if x_bin(N_bits-1) == '1'
        x_offset = x_offset + 1;
    end
    if x_bin(N_bits-2) == '1'
        x_offset = x_offset + 2;
    end
end

%%
function ROM_trans = find_ROM_trans(lz_loc)
    ROM_trans = 2^(lz_loc+1)-2;
end

%%
function f = icdf(x)
    ICDF_x = sqrt(2)./erf(2*x-1);
end        

%%
function get_limits
% Copy/paste from ref: Hierarchical Segmentation Schemes for Function Evaluation
% from Figure 1
% Not working as minimax is not available.

    % Inputs: a, b. d. €, e-max. ulp
    % 0Urpl.f : "0
    a = 0;
    b = 1-2^(-64);
    ulp = 2^(-64);
    x1 = a; x2 = b; m = 1; done = 0; check_x2 = 0 ; prev_x2 = a;
    oscillating = 0;
    while (~done)
        error = minimax(f,d,xl,x2,ulp);
        if (error <= req_error)
            if (x2 == b)
                u(m) = x2;
                done = 1;
            else
                if (oscillating)
                    u(m) = x2;
                    prev_x2 = x2;
                    x1 = x2;
                    x2 = b;
                    m = m+l;
                    oscillating = 0;
                else
                    change_x2 = abs(x2-prev_x2)/2;
                    prev_x2 = x2;
                    if (change_x2 > ulp)
                        x2 = x2 + change_x2;
                    else    
                        x2 = x2 + u1p;
                    end
                end
            end
        else
            change_x2 = abs(x2-prev_x2)/2;
            prev_x2 = x2;
            if (change_x2 > ulp)
                x2 = x2 - change_x2;
            else
                x2 = x2 - ulp;
                if (check_x2 == x2)
                    oscillating = 1;
                else
                    check_x2 = x2;
                end
            end
        end
    end
end