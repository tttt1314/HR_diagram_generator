%% Clear
close all;                                                                 % Clear all image that open by MATLAB
clear all;                                                                 % Clear all pervious parameter store in MATLAB
%% Input region 1
name = 'M93';                                                              % Input of the mane of the cluster
b_sim = FITS.read2sim('m93_b.fits');                      % First input of fits file with B filter (with WCS)
v_sim = FITS.read2sim('m93_v.fits');                      % First input of fits file with V filter (with WCS)
b_image = abs(fitsread('m93_b.fits'));                    % Second input of fits file with B filter (with WCS)
v_image = abs(fitsread('m93_v.fits'));                    % Second input of fits file with V filter (with WCS)
star_list = readtable('m93_cata');                   % Input of star list from DS9
fitting_data_1 = readtable('theorem\387M.csv');                % Input of fitting curve data 1 (first column is B_V second column is V)
fitting_data_2 = readtable('theorem\387M.csv');              % Input of fitting curve data 2 (first column is B_V second column is V)
fit_1 = 'Theory curve from Isochrone and LF Generator';                    % Input of the name of fitting curve 1
fit_2 = 'Theory curve form CMD 3.6 input form';                            % Input of the name of fitting curve 2
distance = 881;                                                            % Input of the cluster (pc)
error_m = 5;                                                               % Input of error permittance for matching between B filter image and V filter image
error_b = 0.005;                                                           % Input of error permittance for matching between B filter image and star list from DS9 (degree)
error_v = 0.005;                                                           % Input of error permittance for matching between V filter image and star list from DS9 (degree)
radius = 6;                                                                % Input of integral radius (pixel)
inner_radius = 7;                                                          % Input of inner radius of correction ring (pixel)
outer_radius = 8;                                                          % Input of inner radius of correction ring (pixel)
G = fspecial('gaussian', 100, 1);                                          % Input of bluring setting
select_b = 3000;                                                           % Input of star selecting lower value bound of B filter image (0~65535)
select_v = 3000;                                                           % Input of star selecting lower value bound of V filter image (0~65535)
choose = 1;                                                                % Input of which fitting line to use (0: none  1:fitting line 1  2:fitting line 2  3:all)
%% Processing of B filter image
b_image_R = b_image-select_b;                                              % Create a star region map matrix
b_image_R(b_image_R < 0) = 0;                                              % Change the value which is smaller than star selecting lower value bound to 0
b_image_R(b_image_R > 0) = 1;                                              % Change the value which is bigger than star selecting lower value bound to 1
figure();                                                                  % Create Figure 1
imagesc(b_image_R);                                                        % Show the star region in the B filter image
title('Figure 1 : Star Region of B Filter Image');
b_image_B = b_image;                                                       % Create a star brightness local maximum finding matrix
b_image_B = imfilter(b_image_B, G);                                        % Blur the star brightness local maximum finding matrix
b_image_B = b_image_B.*b_image_R;                                          % Leave the star only
b_image_B = imregionalmax(b_image_B);                                      % Find the brightness local maximum (star) postions
fprintf('found %d stars in B layer!\n', sum(sum(b_image_B)));              % Show how many star find in B filter image
%set(gcf,'unit','centimeter','position',[0 0 10 10])
drawnow
%% Processing of V filter image
v_image_R = v_image-select_v;                                               % Create a star region map matrix
v_image_R(v_image_R < 0) = 0;                                              % Change the value which is smaller than star selecting lower value bound to 0
v_image_R(v_image_R > 0) = 1;                                              % Change the value which is bigger than star selecting lower value bound to 1
figure();                                                                  % Create Figure 2
imagesc(v_image_R);                                                        % Show the star region in the B filter image
title('Figure 2 : Star Region of V Filter Image');
v_image_B = v_image;                                                       % Create a star brightness local maximum finding matrix
v_image_B = imfilter(v_image_B, G);                                        % Blur the star brightness local maximum finding matrix
v_image_B = v_image_B.*v_image_R;                                          % Leave the star only
v_image_B = imregionalmax(v_image_B);                                      % Find the brightness local maximum (star) postions
fprintf('found %d stars in V layer!\n', sum(sum(v_image_B)));              % Show how many star find in V filter image
%set(gcf,'unit','centimeter','position',[0 10 10 10])
drawnow
%% Star match among B filter image, V filter image, and star list
size_1 = size(b_image_B);                                                                                 % Determine the size of graph (pixel times pixel)
size_2 = size(star_list);                                                                                 % Determine the size of star list
b_dM = [];
v_dM = [];
B = [];
V = [];
for k1 = 1+(error_m+outer_radius) : size_1(1)-(error_m+outer_radius)
    for k2 = 1+(error_m+outer_radius) : size_1(2)-(error_m+outer_radius)
        match = 0;
        continue1 = 0;
        vx = [];
        vy = [];
        b_dm = [];
        v_dm = [];
        if b_image_B(k1,k2) == 1
            fprintf('A star find in B filter image at x=%d, y=%d\n',k1,k2)
            if sum(sum(v_image_B(k1-error_m:k1+error_m, k2-error_m:k2+error_m))) == 0                     % If we don't match a star
                fprintf('No star matches in V filter image\n')
            end
            if sum(sum(v_image_B(k1-error_m:k1+error_m, k2-error_m:k2+error_m))) > 0                      % If we match a star
                for k3 = k1-error_m:k1+error_m
                    for k4 = k2-error_m:k2+error_m
                        if v_image_B(k3,k4) == 1                                                          % If we find the position of the matched star in V filter image
                            vx = [vx,k3];
                            vy = [vy,k4];
                        end
                    end
                end
                if length(vx) > 1                                                                         % If more than one stars are matched
                    fprintf('There are %d stars match in V filter image\n',length(vx));
                elseif length(vx) == 1                                                                        % If only one star is matched
                    fprintf('Only one star match in  V filter image at x=%d, y=%d\n',vx(1),vy(1));
                    b_x = k1;                                                                             % Input the star's x position in B filter image
                    b_y = k2;                                                                             % Input the star's y position in B filter image
                    v_x = vx(1);                                                                          % Input the star's x position in V filter image
                    v_y = vy(1);                                                                          % Input the star's y position in V filter image
                    match = 1;                                                                            % Send a match signal
                end
            end
        end
        if match == 1                                                                                     % If we match the star
            d1 = sum(sum(b_image_R(b_x-outer_radius:b_x+outer_radius,b_y-outer_radius:b_y+outer_radius)));
            d2 = sum(sum(b_image_R(b_x-inner_radius:b_x+inner_radius,b_y-inner_radius:b_y+inner_radius)));
            d3 = sum(sum(v_image_R(v_x-outer_radius:v_x+outer_radius,v_y-outer_radius:v_y+outer_radius)));
            d4 = sum(sum(v_image_R(v_x-inner_radius:v_x+inner_radius,v_y-inner_radius:v_y+inner_radius)));
            if (d1-d2) > 0 || (d3-d4) > 0                                                                 % If there are some star too close to this star
                fprintf('This star is too close to its neighborhood\n')
            end
            if (d1-d2) == 0 && (d3-d4) == 0                                                               % If no star are close to this star
                if b_image(b_x,b_y) == 65535 || v_image(v_x,v_y) == 65535                                 % If the star is too bright
                    fprintf('This star is too bright\n')
                else
                    continue1 = 1;
                end
            end
        end
        if continue1 == 1
            b = (-2.5)*log10(sum(sum(b_image(b_x-radius:b_x+radius,b_y-radius:b_y+radius))) ...
                -(sum(sum(b_image(b_x-outer_radius:b_x+outer_radius,b_y-outer_radius:b_y+outer_radius))) ...
                -sum(sum(b_image(b_x-inner_radius:b_x+inner_radius,b_y-inner_radius:b_y+inner_radius))))* ...
                (2*radius+1)^2/((2*outer_radius+1)^2- (2*inner_radius+1)^2));
            % Determine the star's equipment magnitude in B filter image
            
            v = (-2.5)*log10(sum(sum(v_image(v_x-radius:v_x+radius,v_y-radius:v_y+radius))) ...
                -(sum(sum(v_image(v_x-outer_radius:v_x+outer_radius,v_y-outer_radius:v_y+outer_radius))) ...
                -sum(sum(v_image(v_x-inner_radius:b_x+inner_radius,v_y-inner_radius:v_y+inner_radius))))* ...
                (2*radius+1)^2/((2*outer_radius+1)^2-(2*inner_radius)^2));                                          % Determine the star's equipment magnitude in V filter image
            B = [B,b];
            V = [V,v];
            b_wcs = xy2coo(b_sim(1), b_x, b_y).Cat;                                                       % Find the RA/DEC of the star in B filter image
            v_wcs = xy2coo(v_sim(1), v_x, v_y).Cat;                                                       % Find the RA/DEC of the star in V filter image
            for k5 = 1:size_2(1)
                dB_RA = abs(star_list{k5, 1}-rad2deg(b_wcs(1)));
                dB_DEC = abs(star_list{k5, 2}-rad2deg(b_wcs(2)));
                dV_RA = abs(star_list{k5, 1}-rad2deg(v_wcs(1)));
                dV_DEC = abs(star_list{k5, 2}-rad2deg(v_wcs(2)));
                if dB_RA < error_b && dB_DEC < error_b && dV_RA < error_v && dV_DEC < error_v             % If the star we find match the star on the star list
                    b_dm = [b_dm,star_list{k5, 'Bmag'}-b];                                                       % Calculate dm of B filter image
                    v_dm = [v_dm,star_list{k5, 'Vmag'}-v];                                                       % Calculate dm of V filter image
                end
            end
            if length(b_dm) == 1 && length(v_dm) == 1
                b_dM = [b_dM, b_dm];
                v_dM = [v_dM, v_dm];
                fprintf('With dm contribution\n')
            end
        end   
    end
end
%% Final calculation
B = B+median(b_dM);                                                        % Calculate the actual apparent magnitude of star in B filter image
V = V+median(v_dM);                                                        % Calculate the actual apparent magnitude of star in V filter image
B_V = B-V;                                                                 % Calculate the value of B-V
V = V-5*log10(distance)+5;                                                 % Calculate the absolute magnitude of star in V filter image
fprintf('dM of B: %f\n', median(b_dM));                                    % Show the dm we find for B filter image
fprintf('dM of V: %f\n', median(v_dM));                                    % Show the dm we find for V filter image
B_V_true1 = fitting_data_1{:,1};
V_true1 = fitting_data_1{:,2};
B_V_true2 = fitting_data_2{:,1};
V_true2 = fitting_data_2{:,2};
figure();
HR = axes;
scatter(HR ,B_V, V, 'filled');
hold on
if choose == 0
    title_1 = ['HR diagram of ',name];
    title(title_1,'FontSize',16);
elseif choose == 1
    scatter(HR, B_V_true1, V_true1, 'filled','r');
    title_1 = ['HR diagram of ',name];
    title_2 = ['with ',fit_1];
    title({title_1;title_2},'FontSize',16)
    legend('Observation','Theory')
elseif choose == 2
    scatter(HR, B_V_true2, V_true2, 'filled','r');
    title_1 = ['HR diagram of ',name];
        title_2 = ['with ',fit_2];
    title({title_1;title_2},'FontSize',16)
    legend('Observation','Theory')
elseif choose == 3
    scatter(HR, B_V_true1, V_true1, 'filled','r');
    scatter(HR, B_V_true2, V_true2, 'filled','g');
    title_1 = ['HR diagram of ',name];
    title(title_1,'FontSize',16)
    legend('Observation',fit_1,fit_2)
end
HR.YDir = 'reverse';
xlabel('B-V');
ylabel('V');
set(gcf,'unit','centimeter','position',[20 0 20 20])
