%%
clear all;
%% read fits file
Sb=FITS.read2sim('m37_b.fits');
Sv=FITS.read2sim('m37_v.fits');
table = readtable('m37_b_cata');
gnd_true_table= readtable('theorem/387M.csv'); %first column is B_V second col is V
distance = 1383;%distance of the cluster in pc
%% B-layer processing
image_b = abs(fitsread('m37_b.fits'));
image_b_raw = image_b;
%whos data
%fitsdisp('20200508_M57-001b100.fits');
image_b = log(image_b);
figure();
imagesc(abs(image_b));
% Low pass filtering
h = fspecial('gaussian', 100, 1);
image_b_f = imfilter(image_b, h);
% threshold
% m36: 6.7
image_b_f(image_b_f<7.8) = 0; %MODIFY THIS VALUE!!
% star region
image_b_ff = image_b_f;
image_b_ff(image_b_ff > 0) = 1;
figure();
imagesc(image_b_ff);
title('star region');
image_b_raw = image_b_raw.*image_b_ff;
figure();
imagesc(image_b_raw);
title('Processed Image');
% find local maximum
image_b_max = imregionalmax(image_b_f);
figure();
imagesc(image_b_max);
title('Star peaks');
fprintf('found %d stars in B layer!\n', sum(sum(image_b_max)));
%% V layer processing
image_v = abs(fitsread('m37_v.fits'));
image_v_raw = image_v;
%whos data
%fitsdisp('20200508_M57-001b100.fits');
image_v = log(image_v);
figure();
imagesc(image_v);
% Low pass filtering
h = fspecial('gaussian', 100, 1);
image_v_f = imfilter(image_v, h);
figure();
imagesc(image_v_f);
% threshold
% m36 :7
image_v_f(image_v_f<7.8) = 0; %MODIFY THIS VALUE!!
% star region 
image_v_ff = image_v_f;
image_v_ff(image_v_ff > 0) = 1;
figure();
imagesc(image_v_ff);
title('star region');
image_v_raw = image_v_raw.*image_v_ff;
% find local maximum
image_v_max = imregionalmax(image_v_f);
fprintf('found %d stars in V layer!\n', sum(sum(image_v_max)));

%% Star Matching, star magnitude, dM caculate
[height, width] = size(image_v_max);
result = [];
match = 0;
error =2;
r = 5;
err_b = 0.003; %error of catalog
err_v = 0.003;
sz = size(table);
h = sz(1);
dMb =[];
dMv =[];
for i = 1+r : height -r
    for k = 1+r : width -r
        match = 0;
        if image_v_max(i,k) == 1
            wcs_b = xy2coo(Sb(1), i, k).Cat;
            wcs_v = xy2coo(Sv(1), i, k).Cat;
            xy = coo2xy(Sb(1), wcs_v(1), wcs_v(2)).Cat;
            ii = round(xy(1));
            kk = round(xy(2));
            % find whether B also match
            if(sum(sum(image_b_max(ii-error:ii+error, kk-error:kk+error))) == 1)
                match = 1;
            end
            if match == 1
                %1.1 calculate dMb
                for g = 1:h
                    dV = sqrt((table{g, 1}- rad2deg(wcs_v(1)))^2 + (table{g,2}-rad2deg(wcs_v(2)))^2);
                    dB = sqrt((table{g, 1}- rad2deg(wcs_b(1)))^2 + (table{g,2}-rad2deg(wcs_b(2)))^2); %error distance
                    if(dB < err_b)
                        t = table{g, 'Bmag'}-(-2.5)*log10(sum(sum(image_b_raw(i-r: i+r,k-r:k+r))));
                        dMb = [dMb, t];
                        fprintf('B match at x:%d, y:%d, dM=%d\n', i, k, t);
                    end
                    if(dV < err_v)
                        t = table{g, 'Vmag'}-(-2.5)*log10(sum(sum(image_b_raw(i-r: i+r,k-r:k+r))));
                        dMv = [dMv, t];
                       fprintf('V match at x:%d, y:%d, dM=%d\n', i, k, t);
                    end
                end
               
                %2. search through a small region in B image to find peak 
                peak = 0;

                for m = ii-error: ii+error
                    for n = kk - error: kk+error
                        if  m>0 && m<=2048 && n>0 &&n<=2048&&image_b_max(m, n) == 1 
                            peak = 1;                
                            result = [result; log10(sum(sum(image_b_raw(m-r: m+r,n-r:n+r)))),...
                                log10(sum(sum(image_v_raw(i-r: i+r, k-r:k+r))))];
                            break;
                        end
                    end
                    if peak == 1
                        break;
                    end
                end
                if peak ~= 1
                    fprintf('NO MATCHED PEAK IN B IMAGE!\n');
                end
            end
        end
    end
end
% instrument mag

result = result * -2.5;
%% diagram
result(:,1) = result(:,1) + median(dMb);
result(:,2) = result(:,2) + median(dMv);
fprintf('dM of B: %f\n', median(dMb));
fprintf('dM of V: %f\n', median(dMv));
B_V = result(:,1) - result(:,2);
V = result(:,2) -5*log10(distance)+5;
B_V_true = gnd_true_table{:,1};
V_true = gnd_true_table{:,2};
figure();
HR = axes;
scatter(HR ,B_V, V, 'filled'); % our data
hold;
scatter(HR, B_V_true, V_true, 'filled','r');
HR.YDir = 'reverse';
title('sudo HR diagram');
xlabel('B-V');
ylabel('V');


