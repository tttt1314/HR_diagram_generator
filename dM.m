Sb=FITS.read2sim('m37_b.fits');
table = readtable('m37_b_cata');
gnd_true_tablke= readtable('theorem/1G.csv');
%table_matrix = table{:,:};
%% stars's diagonal is about 30~40 pixels
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
image_b_f(image_b_f<8.6) = 0;
% star region
image_b_ff = image_b_f;
image_b_ff(image_b_ff > 0) = 1;
figure();
imagesc(image_b_ff);
title('star region');
image_b_raw = image_b_raw.*image_b_ff;
% find local maximum
image_b_max = imregionalmax(image_b_f);
fprintf('found %d stars in B layer!\n', sum(sum(image_b_max)));
%% dM calculate
% 1. find stars in image
err = 0.005
r = 5;
target = 'Bmag' % target layer!!
[height, width] = size(image_b_max);
sz = size(table);
h = sz(1);
dM = [];
for m = 1+r : height-r
    for n = 1+r : width-r
        match = 0;
        if(image_b_max(m,n) == 1)
            wcs = xy2coo(Sb(1), m, n).Cat; % get the star ra/dec
            for k = 1:h
                d = sqrt((table{k, 1}- rad2deg(wcs(1)))^2 + (table{k,2}-rad2deg(wcs(2)))^2);
                if(d < err)
                    t = table{k, target}-log(sum(sum(image_b_raw(m-r: m+r,n-r:n+r))));
                    dM = [dM, t];
                    fprintf('match at x:%d, y:%d, dM=%d\n', m, n, t);
                    match = 1;
                    break;
                end
            end
        end
    end
end
dM_b = median(dM); 
fprintf('dM of B: %.2f\n', dM_b);

