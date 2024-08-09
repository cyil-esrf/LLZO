%% VCM
% User-defined shift values for each image
x_shift = -2; % Shift in the x direction for each subsequent image
y_shift = 11;  % Shift in the y direction for each subsequent image

% Initial cropping coordinates
x_start = 500;
y_start = 550;
x_end = 1600;
y_end = 1250;

for i = 0:14
    a = sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/mosa_10x_%02d.mat', i);
    load(a); % Load the .mat file
    
    % Update cropping coordinates based on the shift
    x1 = x_start + (i * x_shift);
    y1 = y_start + (i * y_shift);
    x2 = x_end + (i * x_shift);
    y2 = y_end + (i * y_shift);

    % Crop the vcm with updated coordinates
    vcm = mosa_10x.vcm(y1:y2, x1:x2);
    
    % Apply median filter
    filteredvcm =medfilt2(medfilt2(vcm));

    % Save the filtered image as floating-point TIFF
    tiffFileName = sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/mosa_vcm_%02d.tiff', i);
    t = Tiff(tiffFileName, 'w');
    tagstruct.ImageLength = size(filteredvcm, 1);
    tagstruct.ImageWidth = size(filteredvcm, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;  % Assuming single precision is sufficient
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;  % Floating point
    tagstruct.Compression = Tiff.Compression.None;
    t.setTag(tagstruct);
    t.write(single(filteredvcm)); % Convert to single for saving
    t.close();

    % Create a figure and display the cropped image with specific properties for PNG
    fig = figure('Visible', 'off'); % Invisible figure window
    imagesc(filteredvcm, [-0.1, 0.1]);
    colormap(brewermap([], '*RdBu'));
    axis image off; % Remove axes and axis labels
    set(gca, 'Position', [0, 0, 1, 1], 'Units', 'normalized'); % Remove white space

    % Save the figure without white space, in full resolution, as PNG
    print(fig, sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/mosa_vcm_%02d.png', i), '-dpng', '-r0');
    print(fig, sprintf('/tmp_14_days/yildirim/mosa_h6_vcm_%02d.png', i), '-dpng', '-r0');

    % Close the figure to release resources
    close(fig);
end


%% UCM
% User-defined shift values for each image
x_shift = -2; % Shift in the x direction for each subsequent image
y_shift = 11;  % Shift in the y direction for each subsequent image

% Initial cropping coordinates
x_start = 500;
y_start = 550;
x_end = 1600;
y_end = 1250;

for i = 0:14
    a = sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/mosa_10x_%02d.mat', i);
    load(a); % Load the .mat file
    
    % Update cropping coordinates based on the shift
    x1 = x_start + (i * x_shift);
    y1 = y_start + (i * y_shift);
    x2 = x_end + (i * x_shift);
    y2 = y_end + (i * y_shift);

    % Crop the ucm with updated coordinates
    ucm = mosa_10x.ucm(y1:y2, x1:x2);
    
    % Apply median filter
    filteredUcm =medfilt2(medfilt2(ucm));

    % Save the filtered image as floating-point TIFF
    tiffFileName = sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/mosa_ucm_%02d.tiff', i);
    t = Tiff(tiffFileName, 'w');
    tagstruct.ImageLength = size(filteredUcm, 1);
    tagstruct.ImageWidth = size(filteredUcm, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;  % Assuming single precision is sufficient
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;  % Floating point
    tagstruct.Compression = Tiff.Compression.None;
    t.setTag(tagstruct);
    t.write(single(filteredUcm)); % Convert to single for saving
    t.close();

    % Create a figure and display the cropped image with specific properties for PNG
    fig = figure('Visible', 'off'); % Invisible figure window
    imagesc(filteredUcm, [-0.1, 0.1]);
    colormap(brewermap([], '*RdBu'));
    axis image off; % Remove axes and axis labels
    set(gca, 'Position', [0, 0, 1, 1], 'Units', 'normalized'); % Remove white space

    % Save the figure without white space, in full resolution, as PNG
    print(fig, sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/mosa_ucm_%02d.png', i), '-dpng', '-r0');
    print(fig, sprintf('/tmp_14_days/yildirim/mosa_h6_ucm_%02d.png', i), '-dpng', '-r0');

    % Close the figure to release resources
    close(fig);
end


%%
for i = 0:40
a=sprintf('/tmp_14_days/yildirim/rex_275s/mosa_275s_%02d.mat', i)    
load(a);

figure;imagesc(x,y,medfilt2(mosa_10x.ucm),[-0.01, 0.05]);colormap(brewermap([], '*RdBu'));axis image;axis([20 100 80 160])
%figure;imagesc(x,y,hsv2rgb(mosa_275s.mosaicity));axis image;axis([20 100 80 160])
saveas(gcf,sprintf('/tmp_14_days/yildirim/rex_275s/mosa_ucm_%02d.tif', i))
end
%%
for i = 0:40
%a=sprintf('/tmp_14_days/yildirim/LLZO/mosa/mosa_10x_%02d.mat', i)    
a=sprintf('/tmp_14_days/yildirim/mosa_10x_%02d.mat', i)
load(a);
ucm=sprintf('mosa_10x_%2.ucm', i);

figure;imagesc(x,y,medfilt2(mosa_10x.ucm),[-0.05, 0.05]);colormap(brewermap([], '*RdBu'));axis image;set(gca,'Visible', 'off')
%figure;imagesc(x,y,hsv2rgb(mosa_275s.mosaicity));axis image;axis([20 100 80 160])
saveas(gcf,sprintf('/tmp_14_days/yildirim/mosa_ucm_%02d.png', i));
end


%% Alternative

% User-defined shift values for each image
x_shift = -2; % Shift in the x direction for each subsequent image
y_shift = 15;  % Shift in the y direction for each subsequent image

% Initial cropping coordinates
x_start = 500;
y_start = 400;
x_end = 1600;
y_end = 1500;

for i = 0:13
    a = sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h3_mosa/mosa_10x_%02d.mat', i);
    load(a); % Load the .mat file
    
    % Update cropping coordinates based on the shift
    x1 = x_start + (i * x_shift);
    y1 = y_start + (i * y_shift);
    x2 = x_end + (i * x_shift);
    y2 = y_end + (i * y_shift);

    % Crop the ucm with updated coordinates
    ucm = mosa_10x.ucm(y1:y2, x1:x2);

    % Apply median filter
    filteredUcm = medfilt2(medfilt2(ucm));

    % Apply Gaussian smoothing to smoothen the image
    smoothUcm = imgaussfilt(filteredUcm, 1); % Adjust the sigma value as needed

    % Save the smoothed image as floating-point TIFF
    tiffFileName = sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h5_mosa/mosa_ucm_%02d.tiff', i);
    t = Tiff(tiffFileName, 'w');
    tagstruct.ImageLength = size(smoothUcm, 1);
    tagstruct.ImageWidth = size(smoothUcm, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;  % Assuming single precision is sufficient
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;  % Floating point
    tagstruct.Compression = Tiff.Compression.None;
    t.setTag(tagstruct);
    t.write(single(smoothUcm)); % Convert to single for saving
    t.close();

    % Create a figure and display the cropped image with specific properties for PNG
    fig = figure('Visible', 'off'); % Invisible figure window
    imagesc(smoothUcm, [-0.1, 0.1]);
    colormap(brewermap([], '*RdBu'));
    axis image off; % Remove axes and axis labels
    set(gca, 'Position', [0, 0, 1, 1], 'Units', 'normalized'); % Remove white space

    % Save the figure without white space, in full resolution, as PNG
    print(fig, sprintf('/data/projects/drex/H2IRONMAKE/new_analysis/h5_mosa/mosa_ucm_%02d.png', i), '-dpng', '-r0');
    print(fig, sprintf('/tmp_14_days/yildirim/mosa_h5_ucm_%02d.png', i), '-dpng', '-r0');

    % Close the figure to release resources
    close(fig);
end



%%


for i = 0:40
    a = sprintf('/tmp_14_days/yildirim/matlab_h1/mosa_10x_%02d.mat', i);
    load(a);
    ucm = mosa_10x.ucm(400:1200, 800:1800);  % Crop the ucm to specific pixel limits
    
    % Create a figure and display the cropped image with specific properties
    fig = figure;
    imagesc(medfilt2(ucm), [-0.1, 0.1]);
    colormap(brewermap([], '*RdBu'));
    axis image;
    set(gca, 'Visible', 'off');

    % Set the 'PaperPositionMode' to 'auto' to save in full resolution
    set(gcf, 'PaperPositionMode', 'auto');

    % Save the figure without white space
    saveas(fig, sprintf('/tmp_14_days/yildirim/mosa_ucm_%02d.png', i));
    
    % Close the figure to release resources (optional)
    close(fig);
end

%%
%%


for i = 0:40
    a = sprintf('/tmp_14_days/yildirim/mosa_10x_%02d.mat', i);
    load(a);
    mosa = mosa_10x.mosaicity;%(1200:1700, 600:1400,[]);  % Crop the ucm to specific pixel limits
    
    % Create a figure and display the cropped image with specific properties
    fig = figure;
    imagesc(hsv2rgb(mosa));
    %colormap(brewermap([], '*RdBu'));
    axis image;
    set(gca, 'Visible', 'off');

    % Set the 'PaperPositionMode' to 'auto' to save in full resolution
    set(gcf, 'PaperPositionMode', 'auto');

    % Save the figure without white space
    saveas(fig, sprintf('/tmp_14_days/yildirim/mosa_mosa_%02d.png', i));
    
    % Close the figure to release resources (optional)
    close(fig);
end

%%

figure(1);
subplot(1,5,1);imagesc(x,y,mosa_10x.shape,[1000 8500]);axis image;title '240s';
subplot(1,5,2);imagesc(x,y,mosa_10x.shape,[1000 8500]);axis image;title '275s';
subplot(1,5,3);imagesc(x,y,mosa_10x.shape,[2000 8500]);axis image;title '350s';
subplot(1,5,4);imagesc(x,y,mosa_10x.shape,[1000 8500]);axis image;title '410s'
subplot(1,5,5);imagesc(x,y,mosa_10x.shape,[2000 8500]);axis image;title '475s';

figure(2);
subplot(1,5,1);imagesc(x,y,mosa_10x.ucm,[-0.02 0.03]);axis image;title '240s';
subplot(1,5,2);imagesc(x,y,mosa_10x.ucm,[-0.02 0.03]);axis image;title '275s';
subplot(1,5,3);imagesc(x,y,mosa_10x.ucm,[-0.02 0.03]);axis image;title '350s';
subplot(1,5,4);imagesc(x,y,mosa_10x.ucm,[-0.02 0.03]);axis image;title '410s'
subplot(1,5,5);imagesc(x,y,mosa_10x.ucm,[-0.02 0.03]);axis image;title '475s';

%% new stuff
%mosaFiles =fullfile('/gpfs/easy/data/visitor/ma5408/id06-hxm/restored_dataset/DR700C1minS3take2DR700C1minS3take2posth3DFXM')
mosaFiles =fullfile('/gpfs/easy/data/visitor/ma5408/id06-hxm/restored_dataset/DR700C1minS3take2DR700C1minS3take2posth6DFXM/','mosalayer10x_*');
imageFilesMosa = dir(mosaFiles);
%mkdir /tmp_14_days/yildirim/rex_240s_july_fesi/;
for i = 1: length(imageFilesMosa)
   name = imageFilesMosa(i).name;
                                                                                                                                                                                                                                                                                                                                            
    %fileinfo = {sprintf('/tmp_14_days/yildirim/%s/', name), '*'}
    %fileinfo = {sprintf('/data/visitor/ma5408/id06-hxm/DR_700C_1min_S3_take2/DR_700C_1min_S3_take2_post_h2_DFXM/%s/', name), '*'}
    fileinfo = {sprintf('/gpfs/easy/data/visitor/ma5408/id06-hxm/restored_dataset/DR700C1minS3take2DR700C1minS3take2posth6DFXM/%s/', name), '*'};
    upars = {-0.3,0.3, 23} ; %diffry range in mosaicity PN1; num steps is 10-1
    vpars = {-0.375, 0.375, 10} ; %phi range in mosaicity PN1; num steps is 14
    %fastgrainplot({'C:/mosaicity/', '3Dmap1_*'},{-0.1,0.1,40},{-0.1,0.1,10},bg);
    mosa_10x= fastgrainplot(fileinfo,upars,vpars, bg_pco_zap1s+20);
save(fullfile('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/',sprintf('mosa_10x_%s.mat', name(length(name)-1:length(name)))), 'mosa_10x');
    imwrite(ind2rgb(im2uint8(mat2gray(mosa_10x.shape)), jet(256)), fullfile('/data/projects/drex/H2IRONMAKE/new_analysis/h6_mosa/',sprintf('zsum_layer_%s.png', name(length(name)-1:length(name)))));
    close all

end


%%




%%

%% loop for all layers rocking 


mosaFiles = fullfile('/data/id06-hxm/inhouse/opid06/blc12704/id06/S1_590A/S1_590A_Al_200_ff_saturday/', 'Al_111_rockinglayer_*');
imageFilesMosa = dir(mosaFiles);
mkdir /tmp_14_days/yildirim/Al200/new/ ;%/users/yildirim/Documents/correction/Al_multireflection_27-1-21/Al_200_rockinglayers/;
for i = 1: length(imageFilesMosa)
   name = imageFilesMosa(i).name;

    fileinfo = {sprintf('/data/id06-hxm/inhouse/opid06/blc12704/id06/S1_590A/S1_590A_Al_200_ff_saturday/%s/', name), '*'}

    upars = {-0.06, -0.06, 30} %diffry range in mosaicity PN1; num steps is 10-1
    %vpars = {-0.6, 0.6, 12} %phi range in mosaicity PN1; num steps is 14

    %fastgrainplot({'C:/mosaicity/', '3Dmap1_*'},{-0.1,0.1,40},{-0.1,0.1,10},bg);
   
        
        rocking_200= rocking_limiter(fileinfo,upars,bg+10);

        
     

save(fullfile('/tmp_14_days/yildirim/Al200/new/',sprintf('rocking_200_%s.mat', name(length(name)-1:length(name)))), 'rocking_200');
    imwrite(rocking_200.shape, fullfile('/tmp_14_days/yildirim/Al200/new/',sprintf('Al200_%s.png', name(length(name)-1:length(name)))));
    close all


end
%%
N = 10;          % number of variables
%method 1
for k=1:N
    temp_var = strcat( 'variable_',num2str(k) );
    eval(sprintf('%s = %g',temp_var,k*2));
end
% method 2, more advisable
%%
for k=1:10
    my_field = strcat('v',num2str(k));
    variable.(my_field) = k*2;
end


%%




%%



% names and path
data_dir = ('/data/id06-hxm/inhouse/opid06/blc12704/id06/S1_590A/S1_590A_Al_1-1-1_ff_2xobj/');
grain = 'Al_1-1-1';
mkdir /tmp_14_days/yildirim/Al_rex_analysis/Al_200/;
save_dir = sprintf('/tmp_14_days/yildirim/Al_rex_analysis/Al_200/');
save_path =sprintf('%s/Al_200.h5', save_dir);
%%
% load bg
%background = load('bg.mat')
%bg = background.bg;

% process all mosaicity layers in Ni_s11_farfield_g3. 
% Grain 3 coarse mosaicity scan. Scan are -3/+3 deg in diffry with 
% 40 images, -0.6/+0.6 deg in chi with 12 images. 
mosa_10x = {}
for i = 0:10
    acq_dir = sprintf('%s/Al_1-1-1_ff_2xobj_mosalayer_%02d/', data_dir, i)
    mosa_10x{i + 1} = fastgrainplot({acq_dir, '*'}, {-0.06, 0.06, 30}, {-0.12, 0.12, 10}, bg + 25);
    %save(sprintf('mosa_200_%02d.mat', i),'mosa_200');
    close all;
end
%%
% use a 4D array to store mosaicity
n = size(mosa_10x, 2);
mosaicity_200 = zeros(n, 2160, 2560, 3);
for i = 1:n
    mosaicity_200(i, :, :, :) = hsv2rgb(mosa_10x{i}.mosaicity);
end

% save the results
h5create(save_path, '/mosaicity', size(mosaicity_200), 'ChunkSize', [n, 54, 64, 3], 'Deflate', 9);
h5write(save_path, '/mosaicity', mosaicity_200);

%%

%bg = background.bg;

% process all mosaicity layers in Ni_s11_farfield_g3. 
% Grain 3 coarse mosaicity scan. Scan are -3/+3 deg in diffry with 
% 40 images, -0.6/+0.6 deg in chi with 12 images. 
mosa_10x = {}
for i = 0:40
    acq_dir = sprintf('%s/mosa_10x_%02d/', data_dir, i)
   mosa_10x{i + 1} = fastgrainplot({acq_dir, '*'}, {-0.06, 0.06, 30}, {-0.12, 0.12, 10}, bg + 15);
   save(sprintf('mosa_200_%02d.mat', i),'mosa_200');
    close all;
end
%%
% use a 4D array to store mosaicity
n = size(mosa_10x, 2);
mosaicity_200 = zeros(n, 2160, 2560, 3);
for i = 1:n
    mosaicity_200(i, :, :, :) = hsv2rgb(mosa_10x{i}.mosaicity);
end

% save the results
h5create(save_path, '/mosaicity', size(mosaicity_200), 'ChunkSize', [n, 54, 64, 3], 'Deflate', 9);
h5write(save_path, '/mosaicity', mosaicity_200);

%%
n = size(mosa_10x, 2);
oridist = zeros(n, 31, 11);
for i = 1:n
    oridist_200(i, :, :) = mosa_10x{i}.oridist;
end


%%

[m,n,p] = size(oridist);
[X,Y,Z]= meshgrid(1:n, 1:m,1:p);
isosurface(X,Y,Z,oridist,100);
%%
% process all strain layers in Ni_s11_farfield_g3. 
% Grain 3 coarse strain scan. Scan are -3/+3 deg in the first axis with 
% 40 images, -0.05/+0.05 deg in the second axis with 20 images.
rocking_200 = {}
for i = 0:20
    acq_dir = sprintf('%s/Al_111_mosalayer_%02d/', data_dir, i)
    rocking_200{i + 1} = fastgrainplot({acq_dir, '*'}, {-0.6, 0.6, 30}, {-0.06, 0.06, 15}, bg + 15);
    close all;
end

% compute and plot the strain
%strain= -(sind(stra_G3{5}.vcm/2))/(tand(20.7/2));
%figure;
%imagesc(strain, [-0.001, 0.001]);
%daspect([1 tand(20.5) 1])

% use a 4D array to store strain
n = size(rocking_200, 2);
rocking_200 = zeros(n, 2160, 2560);
for i = 1:n
    rocking_200(i, :, :) = -(sind(rocking_200{i}.vcm/2)) / (tand(20.7/2));
end

% save the results
h5create(save_path, 'strain', size(rocking_200), 'ChunkSize', [n, 54, 64], 'Deflate', 9);
h5write(save_path, 'strain', rocking_200);


%% ROCKING LIMITER
% names and path
data_dir = sprintf('/data/id06-hxm/inhouse/opid06/blc12704/id06/S1_590A/');
grain = 'Al_200'
save_dir = sprintf('/tmp_14_days/yildirim/Al_rex_analysis_2/Al_200_2');
save_path = sprintf('%s/Al_200.h5', save_dir);

%%
rocking_200 = {}
for i = 150:300
    acq_dir = sprintf('%s/Al_111_rockinglayer_%02d/', data_dir, i)
    rocking_200{i + 1} = rocking_limiter({acq_dir, '*'}, {-0.06, 0.06, 30}, bg + 20);
    close all;
end

% compute and plot the strain
%strain= -(sind(stra_G3{5}.vcm/2))/(tand(20.7/2));
%figure;
%imagesc(strain, [-0.001, 0.001]);
%daspect([1 tand(20.5) 1])

% use a 4D array to store strain
n = size(rocking_200, 2);
rocking = zeros(n, 2160, 2560);
for i = 1:n
    rocking(i, :, :) = rocking_200{i}.ucm;
   
end

for k=1:n
     rgbImage = ind2rgb(rocking_200{1,k}.spreadu, jet);
     imwrite(rgbImage, sprintf('rocking_200_fwhm_%02d.png', k));
end
%%

for i=1:11
    figure(i);imagesc(x,y, medfilt2(mosa_111_real{i}.ucm), [-0.02 0.01]);
    colormap(brewermap([], '*Spectral'));
    axis image;colorbar
    
end
%%
for i = 0:40
a=sprintf('/tmp_14_days/yildirim/rex_275s/mosa_275s_%02d.mat', i)    
load(a);

figure;imagesc(x,y,medfilt2(mosa_10x.ucm),[-0.01, 0.05]);colormap(brewermap([], '*RdBu'));axis image;axis([20 100 80 160])
%figure;imagesc(x,y,hsv2rgb(mosa_275s.mosaicity));axis image;axis([20 100 80 160])
saveas(gcf,sprintf('/tmp_14_days/yildirim/rex_275s/mosa_ucm_%02d.tif', i))
end

close all



%%

figure; 

subplot (4,5,1);imagesc(hsv2rgb(mosa_2x_1.mosaicity)); title 'layer 1';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,2);imagesc(hsv2rgb(mosa_2x_2.mosaicity)); title 'layer 2';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,3);imagesc(hsv2rgb(mosa_2x_3.mosaicity)); title 'layer 3';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,4);imagesc(hsv2rgb(mosa_2x_4.mosaicity)); title 'layer 4';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,5);imagesc(hsv2rgb(mosa_2x_5.mosaicity)); title 'layer 5';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,6);imagesc(hsv2rgb(mosa_2x_6.mosaicity)); title 'layer 6';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,7);imagesc(hsv2rgb(mosa_2x_7.mosaicity)); title 'layer 7';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,8);imagesc(hsv2rgb(mosa_2x_8.mosaicity)); title 'layer 8';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,9);imagesc(hsv2rgb(mosa_2x_9.mosaicity)); title 'layer 9';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,10);imagesc(hsv2rgb(mosa_2x_10.mosaicity)); title 'layer 10';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,11);imagesc(hsv2rgb(mosa_2x_11.mosaicity)); title 'layer 11';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,12);imagesc(hsv2rgb(mosa_2x_12.mosaicity)); title 'layer 12';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,13);imagesc(hsv2rgb(mosa_2x_13.mosaicity)); title 'layer 13';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,14);imagesc(hsv2rgb(mosa_2x_14.mosaicity)); title 'layer 14';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,15);imagesc(hsv2rgb(mosa_2x_15.mosaicity)); title 'layer 15';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,16);imagesc(hsv2rgb(mosa_2x_16.mosaicity)); title 'layer 16';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,17);imagesc(hsv2rgb(mosa_2x_17.mosaicity)); title 'layer 17';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,18);imagesc(hsv2rgb(mosa_2x_18.mosaicity)); title 'layer 18';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,19);imagesc(hsv2rgb(mosa_2x_19.mosaicity)); title 'layer 19';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])
subplot (4,5,20);imagesc(hsv2rgb(mosa_2x_20.mosaicity)); title 'layer 20';set(gca,'Xticklabel',[],'Yticklabel', []);axis([700 1900 900 1400])



%%

for i=1:41
    figure(i);imagesc(hsv2rgb(mosa_10x{1,i}.mosaicity));daspect([sind(17.5) 1 1]);set(gca,'Xticklabel',[],'Yticklabel', []);
    %figure;subplot(1,10,5imagesc
end





%%
figure; 

subplot (4,5,1);imagesc(mosa_2x_1.oridist); title 'layer 1';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,2);imagesc(mosa_2x_2.oridist); title 'layer 2';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,3);imagesc(mosa_2x_3.oridist); title 'layer 3';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,4);imagesc(mosa_2x_4.oridist); title 'layer 4';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,5);imagesc(mosa_2x_5.oridist); title 'layer 5';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,6);imagesc(mosa_2x_6.oridist); title 'layer 6';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,7);imagesc(mosa_2x_7.oridist); title 'layer 7';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,8);imagesc(mosa_2x_8.oridist); title 'layer 8';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,9);imagesc(mosa_2x_9.oridist); title 'layer 9';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,10);imagesc(mosa_2x_10.oridist); title 'layer 10';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,11);imagesc(mosa_2x_11.oridist); title 'layer 11';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,12);imagesc(mosa_2x_12.oridist); title 'layer 12';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,13);imagesc(mosa_2x_13.oridist); title 'layer 13';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,14);imagesc(mosa_2x_14.oridist); title 'layer 14';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,15);imagesc(mosa_2x_15.oridist); title 'layer 15';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,16);imagesc(mosa_2x_16.oridist); title 'layer 16';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,17);imagesc(mosa_2x_17.oridist); title 'layer 17';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,18);imagesc(mosa_2x_18.oridist); title 'layer 18';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,19);imagesc(mosa_2x_19.oridist); title 'layer 19';set(gca,'Xticklabel',[],'Yticklabel', []);
subplot (4,5,20);imagesc(mosa_2x_20.oridist); title 'layer 20';set(gca,'Xticklabel',[],'Yticklabel', []);

%%



%%
subplot (4,5,1);contourf(mosa_10x{1,1}.oridist, 10); title 'layer 1';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,2);contourf(mosa_10x{1,2}.oridist, 10); title 'layer 2';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,3);contourf(mosa_10x{1,3}.oridist, 10); title 'layer 3';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,4);contourf(mosa_10x{1,4}.oridist, 10); title 'layer 4';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,5);contourf(mosa_10x{1,5}.oridist, 10); title 'layer 5';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,6);contourf(mosa_10x{1,6}.oridist, 10); title 'layer 6';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,7);contourf(mosa_10x{1,7}.oridist, 10); title 'layer 7';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,8);contourf(mosa_10x{1,8}.oridist, 10); title 'layer 8';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,9);contourf(mosa_10x{1,9}.oridist, 10); title 'layer 9';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,10);contourf(mosa_10x{1,10}.oridist, 10); title 'layer 10';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,11);contourf(mosa_10x{1,11}.oridist, 10); title 'layer 11';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,12);contourf(mosa_10x{1,12}.oridist, 10); title 'layer 12';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,13);contourf(mosa_10x{1,13}.oridist, 10); title 'layer 13';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,14);contourf(mosa_10x{1,14}.oridist, 10); title 'layer 14';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,15);contourf(mosa_10x{1,15}.oridist, 10); title 'layer 15';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,16);contourf(mosa_10x{1,16}.oridist, 10); title 'layer 16';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,17);contourf(mosa_10x{1,17}.oridist, 10); title 'layer 17';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,18);contourf(mosa_10x{1,18}.oridist, 10); title 'layer 18';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,19);contourf(mosa_10x{1,19}.oridist, 10); title 'layer 19';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,20);contourf(mosa_10x{1,20}.oridist, 10); title 'layer 20';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')

%%
subplot (4,5,1);imagesc(mosa_10x{1,1}.oridist); title 'layer 1';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,2);imagesc(mosa_10x{1,2}.oridist); title 'layer 2';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,3);imagesc(mosa_10x{1,3}.oridist); title 'layer 3';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,4);imagesc(mosa_10x{1,4}.oridist); title 'layer 4';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,5);imagesc(mosa_10x{1,5}.oridist); title 'layer 5';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,6);imagesc(mosa_10x{1,6}.oridist); title 'layer 6';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,7);imagesc(mosa_10x{1,7}.oridist); title 'layer 7';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,8);imagesc(mosa_10x{1,8}.oridist); title 'layer 8';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,9);imagesc(mosa_10x{1,9}.oridist); title 'layer 9';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,10);imagesc(mosa_10x{1,10}.oridist); title 'layer 10';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,11);imagesc(mosa_10x{1,11}.oridist); title 'layer 11';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,12);imagesc(mosa_10x{1,12}.oridist); title 'layer 12';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,13);imagesc(mosa_10x{1,13}.oridist); title 'layer 13';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,14);imagesc(mosa_10x{1,14}.oridist); title 'layer 14';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,15);imagesc(mosa_10x{1,15}.oridist); title 'layer 15';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,16);imagesc(mosa_10x{1,16}.oridist); title 'layer 16';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,17);imagesc(mosa_10x{1,17}.oridist); title 'layer 17';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,18);imagesc(mosa_10x{1,18}.oridist); title 'layer 18';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,19);imagesc(mosa_10x{1,19}.oridist); title 'layer 19';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')
subplot (4,5,20);imagesc(mosa_10x{1,20}.oridist); title 'layer 20';set(gca,'Xticklabel',[],'Yticklabel', []);colormap 'jet';set(gca,'ColorScale','log')

%%
h5create(save_path, '/mosaicity', size(mosaicity_200), 'ChunkSize', [n, 54, 64, 3], 'Deflate', 9);
h5write(save_path, '/mosaicity', mosaicity_200);
