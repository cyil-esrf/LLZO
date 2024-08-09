
function [grain] = fastgrainplot(fileinfo,upars,vpars, background)


% This function must be called as eg, with u along theta and v perpendicular:
% grain = fastgrainplot04({'C:/mosaicity/', '3Dmap1_*'},{-0.1,0.1,40},{-0.1,0.1,10},bg);
%
% The background may be a number or a matrix of same size as the images.
% The function generates an output object, grain, of the following
% variables:
%
%   shape - total intensity
%   ucm, vcm - center of mass of u, v for each pixel
%   spreadu, spreadv - fwhm of cm fit for u, v for each pixel
%   uvar, vvar - variance of cm fit (spread = 2.355*sqrt(var))
%   uskew, vskew - skewness of fits (third moment) 
%   ukurt, vkurt - kurtosis of fits (forth moment)
%   oridist - total orientation distribution
%   u, v - list of angles
%   mosaicity - colorcoded image of ucm and vcm combined
%
% Image colorbars are in degrees and optimised by the function opticolor.
% Image axes are given in pixels.
%
% SRA march2017



bg = background;

[umin,umax,uint] = upars{:};    % theta
[vmin,vmax,vint] = vpars{:};    % perpendicular

%input
vstep = ((vmax-vmin)/vint);
ustep = ((umax-umin)/uint);
vlist = vmin:vstep:vmax;      % -0.275:(1/41)*2*0.275:0.275;
ulist = umin:ustep:umax;      %-0.25:(1/90)*2*0.25:0.25;
vnum = length(vlist);
unum = length(ulist);



% make imagelist
tic
[path,prefix] = fileinfo{:};
imglist = dir(strcat(path,prefix,'.edf'));
toc

% initialise
len = vnum*unum
oridist = zeros(len,1);
testimg = edf_read(strcat(path,imglist(end).name));
inttot = zeros(size(testimg));
v1sum = zeros(size(testimg));
v2sum = zeros(size(testimg));
v3sum = zeros(size(testimg));
v4sum = zeros(size(testimg));
u1sum = zeros(size(testimg));
u2sum = zeros(size(testimg));
u3sum = zeros(size(testimg));
u4sum = zeros(size(testimg));




tic
for j=1:len
    img    = edf_read(strcat(path,imglist(j).name))-bg;   % read in file and subtract background %-120 for mar2017
    img( img < 0 ) = 0;        % remove negative points
    imgsum  = sum (img(:));    % total intensity in image
 
    % for grain shape
    inttot = inttot + img;     % total intensity
   
    % for calculating moments
    vv = vlist(floor((j-1)/unum)+1);
    uu = ulist(mod(j-1,unum)+1);
    
    v1sum = v1sum + vv*img;
    v2sum = v2sum + vv.^2*img;
    v3sum = v3sum + vv.^3*img;
    v4sum = v4sum + vv.^4*img;
    
    u1sum = u1sum + uu*img;
    u2sum = u2sum + uu.^2*img;
    u3sum = u3sum + uu.^3*img;
    u4sum = u4sum + uu.^4*img;
   
    % for orientation distribution
    oridist(j) = imgsum;
    
end
toc

% Expectation value (sum(x*I/Itot))
vnorm = v1sum./inttot;
unorm = u1sum./inttot;

% Variance (x^2*I/Itot - mu^2)
vvar  = (v2sum./inttot - vnorm.^2);
uvar  = (u2sum./inttot - unorm.^2);
    % test for negative numbers
    vvar(vvar <= 0)= NaN;
    uvar(uvar <= 0)= NaN;


% FWHM (2 sqrt(2 ln(2))*sigma) (sigma = sqrt(variance))
vfwhm = 2.355*sqrt(vvar);
ufwhm = 2.355*sqrt(uvar);

% Skewness (E[x^3]-3*mu*sigma^2-mu^3)/sigma^3)
vskew = (v3sum./inttot - 3.*vnorm.*vvar - vnorm.^3)./vvar.^(3/2);
uskew = (u3sum./inttot - 3.*unorm.*uvar - unorm.^3)./uvar.^(3/2);

% Kurtosis ((E[x^4] - 4*mu*skew + 6*mu^2*sigma^2-3*mu^4)/sigma^4) 
vkurt = (v4sum./inttot - 4.*vnorm.*v3sum./inttot + 6.*vnorm.^2.*v2sum./inttot - 3*vnorm.^4)./(vvar.^2);
ukurt = (u4sum./inttot - 4.*unorm.*u3sum./inttot + 6.*unorm.^2.*u2sum./inttot - 3*unorm.^4)./(uvar.^2);

% Mosacity

unormnn = unorm(~isnan(unorm));
vnormnn = vnorm(~isnan(vnorm));

sort1 = sort( unormnn(:) );
sort2 = sort( vnormnn(:) );

Imin1 = sort1(round(length(sort1)*0.02));
Imax1 = sort1(round(length(sort1)*0.98));

Imin2 = sort2(round(length(sort2)*0.02));
Imax2 = sort2(round(length(sort2)*0.98));

C(:,:,1) = (unorm-Imin1)/(Imax1-Imin1);
C(:,:,2) = (vnorm-Imin1)/(Imax1-Imin1);
C(:,:,3) = ones(size(bg));
C(isnan(C))=0;
C(C < 0) = 0;
C(C > 1) = 1; 


% make object for output
grain.shape = inttot;
grain.ucm = unorm;
grain.vcm = vnorm;
grain.spreadv = vfwhm;
grain.spreadu = ufwhm;
grain.vvar = vvar;
grain.uvar = uvar;
grain.vskew = vskew;
grain.uskew = uskew;
grain.vkurt = vkurt;
grain.ukurt = ukurt;
grain.oridist = reshape(oridist, [unum,vnum]);
grain.u = ulist;
grain.v = vlist;
grain.mosaicity = C;


figure; 
subplot(2,2,1); imagesc(medfilt2(opticolor(grain.vcm,0.3,0.99)));axis image; axis off;title('Rolling Peak Pos'); colorbar
subplot(2,2,2); imagesc(medfilt2(opticolor(grain.spreadv,0.1,0.99)));axis image; axis off; title('Rolling FWHM'); colorbar
subplot(2,2,3); imagesc(medfilt2(opticolor(grain.vskew,0.08,0.92)));axis image; axis off; title('Rolling Skewness'); colorbar
subplot(2,2,4); imagesc(medfilt2(medfilt2(opticolor(grain.vkurt,0.15,0.85))));axis image; axis off; title('Rolling Kurtosis'); colorbar
figure;
subplot(2,2,1); imagesc(medfilt2(opticolor(grain.ucm,0.02,0.98)));axis image; axis off; title('Rocking Peak Pos'); colorbar
subplot(2,2,2); imagesc(medfilt2(opticolor(grain.spreadu,0.02,0.98)));axis image; axis off; title('Rocking FWHM'); colorbar
subplot(2,2,3); imagesc(medfilt2(opticolor(grain.uskew,0.08,0.92)));axis image; axis off; title('Rocking Skewness'); colorbar
subplot(2,2,4); imagesc(medfilt2(opticolor(grain.ukurt,0.15,0.85)));axis image; axis off; title('Rocking Kurtosis'); colorbar

figure; 
subplot(1,2,1); imagesc(medfilt2(grain.shape)); title('total intensity');axis image; axis off;
subplot(1,2,2); imagesc(medfilt2(grain.oridist)); title('total orientation distribution');axis image; axis off;

figure; imagesc(hsv2rgb(C));axis image; axis off;




end

function C = opticolor(Cin,minc, maxc)

C = Cin;
Cnn = C(~isnan(C));
sortC = sort( Cnn(:) );
Imin = sortC(round(length(sortC)*minc));
Imax = sortC(round(length(sortC)*maxc));
C(C<Imin)=Imin;
C(C>Imax)=Imax;
end


