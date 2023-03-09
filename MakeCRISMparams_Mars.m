%% Program to compute CRISM summary params from input CRISM file and header
%Params from Viviano-Beck, 2014
%Initial inputs from Lulu Pan Denoising script
%Param Code adapted from CRISM CAT IDL code
% Ari Koeppel - 2019
clear
close all;

%Use loop and filelist for multiple files
Dir = char('Y:\CRISM_work\Projects\'); %gets directory
fileID = fopen(join([Dir,'filelist.txt'],""));
list = textscan(fileID,'%s');
fclose('all');
%% loop through files
for g = 1:size(list{1,1},1)
% Find Ratioed cube    
clearvars -except Dir list g
NAME = cell2mat(join([Dir,list{1,1}(g,:),'\'],""));
Files = dir(fullfile(NAME,'*cube.img')); %gets all files
fullFileName = fullfile(NAME, Files.name);
Header = dir(fullfile(NAME,'*cube.hdr')); %gets all header files
fullHeaderName = fullfile(NAME, Header.name);
name = string(Files.name(1:end-4));

%% Read in Data
tic
fprintf(1, 'Now reading %s\n', fullFileName);
HeaderParams = envihdrread(fullHeaderName);
Data = envidataread(fullFileName,HeaderParams);
Data(Data(:,:,:)==65535)=NaN;
Data = single(Data);
wl_char = HeaderParams.wavelength;
wl_char(regexp(wl_char,'[,{}]'))=[' '];
wl_array = str2num(wl_char);
toc

%% Copute Band Params
tic
%%BD1900_2
Bands = 24;
OLINDEX3 = crism_summary_olindex3(Data,wl_array);  %Olivine
R1330=crism_sumutil_single(Data,wl_array,1.330,11); %R1330
%BD1300 = crism_sumutil_band_depth(Data,wl_array,1.080,1.320,1.750,5,15,5); %Fe2+ in Plag
LCPINDEX2 = crism_summary_lcpindex2(Data,wl_array); %Low-Ca Pyroxene                       
HCPINDEX2 = crism_summary_hcpindex2(Data,wl_array); %Hi-Ca Pyroxene
% VAR = crism_summary_var(Data,wl_array); %spectral variance parameter
ISLOPE1=crism_summary_islope1(Data,wl_array); %Spectral Slope
%BD1750_2=crism_sumutil_band_depth(Data,wl_array,1.690, 1.750, 1.815,5,3,5);                     
BD1900_2=crism_sumutil_band_depth(Data,wl_array,1.850,1.930,2.067,5,5,5); %Hydration
BD1900R2=crism_summary_bd1900r2(Data,wl_array); % H2O Continuum removed ratios so it is not longer sensitive to slope effects
BD2100_2=crism_sumutil_band_depth(Data,wl_array, 1.930, 2.132, 2.250,3,5,3); %H2O
BD2165=crism_sumutil_band_depth(Data,wl_array,2.120, 2.165, 2.230,5,3,3); %kaolinite group
BD2190=crism_sumutil_band_depth(Data,wl_array,2.120, 2.185, 2.250,5,3,3); %(Beidellite, Allophane)
D2200=crism_summary_d2200(Data,wl_array); %                                                 
MIN2200=crism_sumutil_band_depth_min(Data,wl_array, 2.120, 2.165, 2.350,2.120, 2.210,2.350,5,3,5,5,3,5); %2.16 and 2.21 micron band depth (DOUB2200): (Kaolinite) 
%BD2210_2=crism_sumutil_band_depth(Data,wl_array,2.165, 2.210, 2.290,5,5,5); %Not Kaolanite
BD2230=crism_sumutil_band_depth(Data,wl_array,2.210, 2.235, 2.252,3,3,3); %Hydroxylated Ferric Sulfate
MIN2250=crism_sumutil_band_depth_min(Data,wl_array,2.165, 2.210, 2.350,2.165, 2.265, 2.350, 5,3,5,5,3,5); %2.21 and 2.26 micron band depth (DOUB2250): (Opal)
BD2250=crism_sumutil_band_depth(Data,wl_array,2.120, 2.245, 2.340,5,7,3); %opal
BD2265=crism_sumutil_band_depth(Data,wl_array,2.210, 2.265, 2.340,5,3,5); %(jarosite, gibbsite)
BD2290=crism_sumutil_band_depth(Data,wl_array,2.250, 2.290, 2.350,5,5,5); %Fe-OH & CO2 Ice
D2300=crism_summary_d2300(Data,wl_array); 
BD2355=crism_sumutil_band_depth(Data,wl_array,2.300, 2.355, 2.450,5,5,5); %CHLORITE              
SINDEX2 = crism_sumutil_band_depth_invert(Data,wl_array, 2.120, 2.290, 2.400,5,7,3); %band depth of the continuum between the 2.1 and 2.4 micron mono- and poly-hydrated sulfate features                                           
BDCARB = crism_summary_bdcarb(Data,wl_array); %Carbonate band overtone
MIN2295_2480 = crism_sumutil_band_depth_min(Data,wl_array, 2.165, 2.295, 2.364,2.364, 2.480, 2.570,5,5,5,5,5,5); %Mg-carbonate overtone band depth (MIN2295_2480): (Mg-Carbonate)             
MIN2345_2537 = crism_sumutil_band_depth_min(Data,wl_array, 2.250, 2.345, 2.430,2.430, 2.537, 2.602,5,5,5,5,5,5); %Fe/Ca-carbonate overtone band depth (MIN2342_2537): (Fe/Ca-Carbonate)
BD2500_2=crism_sumutil_band_depth(Data,wl_array,2.364, 2.480, 2.570,5,5,5); %(Mg-Carbonate)   
%BD3400_2=crism_sumutil_band_depth(Data,wl_array,3.250, 3.420, 3.630,10,15,10); %Broad Carbonate
toc

%% Assign Params to single raster array
tic
ParamData(:,:,1) = OLINDEX3;
ParamData(:,:,2) = SINDEX2;
ParamData(:,:,3) = BDCARB;
ParamData(:,:,4) = LCPINDEX2;
ParamData(:,:,5) = HCPINDEX2;
ParamData(:,:,6) = MIN2295_2480;
ParamData(:,:,7) = ISLOPE1;
ParamData(:,:,8) = MIN2345_2537;
ParamData(:,:,9) = BD1900_2;
ParamData(:,:,10) = BD1900R2;
ParamData(:,:,11) = BD2100_2;
ParamData(:,:,12) = BD2165;
ParamData(:,:,13) = BD2190;
ParamData(:,:,14) = D2200;
ParamData(:,:,15) = MIN2200;
ParamData(:,:,16) = BD2500_2;
ParamData(:,:,17) = BD2230;
ParamData(:,:,18) = MIN2250;
ParamData(:,:,19) = BD2250;
ParamData(:,:,20) = BD2265;
ParamData(:,:,21) = BD2290;
ParamData(:,:,22) = D2300;
ParamData(:,:,23) = BD2355;
ParamData(:,:,24) = R1330;
ParamData(isnan(ParamData(:,:,:)))=65535.0;
toc

%% Write outdata
tic
ParamDataHeader = dir(fullfile(NAME,'*bandmap.hdr'));
ParamDataHeaderName = fullfile(NAME, ParamDataHeader.name);
ParamDataHeaderFormat = envihdrread(ParamDataHeaderName);
ParamDataHeaderFormat.bands = Bands;
ParamDataHeaderFormat.band_names = '{OLINDEX3, SINDEX2, BDCARB, LCPINDEX2, HCPINDEX2, MIN2295_2480, ISLOPE1, MIN2345_2537, BD1900_2, BD1900R2, BD2100_2, BD2165, BD2190, D2200, MIN2200, BD2500_2, BD2230, MIN2250, BD2250, BD2265, BD2290, D2300, BD2355, R1330}';
envihdrwrite(ParamDataHeaderFormat,join([name,'_2014Params.hdr'],""));
envidatawrite(ParamData,join([name,'_2014Params.img'],""),ParamDataHeaderFormat);
toc
end
%% Functions
function BD = crism_sumutil_band_depth(cube,wvt,low,mid,hi,low_width,mid_width,hi_width) 
%extract single bands from the cube, replacing CRISM_NAN with IEEE NAN
    Rlow = crism_sumutil_single(cube,wvt,low,low_width);
    Rmid = crism_sumutil_single(cube,wvt,mid,mid_width);
    Rhi = crism_sumutil_single(cube,wvt,hi,hi_width);

    % determine wavelength values for closest crism channels
    WL = wvt(mro_crism_lookupwv(low,wvt));
    WC = wvt(mro_crism_lookupwv(mid,wvt));
    WH = wvt(mro_crism_lookupwv(hi,wvt));
    a = (WC-WL)./(WH-WL);     %a gets multipled by the longer band
    b = 1.0-a;               %b gets multiplied by the shorter band

    %compute the band depth using precomputed a and b
    BD = 1.0-(Rmid./(b.*Rlow+a.*Rhi));
end

function img = crism_sumutil_band_depth_invert(cube,wvt,low,mid,hi,low_width,mid_width,hi_width)
    %extract single bands from the cube, replacing CRISM_NAN with IEEE NAN
    Rlow = crism_sumutil_single(cube,wvt,low,low_width);
    Rmid = crism_sumutil_single(cube, wvt,mid,mid_width);
    Rhi = crism_sumutil_single(cube, wvt,hi,hi_width);

    %determine wavelength values for closest crism channels
    WL = wvt(mro_crism_lookupwv(low,wvt));
    WC = wvt(mro_crism_lookupwv(mid,wvt));
    WH = wvt(mro_crism_lookupwv(hi, wvt));
    a = (WC-WL)./(WH-WL);     % a gets multipled by the longer band
    b = 1.0-a;               % b gets multiplied by the shorter band

    %compute the band depth using precomputed a and b
    img = 1.0 - (( b.*Rlow + a.* Rhi )./ Rmid );
end

function img = crism_sumutil_band_depth_min(cube,wvt,low1,mid1,hi1,low2,mid2,hi2,low_width1,mid_width1,hi_width1,low_width2,mid_width2,hi_width2)
    img1=crism_sumutil_band_depth(cube, wvt, low1, mid1, hi1,low_width1,mid_width1,hi_width1);
    img2=crism_sumutil_band_depth(cube, wvt, low2, mid2, hi2,low_width2,mid_width2,hi_width2);
    sz = size(cube);
    NX = sz(1);  % spatial detector
    NY = sz(2);  % spatial along track
    img=NaN(NX, NY);
    for i=1:NX
        for j=1:NY           
        	d1=img1(i,j);
            d2=img2(i,j);
            d=[d1,d2];
            img(i,j)=min(d);
        end
    end
end

function single = crism_sumutil_single(cube,wvt,band_wavelength,kernel_width)
    if nargin < 4
        kernel_width = 5;
    end
    R_indx = mro_crism_lookupwv(band_wavelength,wvt);
    half_width = fix(kernel_width/2);
    if (R_indx-half_width)<0
        minidx = 0;  %clamp to 0 if R-half_width is less than 0
        maxidx = R_indx+half_width;
    elseif (R_indx+half_width) > (length(wvt)-1)
    	maxidx = (length(wvt)-1); %clamp to max index if R+half_width is past the end
        minidx = R_indx-half_width;
    else
        minidx = R_indx-half_width;
        maxidx = R_indx+half_width;
    end
%check the wavelength extremes to ensure there are no gaps in wavelength within the kernel
    wavelength_coverage = abs(wvt(minidx)-wvt(maxidx));
    mean_wavelength_per_channel = 6.55e-3; %um/channel
    max_allowable_missing_channels = 2.0; 
    wavelength_gap_tolerance = mean_wavelength_per_channel*(kernel_width + max_allowable_missing_channels);
    if wavelength_coverage > wavelength_gap_tolerance
        fprintf('kernel wavelength coverage exceeds max wavelength gap tolerance. (%.2f um < %.2f um)',wavelength_coverage,wavelength_gap_tolerance);
        fprint('    kernel width = %.2f',kernel_width);
        fprint('    center wavelength = %.2f',wvt(R_indx))
        fprint('    wavelength table index = %.2f',R_indx)
        fprint('    kernel wavelengths: %.2f, %.2f, %.2f, %.2f, %.2f', wvt(minidx:maxidx));
    end
    %create a subsetted wavelength table and find the corresponding index in the subsetted table
    wvt_n=wvt(minidx:maxidx);  
    R_indx_n=mro_crism_lookupwv(band_wavelength,wvt_n);
    %subset the "kernel_width" spectral pixels centered on the base wavelength and
    %filter the subsetted cube using either boxcar or median filtering
    dim = size(cube);
    smthcube=zeros(dim(1),dim(2),kernel_width);
    for i = 1:dim(1)
       for j = 1:dim(2)
           smthcube(i,j,:) = smooth(permute(cube(i,j,minidx:maxidx),[3 2 1]), kernel_width);
       end
    end
    single = smthcube(:,:,R_indx_n);
end

function index = mro_crism_lookupwv(lam, wvlarr)
[val,index] = min(abs(wvlarr-lam));
end

function img = crism_summary_olindex3(cube, wvt)
    R1210 = crism_sumutil_single(cube, wvt, 1.210,7);
    R1250 = crism_sumutil_single(cube, wvt, 1.250,7);
    R1263 = crism_sumutil_single(cube, wvt, 1.263,7);
    R1276 = crism_sumutil_single(cube, wvt, 1.276,7);
    R1330 = crism_sumutil_single(cube, wvt, 1.330,7); 
    R1750 = crism_sumutil_single(cube, wvt, 1.750,7); 
    R1862 = crism_sumutil_single(cube, wvt, 1.862,7);
    
    %identify nearest CRISM wavelength
    W1210 = wvt(mro_crism_lookupwv(1.210,wvt));
    W1250 = wvt(mro_crism_lookupwv(1.250,wvt));
    W1263 = wvt(mro_crism_lookupwv(1.263,wvt));
    W1276 = wvt(mro_crism_lookupwv(1.276,wvt));
    W1330 = wvt(mro_crism_lookupwv(1.330,wvt));
    W1750 = wvt(mro_crism_lookupwv(1.750,wvt));
    W1862 = wvt(mro_crism_lookupwv(1.862,wvt));
    
    %compute the corrected reflectance interpolating 
    slope = (R1862-R1750)./(W1862-W1750);   %slope = ( R2120 - R1690 ) / ( W2120 - W1690 )
    intercept = R1862-slope.*W1862;              %intercept = R2120 - slope * W2120

    %weighted sum of relative differences

    Rc1210 = slope.*W1210+intercept;
    Rc1250 = slope.*W1250+intercept;
    Rc1263 = slope.*W1263+intercept;
    Rc1276 = slope.*W1276+intercept;
    Rc1330 = slope.*W1330+intercept;

    img = (((Rc1210-R1210)./(abs(Rc1210))).*0.1)+(((Rc1250-R1250)./(abs(Rc1250))).*0.1)+...
        (((Rc1263-R1263)./(abs(Rc1263))).*0.2)+(((Rc1276-R1276)./(abs(Rc1276))).*0.2)+(((Rc1330-R1330)./(abs(Rc1330))).* 0.4);  
end

function img = crism_summary_lcpindex2(cube,wvt)
    % extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R1690 = crism_sumutil_single(cube, wvt, 1.690,7);
    R1750 = crism_sumutil_single(cube, wvt, 1.750,7);
    R1810 = crism_sumutil_single(cube, wvt, 1.810,7);
    R1870 = crism_sumutil_single(cube, wvt, 1.870,7);
    R1560 = crism_sumutil_single(cube, wvt, 1.560,7); 
    R2450 = crism_sumutil_single(cube, wvt, 2.450,7);
    
    %identify nearest CRISM wavelength
    W1690 = wvt(mro_crism_lookupwv(1.690,wvt));
    W1750 = wvt(mro_crism_lookupwv(1.750,wvt));
    W1810 = wvt(mro_crism_lookupwv(1.810,wvt));
    W1870 = wvt(mro_crism_lookupwv(1.870,wvt));   
    W1560 = wvt(mro_crism_lookupwv(1.560,wvt));
    W2450 = wvt(mro_crism_lookupwv(2.450,wvt));
            
    %compute the corrected reflectance interpolating 
    slope = (R2450-R1560)./(W2450-W1560);
    intercept = R2450-slope.*W2450;

    %weighted sum of relative differences
    Rc1690 = slope.* W1690 + intercept;
    Rc1750 = slope.* W1750 + intercept;
    Rc1810 = slope.* W1810 + intercept;
    Rc1870 = slope.* W1870 + intercept;
       
    img=((1-(R1690./Rc1690)).*0.2) + ((1-(R1750./Rc1750)).*0.2) + ((1-(R1810./Rc1810)).*0.3) + ((1-(R1870./Rc1870)).*0.3);    
end

function img = crism_summary_hcpindex2(cube,wvt)
    %extract channels from cube replacing CRISM_NAN with IEEE_NaN
    R2120 = crism_sumutil_single(cube, wvt, 2.120,5); 
    R2140 = crism_sumutil_single(cube, wvt, 2.140,7);
    R2230 = crism_sumutil_single(cube, wvt, 2.230,7);
    R2250 = crism_sumutil_single(cube, wvt, 2.250,7);
    R2430 = crism_sumutil_single(cube, wvt, 2.430,7);
    R2460 = crism_sumutil_single(cube, wvt, 2.460,7);
    R1690 = crism_sumutil_single(cube, wvt, 1.690,7);
    %R1810 = crism_sumutil_single(cube, wvt, 1.810,7);
    R2530 = crism_sumutil_single(cube, wvt, 2.530,7);
        
    %identify nearest CRISM wavelength
    W2120 = wvt(mro_crism_lookupwv(2.120,wvt));
    W2140 = wvt(mro_crism_lookupwv(2.140,wvt));
    W2230 = wvt(mro_crism_lookupwv(2.230,wvt));
    W2250 = wvt(mro_crism_lookupwv(2.250,wvt));
    W2430 = wvt(mro_crism_lookupwv(2.430,wvt));
    W2460 = wvt(mro_crism_lookupwv(2.460,wvt));
    W1690 = wvt(mro_crism_lookupwv(1.690,wvt));
    %W1810 = wvt(mro_crism_lookupwv(1.810,wvt));
    W2530 = wvt(mro_crism_lookupwv(2.530,wvt));
        
	%compute the corrected reflectance interpolating 
	%slope = ( R2530 - R1810 ) / ( W2530 - W1810 )         
    slope = ( R2530 - R1690 )./( W2530 - W1690 );      
    intercept = R2530 - slope.* W2530;

	%weighted sum of relative differences
    Rc2120 = slope.* W2120 + intercept;
    Rc2140 = slope.* W2140 + intercept;
    Rc2230 = slope.* W2230 + intercept;
    Rc2250 = slope.* W2250 + intercept;
    Rc2430 = slope.* W2430 + intercept;
    Rc2460 = slope.* W2460 + intercept;

    img=((1-(R2120./Rc2120)).*0.1) + ((1-(R2140./Rc2140)).*0.1) + ((1-(R2230./Rc2230)).*0.15) + ((1-(R2250./Rc2250)).*0.3) + ((1-(R2430./Rc2430)).*0.2) + ((1-(R2460./Rc2460)).*0.15);
end

function img = crism_summary_islope1(cube,wvt)
    %extract individual bands with CRISM_NAN replaced with IEEE NaN
    R1815 = crism_sumutil_single(cube, wvt,1.815,5);
    R2530 = crism_sumutil_single(cube, wvt,2.530,5);
    W1815 = wvt(mro_crism_lookupwv(1.815,wvt));
    W2530 = wvt(mro_crism_lookupwv(2.530,wvt));

    %want in units of reflectance / um
    img = 1000.*(R1815-R2530)./(W2530-W1815);
end

function img = crism_summary_bd1900r2(cube,wvt)
    %extract individual channels, replacing CRISM_NANs with IEEE_NaNs
    R1908 = crism_sumutil_single(cube, wvt, 1.908,1); 
    R1914 = crism_sumutil_single(cube, wvt, 1.914,1); 
    R1921 = crism_sumutil_single(cube, wvt, 1.921,1); 
    R1928 = crism_sumutil_single(cube, wvt, 1.928,1); 
    R1934 = crism_sumutil_single(cube, wvt, 1.934,1); 
    R1941 = crism_sumutil_single(cube, wvt, 1.941,1); 
    R1862 = crism_sumutil_single(cube, wvt, 1.862,1); 
    R1869 = crism_sumutil_single(cube, wvt, 1.869,1); 
    R1875 = crism_sumutil_single(cube, wvt, 1.875,1); 
    R2112 = crism_sumutil_single(cube, wvt, 2.112,1); 
    R2120 = crism_sumutil_single(cube, wvt, 2.120,1); 
    R2126 = crism_sumutil_single(cube, wvt, 2.126,1); 
    R1815 = crism_sumutil_single(cube, wvt, 1.815,5);
    R2132 = crism_sumutil_single(cube, wvt, 2.132,5); 
    
    %retrieve the CRISM wavelengths nearest the requested values
    W1908 = wvt(mro_crism_lookupwv(1.908, wvt)); 
    W1914 = wvt(mro_crism_lookupwv(1.914, wvt));
    W1921 = wvt(mro_crism_lookupwv(1.921, wvt));
    W1928 = wvt(mro_crism_lookupwv(1.928, wvt));
    W1934 = wvt(mro_crism_lookupwv(1.934, wvt));
    W1941 = wvt(mro_crism_lookupwv(1.941, wvt));
    W1862 = wvt(mro_crism_lookupwv(1.862, wvt));
    W1869 = wvt(mro_crism_lookupwv(1.869, wvt));
    W1875 = wvt(mro_crism_lookupwv(1.875, wvt));
    W2112 = wvt(mro_crism_lookupwv(2.112, wvt));
    W2120 = wvt(mro_crism_lookupwv(2.120, wvt));
    W2126 = wvt(mro_crism_lookupwv(2.126, wvt));
    W1815 = wvt(mro_crism_lookupwv(1.815, wvt));    
    W2132 = wvt(mro_crism_lookupwv(2.132, wvt));   
    
    %compute the interpolated continuum values at selected wavelengths between 1815 and 2530
    slope = ( R2132 - R1815 )./( W2132 - W1815 );
    CR1908 = R1815 + slope .* ( W1908 - W1815 );
    CR1914 = R1815 + slope .* ( W1914 - W1815 );
    CR1921 = R1815 + slope .* ( W1921 - W1815 ); 
    CR1928 = R1815 + slope .* ( W1928 - W1815 );
    CR1934 = R1815 + slope .* ( W1934 - W1815 );
    CR1941 = R1815 + slope .* ( W1941 - W1815 ); 
    CR1862 = R1815 + slope .* ( W1862 - W1815 );
    CR1869 = R1815 + slope .* ( W1869 - W1815 );
    CR1875 = R1815 + slope .* ( W1875 - W1815 );   
    CR2112 = R1815 + slope .* ( W2112 - W1815 );
    CR2120 = R1815 + slope .* ( W2120 - W1815 );
    CR2126 = R1815 + slope .* ( W2126 - W1815 );
    
    img= 1.0-((R1908./CR1908+R1914./CR1914+R1921./CR1921+R1928./CR1928+R1934./CR1934+R1941./CR1941)./(R1862./CR1862+R1869./CR1869+R1875./CR1875+R2112./CR2112+R2120./CR2120+R2126./CR2126));
end

function img = crism_summary_d2200(cube, wvt)
    % extract individual channels, replacing CRISM_NANs with IEEE_NaNs
    R1815 = crism_sumutil_single(cube, wvt, 1.815,7);
    R2165 = crism_sumutil_single(cube, wvt, 2.165,5); 
    R2210 = crism_sumutil_single(cube, wvt, 2.210,7);     
    R2230 = crism_sumutil_single(cube, wvt, 2.230,7); 
    R2430 = crism_sumutil_single(cube, wvt, 2.430,7); %2530


    %retrieve the CRISM wavelengths nearest the requested values
    W1815 = (mro_crism_lookupwv(1.815, wvt));
    W2165 = (mro_crism_lookupwv(2.165, wvt));  
    W2210 = (mro_crism_lookupwv(2.210, wvt));
    W2230 = (mro_crism_lookupwv(2.230, wvt));
    W2430 = (mro_crism_lookupwv(2.430, wvt));

    %compute the interpolated continuum values at selected wavelengths between 1815 and 2530
    slope = ( R2430 - R1815 ) ./ ( W2430 - W1815 );
    CR2165 = R1815 + slope .* ( W2165 - W1815 );
    CR2210 = R1815 + slope .* ( W2210 - W1815 );
    CR2230 = R1815 + slope .* ( W2230 - W1815 );
    %compute d2300 with IEEE NaN values in place of CRISM NaN
    img = 1-(((R2210./CR2210)+(R2230./CR2230))./(2.*(R2165./CR2165)));    
end

function img = crism_summary_bdcarb(cube,wvt)
    % extract channels, replacing CRISM_NAN with IEEE NAN
    R2230 = crism_sumutil_single( cube, wvt,2.230,5);
    R2320 = crism_sumutil_single( cube, wvt,2.320,5);
    R2330 = crism_sumutil_single( cube, wvt,2.330,5);
    R2390 = crism_sumutil_single( cube, wvt,2.390,5);
    R2520 = crism_sumutil_single( cube, wvt,2.520,5);
    R2530 = crism_sumutil_single( cube, wvt,2.530,5);
    R2600 = crism_sumutil_single( cube, wvt,2.600,5);

    %identify nearest CRISM wavelengths
    WL1 = wvt(mro_crism_lookupwv(2.230,wvt));
    WC1 = (wvt(mro_crism_lookupwv(2.330,wvt))+wvt(mro_crism_lookupwv(2.320,wvt))).*0.5;
    WH1 =  wvt(mro_crism_lookupwv(2.390,wvt));
    a =  ( WC1 - WL1 )./ ( WH1 - WL1 );  % a gets multipled by the longer (higher wvln)  band
    b = 1.0-a;                          % b gets multiplied by the shorter (lower wvln) band

    WL2 =  wvt(mro_crism_lookupwv(2.390,wvt));
    WC2 = (wvt(mro_crism_lookupwv(2.530,wvt))+wvt(mro_crism_lookupwv(2.520,wvt))).*0.5;
    WH2 =  wvt(mro_crism_lookupwv(2.600,wvt));
    c = ( WC2 - WL2 ) ./ ( WH2 - WL2 );  % c gets multipled by the longer (higher wvln)  band
    d = 1.0-c;                           % d gets multiplied by the shorter (lower wvln) band

    %compute bdcarb
    img = 1.0-(sqrt((((R2320 + R2330).*0.5)./(b.*R2230+a.*R2390)).*(((R2520+R2530).*0.5)./(d.*R2390+c.*R2600))));  %MISTAKE d was accidently multiplied by 2230 instead of 2390  (CEV 4/12)
end

function img = crism_summary_d2300(cube, wvt)
    % extract individual channels, replacing CRISM_NANs with IEEE_NaNs
    R1815 = crism_sumutil_single(cube, wvt, 1.815,5);
    R2120 = crism_sumutil_single(cube, wvt, 2.120,5); 
    R2170 = crism_sumutil_single(cube, wvt, 2.170,5); 
    R2210 = crism_sumutil_single(cube, wvt, 2.210,5); 
    R2290 = crism_sumutil_single(cube, wvt, 2.290,3); 
    R2320 = crism_sumutil_single(cube, wvt, 2.320,3); 
    R2330 = crism_sumutil_single(cube, wvt, 2.330,3); 
    R2530 = crism_sumutil_single(cube, wvt, 2.530,5); %2530


    %retrieve the CRISM wavelengths nearest the requested values
    W1815 = wvt(mro_crism_lookupwv(1.815, wvt));
    W2120 = wvt(mro_crism_lookupwv(2.120, wvt));
    W2170 = wvt(mro_crism_lookupwv(2.170, wvt));
    W2210 = wvt(mro_crism_lookupwv(2.210, wvt));
    W2290 = wvt(mro_crism_lookupwv(2.290, wvt));
    W2320 = wvt(mro_crism_lookupwv(2.320, wvt));
    W2330 = wvt(mro_crism_lookupwv(2.330, wvt));
    W2530 = wvt(mro_crism_lookupwv(2.530, wvt));

    %compute the interpolated continuum values at selected wavelengths between 1815 and 2530
    slope = ( R2530 - R1815 )./ ( W2530 - W1815 );
    CR2120 = R1815 + slope.* ( W2120 - W1815 );
    CR2170 = R1815 + slope.* ( W2170 - W1815 );
    CR2210 = R1815 + slope.* ( W2210 - W1815 );
    CR2290 = R1815 + slope .* ( W2290 - W1815 );
    CR2320 = R1815 + slope .* ( W2320 - W1815 );
    CR2330 = R1815 + slope .* ( W2330 - W1815 );
    %compute d2300 with IEEE NaN values in place of CRISM NaN
    img = 1 - (((R2290./CR2290) + (R2320./CR2320) + (R2330./CR2330))./((R2120./CR2120) + (R2170./CR2170) + (R2210./CR2210))); 
end

% function var = crism_summary_var(cube,wvt)
%     sz = size(cube);
%     NX = sz(1);% spatial detector ( 640 detector columns unless binned )
%     NY = sz(2);%spatial # of frames along track
%     %NB = sz(3);%spectral bands ( # wavelengths )
% 
%     R1021_indx = mro_crism_lookupwv(1.021, wvt);
%     R2253_indx = mro_crism_lookupwv(2.253, wvt);
% %     indices = [R1021_indx, R2253_indx];
% %     wavelength_indx = (max(indices) - min(indices) + 1) + min(indices);
%     wvs=wvt(R1021_indx:R2253_indx);
%     varslopelams = [wvt(mro_crism_lookupwv(1.014,wvt)), wvt(mro_crism_lookupwv(2.287,wvt))];
%     varslopewl   = [mro_crism_lookupwv(1.014,wvt), mro_crism_lookupwv(2.287,wvt)];
%                     
%     varslopers = zeros(NX,NY,length(varslopelams));
% 
%     for k = 1:length(varslopelams)
%         varslopers(:,:,k) = cube(:,:,varslopewl(k));
%     end    
%  
%     %Find the actual reflectances
%     obsrs = NaN(NX,NY,length(wvs));
%     for  k = 1:length(wvs)
%         indx = mro_crism_lookupwv(wvs(k),wvt);
%         obsrs(:,:,k) = cube(:,:,indx);
%     end
% 
%     %Fit a line and find the variance:
%     var = NaN(NX,NY);
%     for j= 1:NY
%         for i=1:NX
%             fit = fitlm(varslopelams, varslopers(i,j,:)); %THIS NEEDS to be FIXED
%             predrs = fit.Coefficients{1,1} + fit.Coefficients{2,1}.*wvs;
%             var(i,j) = sum((predrs - obsrs(i,j,:))^2);
%         end
%     end
% end
