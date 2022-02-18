function dir = drLag4Drev4_hcp( name, TR, vols, PosiMax, THR, FIXED, Sm, range)
%
%	No resampling version - tracking by TR
%	Lag mapping of 4D BOLD data by Toshihiko Aso@RIKEN-BDR, Kobe, Japan
%
%	Images are expected to have been slice-timing corrected, realigned, and spatially normalized to the MNI space.Spatial smoothing by 8 mm FWHM is applied in the script, under the assumption of smooth perfusion lag structure (compared to the neuronal activity). This parameter is subject to change if you assume otherwise, but can result in removal of neurovascular coupling during the "deperfusioning" treatment based on this procedure (Aso 2019; Erdogan 2016 Front Human Neurosci).
%
%	function dir = drLag4D( name, TR, vols, PosiMax, THR, FIXED, range)
%	name: String to be added to the result folder name
%	TR: Repetition time in second.
%	vols: Specify single 4D BOLD data file to process.
%	PosiMax: Determines tracking range IN *TR* -PosiMax - +PosiMax 
%			Note this parameter affects the bandpass filter. 
%	THR: Minimum height of the valid cross-correlogram peak for lag mapping
%	FIXED: Specify tracking method. 
%			0 = recursive tracking, nonzero = fixed seed tracking
%	range: Specify fourth dimension (time) of the 4D data to be used.
%			ex. 1:500 - use first 300 volumes
%			ex.	1:2:500 - "decimate" to see the effect of sampling rate (see Aso 2017; don't forget to double the TR!)
%
%	References
%	Aso, T., et al. (2017) Frontiers in Neuroscience, 11(MAY), 1-13. https://doi.org/10.3389/fnins.2017.00256
%	Satow, T., et al. (2017) Alteration of venous drainage route in idiopathic normal pressure hydrocephalus and normal aging. Frontiers in Aging Neuroscience, 9(NOV), 1-10. https://doi.org/10.3389/fnagi.2017.00387
%	Nishida, S., et al. (2018) Resting-state Functional Magnetic Resonance Imaging Identifies Cerebrovascular Reactivity Impairment in Patients With Arterial Occlusive Diseases: A Pilot Study. Neurosurgery, 0(0), 1-9. https://doi.org/10.1093/neuros/nyy434
%	Aso, T., et al. (2019) Axial variation of deoxyhemoglobin density as a source of the low-frequency time lag structure in blood oxygenation level-dependent signals. PLOS One, https://doi.org/10.1101/658377

global Fdir
Fdir = '/mnt/pub/devel/fsl-5.0.9/bin/'

% default values

if nargin< 8
	range = [];
 if nargin< 7
	Sm = '8';
  if nargin< 6
	FIXED = '0';
	if nargin< 5
		THR = '0.3';
		if nargin< 4
			PosiMax = '7';
			if nargin< 3
				error
			end
		end
	end
  end
 end
else
	range = str2num( range);
end

TR = str2num( TR)
PosiMax = str2num( PosiMax)
THR = str2num( THR)
FIXED = str2num( FIXED)
Sm = str2num( Sm)

save params.mat

setenv('FSLOUTPUTTYPE', 'NIFTI'); % to tell what the output type should be

MaxLag = PosiMax*2;		% lag range in TR
ULfreq = 1/MaxLag/TR*.8;		% upper limit of frequency for uniquely determining the peak

Nseeds = PosiMax*2 + 1;% number of regions/seed timeseries

limit = ceil( PosiMax);
Cmap = jet( Nseeds);

if FIXED==0
	dir = [ pwd '/Lag_rec_' num2str( MaxLag) 'TR_thr' num2str( round( 10*THR)) '_sm' num2str( round( Sm)) '_' name]
else
	dir = [ pwd '/Lag_fix_' num2str( MaxLag) 'TR_thr' num2str( round( 10*THR)) '_sm' num2str( round( Sm)) '_' name]
end

if exist( [ dir '/LagMap.nii'], 'file')
	disp( '..using existing LagMap')
	return
end

Sm = Sm / sqrt( 8*log(2)); % Sm-mm smoothing

disp('Preparing the data...') %-------------------------------------------

S = system( [ Fdir 'fslmaths ' vols ' -Tstd SD.nii']); if S, return, end
S = system( [ Fdir 'fslmaths ' vols ' -Tmean Tmean']); if S, return, end
S = system( [ Fdir 'fslmaths Tmean -thrp 25 Mask']); % Mask 

%S = system( ['fslmaths ' vols ' -sub Tmean -mul 1000 -div SD Data']); if S, return, end
S = system( [ Fdir 'fslmaths ' vols ' -mul 100 -div Tmean Data' ]); if S, error, end

disp('Filtering...')
if Sm>0
	S = system( [ Fdir 'fslmaths Data -nan  -s ' num2str( Sm) ' -bptf ' ...
	num2str( 1/( 0.008*2.35*TR)) ' ' num2str( 1/( ULfreq*2.35*TR)) ' ' pwd '/masked_' num2str( MaxLag) 's.nii']);
	if S, return, end
else
	S = system( [ Fdir 'fslmaths Data -nan -bptf ' ...
	num2str( 1/( 0.008*2.35*TR)) ' ' num2str( 1/( ULfreq*2.35*TR)) ' ' pwd '/masked_' num2str( MaxLag) 's.nii']); 
	if S, return, end
end

% initial seed is from the global cerebral signal
[ P,~,~] = fileparts( mfilename('fullpath'));
%ROIimage = spm_read_vols( spm_vol( drReslice( 'Mask.nii', [ P '/BrainMask_lag.nii'])));
ROIimage =  spm_read_vols( spm_vol( [ P '/BrainMask_lag_subsamp2offc.nii']));
Bmask =  spm_read_vols( spm_vol( 'Mask.nii'));

V = spm_vol( ['masked_' num2str( MaxLag) 's.nii']);
disp('Reading volumes...')
Y = spm_read_vols( V);
system( [ 'rm masked_' num2str( MaxLag) 's.nii']);

if ~isempty( range)
	Y = Y(:,:,:, range);
end

MAX = max( abs( Y),[], 4); % excluding noisy voxels with large amplitude
Y = Y.* repmat( MAX<=4, [1 1 1 size( Y,4)]);

mkdir( dir)
cd( dir)

Ysize = size( Y);

ROIimage(  ROIimage < 0.5) = NaN;
sY = Y.* repmat( ROIimage, [1 1 1 Ysize(4)]);

% reshape to space x time
%if TR~=reso
%	Y = resample( reshape( Y, [ Ysize(1)*Ysize(2)*Ysize(3) Ysize(4)])', double( round( TR*100)), double( round( reso*100)));
%	sY = resample( reshape( sY, [ Ysize(1)*Ysize(2)*Ysize(3) Ysize(4)])', double( round( TR*100)),  double( round( reso*100)));
%else
	Y = reshape( Y, [ Ysize(1)*Ysize(2)*Ysize(3) Ysize(4)])';
	sY = reshape( sY, [ Ysize(1)*Ysize(2)*Ysize(3) Ysize(4)])';
%end

Seed = nanmean( sY,2);
clear sY
save InitSeed.mat Seed

if max( Seed)==0, error, end

disp('Extracting sLFO...') %-------------------------------------------

Lag = NaN * ones( 1, Ysize(1)*Ysize(2)*Ysize(3));
Lag( Bmask ==0) = 100;

Lim = 2;
YY = [];
for Sft = Lim:-1:-Lim
	YY = cat( 3, YY, Y( Lim+Sft+1:end-Lim+Sft, :));
end

disp('Calculating correlation...')
XX = repmat( Seed( Lim+1:end-Lim), [ 1 size( YY,2) size( YY,3)]);
CC = sum( XX.*YY, 1)./( sum( XX.*XX, 1).^.5 .* sum( YY.*YY, 1).^.5);
[ R, I] = max( CC, [], 3);
I( R<THR) = 0;
I( isnan( Bmask)) = 0;
Lag( I==3) = 0; % cross-correlogram peaking at zero: lag=0

Seed0 = nanmean( Y( :, I==3),2);
Seeds = Seed0;
Lim = 1;

figure( 5), clf, hold on
set( gcf, 'position', [10 10 1600 400])
set( gca, 'Color', [.5 .5 .5])
plot( Seed0( 1+limit:end-limit), 'w-')
set( gca, 'Xlim', [ 0 length( Seed0( 1+limit:end-limit))])

% Downstream (I==1)
SeedU = Seed0;
SeedD = Seed0;
Downward = Y;
Upward = Y;

for p=1:limit

	YY = [];
	Downward = [ Downward( 2:end,:); nanmean( Downward,1)];
	for Sft = Lim:-1:-Lim
		YY = cat( 3, YY, Downward( Lim+Sft+1:end-Lim+Sft, :));
	end

	XX = repmat( SeedD( Lim+1:end-Lim), [ 1 size( YY,2) size( YY,3)]);
	CC = sum( XX.*YY, 1)./( sum( XX.*XX, 1).^.5 .* sum( YY.*YY, 1).^.5);

	[ R, I] = max( CC, [], 3);
	I( R<THR) = 0;
    if FIXED==0
    	SeedD = mean( Downward( :,I==2),2); % recursive lag tracking
	end
	I( ~isnan( Lag)) = 0;
	Lag( I==2) = -p;
	Seeds = [ Seeds SeedD(:) ];
	disp( [ sprintf( '%+0.1f', -p*TR) ' sec'])

	plot( Seeds( limit+1-p:end-limit-p,end), 'Color',  Cmap( limit+1-p,:)), drawnow

	YY = [];
	Upward = [ nanmean( Upward,1); Upward( 1:end-1,:)];
	for Sft = Lim:-1:-Lim
		YY = cat( 3, YY, Upward( Lim+Sft+1:end-Lim+Sft, :));
	end

	XX = repmat( SeedU( Lim+1:end-Lim), [ 1 size( YY,2) size( YY,3)]);
	CC = sum( XX.*YY, 1)./( sum( XX.*XX, 1).^.5 .* sum( YY.*YY, 1).^.5);

	[ R, I] = max( CC, [], 3);
	I( R<THR) = 0;
    if FIXED==0
		SeedU = mean( Upward( :,I==2),2);
	end
	I( ~isnan( Lag)) = 0;
	Lag( I==2) = p;
	
	drSaveLag( Lag, V, Ysize, TR)

	Seeds = [ SeedU(:) Seeds ];
	disp( [ sprintf( '%+0.1f', p*TR) ' sec'])
	
	plot( Seeds( limit+1+p:end-limit+p,1), 'Color', Cmap( limit+p,:))
	
end
save Seeds.mat Seeds

drawnow
F = getframe( gca);
imwrite( F.cdata, 'Temporal.png')

Lag( Lag==100) = NaN;

Lag1 = reshape( Lag, Ysize(1:3))*TR;
Vout = V(1); Vout.fname = [ 'LagOrig.nii'];
spm_write_vol( Vout, Lag1);

drErode_Lag( 'LagOrig.nii', Bmask, MaxLag);
!mv eLagOrig.nii LagMap.nii
disp('Finished')
cd ..

!rm Data.nii

setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % modified

close( 5)
return


function out=drInsert( in, opt)
out = [];
for p=1:size(in,1)
	[P N E] = fileparts( in(p,:));
	if length( P)>0
		out = [ out; [ P '/' opt N E]];
	else
		out = [ out; [ opt N E]];
	end
		
end
return

function out = drErode_Lag( in, Brain, MaxLag, Mask)
if ~exist( in, 'file'), out = []; return, end

V = spm_vol( in);
Y = spm_read_vols( V);

if nargin>2
	Y( abs( Y)>=MaxLag) = NaN;
end

V.fname = drInsert( in, 'e');

while length( find( isnan( Y(:))))
	nanmask = isnan( Y);
	try
		newY = cat( 4,  circshift( Y, -1, 1), circshift( Y, 1, 1), ...
				circshift( Y, -1, 2), circshift( Y, 1, 2), ...
				circshift( Y, -1, 3), circshift( Y, 1, 3));
	catch
		newY = cat( 4,  circshift( Y, [-1 0 0]), circshift( Y, [ 1 0 0]), ...
				circshift( Y, [ 0 -1 0]), circshift( Y, [ 0  1 0]), ...
				circshift( Y, [ 0 0 -1]), circshift( Y, [ 0 0  1]));
	end
	newY = nanmean( newY, 4);
	Y( nanmask) = newY( nanmask);
	if isnan( max( Y(:)))
		break
	end
	spm_write_vol( V, Y);
end
Brain( Brain==0) = NaN;
Brain( ~isnan( Brain)) = 1;
Y = Y.*Brain;

if nargin>3
	Mask = spm_read_vols( spm_vol( Mask));
	Y( Mask == 0) = NaN;
	Y( isnan( Mask)) = NaN;
end
spm_write_vol( V, Y);
out = V.fname;

function drSaveLag( Lag, V, Ysize, TR)
Lag( Lag>30) = NaN;

Lag1 = reshape( Lag, Ysize(1:3))*TR;
Vout = V(1); Vout.fname = [ 'LagOrig_temp.nii'];
spm_write_vol( Vout, Lag1);


