function drDeperf_hcp_seed( orig_vols, Lag, TR, section, Nruns) 
% Removing BOLD perfusion structure by Toshihiko Aso
%	
%
%	function dir = drDeperf(  vols, Lag, TR, reso, range)
% 
%
%	vols: fMRI 4D file. Can be relative path assuming lag map folder is in the same folder.
%
%	Lag: Full path to the lag map. Typically
%	'/XXX/XXX/Lag_XXs_thrX_XX/LagMap.nii'. In case you had downsampled the
%	fMRI data to say 4 mm voxel size from 2 mm for example,  
%	for quick and safe lag mapping, 
%	the map must be upsampled (= resliced to the original space) beforehand. 
%	If it was done on SPM, it would be like: '/XXX/XXX/Lag_XXs_thrX_XX/rLagMap.nii'
%
%	TR: Repetition time in second.
%
%	Below are non-compulsory options
%	
%	reso: Lag resolution in second (set to 1 by default)
%
%	range: Specify fourth dimension (time) of the Seeds.mat data to be
%		used. This is only necessary when lag map was created using
%		concatenated runs, but the deperfusioning must be done for each
%		run. Causes	error or unfavorable phase shift depending on the combination of TR and reso.

setenv('FSLOUTPUTTYPE', 'NIFTI_GZ')
global Fdir
Fdir = '/mnt/pub/devel/fsl-5.0.9/bin/'

C = clock;

tDir = [ '/tmp/aso/' mfilename num2str( C( 5)) num2str( C( 6)*1000000)];
PWD = pwd;

Lag = absolute_path( Lag);
cd( fileparts( Lag))

V = spm_vol( Lag);
Lag = spm_read_vols( V);
Mask = 10000*double( Lag<100);
Lag( Lag>100) = nanmin( Lag(:));
Lag = round( Lag/TR);

V.fname = 'Mask.nii';
spm_write_vol( V, Mask);
Mask = absolute_path( 'Mask.nii');

load Seeds.mat
MaxLag =  ( size( Seeds, 2)-1)/2; % in TR

Nvols = length( Seeds)/Nruns;
if Nvols ~= dr4Dlen( orig_vols), error, end

range = ( section-1)*Nvols+1:section*Nvols;
Seeds = Seeds( range, :);

Motodata = [];
LL = ( -MaxLag:MaxLag);
for p = 1:MaxLag
	temp = Seeds( :,1+size( Seeds, 2)-p);
	Motodata = [ Motodata [ zeros( MaxLag-p+1,1); temp( 1:end-MaxLag-1+p)]];
end
temp = Seeds( :,1+MaxLag);
Motodata = [ Motodata temp];
for p = 1:MaxLag
	temp = Seeds( :,1+MaxLag-p);
	Motodata = [ Motodata [ temp( p+1:end); zeros( p,1)]];
end

L = LL;	
	
ULfreq = 1/(2*MaxLag*TR);

mkdir( tDir)
cd( tDir)

disp('Filtering...')
S = system( [ 'niimath ' orig_vols ' -bptf ' ...
	num2str( 1/( 0.008*2.35*TR)) ' -1 ' pwd '/data_lowcut']); 

% - - - - - - - - - - - - - - - - - - - - - - - - - 

V.fname = 'MaskDeperf.nii';
disp('Cleaning images...')

Nseed = size( Seeds, 2);

save( [ PWD '/sLFO.mat'], 'Motodata')

S = system( [ 'niimath ' orig_vols ' -Tmean out.nii.gz' ]);
if S, save all.mat, error, end

figure( 19)
clf, subplot( 2,1,1), drPlotRainbow( Seeds), subplot( 2,1,2), drPlotRainbow( Motodata)
set( gca, 'Position', [ 0 0 1800 600])

for p = 1:Nseed
	Mask = ( Lag==L( p));

	if ~isempty( find( Mask(:)))
		spm_write_vol( V, Mask);
		
		if system( [ 'niimath MaskDeperf -mas out Mask']), error, end

		Reg = Motodata( :,p)./kaiser( Nvols, 2);

		dlmwrite( 'Reg.txt', Reg, 'precision','%.4f')
		
		S = system([ Fdir 'fsl_regfilt --in=data_lowcut '...
				' --out=rfmasked' num2str( p) ' --design=Reg.txt -f "1"' ]);
		if S, save all.mat, error, end
		S = system( [ 'niimath rfmasked' num2str( p) ' -mas Mask rfmasked' num2str( p) ]);
		if S, save all.mat, error, end
		
	end
end

%if system( [ 'cp ' Ref ' out.nii.gz']), error, end
		
disp('Merging images...')
	for p =  1: size( Seeds, 2)
		if exist( ['rfmasked' num2str( p) '.nii.gz'], 'file')
			S = system( [ 'niimath out -add rfmasked' num2str( p) '.nii.gz out.nii.gz']); if S, error, end
		end
		disp( [ 'rfmasked' num2str( p) '.nii.gz'])
	end

[ ~,newvols, ext] = fileparts( orig_vols);

if system( [ 'mv out.nii.gz ' PWD '/dep_' newvols ext]), error, end
%gunzip( [ 'dep_' newvols ext '.gz'], PWD)
cd( PWD)
system( [ 'rm -rf ' tDir])
return

function out = dr4Dlen( in)
global Fdir
if nargin<1
	in = spm_select( 1,'any')
end
[ S, out] = system( [ Fdir 'fslnvols "' in '"']);
out = str2double( out);
if isnan( out), save all.mat, end


function out = absolute_path( in)

[ d, f, e] = fileparts( in);
[ ~, pinfo] = fileattrib( d);

out = fullfile( pinfo.Name, strcat([ f, e]));
