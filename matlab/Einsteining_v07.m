function Einsteining_v07( Runs, TR, Nvols, MaxLag, MinR, Fixed, Sm, onlyLag)
% Pipeline for removal of perfusion lag structure in HCP-style fMRI data
% WITH two-stage scrubbing
%
% All directories named "REST" are treated as source and concatenated
%
%
% Example usage:
% Einsteining_v07( { '/data/subject1/MNINonLinear/Results/rfMRI_REST1_PA';...
% '/data/subject1/MNINonLinear/Results/rfMRI_REST1_PA/rfMRI_REST1_PA.nii.gz'}, '0.72', '1200', '9', '0.2', '0', '8', '0')
%
% Runs = Cell containing target files
% TR = in seconds
% Number of volumes per run
% MaxLag = 8 -> lag mapping up to -8TR / +8TR
% MinR = 0.2 -> Cross-correlogram peak < 0.2 will not be used
% Fixed = 0  -> Recursive (flexible) sLFO instead of Fixed-sLFO algorithm
% Spatial smoothing at 8mm FWHM
% onlyLag = 0 -> Lag mapping + deperfusioning

if ischar( TR), TR = str2double( TR); end
if ischar( Nvols), Nvols = str2double( Nvols); end
if ischar( MaxLag), MaxLag = str2double( MaxLag); end
if ischar( MinR), MinR = str2double( MinR); end
if ischar( Fixed), Fixed = str2double( Fixed); end
if ischar( Sm), Sm = str2double( Sm); end
if ischar( onlyLag), onlyLag = str2double( onlyLag); end

[S R] = system( 'niimath');
if S==0
	fslmath = 'niimath';
else
	fslmath = '/usr/local/fsl-5.0.9/bin/fslmaths'
end

if nargin < 8
	onlyLag = 0;
end

cd( fileparts( Runs{ 1}))
cd ../../../
Sdir = pwd;
cd MNINonLinear/Results

Ref = drDir( [ fileparts( Runs{1}) '/*SBRef.nii.gz']);
if isempty( Ref), error, end

if exist( 'Lag_concat_scrub', 'file')
	!mv Lag_concat_scrub old_Lag_concat_scrub
end
mkdir Lag_concat_scrub
cd Lag_concat_scrub

save Runs.mat Runs Ref

if ~exist( [ 'REST' num2str( length( Runs)) 'run.nii.gz'], 'file')
	zRuns = Runs;
	
	for p= 1:length( Runs)
		[fn ext] = drFilename( Runs{ p});
		fname = [ fn ext];
		fn = fname( 1:end-7);
%		zfile = spm_select( 'FPList', pwd, [ '^mreg_z' fn '.*']);
%		if size( zfile,1) == 0

			setenv( 'FSLOUTPUTTYPE', 'NIFTI')
			disp( ['Reducing resolution Run #' num2str( p)])
			if system( [ fslmath ' ' Runs{ p} ' -subsamp2offc z' fn]), error, end
			
			Y = spm_read_vols( spm_vol( ['z' fn '.nii']));
			Siz = size( Y);
			% DVARS (fix for issue #2): RMS of the frame-to-frame difference of the
			% raw intensity within the brain (voxels that are never zero), in percent
			% of the mean brain signal. A frame (and the one before it) is a spike
			% when DVARS exceeds SpikeThr x median( DVARS).
			SpikeThr = 1.5;
			brain = all( Y~=0, 4) & mean( Y, 4)>0;
			Y2 = reshape( Y, [], Siz(4));
			Y2 = Y2( brain(:), :);
			dvars = [ 0; sqrt( mean( diff( Y2, 1, 2).^2, 1))'] / mean( Y2(:)) * 100;
			clear Y2
			Spike = find( dvars > SpikeThr*median( dvars));
			Spike = [ Spike(:); Spike(:)-1];
			Spike = unique( Spike( Spike>0));

			Sreg = zeros( Siz(4),length( Spike));
			RegN = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19';
			for q=1:length( Spike)
				Sreg( Spike( q),q) = 1;
				RegN = [ RegN ',' num2str( 19+q)]; % fix for issue #2 (was str2double( RegN( end))+q)
			end
			disp( [ 'Reading ' Sdir '/MNINonLinear/Results/' fn '/Movement_Regressors.txt'])
			[ FD MP] = drFD( [ 	Sdir '/MNINonLinear/Results/' fn '/Movement_Regressors.txt']);
			DP = diff( MP(:,1:6));

			out = [ MP(:,1:6) [ zeros( 1,6);DP] [ DP; zeros( 1,6)] FD Sreg];
		%	writematrix( out, ['regout' num2str( p) '.txt'],'Delimiter','tab')
			dlmwrite( ['regout' num2str( p) '.txt'], out, 'delimiter','\t','precision','%.6f')

			setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ')
			S = ['export FSLOUTPUTTYPE=NIFTI_GZ; fsl_regfilt -i z' fn '.nii -d regout' num2str( p) '.txt -f "' RegN '" -v -o mreg_z' fn ];
			if system( S), error, end
			zfile = drDir( [ 'mreg_z' fn '*']);
%		end
		zRuns{ p} = zfile( end,:);
	end
	
	out = drMerge4D( [ 'REST' num2str( length( Runs)) 'run'], TR, zRuns, fslmath)
else
	out = drDir( ['REST' num2str( length( Runs)) '*.nii.gz']);
end

%!rm zrfMRI*.nii

disp( [ 'drLag4Drev7( ''cat' num2str( length( Runs)) ''',' num2str( TR) ',''' out ''',' num2str( MaxLag) ',' num2str( MinR) ',' num2str( Fixed) ',' num2str( Sm) ')'])

Lag = drLag4Drev7( [ 'cat' num2str( length( Runs)) ], num2str( TR), out, num2str( MaxLag), num2str( MinR), '1', num2str( Sm));

if onlyLag, return, end


C = clock;
tDir = [ '/tmp/Ein' mfilename num2str( C( 5)) num2str( C( 6)*1000000)];

gunzip( Ref, tDir)
Ref = drDir( [ tDir '/*SBRef.nii']);
% Ref = spm_select( 'FPList', tDir, '.*SBRef.nii');
if isempty( Ref), error, end
rLag = drReslice_Lag( Ref, [ Lag '/LagMap.nii']);
movefile( rLag, 'rLagMap.nii')
copyfile( 'rLagMap.nii', Lag)
system( [ 'rm -rf ' tDir])

Lag = [ Lag '/rLagMap.nii'];

for r=1:size( Runs, 1)
	Deps = [ drFilename( Runs{ r}) '_dep' ];
	disp( [ 'drDeperf_hcp_seed( ''' Runs{ r} ', ' Lag ', ' num2str( TR) ', ' num2str( r) ', ' num2str( length( Runs)) ', ' fslmath])

	drDeperf_hcp_seed( Runs{ r}, Lag, TR, r, length( Runs), fslmath)
end

cd ..	

for rr=1:length( Runs)
	in = drFilename( Runs{ rr});
	in = in( 1:end-4);
	if ~exist( [ in '_dep'], 'file')
		mkdir( [ in '_dep'])
	end
	cd( [ in '_dep'])
	
	if exist( 'Movement_Regressors.txt', 'file')
		system( [ 'unlink Movement_Regressors.txt']);
	end
	if system( [ 'ln -s ../' in '/Movement_Regressors.txt .']); error, end

	if exist( [ in '_dep_SBRef.nii.gz'], 'file')
		system( [ 'unlink ' in '_dep_SBRef.nii.gz']);
	end
	if system( [ 'ln -s ../' in '/' in '_SBRef.nii.gz ' in '_dep_SBRef.nii.gz']); error, end
	
%	if exist( [ 'dep_' in '.nii.gz'], 'file')
		system( [ 'unlink ' in '_dep.nii.gz']);
%	end
	if system( [ 'ln -s ../Lag_concat_scrub/' in '_dep.nii.gz ' in '_dep.nii.gz']), error, end
	cd ..
end


function out = drFoldername( in, N)
if size( in,1)>1
	out = [];
	for p=1:size( in,1)
		out = [ out; { drFoldername( deblank( in(p,:)))}];
	end
	out = drCellMat( out);
	return
end
if in( end)=='/'
	in = in( 1:end-1);
end
%in = absolute_path( deblank( in));

if nargin==1
	N=1;
end

slash = strfind( in, '/');
if isempty( slash)
	slash = strfind( in, '\');
end
%out = deblank( in( slash( end-1)+1:slash( end)-1));
if N>1
	out = deblank( in( slash( end+1-N)+1:min( end, slash( end+2-N)-1)));
else
	out = deblank( in( slash( end+1-N)+1:end));
end


function out = dr4Dlen( in)
[ S, out] = system( [ 'fslnvols "' in '"']);
if S, error, end
out = str2double( out);
if isnan( out), save all.mat, end


function out = drCellMat( in)
out = [];
if ~iscell( in), out = in; return, end
Nrow = size( in, 1);

Ncol = 0;
for p=1:Nrow
	Ncol = max( Ncol, size( in{p},2));
end

for p=1:Nrow
	out = [ out; [ in{p} repmat( ' ', 1, Ncol-size( in{p},2))]];
end


function [ out ext] = drFilename( in)
[ P out ext ] = fileparts( in);


function out = drMerge4D( name, TR, vols, fslmath)
%  Dr. Merge4D by Toshihiko Aso

StartDir = pwd;

C = clock;
WD = squeeze( deblank( num2str( C)));
WD = [ '/tmp/' WD( ~isspace( WD))];
mkdir( WD)

if isempty( vols), error, end

if ~iscell( vols)
	vols = drMatCell( vols);
end

cd( WD)

for lp = 1:size( vols,1)
	test = drDir( 'temp*.nii');
	if ~isempty( test)
		!rm temp*.nii
	end

	S = system( [ fslmath ' ' vols{ lp} ' -nan Sm']);

	setenv( 'FSLOUTPUTTYPE', 'NIFTI')
	if system( [ 'fslchfiletype NIFTI Sm']), error, end
	
	N = dr4Dlen( 'Sm.nii');
	V = spm_vol( 'Sm.nii,1');
	Y = spm_read_vols( V);
%	W = kbdwin( round( ( N*5/3)/2)*2, 100); % reduce deflection at both ends
%	W = W( round( N/3):round( N/3)+N-1);
	W = [ 0:1/((N-1)*.05):1 ones( 1, ceil( N*.9)) 1:-1/((N-1)*.05):0 0];
	for n=1:N
		V.fname = [ 'temp' sprintf( '%04g', n) '.nii'];
		Y = Y*0 + W( n);
		spm_write_vol( V, Y);
	end
	disp( '...windowing')
	if system( [ 'fslmerge -t Window.nii temp*.nii']), error, end

	disp( '...detrending')
	S = system( [ fslmath ' Sm.nii -Tmean Mean.nii']);
	if S, error, end

	S = system( [ fslmath ' Sm.nii -bptf ' num2str( 1/(0.008*2.35*TR)) ' -1 -mul Window.nii Sm.nii' ]);
	if S, error, end	
	
	if lp==1
		if system( [ 'mv Mean.nii Sum.nii']), error, end
	else
		if system( [ fslmath ' Sum.nii -add Mean.nii Sum.nii']), error, end
	end
	if system( [ 'mv Sm.nii Seg' num2str( lp) '.nii']), error, end		
end
!rm temp*.nii

if system( [ fslmath ' Sum.nii -div ' num2str( lp) ' SumMean.nii']), error, end

if system( [ 'fslmerge -t All.nii Seg*.nii']), error, end

S = system( [ fslmath ' All.nii -add SumMean ' name '.nii']);
if ~isempty( drDir( 'Seg.*.nii'))
	!rm Seg*.nii
end

setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ')
S = system( [ 'fslchfiletype NIFTI_GZ ' name]);

cd( StartDir)
if system( [ 'mv ' WD '/' name '.nii.gz ./' ]), error, end
if system( [ 'rm -rf ' WD ]), error, end

out = [ pwd '/' name '.nii.gz'];

return


function out = drReslice_Lag( ref, source, INT)
spm('defaults', 'FMRI');
if nargin<3
	INT = 0;
end

if iscell( source)
	source = drCellMat( source);
end

if size( source, 1)>1
	for p=1:size( source, 1)
		drReslice_Lag( ref, deblank( source( p,:)), INT)
	end
	return
end

if isstruct( ref)
	ref = ref.fname;
end
ref = absolute_path( ref);
if isstruct( source)
	source = source.fname;
end

source = absolute_path( source);
[~,N,E] = fileparts( source);

v1 = spm_vol( ref);
v1 = v1(1);
v2 = spm_vol( source);
v2 = v2(1);

if isequal( v1.mat, v2.dim)==0

	Y = spm_read_vols( v2);
	tempV = v2;
	tempV.fname = 'temp_reslice_lag.nii';
	Y( isnan( Y)) = 10000;
	spm_write_vol( tempV, Y);

	in{1}.spm.spatial.coreg.write.ref = { ref};
	in{1}.spm.spatial.coreg.write.source = { tempV.fname};
	in{1}.spm.spatial.coreg.write.roptions.interp = INT;
	in{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
	in{1}.spm.spatial.coreg.write.roptions.mask = 0;
	in{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
	matlabbatch = in;
	save Reslice.mat matlabbatch
	spm_jobman( 'run', in)

	out = [ pwd '/r' tempV.fname];

else
	try
		copyfile( source, pwd);
	catch
%		disp('Could not move file')
	end
	out = [ pwd '/' N E];
	
end
	
return


function out = absolute_path( in)

[ d, f, e] = fileparts( in);
[ ~, pinfo] = fileattrib( d);

out = fullfile( pinfo.Name, strcat([ f, e]));

function [ FD in] = drFD( filename)
in = textread(filename);
in = in( :,1:12);
temp = in(:,4:6)*180/pi;
temp = temp * (2*50*pi/360);
FD = [ 0; sum( abs( diff( [in( :, 1:3) temp])), 2)];


function drDeperf_hcp_seed( orig_vols, Lag, TR, section, Nruns, fslmath) 
% Removing BOLD perfusion structure by Toshihiko Aso
%	
%	function dir = drDeperf_hcp_seed(  vols, Lag, TR, section, Nruns, fslmath)
% 
%	vols: original fMRI 4D file names. Can be relative path assuming lag map folder is in the same folder.
%
%	Lag: Full path to the lag map. Typically
%	'/XXX/XXX/Lag_XXs_thrX_XX/LagMap.nii'. In case you had downsampled the
%	fMRI data to say 4 mm voxel size from 2 mm for example,  
%	for quick and safe lag mapping, 
%	the map must be upsampled (= resliced to the original space) beforehand (which is implemented in Einsteining). 
%
%	TR: Repetition time in second.
%	section: Run number 1-Nruns
%	fslmath: 	
%	reso: Lag resolution in second (set to 1 by default)
%
%	range: Specify fourth dimension (time) of the Seeds.mat data to be
%		used. This is only necessary when lag map was created using
%		concatenated runs, but the deperfusioning must be done for each
%		run. Causes	error or unfavorable phase shift depending on the combination of TR and reso.

C = clock;

tDir = [ '/tmp/' mfilename num2str( C( 5)) num2str( C( 6)*1000000)];
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

mkdir( tDir)
cd( tDir)

disp('Filtering...')
S = system( [ fslmath ' ' orig_vols ' -bptf ' ...
	num2str( 1/( 0.008*2.35*TR)) ' -1 ' pwd '/data_lowcut']); 

% - - - - - - - - - - - - - - - - - - - - - - - - - 

V.fname = 'MaskDeperf.nii';
disp('Cleaning images...')

Nseed = size( Seeds, 2);

save( [ PWD '/sLFO.mat'], 'Motodata')

S = system( [ 'export FSLOUTPUTTYPE=NIFTI_GZ;' fslmath ' ' orig_vols ' -Tmean out.nii.gz' ]);
if S, save all.mat, error, end

%figure( 19)
%clf, subplot( 2,1,1), drPlotRainbow( Seeds), subplot( 2,1,2), drPlotRainbow( Motodata)
%set( gca, 'Position', [ 0 0 1800 600])

for p = 1:Nseed
	Mask = ( Lag==L( p));

	if ~isempty( find( Mask(:)))
		spm_write_vol( V, Mask);
		
		if system( [ fslmath ' MaskDeperf -mas out Mask']), error, end

		N = size( Motodata, 1);
%		W = kbdwin( round( ( N*5/3)/2)*2, 100); % reduce deflection at both ends
%		W = W( round( N/3):round( N/3)+N-1);
		
		Reg = Motodata( :,p); % kaiser( Nvols, 2);

		dlmwrite( 'Reg.txt', Reg, 'precision','%.4f')
		
		S = system([ 'export FSLOUTPUTTYPE=NIFTI_GZ; fsl_regfilt --in=data_lowcut '...
				' --out=rfmasked' num2str( p) ' --design=Reg.txt -f "1"' ]);
		if S, save all.mat, error, end
		S = system( [ 'export FSLOUTPUTTYPE=NIFTI_GZ;' fslmath ' rfmasked' num2str( p) ' -mas Mask rfmasked' num2str( p) ]);
		if S, save all.mat, error, end
	end
end

if system( 'cp out.nii.gz Tmean.nii.gz'), error, end

disp('Merging images...')
for p =  1: size( Seeds, 2)
	if exist( ['rfmasked' num2str( p) '.nii.gz'], 'file')
		S = system( [ 'export FSLOUTPUTTYPE=NIFTI_GZ;' fslmath ' out -add rfmasked' num2str( p) ' out.nii.gz']); if S, error, end
	end
end

%N = dr4Dlen( 'out.nii.gz');
%Y = spm_read_vols( V);
%W = max( 0.75, [ 0:1/((N-1)*.05):1 ones( 1, ceil( N*.9)) 1:-1/((N-1)*.05):0 0]);
%for n=1:N
%	V.fname = [ 'temp' sprintf( '%04g', n) '.nii'];
%	Y = Y*0 + W( n);
%	spm_write_vol( V, Y);
%end
%disp( 'Windowing...')
%if system( [ 'fslmerge -t Window.nii temp*.nii']), error, end

[ ~,newvols, ext] = fileparts( orig_vols);
%if system( [ fslmath ' out.nii.gz -sub Tmean.nii.gz -mul Window.nii -add Tmean.nii.gz ' PWD '/' newvols( 1:end-4) '_dep.nii.gz' ]), error, end

if system( [ 'export FSLOUTPUTTYPE=NIFTI_GZ;' fslmath ' out.nii.gz ' PWD '/' newvols( 1:end-4) '_dep.nii.gz' ]), error, end

cd( PWD)
system( [ 'rm -rf ' tDir])
return

function out = drDir( in)
temp = dir( in);
if isempty( temp)
	out = [];
else
	out = [ temp( 1).folder '/' temp( 1).name];
end