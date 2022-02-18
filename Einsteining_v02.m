function Einsteining_v02( Sdir, TR, Nvols, MaxLag, MinR, Fixed, Sm, onlyLag)
% Pipeline for removal of perfusion lag structure in HCP-style fMRI data
%
% All directories named "REST" are treated as source and concatenated
%
% Example: 
% For input data BOLD_REST1_AP, BOLD_REST1_PA, BOLD_REST2_AP, BOLD_REST2_PA 
% output directories will be dep_BOLD_REST1_AP, dep_BOLD_REST1_PA, dep_BOLD_REST2_AP, dep_BOLD_REST2_PA 
%
% Example usage:
% Einsteining_v01( Sdir, 0.80, 375, 9, 0.3, 0, 8, 0)
%
% Sdir: Subject directory
% TR = 0.8s
% Number of volumes per run = 375
% MaxLag = 8 -> lag mapping up to -8TR / +8TR
% MinR = 0.2 -> Cross-correlogram peak < 0.2 will not be used
% Fixed = 0  -> Recursive (flexible) LFO instead of Fixed-sLFO algorithm
% Spatial smoothing at 8mm FWHM
% onlyLag = 0 -> Lag mapping + deperfusioning
%
% Modify following variables:
% Fdir
% 

Sdir = deblank( Sdir);

global Fdir
Fdir = '/mnt/pub/devel/fsl-5.0.9/bin/'

setenv( 'FSLDIR', Fdir( 1:end-4))
%if system( [ 'source /mnt/pub/devel/fsl-5.0.9/etc/fslconf/fsl.sh']), error, end

if nargin < 8
	onlyLag == 0;
end

%try
	cd( Sdir)
	cd MNINonLinear/Results
	disp( Sdir)
	
	[ ~, in] = spm_select( 'FPList', pwd, '.*');
	Runs = [];
	for p=1:size( in, 1)
		if strfind( in( p,:), 'BOLD')*~length( strfind( in( p,:), 'dep_'))
			Rname = drFoldername( in( p,:))
			Vol = spm_select( 'FPList', deblank( in( p,:)), [ '^' Rname '.nii.gz$']);

			if isempty( Vol)
				Vol = spm_select( 'FPList', deblank( in( p,:)), [ '^' Rname '.nii$']);
				if isempty( Vol)
					error( 'no data')
				else
					setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
					if system( [ Fdir 'fslchfiletype NIFTI_GZ ' Vol]), error, end
					Vol = spm_select( 'FPList', deblank( in( p,:)), [ '^' Rname '.nii.gz$']);
					Runs = [ Runs; { Vol}];
					
				end
			else
				dr4Dlen( Vol)
				if dr4Dlen( Vol)==Nvols
					Runs = [ Runs; { Vol}];
				end
			end
		end
	end
	if isempty( Runs), error( 'no data'), end
	Ref = spm_select( 'FPList', fileparts( Runs{1}), '.*SBRef.nii.gz');
%catch
%	return
%end

drCellMat( Runs)

if ~exist( 'Lag_concat', 'file')
	mkdir Lag_concat
end
cd Lag_concat

save Runs.mat Runs Ref

%!rm z*.nii

if ~exist( [ 'REST' num2str( length( Runs)) 'run.nii.gz'], 'file')
	
	for p=1:length( Runs)
		fname = drFilename( Runs{ p});
		zfile = spm_select( 'FPList', pwd, [ '^z' fname( 1:end-7) '.*']);
		if isempty( zfile)
			setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ')
			disp( ['Reducing resolution Run #' num2str( p)])
%			try
%				if system( ['/usr/local/fsl/bin/fslmaths ' Runs{ p} ' -subsamp2offc z' fname( 1:end-7)]), error, end
%			catch				
				if system( [ 'niimath ' Runs{ p} ' -subsamp2offc z' fname( 1:end-7)]), error, end
%			end
	%		Runs{ p} = absolute_path( nii_scale_dims( Runs{ p}, [1 1 1]*.5))
			zfile = spm_select( 'FPList', pwd, [ '^z' fname( 1:end-7) '.*']);
		end
		Runs{ p} = zfile( end,:);
	end
	out = drMerge4D_LagmapGZ_fsl5( [ 'REST' num2str( length( Runs)) 'run'], TR, Runs)
else
	out = spm_select( 'FPList', pwd, ['^REST' num2str( length( Runs)) '.*.nii.gz']);
end

Lag = drLag4Drev4_hcp_niimath( [ 'cat' num2str( length( Runs)) ], num2str( TR), out, num2str( MaxLag), num2str( MinR), num2str( Fixed), num2str( Sm));

if onlyLag, return, end

if exist( [ Lag '/rLagMap.nii'], 'file')
	
else
	C = clock;
	tDir = [ '/tmp/' mfilename num2str( C( 5)) num2str( C( 6)*1000000)];

	gunzip( Ref, tDir)
	Ref = spm_select( 'FPList', tDir, '.*SBRef.nii');
	if isempty( Ref), error, end
	rLag = drReslice_Lag( Ref, [ Lag '/LagMap.nii']);
	movefile( rLag, 'rLagMap.nii')
	copyfile( 'rLagMap.nii', Lag)
	system( [ 'rm -rf ' tDir])
end
Lag = [ Lag '/rLagMap.nii'];

for r=1:size( Runs, 1)
	Deps = [ 'dep_' drFilename( Runs{ r}) ];
	if ~exist( Deps, 'file')
		
	%	drDeperf_hcp_seed_fixit( Deps, Runs{ r}, Lag, TR, r, length( Runs))
		drDeperf_hcp_seed_niimath( Runs{ r}, Lag, TR, r, length( Runs))
	end
	
		if r==length( Runs)
			cd ..
			for rr=1:length( Runs)
				in1 = drFilename( Runs{ rr});
				in = in1( 1:end-7);
				in1 = in1( 1:end-4);
				if exist( [ 'dep_' in1], 'file')
					if system( [ 'rm -r dep_' in1]), error, end
				end
				mkdir( [ 'dep_' in])
				cd( [ 'dep_' in])
	
				system( [ 'rm Movement_Regressors.txt']);
				if system( [ 'ln -s ../' in '/Movement_Regressors.txt .']), error, end;
				system( [ 'rm dep_' in '_SBRef.nii.gz']);
				if system( [ 'ln -s ../' in '/' in '_SBRef.nii.gz dep_' in '_SBRef.nii.gz']), error, end;
				system( [ 'rm dep_' in '.nii.gz']);
				if system( [ 'ln -s ../Lag_concat/dep_' in '.nii.gz dep_' in '.nii.gz']), error, end;
				cd ..
			end
		end

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
global Fdir
if nargin<1
	in = spm_select( 1,'any')
end
[ S, out] = system( [ Fdir 'fslnvols "' in '"']);
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

function [ out] = drFilename( in)

[ P out ext ] = fileparts( in);
out = [ out ext];

function out = drMerge4D_LagmapGZ_fsl5( name, TR, vols)

%  Dr. Merge4D by Toshihiko Aso
global Fdir

StartDir = pwd;

C = clock;
WD = squeeze( deblank( num2str( C)));
WD = [ '/tmp/aso/' WD( ~isspace( WD))];
mkdir( WD)

if isempty( vols), error, end

if ~iscell( vols)
	vols = drMatCell( vols);
end

%if ~strcmp( vols{ 1}( end-2:end), '.gz')
%	for n = 1:length( vols)
%		vols{n} = [vols{n} '.nii'];
%	end
%end

cd( WD)

for lp = 1:size( vols,1)
	test = spm_select( 'List', pwd, '^temp.*.nii');
	if ~isempty( test)
		!rm temp*.nii
	end

	S = system( [ 'niimath ' vols{ lp} ' -nan Sm']);
	setenv( 'FSLOUTPUTTYPE', 'NIFTI')
	setenv( 'FSLDIR', Fdir( 1:end-4))
	if system( [ Fdir 'fslchfiletype NIFTI Sm']), error, end
	
	N = dr4Dlen( 'Sm.nii');
	V = spm_vol( 'Sm.nii,1');
	Y = spm_read_vols( V);
	W = kaiser( N, 2);
	for n=1:N
		V.fname = [ 'temp' sprintf( '%04g', n) '.nii'];
		Y = Y*0 + W( n);
		spm_write_vol( V, Y);
	end
	disp( '...windowing')
	if system( [ Fdir 'fslmerge -t Window.nii temp*.nii']), error, end

	hyojunka( 'Sm.nii', TR)
	
	if system( [ 'niimath Sm.nii -mul Window.nii Sm.nii']), error, end

	if lp==1
		if system( [ 'mv Mean.nii Sum.nii']), error, end
		if system( [ 'mv SD.nii SDall.nii']), error, end
	else
		if system( [ 'niimath Sum.nii -add Mean.nii Sum.nii']), error, end
		if system( [ 'niimath SDall.nii -add SD.nii SDall.nii']), error, end
	end
	if system( [ 'mv Sm.nii Seg' num2str( lp) '.nii']), error, end		
end
!rm temp*.nii

if system( [ 'niimath Sum.nii -div ' num2str( lp) ' SumMean.nii']), error, end
if system( [ 'niimath SDall.nii -div ' num2str( lp) ' SDallMean.nii']), error, end

if system( [ Fdir 'fslmerge -t All.nii Seg*.nii']), error, end

S = system( [ 'niimath All.nii -mul SDallMean -add SumMean ' name '.nii']);
if ~isempty( spm_select( 'List', pwd, '^Seg.*.nii'))
	!rm Seg*.nii
end

setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ')
S = system( [ Fdir 'fslchfiletype NIFTI_GZ ' name]);

cd( StartDir)
if system( [ 'mv ' WD '/' name '.nii.gz ./' ]), error, end
if system( [ 'rm -rf ' WD ]), error, end

out = [ pwd '/' name '.nii.gz'];

return

function hyojunka( in, TR)
disp( '...Hyojunka')
global Fdir

S = system( [ 'niimath ' in ' -Tmean Mean.nii']);
if S, error, end

S = system( [ 'niimath ' in ' -bptf ' num2str( 1/(0.005*2.35*TR)) ' -1 ' in ]);
if S, error, end

S = system( [ 'niimath ' in ' -Tstd SD.nii']);
if S, error, end
S = system( [ 'niimath ' in ' -div SD.nii ' in]);
if S, error, end

return

function out = drReslice_Lag( ref, source, INT)

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

