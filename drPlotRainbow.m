function drPlotRainbow( Y)

Nseed = size( Y,2);
CM = jet( Nseed);

for p=1:Nseed
	plot( Y(:, p), 'Color', CM( p,:)), hold on
end
set( gca, 'Color', [ .4 .4 .4])



