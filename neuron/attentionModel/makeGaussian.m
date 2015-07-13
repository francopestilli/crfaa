function gaussian = makeGaussian(space,center,width,height)
%
% gaussian = makeGaussian(space,center,width,[height])
%
% This is a function creates gaussian centered at "center", 
% over the values defined by vector "space".  
%
% width is the standard deviation of the gaussian
%
% height, if specified is the height of the peak of the gaussian.
% Otherwise, it is scaled to unit volume

plot_figure = 0; % if 1, then plot the two kernels.

gaussian = normpdf(space,center,width); 

if exist('height','var')
  gaussian = height * width * sqrt(2*pi) * gaussian;
end

if plot_figure == 1
    figure; plot (space,gaussian)
end
