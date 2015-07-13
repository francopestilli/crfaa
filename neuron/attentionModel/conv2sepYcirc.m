function result = conv2sepYcirc(im,rowfilt,colfilt)
% conv2SepYcirc: Separable, convolution using upConv, with circular
% convolution for columns and zero-padding for rows.
% 
%      result=conv2sepYcirc(im,rowfilt,colfilt)
%
%      im - input image.
%      rowfilt - 1d filter applied to the rows
%      colfilt - 1d filter applied to the cols
%
% Example: foo=cconv2sep(im,[1 4 6 4 1],[-1 0 1]);
%
% DJH '96

rowfilt=rowfilt(:)';
colfilt=colfilt(:);

tmp = upConv(im,rowfilt,'zero',[1,1]);
result = upConv(tmp,colfilt,'circular',[1,1]);
return

%%% Debug
im=zeros(7,7);
im(4,4)=1;
filter = [1,2,4,2,1];
filter=filter/sum(filter);

res1=conv2sepYcirc(im,filter,filter)
res2=cconv2(im,filter'*filter);
mse(res1,res2)
