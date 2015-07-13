% Creates all figures for the paper

% clear all
% close all

quickcheck = 1;
if quickcheck == 1
  numContrasts = 9;
else
  numContrasts = 25;
end
if quickcheck == 1
  numOrientations = 9;
else
  numOrientations = 25;
end

tic
Figure2A
Figure2B              
Figure3C
Figure3F
Figure4C
Figure4E
Figure5C
Figure6C
Figure7C
toc
