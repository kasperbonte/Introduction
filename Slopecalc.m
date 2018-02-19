clearvars
clc
close all
%nested exercise
load('sierrademml.mat' )
imagesc(sierrademml)

%%
% for loop slope calculation
Size = size(sierrademml)
row = Size (1)
col = Size(2)
Slope = NaN(Size)
r = 30
for i = 2:(Size(1)-1)
    for j = 2:(Size(2)-1)
        Slopex = (sierrademml(i,j+1)-sierrademml(i,j-1))/(2*r);
        Slopey = (sierrademml(i-1,j)-sierrademml(i+1,j))/(2*r);
        Slope(i,j) = sqrt((Slopex)^2+ (Slopey)^2);
    end
end
imagesc(Slope)
%% 
% indexes calculate slope
Size = size(sierrademml)
row = Size (1)
col = Size(2)
Slope = NaN(Size)
r = 30
Slopex = NaN(Size)
Slopey = NaN(Size)
Slope = NaN(Size)
Slopex(2:row-1,2:col-1) = (sierrademml(2:row-1,3:col)-sierrademml(2:row-1,1:col-2))/(2*r)
Slopey(2:row-1,2:col-1) = (sierrademml(1:row-2,2:col-1)-sierrademml(3:row,2:col-1))/(2*r)
Slope(2:row-1,2:col-1) = sqrt((Slopex(2:row-1,2:col-1).^2)+(Slopey(2:row-1,2:col-1).^2))
imagesc(Slope)