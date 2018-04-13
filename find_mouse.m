function [final_image,bwObj,obj_props] = findMouse(BW_track,bkg_image,actual_image,mouse_area)

[rows,cols] = size(bkg_image);

actual_image(BW_track == 0) = 0;
bkg_image(BW_track == 0) = 0;
%%% isolate mouse by subtracting frame from background image
sub_image = abs(bkg_image-actual_image);
%%% Binarize subtracted image
thresh = graythresh(sub_image);
sub_image = im2bw(sub_image,thresh);

final_image = imfill(sub_image,'holes');

%%% Get rid of smaller objects that are not the mouse
final_image = bwareaopen(sub_image,50,4);

%%% Fill in potential mouse objects
se = [1,1,1,1,1,1,1];
se_v = [1;1;1;1;1;1;1];
final_image = imdilate(final_image,se_v);
final_image = imerode(final_image,se_v);
final_image = imdilate(final_image,se);
final_image = imerode(final_image,se);

bwObj = bwconncomp(final_image);
       
%%% Remove objects that are too small
for i=1:bwObj.NumObjects
    if(length(bwObj.PixelIdxList{i}) < mouse_area/3)
        final_image(bwObj.PixelIdxList{i}) = 0;
    end
end

bwObj = bwconncomp(final_image);
obj_props = regionprops(bwObj,'area','centroid','minoraxislength','majoraxislength','eccentricity','solidity','area');


