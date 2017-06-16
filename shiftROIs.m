function shiftCellList = shiftROIs(im1, im2, CellList)
%Shift ROI pixel indexes of CellList from im1 image space to im2 image space.
%This does a rigid shift for the ROIs

%Code for determining image shifts
[ydim,xdim] = size(im2);
imxcorrs = ifftshift(ifft2(fft2(im1).*conj(fft2(im2))));
[y,x] = find(imxcorrs == max(max(imxcorrs)));

%Code for shifting ROIs
centerind = sub2ind([ydim,xdim],ydim/2,xdim/2);
shiftind = sub2ind([ydim,xdim],y,x);
ind_diff = shiftind-centerind;
shiftCellList = CellList;
for iter=1:numel(shiftCellList)
    templist= CellList(iter).PixelIdxList - ind_diff;
    shiftCellList(iter).PixelIdxList = templist(and(templist < ydim*xdim, templist > 0));
    shiftCellList(iter).perimeter = [];
end

end
