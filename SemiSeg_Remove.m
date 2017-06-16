% Pass this function a 3D uint8/16 data struct of a video volume (imArray).
%Allows you to manually select ROIs based on maximum values for each pixel in data.
function [CellList, allRemoved] = SemiSeg_Remove(singleFrame, inputCellList)
    %Start up
    sz = size(singleFrame);
    Nx = sz(2);
    Ny = sz(1);
    NumMinPixels = 5; %Minimum Number of Pixels for ROI
    Min = 0; Max = 5000;
    mainfig = figure;
    H = imshow(singleFrame, [Min Max]); title('Select rois to Remove');
    while(true) %Adjust Contrast
        answer = inputdlg({'Min','Max'}, 'Contrast',1,{num2str(Min),num2str(Max)});
        Min = str2double(answer{1});
        Max = str2double(answer{2});
        figure(mainfig)
        H = imshow(singleFrame, [Min Max]); title('Adjust Contrast');
        button = questdlg('Is this okay?',...
            'Select this Contrast Level?','Yes','No','Yes');
        if strcmp(button,'Yes')
            break
        else
            continue
        end
    end
    hold on
    
    %Select between Kyle & Simon's ROIs
    rois_type = questdlg('Select ROI Type','What Code Generated the ROIs?','Simons','Kyles','Simons');
    switch rois_type
        case 'Simons'
            pixel_list = 'pixel_idx';
        case 'Kyles'
            pixel_list = 'PixelIdxList';
    end
    
    smpos = 1;
    if ~isempty(inputCellList) %If given input CellList
        CellList = inputCellList;
        for idx = 1:numel(CellList)
            binmask = zeros(sz);
            binmask(CellList(idx).(pixel_list)) = 1;
            binmask = bwareaopen(binmask, NumMinPixels);
            binmask = imfill(binmask,'holes');
            b = bwboundaries(binmask);
            if isempty(b)
                fprintf('Empty Pixel Indicies for %f\n',idx)
                smallrois(smpos) = idx;
                smpos = smpos+1;
            else
                outline(idx) = plot(b{1}(:,2),b{1}(:,1),'g');
            end
        end
        CellNum = numel(CellList)+1;
    else
        CellList = struct();
        CellNum = 1;
    end
    
    %Initialize Points for Looping
    pos = 1;
    
    %Start in on ROI Selection to remove
    while(true)
        mb = msgbox('Zoom to desired area');
        zoom on;
        pause
        [clickx, clicky]  = getpts;
        selx = round(clickx);
        sely = round(clicky);
        rmpixels = sub2ind(sz, sely, selx);
        for step = 1:numel(inputCellList)
            temp = intersect(inputCellList(step).(pixel_list), rmpixels);    
            if ~isempty(temp)
                rmrois(pos) = step;
                pos = pos+1;
                delete(outline(step));
            end
        end
        
        %Determine if Done
        button = questdlg('Would you like to remove more rois?','Remove more rois?','Yes','No','Yes');
        if ~strcmp(button,'Yes')
            break
        end
    end
    
    if exist('smallrois','var')
        allRemoved = sort(horzcat(smallrois,rmrois));
    else
        allRemoved=rmrois;
    end
    
    CellList(allRemoved) = [];

end



