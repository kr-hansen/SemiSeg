% Pass this function a 3D or 2D uint8/16 data struct of a video volume (imArray).
%Allows you to manually select ROIs based on maximum values for each pixel in data.
function CellList = SemiSeg(imArray, inputCellList)
    %Start up
    fullsz = size(imArray);
    if numel(fullsz) == 3
        nimages = fullsz(3);
        singleFrame = max(imArray,[],3) - min(imArray,[],3);
        sz = size(singleFrame);
    elseif numel(fullsz) == 2
        nimages = 1;
        singleFrame = imArray;
        sz = size(singleFrame);
    end
    NumMinPixels = 5; %Minimum Number of Pixels for ROI
    Min = 0; Max = 5000;
    mainfig = figure;
    H = imshow(singleFrame, [Min Max]); title('Select rois for this video');
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
    rois_type = questdlg('Select ROI Type','What Code Generated the ROIs?','Simon''s','Kyle''s','Kyle''s');
    switch rois_type
        case 'Simon''s'
            pixel_list = 'pixel_idx';
        case 'Kyle''s'
            pixel_list = 'PixelIdxList';
    end
    
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
            else
                outline = plot(b{1}(:,2),b{1}(:,1),'g');
            end
        end
        CellNum = numel(CellList)+1;
    else
        CellList = struct();
        CellNum = 1;
    end
    
    %Select ROI clicking type
    click_type = questdlg('Select Click Type for ROI Selection','Would you like just Circles or Full ROI Selection?','Circles','Full ROI','Full ROI');
    switch click_type
        case 'Circles' 
            while(true)
                mb = msgbox('Zoom to desired area');
                zoom on;
                pause
                [xList,yList] = getpts;
                for inputs=1:length(xList) %Loop through points
                    %Taken and altered from Mike's selectSingleFrameRois.m borrowed from https://www.mathworks.com/matlabcentral/newsreader/view_thread/146031
                    cx=xList(inputs); cy=yList(inputs); ix=size(singleFrame,1); iy=size(singleFrame,2); r=6;
                    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                    c_mask=((x.^2+y.^2)<=r^2);
                    [b,~,~,~] = bwboundaries(c_mask);
                    circ = plot(b{1}(:,2),b{1}(:,1),'g'); %Plot for visualization
                    %Add to output structure
                    CC = bwconncomp(c_mask); 
                    CellList(CellNum).(pixel_list) = CC.PixelIdxList{1};
                    CellList(CellNum).perimeter = bwperim(c_mask);
                    assignin('base','tempCellList',CellList); %Temporarily Assign tempCellList to workspace
                    CellNum = CellNum+1;
                end
                button = questdlg('Would you  like to add more rois?','More rois?','Yes','No','Yes');
                if ~strcmp(button,'Yes')
                    break
                end
            end

        case 'Full ROI'
            %Start in on ROI Selection
            while(true)
                mb = msgbox('Zoom to desired area');
                zoom on;
                pause
                processtype = questdlg('Process Image?','Type of Processing','Raw','Corr','PCA','Raw');
                switch processtype
                    case 'Raw' %Threshold from Raw Input Image (Only case that works if only single frame is inputted
                        while(true)
                            compFrame = double(singleFrame);%/256; %Hard Coded for uint8 input
                            percentile = .95;
                            threshold = quantile(compFrame(:),percentile);
                            polyselect = roipoly();
                            adjim = polyselect.*compFrame;
                            binmask = zeros(size(adjim));
                            binmask(adjim >= threshold) = 1;
                            binmask = bwareaopen(binmask, NumMinPixels);
                            binmask = imfill(binmask,'holes');
                            b = bwboundaries(binmask);
                            tempim = adjim;
                            if isempty(b)  %Catch if Threshold too High
                                catchthreshold = quantile(tempim(:),0.5);
                                binmask(adjim >= catchthreshold) = 1;
                                binmask = bwareaopen(binmask, NumMinPixels);
                                binmask = imfill(binmask,'holes');
                                b = bwboundaries(binmask);
                            end
                            outline = plot(b{1}(:,2),b{1}(:,1),'g');
                            %Keep this ROI?
                            keepans = questdlg('Keep ROI','Keep Selected ROI?','No','Yes','Adjust Thresh','No');
                            switch keepans
                                case 'Adjust Thresh' %Adjust Threshold & Try again
                                    while(true)
                                        delete(outline)
                                        answer = inputdlg({'Percentile Threshold'},'Contrast',1,{num2str(percentile)});
                                        newpercentile = str2double(answer{1});
                                        newthresh = quantile(compFrame(:),newpercentile);
                                        binmask = zeros(size(adjim));
                                        binmask(adjim >= newthresh) = 1;
                                        binmask = bwareaopen(binmask, NumMinPixels);
                                        binmask = imfill(binmask,'holes');
                                        b = bwboundaries(binmask);
                                        if isempty(b)  %Catch if Threshold too High
                                            catchthreshold = quantile(tempim(:),0.5);
                                            binmask(adjim >= catchthreshold) = 1;
                                            binmask = bwareaopen(binmask, NumMinPixels);
                                            binmask = imfill(binmask,'holes');
                                            b = bwboundaries(binmask);
                                        end
                                        outline = plot(b{1}(:,2),b{1}(:,1),'g');
                                        again = questdlg('Again?','Adjust Threshold Again?','Yes','No','Yes');
                                        if ~strcmp(again,'Yes')
                                            break
                                        end
                                    end
                                    fullmask = binmask;
                                    CC = bwconncomp(fullmask);
                                    CellList(CellNum).(pixel_list) = CC.PixelIdxList{1};
                                    CellList(CellNum).perimeter = bwperim(fullmask);
                                    assignin('base','tempCellList',CellList); %Temporarily Assign tempCellList to workspace
                                    CellNum = CellNum+1;
                                case 'No' %Remove ROI
                                    delete(outline)
                                case 'Yes'
                                    fullmask = binmask;
                                    CC = bwconncomp(fullmask);
                                    CellList(CellNum).(pixel_list) = CC.PixelIdxList{1};
                                    CellList(CellNum).perimeter = bwperim(fullmask);
                                    assignin('base','tempCellList',CellList); %Temporarily Assign tempCellList to workspace
                                    CellNum = CellNum+1;
                            end
                            %Leave this view, or keep this view?
                            exit = questdlg('Select More ROIs from this FOV?','More ROIs from this Zoom?','Yes','No','Yes');
                            if ~strcmp(exit,'Yes')
                                break
                            end
                        end
                    case 'Corr' %Local Area Pixelwise Cross-Correlation
                        %Get Data for Zoomed Area
                        percentile = .95;
                        xl = get(gca,'XLim');
                        yl = get(gca,'YLim');
                        xstart = floor(xl(1));
                        xend = ceil(xl(2));
                        ystart = floor(yl(1));
                        yend = ceil(yl(2));
                        tempim = zeros(yend-ystart,xend-xstart);
                        temptraces = zeros(nimages,9);
                        for xshift = 0:(xend-xstart-1) %Loop Through Pixels
                            for yshift = 0:(yend-ystart-1)
                                curx = xstart+xshift;
                                cury = ystart+yshift;
                                temptraces(:,1) = imArray(cury,curx,:);
                                temptraces(:,2) = imArray(cury-1,curx-1,:);
                                temptraces(:,3) = imArray(cury,curx-1,:);
                                temptraces(:,4) = imArray(cury+1,curx-1,:);
                                temptraces(:,5) = imArray(cury-1,curx,:);
                                temptraces(:,6) = imArray(cury+1,curx,:);
                                temptraces(:,7) = imArray(cury-1,curx+1,:);
                                temptraces(:,8) = imArray(cury,curx+1,:);
                                temptraces(:,9) = imArray(cury+1,curx+1,:);
                                tempcorrcoefs = corrcoef(temptraces);
                                tempim(yshift+1,xshift+1) = mean(tempcorrcoefs(1,2:end));
                            end
                        end
                        tempfig = figure(); imagesc(tempim); axis image; title('Correlation Image'); hold on
                        while(true)
                            threshold = quantile(tempim(:),percentile);
                            polyselect = roipoly();
                            adjim = polyselect.*tempim;
                            binmask = zeros(size(adjim));
                            binmask(adjim >= threshold) = 1;
                            binmask = bwareaopen(binmask, NumMinPixels);
                            binmask = imfill(binmask,'holes');
                            b = bwboundaries(binmask);
                            if isempty(b)  %Catch if Threshold too High
                                catchthreshold = quantile(tempim(:),0.5);
                                binmask(adjim >= catchthreshold) = 1;
                                binmask = bwareaopen(binmask, NumMinPixels);
                                binmask = imfill(binmask,'holes');
                                b = bwboundaries(binmask);
                            end
                            outline = plot(b{1}(:,2),b{1}(:,1),'g');
                            %Keep this ROI?
                            keepans = questdlg('Keep ROI','Keep Selected ROI?','No','Yes','Adjust Thresh','No');
                            switch keepans
                                case 'Adjust Thresh' %Adjust Threshold & Try again
                                    newpercentile = percentile;
                                    while(true)
                                        delete(outline)
                                        answer = inputdlg({'Percentile Threshold'},'Contrast',1,{num2str(newpercentile)});
                                        newpercentile = str2double(answer{1});
                                        newthresh = quantile(tempim(:),newpercentile);
                                        binmask = zeros(size(adjim));
                                        binmask(adjim >= newthresh) = 1;
                                        binmask = bwareaopen(binmask, NumMinPixels);
                                        binmask = imfill(binmask,'holes');
                                        b = bwboundaries(binmask);
                                        if isempty(b)  %Catch if Threshold too High
                                            catchthreshold = quantile(tempim(:),0.5);
                                            binmask(adjim >= catchthreshold) = 1;
                                            binmask = bwareaopen(binmask, NumMinPixels);
                                            binmask = imfill(binmask,'holes');
                                            b = bwboundaries(binmask);
                                        end
                                        outline = plot(b{1}(:,2),b{1}(:,1),'g');
                                        again = questdlg('Again?','Adjust Threshold Again?','Yes','No','Delete ROI','Yes');
                                        if strcmp(again,'No')
                                            break
                                        elseif strcmp(again,'Delete ROI') %Remove ROI
                                            break
                                        end
                                    end
                                    if strcmp(again,'Delete ROI')
                                        delete(outline)
                                    else
                                        fullmask = zeros(sz);
                                        binsz = size(binmask);
                                        fullmask(ystart:ystart+(binsz(1)-1),xstart:xstart+(binsz(2)-1)) = binmask;
                                        CC = bwconncomp(fullmask);
                                        CellList(CellNum).(pixel_list) = CC.PixelIdxList{1};
                                        CellList(CellNum).perimeter = bwperim(fullmask);
                                        assignin('base','tempCellList',CellList); %Temporarily Assign tempCellList to workspace
                                        CellNum = CellNum+1;
                                        %Add Outline to Main Figure
                                        figure(mainfig);
                                        b = bwboundaries(fullmask);
                                        if isempty(b)  %Catch if Threshold too High
                                            catchthreshold = quantile(tempim(:),0.5);
                                            binmask(adjim >= catchthreshold) = 1;
                                            binmask = bwareaopen(binmask, NumMinPixels);
                                            binmask = imfill(binmask,'holes');
                                            b = bwboundaries(binmask);
                                        end
                                        plot(b{1}(:,2),b{1}(:,1),'g');
                                        figure(tempfig);
                                    end
                                case 'No' %Remove ROI
                                    delete(outline)
                                case 'Yes'
                                    fullmask = zeros(sz);
                                    binsz = size(binmask);
                                    fullmask(ystart:ystart+(binsz(1)-1),xstart:xstart+(binsz(2)-1)) = binmask;
                                    CC = bwconncomp(fullmask);
                                    CellList(CellNum).(pixel_list) = CC.PixelIdxList{1};
                                    CellList(CellNum).perimeter = bwperim(fullmask);
                                    assignin('base','tempCellList',CellList); %Temporarily Assign tempCellList to workspace
                                    CellNum = CellNum+1;
                                    %Add outline to main figure
                                    figure(mainfig);
                                    b = bwboundaries(fullmask);
                                    if isempty(b)  %Catch if Threshold too High
                                        catchthreshold = quantile(tempim(:),0.5);
                                        binmask(adjim >= catchthreshold) = 1;
                                        binmask = bwareaopen(binmask, NumMinPixels);
                                        binmask = imfill(binmask,'holes');
                                        b = bwboundaries(binmask);
                                    end
                                    plot(b{1}(:,2),b{1}(:,1),'g');
                                    figure(tempfig);
                            end
                            %Leave this view, or keep this view?
                            exit = questdlg('Select More ROIs from this FOV?','More ROIs from this Zoom?','Yes','No','Yes');
                            if ~strcmp(exit,'Yes')
                                break
                            end
                            figure(tempfig);
                        end
                        close(tempfig);
                    case 'PCA' %Local Area Principle Component Analysis
                        %Get Data for Zoomed Area
                        percentile = .975;
                        xl = get(gca,'XLim');
                        yl = get(gca,'YLim');
                        xstart = floor(xl(1));
                        xend = ceil(xl(2));
                        ystart = floor(yl(1));
                        yend = ceil(yl(2));
                        xdim = xend-xstart+1;
                        ydim = yend-ystart+1;
                        zoomVolume = double(imArray(ystart:yend,xstart:xend,:));
                        zoomVolReshape = squeeze(reshape(zoomVolume,[],xdim*ydim,nimages))';
                        PCs = pca(zoomVolReshape,'NumComponents',3); %3 PCs for RGB Colors
                        unscaledim = reshape(permute(PCs,[3,1,2]),ydim,xdim,3);
                        shiftedim = zeros(size(unscaledim));
                        tempim = zeros(size(unscaledim));
                        for idx = 1:3 %Hard Coded 3 Colors
                            colmin = min(min(unscaledim(:,:,idx)));
                            shiftedim(:,:,idx) = unscaledim(:,:,idx) + abs(colmin);
                            shiftmax = max(max(shiftedim(:,:,idx)));
                            tempim(:,:,idx) = shiftedim(:,:,idx)*(256/shiftmax); %Hard Coded 8-bit
                        end
                        tempim8 = uint8(tempim);
                        %Display Image and select ROIs
                        tempfig = figure(); imagesc(tempim8); axis image; title('PCA Image'); hold on
                        while(true)
                            polyselect = roipoly();
                            adjim = polyselect.*tempim;
                            [y,x,c] = size(adjim);
                            binmask = zeros(y,x);
                            for idx = 1:3 %Loop through 3 Colors independently
                                colim = tempim(:,:,idx);
                                threshold = quantile(colim(:),percentile);
                                binmask(adjim(:,:,idx) >= threshold) = 1;
                            end
                            binmask = bwareaopen(binmask, NumMinPixels);
                            binmask = imfill(binmask,'holes');
                            b = bwboundaries(binmask);
                            if isempty(b) %Catch if Threshold too High
                                binmask = zeros(y,x);
                                for idx = 1:3 %Loop through 3 Colors independently
                                    colim = tempim(:,:,idx);
                                    catchthreshold = quantile(colim(:), 0.5);
                                    binmask(adjim(:,:,idx) >= catchthreshold) = 1;
                                end
                                binmask = bwareaopen(binmask, NumMinPixels);
                                binmask = imfill(binmask,'holes');
                                b = bwboundaries(binmask);
                            end
                            outline = plot(b{1}(:,2),b{1}(:,1),'g');
                            %Keep this ROI?
                            keepans = questdlg('Keep ROI','Keep Selected ROI?','No','Yes','Adjust Thresh','No');
                            switch keepans
                                case 'Adjust Thresh' %Adjust Threshold & Try again
                                    newpercentile = percentile;
                                    while(true)
                                        delete(outline)
                                        answer = inputdlg({'Percentile Threshold'},'Contrast',1,{num2str(newpercentile)});
                                        newpercentile = str2double(answer{1});
                                        binmask = zeros(y,x);
                                        for idx = 1:3 %Loop through 3 Colors independently
                                            colim = tempim(:,:,idx);
                                            newthresh = quantile(colim(:),newpercentile);
                                            binmask(adjim(:,:,idx) >= newthresh) = 1;
                                        end
                                        binmask = bwareaopen(binmask, NumMinPixels);
                                        binmask = imfill(binmask,'holes');
                                        b = bwboundaries(binmask);
                                        if isempty(b) %Catch if Threshold too High
                                            binmask = zeros(y,x);
                                            for idx = 1:3 %Loop through 3 Colors independently
                                                colim = tempim(:,:,idx);
                                                catchthreshold = quantile(colim(:), 0.5);
                                                binmask(adjim(:,:,idx) >= catchthreshold) = 1;
                                            end
                                            binmask = bwareaopen(binmask, NumMinPixels);
                                            binmask = imfill(binmask,'holes');
                                            b = bwboundaries(binmask);
                                        end
                                        outline = plot(b{1}(:,2),b{1}(:,1),'g');
                                        again = questdlg('Again?','Adjust Threshold Again?','Yes','No','Delete ROI','Yes');
                                        if strcmp(again,'No')
                                            break
                                        elseif strcmp(again,'Delete ROI') %Remove ROI
                                            break
                                        end
                                    end
                                    if strcmp(again,'Delete ROI')
                                        delete(outline)
                                    else
                                        fullmask = zeros(sz);
                                        binsz = size(binmask);
                                        fullmask(ystart:ystart+(binsz(1)-1),xstart:xstart+(binsz(2)-1)) = binmask;
                                        CC = bwconncomp(fullmask);
                                        CellList(CellNum).(pixel_list) = CC.PixelIdxList{1};
                                        CellList(CellNum).perimeter = bwperim(fullmask);
                                        assignin('base','tempCellList',CellList); %Temporarily Assign tempCellList to workspace
                                        CellNum = CellNum+1;
                                        %Add Outline to Main Figure
                                        figure(mainfig);
                                        b = bwboundaries(fullmask);
                                        if isempty(b) %Catch if Threshold too High
                                            binmask = zeros(y,x);
                                            for idx = 1:3 %Loop through 3 Colors independently
                                                colim = tempim(:,:,idx);
                                                catchthreshold = quantile(colim(:), 0.5);
                                                binmask(adjim(:,:,idx) >= catchthreshold) = 1;
                                            end
                                            binmask = bwareaopen(binmask, NumMinPixels);
                                            binmask = imfill(binmask,'holes');
                                            b = bwboundaries(binmask);
                                        end
                                        plot(b{1}(:,2),b{1}(:,1),'g');
                                        figure(tempfig);
                                    end
                                case 'No' %Remove ROI
                                    delete(outline)
                                case 'Yes'
                                    fullmask = zeros(sz);
                                    binsz = size(binmask);
                                    fullmask(ystart:ystart+(binsz(1)-1),xstart:xstart+(binsz(2)-1)) = binmask;
                                    CC = bwconncomp(fullmask);
                                    CellList(CellNum).(pixel_list) = CC.PixelIdxList{1};
                                    CellList(CellNum).perimeter = bwperim(fullmask);
                                    assignin('base','tempCellList',CellList); %Temporarily Assign tempCellList to workspace
                                    CellNum = CellNum+1;
                                    %Add outline to main figure
                                    figure(mainfig);
                                    b = bwboundaries(fullmask);
                                    if isempty(b) %Catch if Threshold too High
                                        binmask = zeros(y,x);
                                        for idx = 1:3 %Loop through 3 Colors independently
                                            colim = tempim(:,:,idx);
                                            catchthreshold = quantile(colim(:), 0.5);
                                            binmask(adjim(:,:,idx) >= catchthreshold) = 1;
                                        end
                                        binmask = bwareaopen(binmask, NumMinPixels);
                                        binmask = imfill(binmask,'holes');
                                        b = bwboundaries(binmask);
                                    end
                                    plot(b{1}(:,2),b{1}(:,1),'g');
                                    figure(tempfig);
                            end
                            %Leave this view, or keep this view?
                            exit = questdlg('Select More ROIs from this FOV?','More ROIs from this Zoom?','Yes','No','Yes');
                            if ~strcmp(exit,'Yes')
                                break
                            end
                            figure(tempfig);
                        end
                        close(tempfig);
                end %End Switch-Case for Raw/Corr/PCA
                button = questdlg('Would you  like to add more rois?','More rois?','Yes','No','Yes');
                if ~strcmp(button,'Yes')
                    break
                end
            end
    end

end



