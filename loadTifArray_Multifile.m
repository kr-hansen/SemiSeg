function imArray = loadTifArray_Multifile(num_frame_select)

%Load Tif Video File into 3D Array and subsample it by num_frame_select

[fname,fdir] = uigetfile('*.tif','MultiSelect','on');
switch class(fname)
    case 'char'
        filename{1} = fullfile(fdir,fname);
    case 'cell'
        filename = cell(numel(fname),1);
        for n = 1:numel(fname)
            filename{n} = fullfile(fdir,fname{n});
        end
end

%Determine Starting Array Extent
nfiles = numel(filename);
info = imfinfo(filename{1});
nimages = numel(info);
endinfo = imfinfo(filename{end});
endimages = numel(endinfo);

%Downsampling Calculations
if isempty(num_frame_select) %Set default for number of frames to select
    num_frame_select = 5; %Select every xth frame
end
max_frame_per_vid = ceil(nimages/num_frame_select);
select_frames = nan(max_frame_per_vid, nfiles);
frame_val = 1;
last_samp_frame = floor(nimages/num_frame_select);
for idx = 1:(nfiles-1)
    for step = 1:(last_samp_frame) %Removed +1 because of error when num_frame_select=1.  May need +1 if not taking all frames
        select_frames(step,idx) = frame_val;
        frame_val = frame_val + num_frame_select;
    end
    remain = nimages - (frame_val); %- num_frame_select to take from frame_val %Subtract extra num frames
    frame_val = num_frame_select - remain;
    last_samp_frame = floor((nimages-frame_val)/num_frame_select);
end
last_samp_frame = floor(endimages/num_frame_select);
for step = 1:(last_samp_frame) %Removed +1 because of error.  Shouldn't need on last file?
    select_frames(step,nfiles) = frame_val;
    frame_val = frame_val + num_frame_select;
end

imArray = zeros(info(1).Height, info(1).Width, sum(sum(~isnan(select_frames))), 'uint16');
iter = 0;

%Open Tiffs and write into array
for fnum = 1:nfiles
    TifLink = Tiff(filename{fnum}, 'r');
    iter_frames = select_frames(:,fnum);
    n_subsamp = sum(~isnan(iter_frames));
    for idx=1:n_subsamp
        TifLink.setDirectory(iter_frames(idx));
        imArray(:,:,iter+idx) = TifLink.read();
    end
    TifLink.close();
    iter = iter+n_subsamp;
end

end