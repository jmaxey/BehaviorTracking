function [output_vid,centroids] = track_mouse_position_claudia(vid_name,output_vid_name)

%%% Inputs: 
%%%    vid_name: string with path to the behavior video
%%%    output_vid_name: path and name of the file to save (should be avi)
%%% Outputs:
%%%    output_vid: 4-D structure containing the output video [rows cols color frame]
%%%    centroids: 2-D matrix containing the centroid location of the mouse
%%%    [x y]

fprintf('Reading in behavior video...\n');
behavior_vid = VideoReader(vid_name);
framerate = behavior_vid.FrameRate;
duration = behavior_vid.Duration;
writerObj = VideoWriter(output_vid_name);
writerObj.FrameRate = framerate;
open(writerObj);
num_frames = framerate*duration;

%%% Create average background image
bkg_vid = read(behavior_vid,[1 min(num_frames,5000)]);
bkg_vid = squeeze(bkg_vid(:,:,1,:));
bkg_image = mean(bkg_vid,3);
bkg_image = uint8(bkg_image);
bkg_thresh = graythresh(bkg_image);
bkg_image_bw = im2bw(bkg_image,bkg_thresh);
clear bkg_vid

%%% Get the area of the mouse
ref_image = read(behavior_vid,10);
[rows,cols,colors] = size(ref_image);
imagesc(ref_image); colormap('gray');
title('Draw outline around mouse');
drawnow
BW_mouse = roipoly(ref_image);
props = regionprops(BW_mouse,'Centroid');
mouse_area = length(find(BW_mouse(:) == 1));

%%% Get the outline of the box
imagesc(ref_image); colormap('gray');
title('Draw a box around the track');
h_rect = imrect;
mask_rect = round(getPosition(h_rect));
startRow = mask_rect(2);
stopRow = mask_rect(2)+mask_rect(4);
startCol = mask_rect(1);
stopCol = mask_rect(1)+mask_rect(3);
BW_track = zeros(size(ref_image));
BW_track(startRow:stopRow,startCol:stopCol) = 1;

width_thresh = rows/2;

centroids = zeros(num_frames,2);
mouse_start = 0;

start_frame = 1;
stop_frame = num_frames;

fig = figure(1);
for f=start_frame:stop_frame
    
    disp(f);
    
    image = read(behavior_vid,f);
    
    if(f == start_frame)
        cla 
        imagesc(image); colormap('gray');
        mouse_start = input('Is mouse visible? (Yes (1) or No (0))');
        if(mouse_start == 1)
            loc = ginput(1);
            centroids(f,:) = loc;
        end
        final_image = zeros(size(image,1),size(image,2));
    else
        [final_image_m,bwObj,s] = find_mouse(BW_track,bkg_image,image,mouse_area);
        
        cent_vector = [s.Centroid];
        area_vector = [s.Area];
        min_length_vector = [s.MinorAxisLength];
        extent_vector = [s.Solidity];
        max_length_vector = [s.MajorAxisLength];
        
        if isempty(area_vector) %%% If no objects are found, get user input
            final_image  = final_image_m;
            if(mouse_start == 0)
                centroids(f,:) = c_old;
            else
                imagesc(image);
                colormap('gray');
                loc = ginput(1);
                centroids(f,:) = loc;
            end
        else
            if(f ~= start_frame)
                prev_cent = centroids(max(f-1,1),:);
            end
            count = 0;
            dist = zeros(1,length(cent_vector)/2); %%% Find distance of each object to previous centroid
            for c = 1:2:length(cent_vector)
                count = count+1;
                dist(count) = sqrt((prev_cent(1)-cent_vector(c))^2 + (prev_cent(2)-cent_vector(c+1))^2);
            end
            [ext_idx] = find(extent_vector > 0.975);
            area_vector(ext_idx) = 0;
            [len_idx] = find(abs(min_length_vector-max_length_vector) < 1); %%% Get rid of square objs
            area_vector(len_idx) = 0;
            
            if(mouse_start == 1)
                [id_dist] = find(dist > width_thresh); %%% Get rid of objs too far away from previous mouse location
                area_vector(id_dist) = 0;
            elseif(mouse_start == 0)
                [id_dist] = find(dist > width_thresh*1.5);
                area_vector(id_dist) = 0;
            end
            
            [min_val,id] = min(abs(area_vector-mouse_area)); %%% Find reminaing obj closest in size to the mouse
     
            final_image = zeros(size(final_image_m));
            final_image(bwObj.PixelIdxList{id}) = 1;
            
            clear ext_idx len_idx id_dist
            
            if(min_val == mouse_area) %%% If no objects met above criteria (area = 0)   
                if(mouse_start == 0)
                    centroids(f,:) = c_old;
                    final_image = zeros(size(final_image_m));
                else
                    mouse_start = 1;
                    loc = c_old;
                    centroids(f,:) = loc;
                    center_ind = sub2ind(size(final_image_m),round(loc(2)),round(loc(1)));
                    obj_intersect = 0;
                    %%% Find object that intersects the previous centroid
                    %%% location and is greater than 1/3 of the mouse's
                    %%% area
                    for v = 1:bwObj.NumObjects
                        idx_list = bwObj.PixelIdxList{v};
                        obj_inter = intersect(idx_list,center_ind);
                        if(~isempty(obj_inter) && length(idx_list) > mouse_area/3)
                            obj_intersect = v;
                        end
                    end
                    
                    if(obj_intersect == 0)
                        %%% If not object meets the above criteria, find the
                        %%% object closest to the last centroid and ensure that
                        %%% it is within a wider distance threshold and
                        %%% greater than 1/3 of the mouse's area
                        [~,new_idx] = min(dist);
                        area_vector = [s.Area];
                        if((dist(new_idx) <= width_thresh*1.25) && ...
                                (area_vector(new_idx) > mouse_area/3) && (f > 50))
                            final_image = zeros(size(final_image));
                            final_image(bwObj.PixelIdxList{new_idx}) = 1;
                        else
                            %%% If all of the above fail, get user input
                            final_image = zeros(size(final_image));
                            imagesc(image);
                            colormap('gray');
                            loc = ginput(1);
                            centroids(frame_idx,:) = loc;
                            if(on_plat(frame_idx-1) == 1)
                                on_plat(frame_idx) = 1;
                            end
                        end
                        
                    else
                        final_image = zeros(size(final_image));
                        final_image(bwObj.PixelIdxList{obj_intersect}) = 1;
                    end
                end
                c_old = centroids(f,:);
            elseif(min_val > mouse_area*1.5) 
                %%% If the minimum object area found is greater 1.5 times
                %%% the actual mouse area, get user input
                imagesc(image);
                colormap('gray');
                loc = ginput(1);
                centroids(f,:) = loc;
                final_image = zeros(size(image));
            else
                if(mouse_start == 0 && (min_val > mouse_area*.5))
                    centroids(f,:) = c_old;
                    final_image = zeros(size(image));
                else
                    mouse_start = 1;
                    bwObj = bwconncomp(final_image);
                    %%% If multiple objects pass the above criteria, choose
                    %%% object with area closest to that of the mouse
                    if(bwObj.NumObjects > 1)
                        props = regionprops(bwObj,'Area');
                        area_vector = [props.Area];
                        [min_val,id] = min(abs(area_vector-mouse_area));
                        final_image = zeros(size(final_image));
                        final_image(bwObj.PixelIdxList{id}) = 1;
                    end
                    
                    clear s
                    s = regionprops(final_image,'Centroid');
                    
                    c_new = s.Centroid(1:2);
                    centroids(f,:) = c_new;
                    final_image = zeros(size(final_image));
                    final_image(bwObj.PixelIdxList{1}) = 1;
                    
                    c_old = c_new;
                end
            end
        end
    end
    
    %%% Display results
    m_row = round(centroids(f,2));
    m_col = round(centroids(f,1));
    mouse_center = sub2ind([rows cols],m_row,m_col);
    mouse_ind = find(final_image == 1);
    
    if(isa(image,'uint8'))
        type = 2^8;
        type_min = 0;
    elseif(isa(image,'uint16'))
        type = 2^16;
        type_min = 0;
    elseif(isa(image,'single'))
        type = realmax('single');
        type_min = realmin('single');
    elseif(isa(image,'double'))
        type = realmax('double');
        type_min = realmin('double');
    end
    
    red = image;
    green = image;
    blue = image;
    
    blue(mouse_ind) = type;
    red(mouse_ind) = type_min;
    green(mouse_ind) = type_min;
    
    out_image = cat(3, red, green, blue);
    
    cla
    subplot(121)
    imagesc(out_image);
    subplot(122)
    imagesc(image); hold on; plot(centroids(f,1),centroids(f,2),'*r'); 
    if(f>start_frame+50)
        plot(centroids(f-50:f,1),centroids(f-50:f,2),'r');
    else
        plot(centroids(start_frame:f,1),centroids(start_frame:f,2),'r');
    end
    
    axis tight;
    
    output_vid(:,:,:,f) = out_image;
    
    data(1,:) = getframe(fig);
    im = data.cdata;
    writeVideo(writerObj,im);
    clear data im
    
    gate = 0;
    clear red blue green image len_idx area_vector...
        min_length_vector ecc_vector extent_vector max_length_vector time_idx...
        s bwObj final_image dist out_image
end

close(writerObj);

