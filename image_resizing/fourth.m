img = imread('image14.jpg');
[h, w, ~] = size(img); % not required , remove later

subplot(1, 2, 1);
imshow(img);
title('Original Image');

seams_to_remove = size(img, 2) - 1; 

writerObj = VideoWriter('seam_remove_video.mp4', 'MPEG-4');
writerObj.Quality = 90;
writerObj.FrameRate = 30;

open(writerObj);

writeVideo(writerObj, img);

for i = 1:seams_to_remove
    
    energy = energy_function(img);
    seam = find_vertical_seam(energy);
    img = remove_vertical_seam(img, seam);

    img_padded = padarray(img, [0, i], 0, 'post');

    writeVideo(writerObj, img_padded);
   
    img_seam = img_padded;
    for j = 1:size(seam, 1)
        img_seam(j, seam(j), 1) = 255;
        img_seam(j, seam(j), 2) = 0;
        img_seam(j, seam(j), 3) = 0;
    end
   
    writeVideo(writerObj, img_seam);
    
    %testing purpose 
    if i == 100
    subplot(1, 2, 2);
    imshow(img);
    title('for testing -Seam-Carved Image without padding after removing 100 seams');

    end

end

close(writerObj);


