%% Resizes a video to the desired dimensions and frame rate
%% CHANGE ME! ⌄⌄⌄
ExcelRows = 100;
changeFrameRate = 1;    % Toggle 0 or 1
outputFrameRate = 15;
shrinkFactor = 4;       % Shrink by a factor of 4 in both directions.
tableName = "Vince PatSep Trial Data 06-14-23.xlsx";
mainDir = "P:\Vince Mouse Tracking";
%% CHANGE ME! ^^^

MatlabRows = ExcelRows - 1;     % Offset by 1 row since 1st row becomes the table header in Matlab

cd(mainDir)
fullExcel = readtable(tableName);
ExcelArr = fullExcel(MatlabRows,:);

filePath = ExcelArr.VideoPath;

fontSize = 22;

for n = 1:length(filePath)
    inputVideo = VideoReader(filePath{n});

    %% Determine how many frames there are.
    numberOfFrames = inputVideo.NumFrames;
    inputVideoRows = inputVideo.Height;
    inputVideoColumns = inputVideo.Width;

       %% Create a VideoWriter object to write the video out to a new, different file.
    outputFullFileName = strrep(filePath{n}, '.mp4', '_Red.avi');
    outputVideo = VideoWriter(outputFullFileName,"Motion JPEG AVI");

    %% Specify the frame rate if desired (changeFrameRate = 0 for off, 1 for on)
    if changeFrameRate == 1
        outputVideo.FrameRate = outputFrameRate;
        frameRatio = inputVideo.FrameRate/outputFrameRate;
    else
        outputVideo.FrameRate = inputVideo.FrameRate;
        frameRatio = 1;
    end

    %% Specify the output video size.
    outputVideoRows = round(inputVideoRows / shrinkFactor);
    outputVideoColumns = round(inputVideoColumns / shrinkFactor);

    %% Prepare variables and a figure to show the images in the upper half of the screen.
    open(outputVideo);
    numberOfFramesWritten = 0;

    figure;
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

    %% Loop through the movie, writing all frames out.
    for frame = round(1 : frameRatio : numberOfFrames)
        thisInputFrame = read(inputVideo, frame);
        image(thisInputFrame);
        axis on;
        axis image;
        caption = sprintf('Frame %4d of %d.', frame, numberOfFrames);
        title(caption, 'FontSize', fontSize);
        drawnow;

        outputFrame = imresize(thisInputFrame, [outputVideoRows, outputVideoColumns]);  % Resize the image
        writeVideo(outputVideo, outputFrame);                                           % Write this new, resized frame out to the new video file.

        progressIndication = sprintf('Processed frame %4d of %d.', frame, numberOfFrames);
        numberOfFramesWritten = numberOfFramesWritten + 1;
    end

    %% Close the output video object (You don't need to close the input video reader)
    close(outputVideo);

    %% Show input frame and output frame side by side:
    subplot(1, 2, 1);
    image(thisInputFrame);
    axis on;
    axis image;
    caption = sprintf('Input Frame: %d rows by %d columns', inputVideoRows, inputVideoColumns);
    title(caption, 'FontSize', fontSize);
    
    %% Show the smaller output frame.
    subplot(1, 2, 2);
    image(outputFrame);
    axis on;
    axis image;
    caption = sprintf('Output Frame: %d rows by %d columns.', outputVideoRows, outputVideoColumns);
    title(caption, 'FontSize', fontSize);
    
    finishedMessage = sprintf('Done!  It processed %d frames of\n"%s" and created output video \n', numberOfFramesWritten, outputFullFileName);
    fprintf('%s\n', finishedMessage); % Write to command window.
    close all
end

disp("All files run :)");
