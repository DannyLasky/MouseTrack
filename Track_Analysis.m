function Track_Analysis(fullExcel, ExcelArr, MatlabRows, geom, closeThresh, reducedVid, tableName, mainDir, versionNum, onMac, jumpLimit, wallLength, manual)

% Works with Track_Master to track a mouse in a 2D arena
% Performs the frame-by-frame computational analysis
% Last updated 7/2023, Danny Lasky

%% Begin computational loop
for n = 1:length(MatlabRows)
    % n = 1;    
    tic

    %% Convert video path to match Mac
    if onMac == 1
        tempConvert = char(ExcelArr{n,'VideoPath'});
        tempConvert = strrep(tempConvert, 'P:', '/Volumes/mathewjones');
        videoPathMac = strrep(tempConvert, '\', '/');
        tempObj = videoPathMac;
    elseif onMac == 0
        tempConvert = char(ExcelArr{n,'VideoPath'});
        tempConvert = strrep(tempConvert, '/Volumes/mathewjones', 'P:');
        videoPathMac = strrep(tempConvert, '/', '\');
        tempObj = videoPathMac;
    else
	    tempObj = char(ExcelArr{n,'VideoPath'});
    end

    %% Only need geom parameters for the video currently being run
    geomRun.Edges    = geom.Edges{n};
    geomRun.Ctr      = geom.Ctr{n};
    geomRun.Perim    = geom.Perim{n};
    geomRun.DspEdges = geom.DspEdges{n};
    geomRun.DspPerim = geom.DspPerim{n};
    geomRun.MouseCtr = geom.MouseCtr{n};
    geomRun.ObjPos   = geom.ObjPos{n};

    %% Reduce geom parameters to match reduced video size (only need the geometry variables changed if the coordinates were selected on the unreduced video)
    if reducedVid == 1
        tempObj         = strrep(tempObj,'.mp4','_Red.avi');
        % geomRun.Edges     = geomRun.Edges / 4;
        % geomRun.Ctr       = geomRun.Ctr / 4;
        % geomRun.Perim     = geomRun.Perim / 4;
        % geomRun.DspEdges  = geomRun.DspEdges / 4;
        % geomRun.DspPerim  = geomRun.DspPerim / 4;
        % geomRun.MouseCtr  = geomRun.MouseCtr / 4;
        % geomRun.ObjPos    = geomRun.ObjPos / 4;
    end

    %% Read in video and Excel data
	vidObj          = VideoReader(tempObj);
	vid.Width   	= vidObj.Width;
	vid.Height  	= vidObj.Height;
	vid.ID          = "Coh " + ExcelArr{n,'Cohort'} + " Animal " + ExcelArr{n,'AnimalID'} + " Wk " + ExcelArr{n,'Week'} + " Phase " + ExcelArr{n,'Phase'};
	vid.FullPath    = char(ExcelArr{n,'VideoPath'});
	trialPhaseNum   = ExcelArr{n,'Phase'};
	objectMovedNum  = ExcelArr{n,'ObjectMoved'};
	distMoved       = ExcelArr{n,'DistanceObjectMoved'};
    videoOrient     = ExcelArr{n, 'VideoOrientation'};

    vid.FPS = vidObj.FrameRate;

    if isfield(manual, 'FramesRange')
        startFrame = manual.FramesRange(1);
        endFrame = manual.FramesRange(2);
    else
        startFrame = ExcelArr{n,'StartFrame'};
        endFrame = ExcelArr{n,'EndFrame'};
    end

 	initialFrame = read(vidObj, startFrame);

    bckgStartFrame  = ExcelArr{n,'StartFrame'};    % Always use true start frame for background generation so it doesn't differ by video run length
    bckgEndFrame    = ExcelArr{n,'EndFrame'};      % Always use true end frame for background generation so it doesn't differ by video run length

    %% If phase 2, figure out which object will be moved
    if trialPhaseNum == 2
        toBeMoved = fullExcel{MatlabRows(n)+1, 'ObjectMoved'};
    end

    %% Create a pixels to cm ratio. Assumes that the inner arena walls are 63 cm
    wallX = geomRun.Edges(:,1);
    wallY = geomRun.Edges(:,2);
    wall1 = norm([wallX(1); wallY(1)] - [wallX(2); wallY(2)]);
    wall2 = norm([wallX(2); wallY(2)] - [wallX(3); wallY(3)]);
    wall3 = norm([wallX(3); wallY(3)] - [wallX(4); wallY(4)]);
    wall4 = norm([wallX(4); wallY(4)] - [wallX(1); wallY(1)]);
    wallAvg = mean([wall1, wall2, wall3, wall4]);
    pixelsInCm = wallAvg/wallLength;
    toyRadPixels = round(closeThresh * pixelsInCm);

    %% Create a background image using frames from the start of the video and interspersed throughout true video run length
	initialFrames = vid.FPS;        % Number of frames in first second of video
	interFrames = floor((bckgEndFrame - bckgStartFrame)/vid.FPS);
    
	initialArr  = zeros(vid.Height, vid.Width, round(initialFrames));
	interArr    = zeros(vid.Height, vid.Width, interFrames);
    
    for i = 1:initialFrames
        imgTemp = readFrame(vidObj);
        imgTemp = rgb2gray(imgTemp);
        initialArr(:,:,i) = imgTemp;
    end
    
    for i = 1:interFrames
        imgTemp = read(vidObj, bckgStartFrame + i * round(vid.FPS));
        imgTemp = rgb2gray(imgTemp);
        interArr(:,:,i) = imgTemp;
    end
    
	initialImage = (mean(initialArr,3));
	interImage   = (mean(interArr,3));
    
	tempImage = cat(3, initialImage, interImage, interImage, interImage);   % 3/4ths weighted for inter frames
    
	bckgImage = (mean(tempImage,3));
	bckgImage = cast(bckgImage,'uint8');

    bckgImageAdj = imadjust(bckgImage);

	imgHeat = read(vidObj, startFrame);
    clear initialImage initialArr interImage interArr tempImage

    %% Read in 500 frames, convert to grayscale, delete initial read in, repeat
	numFrames = endFrame - startFrame + 1;
	numCycles = ceil(numFrames/500);
	disp(['Reading ' num2str(numFrames) ' frames — This may take ' num2str(0.005 * numFrames) ' secs'])
	fullVid = zeros(270, 480);
	fullVid = uint8(fullVid);
    
	tic
    for m = 1:numCycles
        startCycle  = startFrame + 500*(m-1);
        endCycle    = startFrame + 500*m - 1;
    
        if endCycle >  endFrame
            endCycle = endFrame;
        end
    
        tempVid = read(vidObj, [startCycle endCycle]);
        sizeVid = size(tempVid);
        grayVid = zeros(sizeVid(1), sizeVid(2), sizeVid(4));
    
        for p = 1:size(tempVid, 4)
            tempFrame = tempVid(:,:,:,p);
            grayVid(:,:,p) = rgb2gray(tempFrame);
        end
    
        grayVid = uint8(grayVid);
        fullVid = cat(3,fullVid,grayVid);
    end
    
	fullVid(:,:,1) = [];
 	clear tempVid grayVid
	toc

    %% Establish tracking variables
	runCount = 1;
    sTime = tic;
    eTime = zeros(numFrames,1);
    xPath = zeros(numFrames, 1);
	yPath =	zeros(numFrames, 1);
    blobArea = zeros(numFrames, 1);
    screencapData = struct('cdata', cell(numFrames, 1), 'colormap', []);

    %% Make sure that corners always match correctly, regardless of selection order or video orientation
    % Remember that y-axis counts downward for some reason!
    [~, topLeft] = min(sum(geomRun.Edges, 2));
    [~, botRight] = max(sum(geomRun.Edges, 2));
 
    geomIndices = 1:4;
    geomIndices(geomIndices == topLeft | geomIndices == botRight) = [];

    if geomRun.Edges(geomIndices(1), 1) > geomRun.Edges(geomIndices(2), 1)
        topRight = geomIndices(1);
        botLeft = geomIndices(2);
    elseif geomRun.Edges(geomIndices(1), 1) < geomRun.Edges(geomIndices(2), 1)
        topRight = geomIndices(2);
        botLeft = geomIndices(1);
    end

    if contains(videoOrient, 'Side')
        geomRun.Corner1 = [geomRun.Edges(botLeft,1), geomRun.Edges(botLeft,2), toyRadPixels];
        geomRun.Corner2 = [geomRun.Edges(topLeft,1), geomRun.Edges(topLeft,2), toyRadPixels];
        geomRun.Corner3 = [geomRun.Edges(topRight,1), geomRun.Edges(topRight,2), toyRadPixels];
        geomRun.Corner4 = [geomRun.Edges(botRight,1), geomRun.Edges(botRight,2), toyRadPixels];
        xEdgeMin = mean([geomRun.Corner1(1), geomRun.Corner2(1)]);
        xEdgeMax = mean([geomRun.Corner3(1), geomRun.Corner4(1)]);
        yEdgeMin = mean([geomRun.Corner2(2), geomRun.Corner3(2)]);
        yEdgeMax = mean([geomRun.Corner1(2), geomRun.Corner4(2)]);
    elseif contains(videoOrient, 'Up')
        geomRun.Corner1 = [geomRun.Edges(topLeft,1), geomRun.Edges(topLeft,2), toyRadPixels];
        geomRun.Corner2 = [geomRun.Edges(topRight,1), geomRun.Edges(topRight,2), toyRadPixels];
        geomRun.Corner3 = [geomRun.Edges(botRight,1), geomRun.Edges(botRight,2), toyRadPixels];
        geomRun.Corner4 = [geomRun.Edges(botLeft,1), geomRun.Edges(botLeft,2), toyRadPixels];
        xEdgeMin = mean([geomRun.Corner1(1), geomRun.Corner4(1)]);
        xEdgeMax = mean([geomRun.Corner2(1), geomRun.Corner3(1)]);
        yEdgeMin = mean([geomRun.Corner1(2), geomRun.Corner2(2)]);
        yEdgeMax = mean([geomRun.Corner3(2), geomRun.Corner4(2)]);
    end

    %% Prepare object circle displays
    if trialPhaseNum == 1
        initialFrameFig = insertShape(initialFrame, "circle", [geomRun.Corner1; geomRun.Corner2; geomRun.Corner3; geomRun.Corner4], ...
            Color = ["green", "green", "green", "green"]);
    else
        geomRun.objCircles1 = [geomRun.ObjPos(1,1), geomRun.ObjPos(1,2), toyRadPixels];
        geomRun.objCircles2 = [geomRun.ObjPos(2,1), geomRun.ObjPos(2,2), toyRadPixels];
        initialFrameFig = insertShape(initialFrame, "circle", [geomRun.objCircles1; geomRun.objCircles2; geomRun.Corner1; geomRun.Corner2; ...
            geomRun.Corner3; geomRun.Corner4], ...
            Color = ["yellow", "yellow", "green", "green", "green", "green"]);
    end

    %% Prepare frame by frame figure
  	figure;
	% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.025 0.025 0.475 0.95], 'Visible', 'on');
    set(gcf, 'Units', 'Inches', 'Position', [1 1 8 10], 'Visible', 'on');
    ax1 = subplot(2,1,1);
        im = image(initialFrameFig); 
 	    axis image
 	    hold on
        plot(geomRun.DspEdges(:,1), geomRun.DspEdges(:,2), 'y-');
	    plot(geomRun.DspPerim(:,1), geomRun.DspPerim(:,2), 'c-');
	    plot(geomRun.Ctr(1), geomRun.Ctr(2), 'm+', 'MarkerSize', 12);
	    pl1 = plot(0, 0, 'ro', 'MarkerFaceColor', 'none', 'MarkerSize', 16);
        pl2 = plot(0, 0, 'r.-', 'MarkerSize', 8, 'LineWidth', 0.5);
        ax1.FontSize = 12;
    ax2 = subplot(2,1,2);
  	    im2 = image(initialFrameFig);
  	    axis image 
        hold on
	    pl3 = plot(0, 0, 'ro', 'MarkerFaceColor', 'none', 'MarkerSize', 16);
        pl4 = plot(0, 0, 'mo', 'MarkerFaceColor', 'none', 'MarkerSize', jumpLimit * 2);
        ax2.FontSize = 12;

    %% Display user feedback and times as frames are run through
    for m = startFrame:endFrame
        if mod(m, 100) == 0
            eTime(m) = toc(sTime);
            disp(strcat("Analyzing frame ", num2str(m), " of frames ", num2str(startFrame), "-", num2str(endFrame), " from ", vid.ID))
            disp(['Elapsed time ' num2str(eTime(m)) ' sec'])   % Display output for user
        end
    
        frameCount = m - startFrame + 1;

        if runCount == 101
            runCount = 1; 
        end

        %% Begin computations
        thisFrame = fullVid(:,:,runCount);                                          % Frame for tracking 
        arenaMask = roipoly(thisFrame, geomRun.Edges(:,1), geomRun.Edges(:,2));     % Restrict tracking to within maze
    
        thisFrameAdj = imadjust(thisFrame);

        diffImg = bckgImageAdj - thisFrameAdj;        % Subtract from background

        gaussImg = imgaussfilt(diffImg, 2);

        %% Make the image binary
        binaryImg = gaussImg;
        binaryImg(arenaMask == 0) = 0;               % Remove values outside of arena

        cutoff = prctile(binaryImg(:), 99.7);
        if cutoff == 0
            cutoff = 1;
        end

        binaryImg(binaryImg < cutoff)  = 0;
        binaryImg(binaryImg >= cutoff) = 1;
        outsideJump = binaryImg;

        %% Remove values outside of jump limit
        if frameCount > 1
            prevX = xPath(frameCount - 1);
            prevY = yPath(frameCount - 1);
        elseif isfield(manual, 'FramesRange')
            if frameCount == 1
                prevX = manual.FrameOneMouse(1);
                prevY = manual.FrameOneMouse(2);
            end
        end

        if jumpLimit ~= 0
            if frameCount == 1
                jumpCenter = geomRun.MouseCtr;
            else
                jumpCenter = [prevX, prevY];
            end

            if isfield(manual, 'Center')
                if any(manual.Center(:,1) == frameCount)
                    jumpCenter = [manual.Center(manual.Center(:,1) == frameCount, 2), manual.Center(manual.Center(:,1) == frameCount, 3)];
                end
            end
 
            jumpLimitFrame = jumpLimit;
            if isfield(manual, 'Radius')
                if any(manual.Radius(:,1) == frameCount)
                    jumpLimitFrame = manual.Radius(manual.Radius(:,1) == frameCount, 2);
                end
            end
            
            mouseMask = createCirclesMask(thisFrame, jumpCenter, jumpLimitFrame);
            binaryImg(mouseMask == 0) = 0;                                             % Remove values outside of arena
            outsideJump(mouseMask == 1) = 0;   
            % figure; imshowpair(thisFrame, mouseMask)
        end

        %% Simple figure showing filtering process, might be good for a poster
        % figure
        % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
        % tiledlayout(3, 3, 'TileSpacing', 'compact')
        % nexttile
        %     imshow(bckgImage)
        %     title("Background Image")
        % nexttile
        %     bckgAdj = imadjust(bckgImage);
        %     imshow(bckgAdj)
        %     title("Background Image Full Range")
        % nexttile
        %     imshow(thisFrame)
        %     title("Current Frame")
        % nexttile
        %     imshow(thisFrameAdj)
        %     title("Current Frame Full Range")
        % nexttile
        %     imshow(diffImg)
        %     title("Current Frame Subtracted")
        % nexttile
        %     imshow(gaussImg)
        %     title("Current Frame Filtered")
        % nexttile
        %     imshow(binaryImg, [])
        %     title("Current Frame Binary")

        %% Find biggest blob
        blobs = bwconncomp(binaryImg);
        regionProps = regionprops(blobs, binaryImg, {'Area', 'Centroid'});
        blobAreas = cat(1, regionProps(:).Area);
        bigBlobIdx = find(blobAreas == max(blobAreas));

        %% Select biggest blob, if a tie, choose the one closer to the previous centroid  
        blobDist = zeros(length(bigBlobIdx), 1);

        if length(bigBlobIdx) > 1
            for p = 1:length(bigBlobIdx)
                blobDist(p) = abs(prevX - regionProps(bigBlobIdx(p)).Centroid(1)) + ...
                    abs(prevY - regionProps(bigBlobIdx(p)).Centroid(2));
            end
            [~, minDistIdx] = min(blobDist(:));
            bigBlobIdx = bigBlobIdx(minDistIdx);
        end

        %% Place on blob, if no blob in range place in previous spot
        if length(bigBlobIdx) == 1
            xPath(frameCount) = regionProps(bigBlobIdx).Centroid(1);
            yPath(frameCount) = regionProps(bigBlobIdx).Centroid(2); 
            blobArea(frameCount)  = regionProps(bigBlobIdx).Area;
        elseif isempty(bigBlobIdx)
            xPath(frameCount) = prevX;
            yPath(frameCount) = prevY;
            blobArea(frameCount)  = 0;
        end
    
        %% Color selected blob green, other blobs yellow, and background blue
        blobImg = binaryImg;

        smallBlobIdx = 1:blobs.NumObjects;
        smallBlobIdx(smallBlobIdx == bigBlobIdx) = [];
        for p = smallBlobIdx
            blobImg(blobs.PixelIdxList{p}) = 230;
        end
        
        if length(bigBlobIdx) == 1
            blobImg(blobs.PixelIdxList{bigBlobIdx}) = 170;
        end
    
        %% Get blobs outside of mouse jump limit and color differently for troubleshooting purposes
        outsideBlobs = bwconncomp(outsideJump);
        if outsideBlobs.NumObjects > 0
            for p = 1:outsideBlobs.NumObjects
                blobImg(outsideBlobs.PixelIdxList{p}) = 60;
            end
        end

        %% Apply manually set path
        if isfield(manual, 'Path')
            if any(manual.Path(:,1) == frameCount)
                xPath(frameCount) = manual.Path(manual.Path(:,1) == frameCount, 2);
                yPath(frameCount) = manual.Path(manual.Path(:,1) == frameCount, 3);
            end
        end

        %% Display image, estimated location, and troubleshooting text
        if trialPhaseNum == 1
            thisFrameFig = insertShape(thisFrame, "circle", [geomRun.Corner1; geomRun.Corner2; geomRun.Corner3; geomRun.Corner4], ...
                Color = ["green", "green", "green", "green"]);
        else
            thisFrameFig = insertShape(thisFrame, "circle", [geomRun.objCircles1; geomRun.objCircles2; geomRun.Corner1; geomRun.Corner2; ...
                geomRun.Corner3; geomRun.Corner4], ...
                Color = ["yellow", "yellow", "green", "green", "green", "green"]);
        end
        axes(ax1)
            colormap(ax1, 'gray')
            set(im, 'cdata', thisFrameFig);
            set(pl1, 'xdata', xPath(frameCount), 'ydata', yPath(frameCount));
            set(pl2, 'xdata', xPath(1:frameCount), 'ydata', yPath(1:frameCount));
            title({vid.ID; "Frame #" + m + " (" + frameCount + "/" + numFrames + ")"}, 'Interpreter', 'none', 'FontSize', 14)
            drawnow
        axes(ax2)
            set(im2, 'cdata', blobImg)
            colormap(ax2, jet(256))
            set(pl3, 'xdata', xPath(frameCount), 'ydata', yPath(frameCount))
            set(pl4, 'xdata', jumpCenter(1), 'ydata', jumpCenter(2), 'MarkerSize', jumpLimitFrame * 1.9)      % Approximate size
            title({"Frame #" + m + " (" + frameCount + "/" + numFrames + ")"; "Area: " + blobArea(frameCount)}, 'FontSize', 14);

        %% Save screen frame for video
        screencapData(frameCount) = getframe(gcf);

        if mod(frameCount, 100) == 0
            fullVid(:,:,1:100) = [];
        end

        runCount = runCount + 1;
    end
    
    clear fullVid
    close all
    
    %% Smooth mouse trajectory and convert to cm
 	z = smoothn({xPath, yPath}, 500);
    smoothX = z{1};
    smoothY = z{2};

    smoothXCm = z{1} / pixelsInCm;
    smoothYCm = z{2} / pixelsInCm;
    
    %% Calculate distance traveled and velocity
    shiftSmoothXCm = circshift(smoothXCm, -1);
    distXSmoothCm = abs(smoothXCm - shiftSmoothXCm);
    distXSmoothCm(end) = [];
    
    shiftSmoothYCm = circshift(smoothYCm, -1);
    distYSmoothCm = abs(smoothYCm - shiftSmoothYCm);
    distYSmoothCm(end) = [];
    
    distTravCm = sqrt(distXSmoothCm.^2 + distYSmoothCm.^2);
    distTravCmSum = sum(distTravCm);
    
    mouseVelCm = distTravCm/(1/vid.FPS);
    mouseVelCmAve = distTravCmSum/(numFrames/vid.FPS);

    %% Find distance and velocity in bins for figure
    figBins = 20;
    framesPerBin = floor((frameCount - 1)/figBins);
    frameAxis = 1:framesPerBin:figBins*framesPerBin;
    timeAxis = frameAxis / vid.FPS;
    distSum = sum(reshape(distTravCm(1:framesPerBin * figBins),framesPerBin,[]));
    velAvg = mean(reshape(mouseVelCm(1:framesPerBin * figBins),framesPerBin,[]));

    %% Calculate time in center and near walls
    insidePerim = inpolygon(smoothX, smoothY, geomRun.Perim(:,1), geomRun.Perim(:,2));
    pctNearWall = (1-length(find(insidePerim)) / length(smoothX)) * 100;
    pctNearCtr  = 100 - pctNearWall;

    %% Calculate distance from corners for figures
    distFromC1 = sqrt((smoothXCm - geomRun.Corner1(1)/pixelsInCm).^2 + (smoothYCm - geomRun.Corner1(2)/pixelsInCm).^2);
    distFromC2 = sqrt((smoothXCm - geomRun.Corner2(1)/pixelsInCm).^2 + (smoothYCm - geomRun.Corner2(2)/pixelsInCm).^2);
    distFromC3 = sqrt((smoothXCm - geomRun.Corner3(1)/pixelsInCm).^2 + (smoothYCm - geomRun.Corner3(2)/pixelsInCm).^2);
    distFromC4 = sqrt((smoothXCm - geomRun.Corner4(1)/pixelsInCm).^2 + (smoothYCm - geomRun.Corner4(2)/pixelsInCm).^2);

    probC1 = (length(find(distFromC1 <= closeThresh)) / numFrames) * 100;
    probC2 = (length(find(distFromC2 <= closeThresh)) / numFrames) * 100;
    probC3 = (length(find(distFromC3 <= closeThresh)) / numFrames) * 100;
    probC4 = (length(find(distFromC4 <= closeThresh)) / numFrames) * 100;

    %% Calculate time in corners and not in corners
    nearCorner = probC1 + probC2 + probC3 + probC4;
    notCorner = 100 - nearCorner;

    %% Calculate average distance from center and each toy in cm and prepare for figures
    geomRun.CtrCm = geomRun.Ctr / pixelsInCm;
    distFromCtr = sqrt((smoothXCm - geomRun.CtrCm(1)).^2 + (smoothYCm - geomRun.CtrCm(2)).^2);

    binMax = floor(wallLength) - 1;
    [ctrValues, distBins] = histcounts(distFromCtr, 0:2:binMax);
    distBins = distBins(1:end-1);

    distFromCtrPDF = ctrValues / sum(ctrValues);
    distFromCtrCDF = cumsum(distFromCtrPDF);
    distFromCtrAve = mean(distFromCtr);

    if trialPhaseNum ~= 1
        geomRun.ObjPosCm = geomRun.ObjPos / pixelsInCm;
        distFromA = sqrt((smoothXCm - geomRun.ObjPosCm(1,1)).^2 + (smoothYCm - geomRun.ObjPosCm(1,2)).^2);
        distFromB = sqrt((smoothXCm - geomRun.ObjPosCm(2,1)).^2 + (smoothYCm - geomRun.ObjPosCm(2,2)).^2);
        histAValues = histcounts(distFromA, 0:2:binMax);
        histBValues = histcounts(distFromB, 0:2:binMax);
        distFromAPDF = histAValues / sum(histAValues);
        distFromBPDF = histBValues / sum(histBValues);
        distFromACDF = cumsum(distFromAPDF);
        distFromBCDF = cumsum(distFromBPDF); 
        probCloseA = (length(find(distFromA <= closeThresh)) / numFrames) * 100;
        probCloseB = (length(find(distFromB <= closeThresh)) / numFrames) * 100;
        distFromAAve = mean(distFromA);
        distFromBAve = mean(distFromB);
    end

    %% Find time spent near to be moved/unmoved objects (phase 2) or moved/unmoved objects (phase 3)
    if trialPhaseNum == 2
        if contains(toBeMoved, 'A')
            probCloseMoved      = probCloseA;
            probCloseUnmoved    = probCloseB;
        elseif contains(toBeMoved, 'B')
            probCloseMoved      = probCloseB;
            probCloseUnmoved    = probCloseA;
        end
    elseif trialPhaseNum == 3
        if contains(objectMovedNum, 'A')
            probCloseMoved      = probCloseA;
            probCloseUnmoved    = probCloseB;
        elseif contains(objectMovedNum, 'B')
            probCloseMoved      = probCloseB;
            probCloseUnmoved    = probCloseA;
        end
    end

    %% Create mouse trajectory figure
    colorVec = linspace(1, 0, numFrames)';
    colorMap = [ones(numFrames, 1) , colorVec, colorVec];
    if trialPhaseNum == 1
        imgHeatFig = insertShape(imgHeat, "circle", [geomRun.Corner1; geomRun.Corner2; geomRun.Corner3; geomRun.Corner4], ...
            Color = ["green", "green", "green", "green"]);
    else
       imgHeatFig = insertShape(imgHeat, "circle", [geomRun.objCircles1; geomRun.objCircles2; geomRun.Corner1; geomRun.Corner2; ...
           geomRun.Corner3; geomRun.Corner4], Color = ["yellow", "yellow", "green", "green", "green", "green"]);
    end
    f2 = figure;
        % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
        set(gcf, 'Units', 'Inches', 'Position', [1 1 14 8]);
        ax = axes('pos', [0.1 0.1 0.8 0.8]);
        image(imgHeatFig);
        axis image
        hold on
        plot(geomRun.DspEdges(:,1), geomRun.DspEdges(:,2), 'y-', 'LineWidth', 1.5);
        plot(geomRun.DspPerim(:,1), geomRun.DspPerim(:,2), 'c-', 'LineWidth', 1.5);
        plot(geomRun.Ctr(1), geomRun.Ctr(2), 'm+', 'MarkerSize', 12, 'LineWidth', 1.5);  
        for k = 1:numFrames-1
            plot([smoothX(k), smoothX(k+1)], [smoothY(k), smoothY(k+1)], '.-', 'Color', colorMap(k,:), 'LineWidth', 1.5, 'MarkerSize', 14)
        end
        ax.Colormap = colorMap;
        c = colorbar;
        ax.FontSize = 14;
        title(vid.ID, 'FontSize', 22', 'Interpreter', 'none')
        c.Label.String = "Duration (%)";
        c.FontSize = 14;
        c.Label.FontSize = 18;
        c.Ticks = c.Ticks * 100;        % Convert to percent run time
        clim([0 100])
        axis off

    %% Heat map figure
    xEdges = linspace(xEdgeMin, xEdgeMax, 10);
    yEdges = linspace(yEdgeMin, yEdgeMax, 10);

    f3 = figure;
        % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
        set(gcf, 'Units', 'Inches', 'Position', [1 1 14 8]);
        ax = axes('pos', [0.1 0.1 0.8 0.8]);
        image(imgHeatFig);
        axis image
        hold on
        plot(geomRun.DspEdges(:,1), geomRun.DspEdges(:,2), 'y-', 'LineWidth', 1.5);
        plot(geomRun.DspPerim(:,1), geomRun.DspPerim(:,2), 'c-', 'LineWidth', 1.5);
        plot(geomRun.Ctr(1), geomRun.Ctr(2), 'm+', 'MarkerSize', 12, 'LineWidth', 1.5); 
        hist2 = histogram2(smoothX, smoothY, 9, 'DisplayStyle', 'tile', 'FaceAlpha', 0.5, 'XBinEdges', xEdges, 'YBinEdges', yEdges);
        
        % Convert histogram from frames to percent time
        histCountsFrames    = hist2.BinCounts;
        histCountsSec       = histCountsFrames/vid.FPS;
        histCountsPerTime   = histCountsSec/sum(histCountsSec, 'all') * 100;
        hist2.BinCounts     = histCountsPerTime;

        ax.FontSize = 14;
        title(vid.ID, 'FontSize', 22', 'Interpreter', 'none')
        c = colorbar;
        c.Label.String = "Time (%)";
        c.FontSize = 14;
        c.Label.FontSize = 18;

    %% Create stats summary figure
    f4 = figure;
    % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
    set(gcf, 'Units', 'Inches', 'Position', [1 1 15 9]);
    t = tiledlayout(2,3);
    title(t, vid.ID, 'FontSize', 24', 'Interpreter', 'none')
    t1 = nexttile(1);
        yyaxis left
            plot(timeAxis, cumsum(distSum), 'LineWidth', 2, 'Color', '#21B3DE')
            ylabel('Total Distance (cm)', 'FontSize', 14)
            ax = gca;
            ax.FontSize = 12;
            ax.YAxis(1).Color = '#21B3DE';
        yyaxis right
            hold on
            plot(timeAxis, velAvg, 'LineWidth', 2, 'Color', '#EA2B16')
            ylabel('Average Velocity (cm/s)', 'FontSize', 14)
            ax = gca;
            ax.FontSize = 12;
            ax.YAxis(2).Color = '#EA2B16';
        xlim([0 timeAxis(end)])
        title('Mouse Distance and Velocity', 'FontSize', 16)
        xlabel('Time (sec)', 'FontSize', 14)
    t2 = nexttile(2);
        hold on; box on
        bar(distBins, distFromCtrPDF, 'FaceColor', '#21B3DE')
        plot(distBins, distFromCtrCDF * max(distFromCtrPDF), 'LineWidth', 2, 'Color', 'k')
        xlim([0 round(binMax * 0.70)])
        title('Mouse Position From Center', 'FontSize', 16)
        xlabel('Distance From Center (cm)', 'FontSize', 14)
        ylabel('Probablity', 'FontSize', 14)
        t2.XAxis.FontSize = 12;
        t2.YAxis.FontSize = 12;
    if trialPhaseNum ~= 1
    t3 = nexttile(3);
        hold on; box on
        b = bar(distBins, [distFromAPDF', distFromBPDF'], 1.2);
        b(1).FaceColor = '#21B3DE';
        b(2).FaceColor = '#EA2B16';
        mx = max([distFromAPDF(:); distFromBPDF(:)]);
        yLimit = ylim;
        if contains(objectMovedNum, 'A')
            plot(distBins, distFromACDF' *mx, '-', 'LineWidth', 2, 'Color', [0.23 0.23 0.23]);
            plot(distBins, distFromACDF' *mx, '-.', 'LineWidth', 2, 'Color', '#21B3DE');
            plot(distBins, distFromBCDF' *mx, '-', 'LineWidth', 2, 'Color', '#EA2B16');
            hatchfill2(b(1), 'single', 'HatchAngle', 30, 'hatchcolor', [0.23 0.23 0.23]);
            [~, object_h, ~, ~] = legendflex(b, {strcat('Obj A (', distMoved{1}, ' cm)'), 'Obj A'}, 'FontSize', 12, 'Anchor', {'nw' 'nw'}, 'buffer', [5 -5], 'Box', 'off');
            hatchfill2(object_h(3), 'single', 'HatchAngle', 210, 'hatchcolor', [0.23 0.23 0.23], 'HatchDensity', 10);
        elseif contains(objectMovedNum, 'B')
            plot(distBins, distFromACDF' *mx, '-', 'LineWidth', 2, 'Color', '#21B3DE');
            plot(distBins, distFromBCDF' *mx, '-', 'LineWidth', 2, 'Color', [0.23 0.23 0.23]);
            plot(distBins, distFromBCDF' *mx, '-.', 'LineWidth', 2, 'Color', '#EA2B16');
            hatchfill2(b(2), 'single', 'HatchAngle', 30, 'hatchcolor', [0.23 0.23 0.23]);
            [~, object_h, ~, ~] = legendflex(b, {'Obj A', strcat('Obj B (', distMoved{1}, ' cm)')}, 'FontSize', 12, 'Anchor', {'nw' 'nw'}, 'buffer', [5 -5], 'Box', 'off');
            hatchfill2(object_h(4), 'single', 'HatchAngle', 210, 'hatchcolor', [0.23 0.23 0.23], 'HatchDensity', 10);
        elseif trialPhaseNum == 2
            legendflex(b, {'Obj A', 'Obj B'}, 'FontSize', 12, 'Anchor', {'nw' 'nw'}, 'buffer', [5 -5], 'Box', 'off');
        end
        xlim([0 binMax])
        t3.YLim = [yLimit(1), yLimit(2)];
        title('Mouse Position From Objects', 'FontSize', 16)
        xlabel('Distance From Objects (cm)', 'FontSize', 34)
        ylabel('Probablity', 'FontSize', 14)
        t3.XAxis.FontSize = 12;
        t3.YAxis.FontSize = 12;
    end
    t4 = nexttile(4);
        bar([pctNearWall  pctNearCtr], 'FaceColor', '#21B3DE')
        yLimits = ylim;
        if pctNearWall > yLimits(2) * 0.1
            text(1, 0.5 * pctNearWall, num2str(pctNearWall, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(1, pctNearCtr, num2str(pctNearWall, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        if pctNearCtr > yLimits(2) * 0.1
            text(2, 0.5 * pctNearCtr, num2str(pctNearCtr, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(2, pctNearCtr, num2str(pctNearCtr, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        xlim([0 3])
        title('Time Near Walls', 'FontSize', 16)
        xticklabels({'Near Wall', 'Open Field'})
        ylabel('% Time', 'FontSize', 14)
        t4.XAxis.FontSize = 14;
        t4.YAxis.FontSize = 12;
    t5 = nexttile(5);
        bar([probC1 probC2 probC4 probC3], 'FaceColor', '#21B3DE')
        yLimits = ylim;
        if probC1 > yLimits(2) * 0.1
            text(1, 0.5 * probC1, num2str(probC1, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(1, probC1, num2str(probC1, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        if probC2 > yLimits(2) * 0.1
            text(2, 0.5 * probC2, num2str(probC2, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(2, probC2, num2str(probC2, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        if probC4 > yLimits(2) * 0.1
            text(3, 0.5 * probC4, num2str(probC4, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(3, probC4, num2str(probC4, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        if probC3 > yLimits(2) * 0.1
            text(4, 0.5 * probC3, num2str(probC3, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(4, probC3, num2str(probC3, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        xlim([0 5])
        if probC1 == 0 && probC2 == 0 && probC3 == 0 && probC4 == 0
            ylim([0 1])
        end
        title('Time Near Corners', 'FontSize', 16)
        xticklabels({'Front Left', 'Front Right', 'Back Left', 'Back Right'})
        ylabel('% Time', 'FontSize', 14)
        t5.XAxis.FontSize = 14;
        t5.YAxis.FontSize = 12;
    if trialPhaseNum ~= 1
    t6 = nexttile(6);
        b1 = bar(1, probCloseA, 'FaceColor', '#21B3DE');
        hold on
        b2 = bar(2, probCloseB, 'FaceColor', '#EA2B16');
        yLimits = ylim;
        if probCloseA > yLimits(2) * 0.1
            text(1, 0.5 * probCloseA, num2str(probCloseA, '%.1f'), 'Color','k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(1, probCloseA, num2str(probCloseA, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        if probCloseB > yLimits(2) * 0.1
            text(2, 0.5 * probCloseB, num2str(probCloseB, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center')
        else
            text(2, probCloseB, num2str(probCloseB, '%.1f'), 'Color', 'k', 'FontSize', 18, 'Horiz', 'center', 'Vert', 'bottom')
        end
        xlim([0 3])
        if probCloseA == 0 && probCloseB == 0
            ylim([0 1])
        end
        yLimit = ylim;
        title('Time Near Objects', 'FontSize', 16)
        xticks([1 2])
        xticklabels({'Object A', 'Object B'})
        ylabel('% Time', 'FontSize', 14)
        t6.XAxis.FontSize = 14;
        t6.YAxis.FontSize = 12;
        if contains(objectMovedNum, 'A')
            xlabel("Moved " + distMoved{1} + " cm", 'Color', '#21B3DE')
            if probCloseA ~= 0
                hatchfill2(b1, 'single', 'HatchAngle', 30, 'hatchcolor', [0.23 0.23 0.23]);
            end
        elseif contains(objectMovedNum, 'B')
            xlabel("Moved " + distMoved{1} + " cm", 'Color', '#EA2B16')
            if probCloseB ~= 0
                hatchfill2(b2, 'single', 'HatchAngle', 30, 'hatchcolor', [0.23 0.23 0.23]);
            end
        end
        t6.YLim = [yLimit(1), yLimit(2)];
    end

    %% Prepare structure variables for output
    vid.PixelsInCm = pixelsInCm;

    taskParams = struct('TrialPhaseNum', trialPhaseNum, 'ObjectMovedNum', objectMovedNum, 'DistMoved', distMoved, 'StartFrame', startFrame, ...
    'EndFrame', endFrame, 'NumFrames', numFrames, 'CloseThresh', closeThresh, 'JumpLimit', jumpLimit);

    traj = struct('SmoothX', smoothX, 'SmoothY', smoothY, 'SmoothXCm', smoothXCm, 'SmoothYCm', smoothYCm, 'DistTravCm', distTravCm, ...
        'DistTravCmSum', distTravCmSum, 'MouseVelCm', mouseVelCm, 'MouseVelCmAve', mouseVelCmAve, 'TimeAxis', timeAxis, 'DistSum', distSum, 'VelAvg', velAvg);

    stats = struct('PctNearWall', pctNearWall, 'PctNearCtr', pctNearCtr, 'DistFromCtr', distFromCtr, 'DistFromCtrPDF', distFromCtrPDF, ...
        'DistFromCtrCDF', distFromCtrPDF, 'DistFromCtrAve', distFromCtrAve, 'HistCountsFrames', histCountsFrames, 'HistCountsSec', histCountsSec, ...
        'HistCountsPerTime', histCountsPerTime, 'PctCorner1', probC1, 'PctCorner2', probC2, 'PctCorner3', probC3, 'PctCorner4', probC4, 'PctCorner', nearCorner, ...
        'PctNotCorner', notCorner);

    if trialPhaseNum ~= 1
        stats.DistFromA         = distFromA;
        stats.DistFromB         = distFromB;
        stats.distFromAPDF      = distFromAPDF;
        stats.distFromBPDF      = distFromBPDF;
        stats.distFromACDF      = distFromACDF;
        stats.distFromBCDF      = distFromBCDF;
        stats.probCloseA        = probCloseA;
        stats.probCloseB        = probCloseB;
        stats.probCloseMoved    = probCloseMoved;
        stats.probCloseUnmoved  = probCloseUnmoved;
        stats.distFromAAve      = distFromAAve;
        stats.distFromBAve      = distFromBAve;
    end

    %% Save figures, Matlab summary, and to Excel sheet
    outputDir = fullfile(mainDir, "Script Output 06-21-23", "Cohort " + ExcelArr{n,'Cohort'}, "Animal " + ExcelArr{n,'AnimalID'});
    mkdir(outputDir);
    cd(outputDir);

    disp('Writing figures to disk')
        print(f2, vid.ID + " Traj.png", '-dpng')
        print(f3, vid.ID + " HeatMap.png", '-dpng')
        print(f4, vid.ID + " Stats.png", '-dpng')

    disp('Writing screencap video to disk')
        screencapName = vid.ID + " Screencap.mp4";
        screencapObj = VideoWriter(screencapName, 'MPEG-4');
        screencapObj.Quality = 100;
        open(screencapObj);
        writeVideo(screencapObj,screencapData); 
        close(screencapObj);

    disp('Writing summary matfile to disk')
        save(strcat(vid.ID, " Values.mat"), 'geomRun', 'vid', 'taskParams', 'traj', 'stats')

    if isfield(manual, 'FramesRange') == 0
        disp('Writing summary matfile to disk')
            fullExcel{MatlabRows(n), 'RunStatus'}       = {versionNum + ", " + string(datetime('today', 'Format', 'MM/dd/uuuu'))};
            fullExcel{MatlabRows(n), 'PctTimeWalls'}    = pctNearWall;
            fullExcel{MatlabRows(n), 'PctTimeMiddle'}   = pctNearCtr;
            fullExcel{MatlabRows(n), 'AvgDistCenter'}   = distFromCtrAve;
            fullExcel{MatlabRows(n), 'DistMovedMouse'}  = distTravCmSum;
            fullExcel{MatlabRows(n), 'AvgVelocity'}     = mouseVelCmAve;
            fullExcel{MatlabRows(n), 'PctCorner'}       = nearCorner;
            fullExcel{MatlabRows(n), 'PctNotCorner'}    = notCorner;
            fullExcel{MatlabRows(n), 'PctTimeC1'}       = probC1;
            fullExcel{MatlabRows(n), 'PctTimeC2'}       = probC2;
            fullExcel{MatlabRows(n), 'PctTimeC3'}       = probC3;
            fullExcel{MatlabRows(n), 'PctTimeC4'}       = probC4;
            fullExcel{MatlabRows(n), 'CloseThresh'}     = closeThresh;
            fullExcel{MatlabRows(n), 'JumpLimit'}       = jumpLimit;

            if fullExcel{MatlabRows(n),'Phase'} >= 2
                fullExcel{MatlabRows(n), 'PctTimeObjAVideo'}    = probCloseA;
                fullExcel{MatlabRows(n), 'PctTimeObjBVideo'}    = probCloseB;
                fullExcel{MatlabRows(n), 'TimeMovedScript'}     = probCloseMoved;
                fullExcel{MatlabRows(n), 'TimeUnmovedScript'}   = probCloseUnmoved;
                fullExcel{MatlabRows(n), 'AgeDistObjA'}         = distFromAAve;
                fullExcel{MatlabRows(n), 'AvgDistObjB'}         = distFromBAve;
            end
    end

    cd(mainDir)
    writetable(fullExcel, tableName)
    
	disp('File completed')
    close all
	toc
end
