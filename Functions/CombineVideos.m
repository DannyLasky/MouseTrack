% Set up parameters
Vid1 = VideoReader('10-27-21_DF221_3_Part1.mp4'); % Video 1
Vid2 = VideoReader('10-27-21_DF221_3_Part2.mp4'); % Video 2

v = VideoWriter('MergedVideo','MPEG-4'); % Create new video file
outputFrameRate = (Vid1.FrameRate * Vid1.Duration + Vid2.FrameRate * Vid2.Duration) / (Vid1.Duration + Vid2.Duration);
v.FrameRate = outputFrameRate;
open(v)

% Iterate on all frames in video 1 and write one frame at a time
while hasFrame(Vid1) 
    Video1 = readFrame(Vid1); % read each frame
    writeVideo(v,Video1) % write each frame
end
% Iterate again in video 2
while hasFrame(Vid2)
    Video2 = readFrame(Vid2);
    writeVideo(v,Video2)
end
close(v)