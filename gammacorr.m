xRes=1920;yRes=1080;
DScreen=20;%in cm
VertCRTSize=30;%in cm
VertScreenDimDeg=atan(VertCRTSize/DScreen)*360/(2*pi);
PixperDeg=yRes/VertScreenDimDeg;
SensorSize=40;

RFPosxDeg=20;
RFPosyDeg=10;
RFSizeDeg=15;

RFRadiusPix=ceil(RFSizeDeg*PixperDeg/2);

x0=floor(xRes/2 + RFPosxDeg*PixperDeg - RFSizeDeg*PixperDeg/2);
y0=floor(yRes/2 - RFPosyDeg*PixperDeg - RFSizeDeg*PixperDeg/2);


AssertOpenGL;
screens=Screen('Screens');
screenNumber=max(screens);
frameRate=Screen('FrameRate',screenNumber);
if(frameRate==0)  %if MacOSX does not know the frame rate the 'FrameRate' will return 0.
    frameRate=60;
end
white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);
gray=(white+black)/2;
if round(gray)==white
    gray=black;
end
inc=white-gray;
Screen('Preference', 'VBLTimestampingMode', -1);
Screen('Preference','SkipSyncTests', 0);
[w,rect]=Screen('OpenWindow',screenNumber);%, 0,[],16,2);
% load('GammaTable.mat');
% CT=(ones(3,1)*correctedTable(:,2)')'/255
% Screen('LoadNormalizedGammaTable',w, CT);

L=0:32:256; L(9) = 255;
for i=4:7
    bg=ones(yRes,xRes)*L(i);%bg(1:SensorSize,1:SensorSize)=black;
    BG=Screen('MakeTexture', w, bg);
    Screen('DrawTexture',w, BG);
    Screen('Flip', w);
    WaitSecs(15);
end



Screen('CloseAll');
Priority(0);