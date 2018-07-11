function [speed, ta] = get_wheelspeed(wheeldata)

%find when rising or falling edge is happening
d = diff(wheeldata);
% only look at rising edges
d(find(d==-1)) = 0;
% try 1s windows
for i = 1:floor(length(d)/30) % because it's sampled at 30kHz
    speed(i) = sum(d((i-1)*30+1:i*30));
end
ta = 1:length(speed);


% for i = 1:floor((length(d)-30000)/3000)
%     beg = (i-1)*3000+1;
%     s(i) = numel(find(d(beg:beg+30000)));
% end
