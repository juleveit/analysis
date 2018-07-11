function stimons = get_timestamps(trigsig)

times = find(trigsig);
offs = find(diff(times)>500);

% for new data use this version
if ~isempty(offs)
    stimons = times(offs+1);
    stimons = [times(1);stimons];
else
    stimons = [];
end

% % for the contaminated data from the first recording day run this instead
% i = 1;
% stimons = [];
% while i < length(times)
%     if length(intersect(times,times(i)+1:times(i)+5)) == 5 % all following five samples were also high
%         stimons = [stimons, times(i)];
%         while find(times == times(i)+1)
%             i = i+1;
%         end
%     else
%         i = i+1;
%     end
% end
