GM_WM_filled = false(size(GM_WM));
for n = 1:size(GM_WM,3)
    if any(any(GM_WM(:,:,n)))
        for m = 1:size(GM_WM,2)
            if any(GM_WM(:,m,n))
                GM_WM_filled(find(GM_WM(:,m,n),1,'first'):find(GM_WM(:,m,n),1,'last'),m,n) = true;
            end
        end
    end
end
% This moves through all planes in the 'z' direction, and within each
% searches for columns with at least one point filled. It then fills the
% gap between the first and last point found in that column so that all
% space is filled within the volume. This is called 'GM_WM_filled'
    