function [y_dst] = Lagrange(x_scr,y_scr,x_dst, pointcount)
% let's assume both sects are regular
% and x_coarse(1) < x_fine(1)
h_dst = x_dst(2) -x_dst(1);
h_scr = x_scr(2) - x_scr(1);
if h_dst > h_scr
    % from fine mesh to coarse
    rel = round(h_dst/h_scr);
    scheme = 0:pointcount-1;
    shift = round(rel/2);
    cond = 1:length(x_dst);
else
    % from coarse mesh to fine
    scheme = 1:pointcount;
    shift = 0;
    rel = 1;
    cond = 1:length(x_dst);
end

y_dst = zeros(1,length(x_dst));

for m = cond
    if h_dst > h_scr
        if m~=1 && shift + rel < length(x_scr)
            shift = shift + rel;    
        end
        if shift + rel > length(x_scr)
            scheme = (1-pointcount):0;    
        end
    else
        if x_scr(shift+2) < x_dst(m) && shift + pointcount < length(x_scr)
            shift = shift + rel;
        end
    end

    
    D = ones(1,pointcount);
    for i = 1:pointcount
        for j = 1:pointcount
            if(i~=j)
                D(i) =D(i) * (x_scr(scheme(i) + shift) - x_scr(scheme(j) + shift));
            end %if (i~=j)
        end % for j
    end %for i
%%
    for i = 1:pointcount
        multiplication = y_scr(scheme(i) + shift);
        for j = 1:pointcount
            if i~=j
                multiplication = multiplication * (x_dst(m) - x_scr(scheme(j) + shift));
            end % if (i~=j)
        end %for j
        y_dst(m) = y_dst(m) + multiplication/D(i);
    end % for i 
end %for m

end