function th=EWT_maxcheckplateau(L)

%find the max
[m,th]=max(L);

if isempty(th)
    %take care of the case there is a single mode in the length histogram
    th=1;
else
    if (th<length(L)) && (L(th+1)==m)  %there is a plateau
        thb=th;
        while (L(thb+1)==m) && (thb<length(L)-1) %search the last index of the plateau
            thb=thb+1;
        end
        th=floor((th+thb)/2);
    end
end


