function [ Ne ] = get_intp_data( alt,lat,Ne,intp )
    intp.intrp = scatteredInterpolant(alt(:),lat(:),Ne(:));
    Ne = zeros(size(intp.alt));
    S = size(intp.alt);
    for i = 1:S(1)
        for j = 1:S(2)
            if intp.alt(j,i) >= intp.threshold
                ne = intp.intrp(intp.alt(j,i),intp.lat(i,j));
                if ne < 0
                    ne = 0;
                end
                Ne(i,j) = ne;
            else
                Ne(i,j) = 0.0;
            end
        end
    end
end