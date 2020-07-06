function [ D ] = get_GC_distance(lat,lon,alt)

    D = zeros(size(lat));
    for j = 2:length(D)
        D(j) = deg2km(lat(j)-lat(1));
    end

end