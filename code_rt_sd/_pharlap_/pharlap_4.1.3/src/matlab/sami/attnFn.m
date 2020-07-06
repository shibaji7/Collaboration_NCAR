function [ ne,attn ] = attnFn( dene,lat,alt,isAtn )
    prompt = 'Center latitude of eclipse (deg):';
    d_lat = input(prompt);
    % d_lat = 0;
    
    S = size(lat);
    
    attn = 0.8/30*abs(lat-d_lat)+0.2;
    for i=1:S(1)
        for j=1:S(2)
            if (attn(i,j))>1
                attn(i,j)=1;
            end
        end 
    end
    if isAtn
        ne = dene.*attn;
    else
        ne = dene;
    end
end
