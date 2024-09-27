function [tfl] = func_InfiltrationWidth(vals, xmesh)
        
        for k = 1:length(vals)
            if(vals(1,k) >= (0.80 * max(vals)))
                ind80 = k;
            end
        end

        for k = length(vals):-1:1
            if(vals(1,k) <= (0.02 * max(vals)))
                ind02 = k;
            end
        end
        
        tfl = xmesh(1,ind02) - xmesh(1,ind80);
    end
