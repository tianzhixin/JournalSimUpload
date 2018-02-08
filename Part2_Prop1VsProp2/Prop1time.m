
function [suc] = Prop1time(Q,no_region,T,A_adj)
z = tf('z',T);
P = inv(z*eye(no_region) - Q);
A_adj_nd = A_adj - eye(no_region);
suc = 1;


A_adj_nd_triu = triu(A_adj_nd);
for Ai = 1:no_region
        for Aj = 1:no_region
            if A_adj_nd_triu(Ai,Aj) == 1
                H = tf(P.Numerator{Aj,Ai},P.Numerator{Ai,Ai},60);
                [ninf,fpeak] = norm(H,inf,0.0000001);
                if ninf >= 1
%                     [ninf root Ai Aj]
                    suc = 0
                    break;
                end

                H = tf(P.Numerator{Ai,Aj},P.Numerator{Aj,Aj},60);
                [ninf,fpeak] = norm(H,inf,0.0000001);
                if ninf >= 1
%                     [ninf root Aj Ai]
                    suc = 0
                    break;
                end
            end
        end
        
        if suc == 0
            break;
        end
end

end
