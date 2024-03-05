function mpgrid(B::SMatrix{E, D, T, ED}, Nkpts) where {D, E, T, ED}
    
    N = prod(Nkpts)
    kpts = SVector{E, T}[]
    
    for ns in CartesianIndices(Nkpts)
            
        k = @SVector zeros(T, E)
        
        for (i, ni) in enumerate(Tuple(ns))
            k += (2*ni - Nkpts[i] - 1)/(2*Nkpts[i])*B[:, i]
        end
        
        push!(kpts, k)
        
    end
    
    return kpts
    
end

function fractional(x)
    return x - ceil(x-0.5)
end

function modlat(x, A, B)
    u = transpose(B)*x/(2pi)
    return A*fractional.(u)
end

function volumeUC(A::SMatrix{E, 2, Tlat, ED}) where {E, Tlat, ED}
    
    return abs(A[1, 1]*A[2, 2] - A[2, 1]*A[1, 2])
    
end

function volumeUC(A::SMatrix{3, 3, Tlat, ED}) where {Tlat, ED}
    
    return abs(det(A))
    
end

function inradius(a1, a2)
    
    len = zero(a1)
    if norm(a1) > norm(a2)
        len = norm(a1)
    else
        len = norm(a2)
    end
    
    area = abs(a1[1]*a2[2] - a1[2]*a2[1])
    
    return area/len/2
    
end