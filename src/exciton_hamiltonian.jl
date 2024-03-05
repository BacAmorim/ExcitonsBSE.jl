struct BSEproblem{E, D, ED, Tlat, Tenergies, Tstates, Tint, W<:Function, V<:Function}
    Q::SVector{E, Tlat}
    real_lattice::SMatrix{E, D, Tlat, ED}
    reciprocal_lattice::SMatrix{E, D, Tlat, ED}
    unitcell_volume::Tlat
    sites::Vector{SVector{E, Tlat}}
    Nkpts::NTuple{D, Int}
    kpts::Vector{SVector{E, Tlat}}
    interaction_bare::V
    interaction_screened::W
    interaction_IR_regularization::Tint
    cutoff::Tlat
    Gstar::Vector{SVector{E, Tlat}}
    energies_valence::Array{Tenergies, 2}
    states_valence::Array{Tstates, 3}
    energies_conduction::Array{Tenergies, 2}
    states_conduction::Array{Tstates, 3}
end

function exciton_problem(qh::Quantica.Hamiltonian{Tlat, E, 2, Th}, Q, fermienergy, Nkpts, Vbare, Wscreened, cutoff) where {Tlat, E, Th}
    
    # extract lattice basis
    A = SMatrix{E, 2, Tlat, E*2}(qh.lattice.bravais.matrix)
    # unit cell volume
    UCvol = volumeUC(A)
    # extract reciprocal lattice basis
    B = transpose(2pi*pinv(A))
    # extract site positions
    positions = sites(lattice(qh))
    
    # span umklap vectors
    gcut = maximum([norm(B[:, n]) for n in axes(B, 2)]) + cutoff 
    Gstar = sites(supercell(lattice(sublat(@SVector zeros(Tlat, E)); bravais = B), region = g -> norm(g) <= gcut))
    
    # generate Monkhorst-Pack grid of k-points
    kpts = mpgrid(B, Nkpts)
    Nk = length(kpts)
    
    # determine number of occupied and empty bands
    dimh = size(qh[()], 1)
    E0, _ = spectrum(qh, @SVector zeros(2))
    Nvbands = count(e -> e < fermienergy, real(E0))
    Ncbands = dimh - Nvbands
    
    # calculate band structures
    Ekv = zeros(Nvbands, Nk)
    Ukv = zeros(Th, dimh, Nvbands, Nk)
    
    EkQc = zeros(Ncbands, Nk)
    UkQc = zeros(Th, dimh, Ncbands, Nk)
    
    for (ik, k) in enumerate(kpts)
        ek, psik = spectrum(qh, A'k)
        
        kQ = modlat(k + Q, B, A)
        ekQ, psikQ = spectrum(qh, A'kQ)

        for n in 1:dimh
            # valence band states
            if n <= Nvbands
                Ekv[n, ik] = real(ek[n])
                Ukv[:, n, ik] = psik[:, n]
            # conduction band states
            else 
                EkQc[n-Nvbands, ik] = real(ekQ[n])
                UkQc[:, n-Nvbands, ik] = psikQ[:, n]
            end
        end
        
    end
    
    W0 = regularize_IR(Wscreened, B[:, 1], B[:, 2], Nkpts[1], Nkpts[2], UCvol)
    
    return BSEproblem(Q, A, B, UCvol, positions, Nkpts, kpts, Vbare, Wscreened, W0, cutoff, Gstar, Ekv, Ukv, EkQc, UkQc)
    
end



function regularize_IR(W, b1, b2, Nkpts1, Nkpts2, UCvol)
    
    k0 = inradius(b1/Nkpts1, b2/Nkpts2)
    Wk0 = W(k0)
    
    Vol = UCvol*Nkpts1*Nkpts2
    
    function Wcut(u)
        
        k = u[1]*b1/Nkpts1 + u[2]*b2/Nkpts2 
        
        if norm(k) > k0
            return W(norm(k))
        else
            return Wk0
        end
        
    end
    
    W1, error1 = hcubature(Wcut, [-0.5, -0.5], [0.5, 0.5])

    W2, error2 = quadgk(k -> k*W(k)/(2pi), 0, k0)
    
    W3 = k0^2*Wk0/(4*pi)
    
    return W1 + Vol*(W2 + W3)
    
end         

function exciton_hamiltonian_tda!(h, Q, A, B, UCvol, positions, kpts, Vbare, Wscreened, W0, cutoff, Gstar, Ev, Uv, EQc, UQc)

    Nks = length(kpts)
    Nvbands = size(Ev, 1)
    Ncbands = size(EQc, 1)

    fill!(h, zero(eltype(h)))

    for jk in 1:Nks, cj in 1:Ncbands, vj in 1:Nvbands
        
        Jkcv = vj + (cj-1)*Nvbands + (jk-1)*Nvbands*Ncbands

        # non-interacting part
        h[Jkcv, Jkcv] = EQc[cj, jk] - Ev[vj, jk]
        

        for ik in 1:Nks, ci in 1:Ncbands, vi in 1:Nvbands

            Ikcv = vi + (ci-1)*Nvbands + (ik-1)*Nvbands*Ncbands

            ukQci = view(UQc, :, ci, ik)
            ukQcj = view(UQc, :, cj, jk)
            
            ukvi = view(Uv, :, vi, ik)
            ukvj = view(Uv, :, vj, jk)

            # Fock contribution
            Wij = zero(complex(eltype(Uv)))
            P = kpts[ik] - kpts[jk]
            p = modlat(kpts[ik] - kpts[jk], B, A)
            @show P, p, ik, jk
            # Hartree contribution
            Vij = zero(complex(eltype(Uv)))

            for G in Gstar

                # accumulate Fock
                if norm(p + G) < cutoff

                    rho_cicj = rho(ukQci, ukQcj, p + G, positions)
                    rho_vjvi = rho(ukvj, ukvi, -p - G, positions)

                    if p == zero(p) && G == zero(G)
                        Wij += W0*rho_cicj*rho_vjvi/(Nks*UCvol)
                    else
                        Wij += Wscreened(norm(p + G))*rho_cicj*rho_vjvi/(Nks*UCvol)
                    end

                end

                # accumulate Hartree
                if norm(Q + G) < cutoff

                    if Q != zero(Q) && G != zero(G)
                    
                        rho_civi = rho(ukQci, ukvi, Q + G, positions)
                        rho_vjcj = rho(ukvj, ukQcj, -Q - G, positions)

                        Vij += Vbare(norm(Q + G))*rho_civi*rho_vjcj/(Nks*UCvol)
                    end
                    
                end

            end

            h[Ikcv, Jkcv] += Wij - 2*Vij

        end

    end

    return h

end

function exciton_hamiltonian_tda!(h, bse::BSEproblem)
    
    # umpack
    Q = bse.Q
    A = bse.real_lattice
    B = bse.reciprocal_lattice
    UCvol = bse.unitcell_volume
    positions = bse.sites
    Nkpts = bse.Nkpts
    kpts = bse.kpts
    Vbare = bse.interaction_bare
    Wscreened = bse.interaction_screened
    W0 = bse.interaction_IR_regularization
    cutoff = bse.cutoff
    Gstar = bse.Gstar
    Ekv = bse.energies_valence
    Ukv = bse.states_valence
    EkQc = bse.energies_conduction
    UkQc = bse.states_conduction
    
    # build hamiltonian
    
    return exciton_hamiltonian_tda!(h, Q, A, B, UCvol, positions, kpts, Vbare, Wscreened, W0, cutoff, Gstar, Ekv, Ukv, EkQc, UkQc)
    
end

function exciton_hamiltonian_tda(bse::BSEproblem)
    
    Nvbands = size(bse.energies_valence, 1)
    Ncbands = size(bse.energies_conduction, 1)
    
    Nks = length(bse.kpts)
    
    N = Nks*Nvbands*Ncbands
    
    h = zeros(eltype(bse.states_valence), N, N)
    
    
    return exciton_hamiltonian_tda!(h, bse)
    
end

