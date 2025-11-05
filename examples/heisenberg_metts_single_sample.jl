using ITensors
using ITensorMPS
using METTS


function my_measurements(psi0::MPS, H::MPO)

    energy = ITensors.inner(psi0', H, psi0)
    vne = entropy_von_neumann(psi0, length(psi0) ÷ 2)

    ### collect in dictionary
    measurements = Dict{AbstractString,Real}()

    measurements["energy"] = energy
    measurements["vne"] = vne
    return measurements
end

let

    N = 16
    J = 1.0

    ### define Hamiltonian
    sites = siteinds("S=1/2", N; conserve_qns=true)

    # Build Hamiltonian
    terms = OpSum()
    for i in 1:(N-1)
        terms += J / 2, "S+", i, "S-", i + 1
        terms += J / 2, "S-", i, "S+", i + 1
        terms += J, "Sz", i, "Sz", i + 1
    end

    H = MPO(terms, sites)


    # start a product state:
    states = [isodd(j) == 1 ? "Z+" : "Z-" for j in 1:N]
    psi0 = MPS(sites, states)


    results, psi0 = metts_single_temperature(my_measurements, psi0, 2.0, H, 0.5, 1E-6, 512, 0.1, 1)

end