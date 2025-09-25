module SeqGen

export normal_SeqGen, golden_SeqGen, silver_SeqGen, thue_morse_SeqGen, cut_to_length, golden_LengthCalc, silver_LengthCalc, thue_morse_LengthCalc, plastic_LengthCalc

function normal_SeqGen(N::Int)
    return fill(1, N)
end

function golden_SeqGen(N::Int)
    sequence = "A"
    while length(sequence) < N
        sequence = replace(sequence, "A" => "AB", "B" => "A")
    end
    # println(sequence)
    
    number_sequence = [ch == 'A' ? 1 : 2 for ch in sequence]

    golden_phi = (1 + sqrt(5)) / 2
    tol = 1e-2
    ratio, is_valid = check_PV_ratio(number_sequence, golden_phi, tol)
    if !is_valid
        @warn "Ratio $ratio is outside the tolerance of φ ≈ $golden_phi (tol=$tol)."
    end
    
    return number_sequence
end

function silver_SeqGen(N::Int)
    sequence = "A"
    while length(sequence) < N
        sequence = replace(sequence, "A" => "BAA", "B" => "A")
    end
    # println(sequence)
    
    number_sequence = [ch == 'A' ? 1 : 2 for ch in sequence]

    silver_phi = 1 + sqrt(2)
    tol = 1e-2
    ratio, is_valid = check_PV_ratio(number_sequence, silver_phi, tol)
    if !is_valid
        @warn "Ratio $ratio is outside the tolerance of φ ≈ $silver_phi (tol=$tol)."
    end
    
    return number_sequence
end

function thue_morse_SeqGen(n::Int)
    if n <= 0 || (n & (n - 1)) != 0
        throw(ArgumentError("n must be a power of 2."))
    end
    
    sequence = Int[]
    for i in 0:n-1
        num_ones = count(c -> c == '1', string(i, base=2))
        push!(sequence, num_ones % 2 + 1)
    end

    thue_morese_phi = 1.0
    tol = 1e-2
    ratio, is_valid = check_PV_ratio(sequence, thue_morese_phi, tol)
    if !is_valid
        @warn "Ratio $ratio is outside the tolerance of φ ≈ $thue_morese_phi (tol=$tol)."
    end
    
    return sequence
end

function plastic_SeqGen(N::Int)
    sequence = "A"
    while length(sequence) < N
        sequence = replace(sequence, "A" => "B", "B" => "AC", "C" => "A") # Van der Laan word generation
        # sequence = replace(sequence, "A" => "AB", "B" => "AC", "C" => "A") # ChatGPT
    end
    
    number_sequence = [ch == 'A' ? 1 : ch == 'B' ? 2 : 3 for ch in sequence]

    plastic_rho = 1.3247179572447460259609088
    tol = 1e-2
    ratio_AB, ratio_BC, AB_is_valid, BC_is_valid = check_plastic_ratio(number_sequence, plastic_rho, tol)
    if !AB_is_valid
        error("Ratio $ratio_AB is outside the tolerance of ρ ≈ $plastic_rho (tol=$tol).")
    end
    if !BC_is_valid
        error("Ratio $ratio_BC is outside the tolerance of ρ ≈ $plastic_rho (tol=$tol).")
    end

    # if !AB_is_valid || !BC_is_valid
    #     error("Ratio $(ifelse(!AB_is_valid, ratio_AB, ratio_BC)) is outside the tolerance of ρ ≈ $plastic_rho (tol=$tol).")
    # end
    
    return number_sequence
end

function check_PV_ratio(number_sequence::Vector{Int}, phi::Float64, tol::Float64)
    count_A = count(x -> x == 1, number_sequence)
    count_B = count(x -> x == 2, number_sequence)
    
    ratio = count_A / count_B
    bool = abs(ratio - phi) ≤ tol
    
    return ratio, bool
end

function check_plastic_ratio(number_sequence::Vector{Int}, rho::Float64, tol::Float64)
    count_A = count(x -> x == 1, number_sequence)
    count_B = count(x -> x == 2, number_sequence)
    count_C = count(x -> x == 3, number_sequence)

    ratio_AB = count_A / count_B
    ratio_BC = count_B / count_C

    bool_AB = abs(ratio_AB - rho) ≤ tol
    bool_BC = abs(ratio_AB - rho) ≤ tol

    return ratio_AB, ratio_BC, bool_AB, bool_BC
end

function cut_to_length(sequence::Vector{Int}, sequence_length::Int)
    cut_sequence = Vector{Int}(undef, sequence_length)

    for i in 1:sequence_length
        cut_sequence[i] = sequence[i]
    end

    return cut_sequence
end

function golden_LengthCalc(exp_range::Vector{Int})
    phi = (1 + sqrt(5) ) / 2
    phi_conj = (1 - sqrt(5) ) / 2
    sequence_lengths = Int[]

    for g in exp_range
        length = real((phi^(g+2) - phi_conj^(g+2)) / sqrt(5))
        length = floor.(Int, length)
        push!(sequence_lengths, length)
    end

    return sequence_lengths
end

function silver_LengthCalc(n_range::Vector{Int})
    lengths = Int[]
    for n in n_range
        if n == 0
            push!(lengths, 1)
        elseif n == 1
            push!(lengths, 3)
        else
            L_prev_2, L_prev_1 = 1, 3
            for _ in 2:n
                L_current = 2 * L_prev_1 + L_prev_2
                L_prev_2, L_prev_1 = L_prev_1, L_current
            end
            push!(lengths, L_prev_1)
        end
    end
    return lengths
end

function thue_morse_LengthCalc(n_range::Vector{Int})
    lengths = Int[]
    for n in n_range
        length = 2^n
        push!(lengths, length)
    end
    return lengths
end

function plastic_LengthCalc(x_range::Vector{Int})
    lengths = Int[]
    for n in x_range
        sequence = "A"
        for _ in 1:n
            sequence = replace(sequence, "A" => "B", "B" => "AC", "C" => "A") # Van der Laan word generation
        end
        # println(typeof(sequence))
        number_sequence = [ch == 'A' ? 1 : ch == 'B' ? 2 : 3 for ch in sequence]
        length = size(number_sequence,1)
        push!(lengths, length)
    end
    return lengths
end




# # # Example usage

# # # Create pairing parameter sequence for N sites
# sequence_length = 1000

# # Normal Crystal
# N_normal = 1000
# normal_sequence = normal_crystal_gen(N_normal)
# normal_sequence = cut_to_length(normal_sequence, sequence_length)

# # Golden Ratio
# N_gold = 1000
# golden_sequence = golden_ratio_sequence_gen(N_gold)
# golden_sequence = cut_to_length(golden_sequence, sequence_length)

# # Silver Ratio
# N_silver = 1000
# silver_sequence = silver_ratio_sequence_gen(N_silver)
# silver_sequence = cut_to_length(silver_sequence, sequence_length)

# # Thue-Morse 
# N_thue_morse = 1024
# thue_morse_sequence = thue_morse_sequence_gen(N_thue_morse)
# thue_morse_sequence = cut_to_length(thue_morse_sequence, sequence_length)

# # Plastic Ratio
# N_plastic = 1000
# plastic_sequence = plastic_ratio_sequence_gen(N_plastic)
# plastic_sequence = cut_to_length(plastic_sequence, sequence_length)

# println("finished")

end # end module