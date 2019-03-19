# Constucts a lookup table for the different household configurations ready for household ODE solver

function CreateStatesTupleExpression(n::Int64,L1::Int64)
    str_start = "StatesTuple = Base.product("
    str = "0:$(L1)"
    for i = 2:n
      str = string("0:$(L1),",str)
    end
    str_end = ")"
    return parse(string(str_start,str,str_end))
end
ex = CreateStatesTupleExpression(6,MaxHouseholdSize)
eval(ex)

StatesTuple = collect(StatesTuple)
StatesTuple = StatesTuple[:]


ConstructStatesMatrix = function(StatesTuple::Vector{NTuple{6,Int64}},L1::Int64)
    y = map(sum,StatesTuple)
    y = (y.<= L1 ).*(y.>0)
    StatesTuple = StatesTuple[y]
    d1 = length(StatesTuple)
    d2 = length(StatesTuple[1])
    N_O1 = [StatesTuple[i][4]+StatesTuple[i][5]+StatesTuple[i][6] for i = 1:d1]
    N_U1 = [StatesTuple[i][1]+StatesTuple[i][2]+StatesTuple[i][3] for i = 1:d1]
    y = (N_O1 .<= MaxNumberOfO1s).*(N_U1 .<= MaxNumberOfU1s).*(N_O1 .> 0)
    StatesTuple = StatesTuple[y]

    d1 = length(StatesTuple)
    d2 = length(StatesTuple[1])
    States = zeros(Int64,d1,d2)
    for i = 1:d1
      for j = 1:d2
        States[i,j] = StatesTuple[i][j]
      end
    end
    return States,d1
end
States,d1 = ConstructStatesMatrix(StatesTuple,MaxHouseholdSize)
StatesTuple = 0

# Construct vectors of household state properties
N_vect = sum(States,2)
N1_vec = sum(States[:,1:3],2)
N2_vec = sum(States[:,4:6],2)
U1_state = N1_vec .> 0
