using AlgebraicPetri

### Predator Prey
# B -> []  
# [] -> A 
# A + 2B <-> 3B 

predatorPrey = LabelledPetriNet([:A, :B], :death=>(:B=>()), :birth=>(()=>:A), :predation=>((:A,:B,:B)=>(:B,:B,:B)))
 
### Replicator
# 2F <-> B_2 
# B_2 + F <-> B_3 
# B_3 + F <-> B_4 
# 2F <-> A_2 
# A_2 + F <-> A_3 
# A_3 + F <-> A_4


### Brusselator
# A -> X
# 2X + Y -> 3X
# B + X -> Y + D
# X -> E

brusselator = LabelledPetriNet([:A, :B, :D, :X, :Y, :E], :a=>(:A=>:X), :b=>((:X, :X, :Y)=>(:X, :X, :X)), :c=>((:B, :X)=>(:Y,:D)), :d=>(:X=>:E))

### Substrate-depletive Clock
# A -> C
# B + C -> prod

clock = LabelledPetriNet([:A, :B, :C], :a=>(:A=>:C), :b=>((:B, :C)=>()))

simplest = LabelledPetriNet([:A, :B], :fwd=>(:A=>:B), :rev=>(:B=>:A))
