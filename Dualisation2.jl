# PROJET ECMA : Ahmed Belmabrouk

using JuMP

using CPLEX

# Lire le fichier des instances.
function  LireInstance(nom_instance :: String ) 

# Si le fichier "instance.gr" existe
if isfile("C:/Users/ahmed/Desktop/ECMA Belmabrouk/Instances_ECMA/" * nom_instance )
    # L’ouvrir
    myFile = open("C:/Users/ahmed/Desktop/ECMA Belmabrouk/Instances_ECMA/" * nom_instance, "r")
    # Lire toutes les lignes d’un fichier
    data = readlines(myFile) #Retourne un tableau de String
    # n nombre de sommets 
    n_as_string=split(data[1]," ")[end]
    n=parse(Int, n_as_string)
    # s numero de ville de depart 
    s_as_string=split(data[2]," ")[end]
    s=parse(Int, s_as_string)
    # t numero de ville d'arrivée 
    t_as_string=split(data[3]," ")[end]
    t=parse(Int, t_as_string)
    # S poids max du chemin 
    S_as_string=split(data[4]," ")[end]
    S=parse(Int, S_as_string)
    # d1 l'augmentation totale de la durée des arcs est limitée par le parametre d1 
    d1_as_string=split(data[5]," ")[end]
    d1=parse(Int, d1_as_string)
    # d2 l'augmentation totale des poids du graphe est limitée par le parametre d2 
    d2_as_string=split(data[6]," ")[end]
    d2=parse(Int, d2_as_string)
    # p  vecteur des poids 
    pstr=data[7][6:1:length(data[7])-1]
    p_as_string=split(pstr,", ")
    p=[]
    for elt in p_as_string 
        append!(p,parse(Int64, elt))
    end
    # ph  vecteur des poids ph
    phstr=data[8][7:1:length(data[8])-1]
    ph_as_string=split(phstr,", ")
    ph=[]
    for elt in ph_as_string 
        append!(ph,parse(Int64, elt))
    end

    # Mat matrice des disctances Dij , dij
    Mat=[]
    for line in data[10:1:length(data)]
        line_as_string=split(line[1:1:length(line)-1]," ")
        coef=[]
        append!(coef,parse(Int64, line_as_string[1]))
        append!(coef,parse(Int64, line_as_string[2]))
        append!(coef,parse(Int64, line_as_string[3]))
        append!(coef,parse(Float64, line_as_string[4]))
        append!(Mat,[coef])
    end
    # Fermer le fichier
    close(myFile)
    end
return (n,s,t,S,d1,d2,p,ph,Mat)
end
(n,s,t,S,d1,d2,p,ph,Mat)=LireInstance("20_USA-road-d.NY.gr")
"""
println(n) 
println(s)
println(t)
println(S)
println(d1)
println(d2)
println(p)
println(ph)
println(Mat)

"""
#***************Dualisation****************************
m1 = Model(CPLEX.Optimizer)
#Définition des variables
@variable(m1, x[i in 1:length(Mat)], Bin)
@variable(m1, z[i in 1:length(Mat)] >=0 )
@variable(m1, t1 >=0 )
@variable(m1, t2 >=0 )
@variable(m1, y[i in 1:n], Bin)
@variable(m1, z1[i in 1:n] >=0 )
#Définition des contraintes

@constraint(m1, (sum(x[i] for i in 1:length(Mat) if Mat[i][1]==s)==1))
@constraint(m1, (sum(x[i] for i in 1:length(Mat) if Mat[i][2]==s)==0))
@constraint(m1, (sum(x[i] for i in 1:length(Mat) if Mat[i][2]==t)==1))
@constraint(m1, (sum(x[i] for i in 1:length(Mat) if Mat[i][1]==t)==0))
@constraint(m1,[v in 1:n ; (v≠t && v≠s)],(sum(x[i] for i in 1:length(Mat) if Mat[i][2]==v) )==(sum(x[i] for i in 1:length(Mat) if Mat[i][1]==v) ))
@constraint(m1,[v in 1:n ; v≠t], y[v]==(sum(x[i] for i in 1:length(Mat) if Mat[i][1]==v) ))
@constraint(m1, y[t]==(sum(x[i] for i in 1:length(Mat) if Mat[i][2]==t) ))
@constraint(m1, (sum(p[i] * y[i] + 2 * z1[i] for i in  1:n) + t2 * d2 ) <= S)
@constraint(m1,[i in 1:length(Mat)],(t1 + z[i]) >= Mat[i][3] * x[i])
@constraint(m1,[i in 1:n],(t2 + z1[i]) >= ph[i] * y[i])
#Définition de l'objectif
@objective(m1, Min, d1*t1+sum(Mat[i][3]*x[i] + Mat[i][4]*z[i] for i in 1:length(Mat)))


optimize!(m1)
vX = value.(x)
for i in 1: length(Mat)
    if vX[i]==1
        println(Mat[i][1], " ",Mat[i][2])
    end
end
print("objective value:  ", objective_value(m1))

