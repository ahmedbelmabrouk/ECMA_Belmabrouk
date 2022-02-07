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
(n,s,t,S,d1,d2,p,ph,Mat)=LireInstance("20_USA-road-d.BAY.gr")
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
#***************Plan_Coupants****************************
#Probleme maitre
m = Model(CPLEX.Optimizer)

@variable(m, x[i in 1:length(Mat)] , Bin)
@variable(m, y[i in 1:n], Bin)
@variable(m, z)
#initialisation
d_1 =[]
for j in 1:length(Mat)
    append!(d_1,Mat[j][3])
end
p2=p
#Définition des contraintes
@constraint(m, (sum(d_1[i] * x[i] for i in 1:length(Mat))) <= z )
@constraint(m, (sum(x[i] for i in 1:length(Mat) if Mat[i][1] == s) )==1)
@constraint(m, (sum(x[i] for i in 1:length(Mat) if Mat[i][2]==t) )==1)
@constraint(m,[v in 1:n ; (v≠s & v≠t)],(sum(x[i] for i in 1:length(Mat) if Mat[i][2]==v) )==(sum(x[i] for i in 1:length(Mat) if Mat[i][1]==v) ))
@constraint(m,[v in 1:n ; v≠t], y[v]==(sum(x[i] for i in 1:length(Mat) if Mat[i][1]==v) ))
@constraint(m, y[t]==(sum(x[i] for i in 1:length(Mat) if Mat[i][2]==t) ))
@constraint(m, (sum(p2[i] * y[i]  for i in  1:n) <= S))

#Définition de l'objectif
@objective(m, Min, z)
optimize!(m)

#**********sous probleme 1************
m1 = Model(CPLEX.Optimizer)

@variable(m1, sigma1[i in 1:length(Mat)] >=0)
@constraint(m1, (sum(sigma1[i] for i in 1:length(Mat))) <=d1)
@constraint(m1 , [i in 1:length(Mat)] , sigma1[i] <= Mat[i][4] )
x_cour=value.(x)
@objective(m1, Max,sum(Mat[i][3]*(1+sigma1[i])*x_cour[i] for i in 1:length(Mat)))


#***********sous problme 2************
m2 = Model(CPLEX.Optimizer)

@variable(m2, sigma2[i in 1:n] >=0)
@constraint(m2, (sum(sigma2[i] for i in 1:n) <= d2))
@constraint(m2 , [i in 1:n] , sigma2[i] <= 2 )
y_cour=value.(y)
@objective(m2, Max,sum((p[i]+ph[i]*sigma2[i])*y_cour[i] for i in 1:n))



#boucle 
while true 
    optimize!(m)
    optimize!(m1)
    vsigma1=value.(sigma1)
    optimize!(m2)
    vsigma2=value.(sigma2)
    if ((objective_value(m) == objective_value(m1)) & (objective_value(m2) <= S ))
        """
        vX = value.(x)
        for i in 1: length(Mat)
            if vX[i]==1
                println(Mat[i][1],Mat[i][2])
            end
        end
        """
        break
    else #Mise a jour de d_1
        for i in 1:length(Mat)
            d_1[i]=Mat[i][3]*(1 + vsigma1[i])
        end
        for i in 1:n
            p2[i]=ph[i]+vsigma2[i]
        end
    end
end