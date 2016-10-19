include("jTeachOPTmain.jl");
include("branchandbound.jl");

ukp = instance(7, zeros(7), zeros(7), 0)

ukp.c[1] = 70;
ukp.c[2] = 20;
ukp.c[3] = 39;
ukp.c[4] = 37;
ukp.c[5] = 7;
ukp.c[6] = 5;
ukp.c[7] = 10;

ukp.w[1] = 31;
ukp.w[2] = 10;
ukp.w[3] = 20;
ukp.w[4] = 19;
ukp.w[5] = 4;
ukp.w[6] = 3;
ukp.w[7] = 6;

ukp.W = 50;

#solution gloutone
sol = computeGreedySolutionUKP(ukp);

#affichage solution gloutone
print("solution gloutone : \n valeur :");
print(sol.z);
print("\n variable a 1 : ");
for i=1:length(sol.v1)
	print(sol.v1[i]);
	print(" - ");
end

#borne duale
print("\nborne duale : ");
sGreedy = deepcopy(sol);
print(computeLinearRelaxationUKP(ukp));

#amelioration heuristique
prec = sol.z - 1;

while (sol.z > prec)
	sol = heuExplore(sol);
	prec = sol.z;
end	

#affichage solution ameliore
print("\nsolution ameliore : \n valeur :");
print(sol.z);
print("\n variable a 1 : ");
for i=1:length(sol.v1)
	print(sol.v1[i]);
	print(" - ");
end

#solution exacte par branch and bound
sol = solution(zeros(Int64, ukp.n), [], [], 0, 0, 0);
sol = branchandbound(ukp, sol, 0, sol);
print("\nsolution exacte : \n valeur :");
print(sol.z);
for i=1:length(sol.v1)
	print(sol.v1[i]);
	print(" - ");
end