include("jTeachOPTmain.jl");
include("branchandbound.jl");

ukp = instance(5, zeros(5), zeros(5), 0)

ukp.c[1] = 15;
ukp.c[2] = 7;
ukp.c[3] = 14;
ukp.c[4] = 18;
ukp.c[5] = 17;

ukp.w[1] = 12;
ukp.w[2] = 7;
ukp.w[3] = 15;
ukp.w[4] = 24;
ukp.w[5] = 23;

ukp.W = 53;

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