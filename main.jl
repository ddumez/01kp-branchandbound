include("jTeachOPTmain.jl");
include("branchandbound.jl");

function heuExplore(s::solution)
  essais = 10 # nombre de voisins construits
  for essai = 1:essais
    sVoisin = deepcopy(swap(s))
    if sVoisin.z > s.z
      s = deepcopy(sVoisin)
    end
  end
  return s;
end

function trier(ukp::instance)
	for i=1:ukp.n
		for j=1:ukp.n-1
			if (ukp.c[j]/ukp.w[j] <= ukp.c[j+1]/ukp.w[j+1])
				tmp = ukp.c[j]
				ukp.c[j] = ukp.c[j+1]
				ukp.c[j+1] = tmp
				tmp = ukp.w[j]
				ukp.w[j] = ukp.w[j+1]
				ukp.w[j+1] = tmp
			end
		end
	end
	return ukp;
end

ukp = instance(7, zeros(7), zeros(7), 7)

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
print("\n objet utilise : ");
for i=1:ukp.n-1
	print(sol.x[i]);
	print(" - ");
end
print(sol.x[ukp.n],"\n \n");

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
print("\n objet utilise : ");
for i=1:ukp.n-1
	print(sol.x[i]);
	print(" - ");
end
print(sol.x[ukp.n],"\n \n");

#solution exacte par branch and bound
sol = solution(zeros(Int64, ukp.n), [], [], 0, ukp.W, 0);
ukp = trier(ukp);
sol = branchandbound(ukp, sol, 1, solution(zeros(Int64, ukp.n), [], [], 0, ukp.W, 0));
print("\nsolution exacte : \n valeur :");
print(sol.z);
print("\n objet utilise : ");
for i=1:ukp.n-1
	print(sol.x[i]);
	print(" - ");
end
print(sol.x[ukp.n],"\n \n");


#experimentations numeriques
for i=1:8
	f = open("./instances/$i/p0$(i)_c.txt");
	ukp.W = parse(Int, readline(f));
	ukp.n = parse(Int, readline(f));
	close(f)
	ukp = instance(ukp.n, zeros(ukp.n), zeros(ukp.n), ukp.W)
	f = open("./instances/$i/p0$(i)_w.txt");
	f2 = open("./instances/$i/p0$(i)_p.txt");
	f3 = open("./instances/$i/p0$(i)_s.txt");
	opt = []
	for k=1:ukp.n
		ukp.w[k] = parse(Int64, readline(f));
		ukp.c[k] = parse(Int64, readline(f2));
		opt = push!(opt,parse(Int64, readline(f3)));
	end
	close(f);
	close(f2);
	close(f3);
	valopt = 0;
	for k = 1:ukp.n
		valopt = valopt + ukp.c[k] * opt[k];
	end
	ukp = trier(ukp);
	sol = solution(zeros(Int64, ukp.n), [], [], 0, ukp.W, 0);
	sol = @time branchandbound(ukp, sol, 1, solution(zeros(Int64, ukp.n), [], [], 0, ukp.W, 0));
	sol2 = solution(zeros(Int64, ukp.n), [], [], 0, ukp.W, 0);
	sol2 = @time branchandbound2(ukp, sol2, 1, solution(zeros(Int64, ukp.n), [], [], 0, ukp.W, 0));
	print("probleme $(i) (",valopt,") : ",sol.z," : ",sol.x," ; ou ",sol2.z," : ",sol2.x);
	print("\n");
end
