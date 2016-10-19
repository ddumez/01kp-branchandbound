function completeGreedy(ukp, sol, k)
    # ---
    # Calcule la solution gloutonne avec pour utilite : u(i) = c(i)/p(i)
    sGreedy = sol;
    sommeW = ukp.w - sol.r;

    for i = (k+1):ukp.n
      if (sommeW + ukp.w[i] <= ukp.W)
        sGreedy.x[i] = 1
        sommeW = sommeW + ukp.w[i]
      end
    end

    sGreedy.z = fctUKP(ukp.c, sGreedy.x)
    sGreedy.r = ukp.W - dot(ukp.w, sGreedy.x)
    sGreedy.somX=sum(sGreedy.x) # somme des x_i = 1 
    return sGreedy

end

function completeLinear(ukp, sol, k)
    # identify the last item integrally selected
    sommeW = ukp.w - sol.r ; s = k+1;
    while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
      sommeW = sommeW + ukp.w[s]
      s = s + 1
    end
    s = s - 1
                
    # compute the upper bound
    if (ukp.W - dot(ukp.w[1:s], sGreedy.x[1:s])) > 0
      # contrainte non saturee => ajout de la partie fractionnaire de l'item bloquant
      zRelax = sum(ukp.c[1:s]) + (ukp.W - dot(ukp.w[1:s], sGreedy.x[1:s])) * (ukp.c[s+1] ./ ukp.w[s+1])
    else
      # contrainte saturee => rien a faire
      zRelax = sum(ukp.c[1:s])
    end
    return zRelax
end

# instance , solution en cour , indice de la variable considere , meilleur solution actuele
function branchandbound(ukp, sol, k, best)
	print("\nsolution : valeur :");
  print(sol.z);
  for i=1:ukp.n
    print(sol.x[i]);
    print(" - ");
  end
  print("\n");

  if(k == ukp.n)
    if (sol.z > ukp.z)
      best = sol;
    end 
    #backtrack
    while (k>0) && (sol[k] == 0)
      k = k-1;
    end
    if (0 == k)
      #plus de backtrack possible
      return best;
    else
      sol[k] = 0;
      return branchandbound(ukp, sol, k, best);
    end
  else 
    if(completeLinear(ukp, sol, k) > best.z)
      return branchandbound(ukp, completeGreedy(ukp, sol, k) , k, best);
    else
      #backtrack
      while (k>0) && (sol[k] == 0)
        k = k-1;
      end
      sol[k] = 0;
      return branchandbound(ukp, sol, k, best);
    end
  end
end