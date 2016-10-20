function completeGreedy(ukp, sol, k)
    # ---
    # Calcule la solution gloutonne avec pour utilite : u(i) = c(i)/p(i)
    sommeW = 0;
    for i=1:k
      sommeW = sommeW + sol.x[i]*ukp.w[i];
    end

    for i = (k+1):ukp.n
      if (sommeW + ukp.w[i] <= ukp.W)
        sol.x[i] = 1
        sommeW = sommeW + ukp.w[i]
      end
    end

    sol.z = fctUKP(ukp.c, sol.x)
    sol.r = ukp.W - dot(ukp.w, sol.x)
    sol.somX=sum(sol.x) # somme des x_i = 1
    return sol;
end

function completeLinear(ukp, sol, k)
    # identify the last item integrally selected
    sommeW = 0;
    for i=1:k
      sommeW = sommeW + sol.x[i]*ukp.w[i];
    end
    s = k+1;
    while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
      sommeW = sommeW + ukp.w[s]
      s = s + 1
    end
    s = s - 1;
    # compute the upper bound
    zRelax = sum(ukp.c[1:s]);
    if (s < ukp.n) && ((ukp.W - dot(ukp.w[1:s], sGreedy.x[1:s])) > 0)
      # contrainte non saturee => ajout de la partie fractionnaire de l'item bloquant
      zRelax = zRelax + (ukp.W - dot(ukp.w[1:s], sGreedy.x[1:s])) * (ukp.c[s+1] ./ ukp.w[s+1]);
    end
    return zRelax;
end

# instance , solution en cour , indice de la variable considere , meilleur solution actuele
function branchandbound(ukp::instance, sol::solution, k, best::solution)
  if(k == ukp.n)
    if (sol.z > best.z)
      best.z = sol.z;
      best.z = sol.z;
      best.r = sol.r;
      for i=1:ukp.n
        best.x[i] = sol.x[i];
      end
    end
    #test si on a la meilleur solution possible
    if(floor(computeLinearRelaxationUKP(ukp)) == best.z)
      return best;
    end

    #backtrack
    while (k>0) && (sol.x[k] == 0)
      k = k-1;
    end
    if (0 == k)
      #plus de backtrack possible
      return best;
    else
      sol.x[k] = 0;
      sol.z = sol.z - ukp.c[k];
      sol.r = sol.r + ukp.w[k];
      return branchandbound(ukp, sol, k, best);
    end
  else
    if(completeLinear(ukp, sol, k) > best.z)
      return branchandbound(ukp, completeGreedy(ukp, sol, k) , ukp.n, best);
    else
      #backtrack
      while (k>0) && (sol[k] == 0)
        k = k-1;
      end
      sol[k] = 0;
      sol.z = sol.z - ukp.c[k];
      sol.r = sol.r + ukp.w[k];
      return branchandbound(ukp, sol, k, best);
    end
  end
end
