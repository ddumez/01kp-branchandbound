include("relaxation.jl")

# instance , solution en cour , indice de la variable considere , meilleur solution actuele
function branchandbound(ukp::instance, sol::solution, k, best::solution, borne)
  if(k == ukp.n+1)
    #la solution est complete
    if (sol.z > best.z)
      best.z = sol.z;
      best.r = sol.r;
      for i=1:ukp.n
        best.x[i] = sol.x[i];
      end
    end
    return best;
  else
    #on essaye d'ajoute un objet de maniere greedy
    if(sol.r >= ukp.w[k])
      #on ajoute l'objet
      sol.x[k] = 1;
      sol.r = sol.r - ukp.w[k];
      sol.z = sol.z + ukp.c[k];
      if(borne(ukp, sol, k-1) > best.z)
        best = branchandbound(ukp, sol, k+1, best, borne);
      end
      sol.r = sol.r + ukp.w[k];
      sol.z = sol.z - ukp.c[k];
    end
    #on n'ajoute pas l'objet
    sol.x[k] = 0;
    if(borne(ukp, sol, k) > best.z)
      best = branchandbound(ukp, sol, k+1, best, borne);
    end
    return best;
  end
end
