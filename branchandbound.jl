function completeLinear(ukp, sol, k)
    # identify the last item integrally selected
    sommeW = ukp.W - sol.r
    zRelax = sol.z
    
    s = k+1;
    while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
      sommeW = sommeW + ukp.w[s]
      zRelax = zRelax + ukp.c[s]
      s = s + 1
    end

    # compute the upper bound
    if (s <= ukp.n) && ((ukp.W - sommeW) > 0)
      # contrainte non saturee => ajout de la partie fractionnaire de l'item bloquant
      zRelax = zRelax + floor((ukp.W - sommeW) * (ukp.c[s] / ukp.w[s]));
    end
    return zRelax;
end

# instance , solution en cour , indice de la variable considere , meilleur solution actuele
function branchandbound(ukp::instance, sol::solution, k, best::solution)
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
      if(completeLinear(ukp, sol, k-1) > best.z)
        best = branchandbound(ukp, sol, k+1, best);
      end
      sol.r = sol.r + ukp.w[k];
      sol.z = sol.z - ukp.c[k];
    end
    #on n'ajoute pas l'objet
    sol.x[k] = 0;
    if(completeLinear(ukp, sol, k-1) > best.z)
      best = branchandbound(ukp, sol, k+1, best);
    end
    return best;
  end
end

function completeMetT(ukp, sol, k)
  sommeW = ukp.W - sol.r; zRelax = sol.z;
  s = k+1;
  while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
    sommeW = sommeW + ukp.w[s]
    zRelax = zRelax + ukp.c[s]
    s = s + 1
  end

  reste = ukp.W - sommeW;
  if(s <= ukp.n)
      zRelax = zRelax + floor(reste / ukp.w[s])*ukp.c[s]
      reste = reste - floor(reste / ukp.w[s])*ukp.w[s]

    if(s < ukp.n)
      u0 = zRelax + floor( reste * (ukp.c[s+1] / ukp.w[s+1]) );
    else
      u0 = zRelax
    end
    if(s>1)
      u1 = zRelax + floor( ukp.c[s] - (ukp.w[s]-reste)*(ukp.c[s-1] / ukp.w[s-1]) );
    else
      u1 = zRelax
    end

    return max(u0,u1);
  else
    return zRelax
  end
end

function branchandbound2(ukp::instance, sol::solution, k, best::solution)
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
      if(completeMetT(ukp, sol, k-1) > best.z)
        best = branchandbound2(ukp, sol, k+1, best);
      end
      sol.r = sol.r + ukp.w[k];
      sol.z = sol.z - ukp.c[k];
    end
    #on n'ajoute pas l'objet
      sol.x[k] = 0;
      if(completeMetT(ukp, sol, k-1) > best.z)
        best = branchandbound2(ukp, sol, k+1, best);
      end
      return best;
  end
end
