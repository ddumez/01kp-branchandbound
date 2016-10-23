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
    if (s <= ukp.n)
      # contrainte non saturee => ajout de la partie fractionnaire de l'item bloquant
      zRelax = zRelax + floor((ukp.W - sommeW) * (ukp.c[s] / ukp.w[s]));
    end
    return zRelax;
end

function completeMetT(ukp, sol, k)
  sommeW = ukp.W - sol.r;
  zRelax = sol.z;
  s = k+1;
  while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
    sommeW = sommeW + ukp.w[s]
    zRelax = zRelax + ukp.c[s]
    s = s + 1
  end

  reste = ukp.W - sommeW;
  if(s <= ukp.n)
      zRelax = zRelax + floor(reste / ukp.w[s]) * ukp.c[s]
      reste = reste - floor(reste / ukp.w[s])*ukp.w[s]

    if(s < ukp.n)
      u0 = zRelax + floor( reste * (ukp.c[s+1] / ukp.w[s+1]) );
    else
      u0 = 0
    end
    if(s>1)
      u1 = zRelax + floor( ukp.c[s] - (ukp.w[s]-reste)*(ukp.c[s-1] / ukp.w[s-1]) );
    else
      u1 = 0
    end

    return max(u0,u1);
  else
    return zRelax
  end
end
