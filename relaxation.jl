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
            u0 = floor( reste * (ukp.c[s+1] / ukp.w[s+1]) );
        else
            u0 = 0
        end
        if(s>1)
            u1 = floor( ukp.c[s] - (ukp.w[s]-reste)*(ukp.c[s-1] / ukp.w[s-1]) );
        else
            u1 = 0
        end

        return zRelax + max(u0,u1);
    else
        return zRelax
    end
end

function completeMetT2(ukp, sol, k)
    sommeW = ukp.W - sol.r;
    zRelax = sol.z;
    s = k+1;
    while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
        sommeW = sommeW + ukp.w[s]
        zRelax = zRelax + ukp.c[s]
        s = s + 1
    end
    reste = ukp.W - sommeW;

    if (s <= ukp.n)
        sigma0s = sigma0(ukp, s)
        sigma1s = sigma1(ukp, s)
        if (sigma0s >= 1) && (sigma0s <= ukp.n)
            u0 = 0
            reste = ukp.W;
            for i = 1:(sigma0s-1)
                if (i != sigma0s)
                    u0 = u0 + ukp.c[i];
                    reste = reste - ukp.w[i];
                end
            end
            u0 = u0 + floor((reste * ukp.c[sigma0s])/ukp.w[sigma0s]);
        else
            u0 = 0;
        end

        if (sigma0s >= 1) && (sigma0s <= ukp.n)
            u1 = ukp.c[s]
            reste = ukp.W - ukp.w[s]
            for i=1:(sigma1s-1)
                u1 = u1 + ukp.c[i]
                reste = reste - ukp.w[i]
            end
            u1 = u1 + floor( (reste * ukp.c[sigma1s])/ukp.w[sigma1s])
        else
            u1 = 0
        end
        return max(u0,u1)
    else
        return zRelax;
    end
end

function sigma1(ukp, j)
    sommeW = ukp.w[j];
    i = 1;
    while (i<= ukp.n) && (sommeW + ukp.w[i] <= ukp.W)
        sommeW = sommeW + ukp.w[i]
        i = i+1
    end
    return i;
end

function sigma0(ukp, j)
    sommeW = 0
    i = 1;
    while (i<=ukp.n) && (sommeW + ukp.w[i] <= ukp.W)
        if (j != i)
            sommeW = sommeW + ukp.w[i]
        end
        i = i+1
    end
    return i;
end
