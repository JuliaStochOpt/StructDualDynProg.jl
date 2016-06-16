            stats.setxtime += @mytime setchildx(parent, i, psol.x)
            stats.nsetx += 1
            stats.solvertime += @mytime childsol = loadAndSolve(child)
            stats.nsolved += 1
            childsolved[i] = true

            coef = childsol.πT
            rhs = childsol.πh + childsol.σd
            if childsol.status == :Infeasible
              infeasibility_detected = true
              feasible = false
              nnewfcuts += 1
              stats.fcutstime += @mytime pushfeasibilitycut!(child, coef, rhs, parent)
              break
            else
              rhs += childsol.ρe
              childocuts[i] = (coef, rhs)
            end
            if npaths[i] == :All || npaths[i] > 0
              push!(curpathss, (child, childsol, z + childsol.objvalx, prob * parent.proba[i], npaths[i]))
            end
          end
        end
        if feasible
          if parent.nlds.cutmode == :MultiCut
            for i in 1:length(parent.children)
              if childsolved[i]
                a, β = childocuts[i][1], childocuts[i][2]
                if mylt(β - dot(a, psol.x), .0, ztol)
                  error("The objectives are supposed to be nonnegative")
                end
                if mylt(psol.θ[i], β - dot(a, psol.x), ztol)
                  stats.ocutstime += @mytime pushoptimalitycutforparent!(parent.children[i], a, β, parent)
                  nnewocuts += 1
                end
              end
            end
          elseif parent.nlds.cutmode == :AveragedCut
            if !isempty(parent.children)
              if isnull(parent.childT)
                a = sum(map(i->childocuts[i][1]*parent.proba[i], 1:length(parent.children)))
              else
                a = sum(map(i->get(parent.childT)[i]'*childocuts[i][1]*parent.proba[i], 1:length(parent.children)))
              end
              β = sum(map(i->childocuts[i][2]*parent.proba[i], 1:length(parent.children)))
              if mylt(β - dot(a, psol.x), .0, ztol)
                error("The objectives are supposed to be nonnegative")
              end
              if mylt(psol.θ[1], β - dot(a, psol.x), ztol)
                stats.ocutstime += @mytime pushoptimalitycut!(parent, a, β, parent)
                nnewocuts += 1
              end
            end
