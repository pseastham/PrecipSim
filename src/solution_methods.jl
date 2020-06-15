function computeStepByCrankNicholson(A,B,u,dt,f1,f2)
    LHS = B + 0.5*dt*A
    RHS = (B - 0.5*dt*A)*u + 0.5*dt*(f1+f2)

    return LHS\RHS
end

function computeStepByBackwardEuler(A,B,u,dt,f)
    LHS = B + dt*A
    RHS = B*u + dt*f

    return LHS\RHS
end

function computeStepByForwardEuler(A,B,u,dt,f)
    LHS = B 
    RHS = (B-dt*A)*u+ dt*f

    return LHS\RHS
end

function initializeChemical(sMesh,CHEMICAL::Symbol)
    ψinit = zeros(length(sMesh.xy))

    if CHEMICAL == :ψA
        ENTRANCE = :lower
    elseif CHEMICAL == :ψB
        ENTRANCE = :upper
    elseif CHEMICAL == :ψC
        ENTRANCE = :none
    end

    for i=1:length(sMesh.xy)
        ψinit[i] = scalarBCf_ydom(sMesh.xy[i].x,sMesh.xy[i].y,ENTRANCE)
    end

    return ψinit
end

function updateChemical(M,Stiff,u,p;solver=:backwarddeuler)
    A = Stiff.Op
    B = M
    uold = u
    dt = p.Δt
    f1 = Stiff.rhs
    f2 = copy(f1)   # change this for general case w/ reaction

    if solver==:forwardeuler
        return computeStepByForwardEuler(A,B,uold,dt,f1)
    elseif solver==:backwardeuler
        return computeStepByBackwardEuler(A,B,uold,dt,f1)
    elseif solver==:cranknicholson
        # NOTE: forcing will change with time in general, so this will need to be modified
        return computeStepByCrankNicholson(A,B,uold,dt,f1,f2)
    else
        error("solver $solver not recognized")
    end
end

function computeReaction!(p::precipParam,soln::precipSolution)
    N = length(soln.ψA.u)

    for i=1:N
        ψA = soln.ψA.u[i]
        ψB = soln.ψB.u[i]
        θs = soln.θs.u[i]
        a,b,c,d = QSSRstep(ψA,ψB,θs,p.Ma,p.Mb,p.ρf,p.ρs,p.a,p.b)
        soln.ψA.u[i] = a
        soln.ψB.u[i] = b
        soln.θs.u[i] = c
        soln.θf.u[i] = d
    end
    
    nothing
end

function QSSRstep(ψA::T,ψB::T,θs::T,Ma::T,Mb::T,ρf::T,ρs::T,a,b) where T<:Real
    θf = one(T) - θs
    J = (ρs-ρf)*θs + (Ma*ψA + Mb*ψB)*θf
    if isAlimitingReagent(ψA,ψB,a,b)
        ψAnew = zero(T)
        ϕB = (ψB - b*ψA/a)*θf
        θsnew = (J - Mb*ϕB)/(ρs-ρf) 
        ψBnew = ϕB/(1-θsnew)
    else
        ψBnew = zero(T)
        ϕA = (ψA - a*ψB/b)*θf
        θsnew = (J - Ma*ϕA)/(ρs-ρf) 
        ψAnew = ϕA/(1-θsnew)
    end
    θfnew = one(T)-θsnew
    return ψAnew,ψBnew,θsnew,θfnew
end

isAlimitingReagent(ψA,ψB,a,b) = (ψA/a <= ψB/b ? true : false)
isBlimitingReagent(ψA,ψB,a,b) = !(isAlimitingReagent(ψA,ψB,a,b))

function QSSRstep_test()
    @testset "QSSR reaction" begin
        a,b,Ma,Mb,Mc,ρs,ρf = defineChemistry()
        θfT = 0.738; θsT = 1.0 - θfT
        ψAT = 0.5*1e-2
        ψBT = 0.5*1e-2
        ψCT = 0.0
        ψA = 0.0
        ψB = 0.0
        θf = 0.0

        @testset "A limiting reagent" begin
            a=2; b=1
            initialMass = computeDensity(ψAT,ψBT,ψCT,θfT,ρf,ρs,Ma,Mb,Mc)
            ψA,ψB,θs,θf = QSSRstep(ψAT,ψBT,θsT,Ma,Mb,ρf,ρs,a,b)
            finalMass = computeDensity(ψA,ψB,ψCT,θf,ρf,ρs,Ma,Mb,Mc)
            @test isAlimitingReagent(ψAT,ψBT,a,b)
            @test abs(finalMass - initialMass)<3*eps()
            @test !(θfT < 0.0)
        end

        a,b,Ma,Mb,Mc,ρs,ρf = defineChemistry()
        θf = 1.0; θs = 1.0-θf
        ψA = 0.5*1e-2
        ψB = 0.5*1e-2
        ψC = 0.0

        @testset "B limiting reagent" begin
            a=1; b=2

            initialMass = computeDensity(ψA,ψB,ψC,θf,ρf,ρs,Ma,Mb,Mc)
            ψAnew,ψBnew,θsnew,θfnew = QSSRstep(ψA,ψB,θs,Ma,Mb,ρf,ρs,a,b)
            finalMass = computeDensity(ψA,ψB,ψC,θf,ρf,ρs,Ma,Mb,Mc)

            @test isBlimitingReagent(ψA,ψB,a,b)
            @test abs(finalMass - initialMass)<3*eps()
            @test !(θfnew < 0.0)
        end
    end

    nothing
end

# assumes that ψC==0
function computeTotalMass(p::precipParam,mesh::eFEM.ScalarMesh,soln::precipSolution)
    N = length(soln.ψA.u)
    density = zeros(N)
    # compute mass density function
    for i=1:N
        ψA = soln.ψA.u[i]
        ψB = soln.ψB.u[i]
        ψC = soln.ψC.u[i]
        θf = soln.θf.u[i]
        θs = soln.θs.u[i]
        density[i] = computeDensity(ψA,ψB,ψC,θf,p.ρf,p.ρs,p.Ma,p.Mb,p.Mc)
    end

    return eFEM.DomainQuad(mesh,density)
end

function computeDensity(ψA::T,ψB::T,ψC::T,θf::T,ρf::T,ρs::T,Ma::T,Mb::T,Mc::T) where T
    return ρf*θf + ρs*(1-θf) + Ma*ψA*θf + Mb*ψB*θf + Mc*ψC*θf
end

function defineChemistry()
    ρs = 4.1        # density (g/cm^3) of Ni(OH)2
    ρf = 1.0        # density (g/cm^3) of water
    Ma = 58.693     # molar mass (g/mol) of Ni2+
    Mb = 17.008     # molar mass (g/mol) of OH-
    Mc = 92.724     # molar mass (g/mol) of (Ni)OH_2
    a = 1
    b = 2
    c = 1

    return a,b,c,Ma,Mb,Mc,ρs,ρf
end

function checkValidChemistry(Ma,Mb,Mc,a,b,c)
    lhs = a*Ma + b*Mb
    rhs = c*Mc
    if abs(lhs - rhs) > 2*eps()
        return false
    else
        return true
    end
end

function computePureAqeuousReaction(ψAn::T,ψBn::T,ψCn::T,Ma::T,Mb::T,Mc::T,a::Int,b::Int,c::Int) where T<:Real
    ψnp1 = zero(typeof(ψAn))
    if isAlimitingReagent(ψAn,ψBn,a,b)
        ψBnp1 = ψBn - (b/a)*ψAn
        ψnp1 = (Ma*ψAn + Mb*ψBn + Mc*ψCn - Mb*ψBnp1)/Mc
    elseif isBlimitingReagent(ψAn,ψBn,a,b)
        ψAnp1 = ψAn - (a/b)*ψBn
        ψnp1 = (Ma*ψAn + Mb*ψBn + Mc*ψCn - Ma*ψAnp1)/Mc
    end
    return ψnp1
end

function checkψCExceedsThreshold(ψC,ψCthreshold)
    return (ψC > ψCthreshold) ? true : false
end