function run_simulation(p::PrecipParam)
    timer, mesh, soln, tsoln, prob = initialize_data(p)

    # initial brinkman
    mpbparam = BrinkmanMPParam(1.0,mpbFriction(soln.θf.u),soln.θf.u)
    sol = solve(prob.fluid,mesh.fMesh,mpbparam)
    tsoln.fluid.u = sol.u
    tsoln.fluid.v = sol.v
    tsoln.fluid.p = sol.p

    intializeTimeSteppingTimer!(timer)
    doPrinting!(timer,p,mesh,tsoln)
    while timer.simulationTime < timer.finalTime
        # reaction
        computeReaction!(p,tsoln)
        soln.θs.u = tsoln.θs.u
        soln.θf.u = tsoln.θf.u

        θfn = normalize_thetaf(tsoln.θf.u)
        peAnormalized = p.PeA./(θfn .+ 0.0 .*(1 .- θfn))        # 0.0 -> 0.5
        peBnormalized = p.PeB./(θfn .+ 0.0 .*(1 .- θfn))        # 0.0 -> 0.25
        ψAParam = AdvDiffParam(tsoln.fluid.u,tsoln.fluid.v,peAnormalized)
        ψBParam = AdvDiffParam(tsoln.fluid.u,tsoln.fluid.v,peBnormalized)
        ψCParam = AdvDiffParam(tsoln.fluid.u,tsoln.fluid.v,p.PeC./θfn)

        # advection-diffusion
        M = generateMassMatrixWScalar(mesh.sMesh,θfn)

        StiffA = GenerateSystem(mesh.sMesh,prob.ψA,ψAParam)
        ApplyBC!(StiffA,mesh.sMesh,prob.ψA,ψAParam,prob.ψA.OperatorType)

        StiffB = GenerateSystem(mesh.sMesh,prob.ψB,ψBParam)
        ApplyBC!(StiffB,mesh.sMesh,prob.ψB,ψBParam,prob.ψB.OperatorType)

        StiffC = GenerateSystem(mesh.sMesh,prob.ψC,ψCParam)
        ApplyBC!(StiffC,mesh.sMesh,prob.ψC,ψCParam,prob.ψC.OperatorType)

        if timer.stepIndex < 5
            soln.ψA.u = updateChemical(M,StiffA,tsoln.ψA.u,p;solver=:backwardeuler)
            soln.ψB.u = updateChemical(M,StiffB,tsoln.ψB.u,p;solver=:backwardeuler)
            soln.ψC.u = updateChemical(M,StiffC,tsoln.ψC.u,p;solver=:backwardeuler)
        else
            soln.ψA.u = updateChemical(M,StiffA,tsoln.ψA.u,p;solver=p.integrator)
            soln.ψB.u = updateChemical(M,StiffB,tsoln.ψB.u,p;solver=p.integrator)
            soln.ψC.u = updateChemical(M,StiffC,tsoln.ψC.u,p;solver=p.integrator)
        end

        scrub_scalars!(soln.ψA.u,soln.ψB.u,soln.ψC.u)

        mpbparam = BrinkmanMPParam(1.0,mpbFriction(θfn),θfn)
        sol = solve(prob.fluid,mesh.fMesh,mpbparam)
        soln.fluid.u = sol.u
        soln.fluid.v = sol.v
        soln.fluid.p = sol.p

        tsoln = deepcopy(soln)
        if isPrintIndex(timer)
            doPrinting!(timer,p,mesh,soln)
        end
        updateSimulationTime!(timer)
        updateStepIndex!(timer)
    end
    printSimulationFinishedMessage(timer)
    nothing
end

function scrub_scalars!(u1,u2,u3)
    scrub_scalar!(u1)
    scrub_scalar!(u2)
    scrub_scalar!(u3)
    nothing
end
function scrub_scalar!(u)
    for i=1:length(u); if u[i] < 0; u[i] = 0.0; end; end
    nothing
end

function initialize_data(p::precipParam)
    timer = initializeTimer(p)
    mesh = initialize_meshes(p)
    soln = initializeSolution(p)
    tsoln = initializeSolution(p)
    prob = defineProblems(p,mesh)

    return timer, mesh, soln, tsoln, prob
end

function initialize_meshes(p::precipParam)
    sMesh,fMesh = loadMesh(p)

    return precipMesh(sMesh,fMesh)
    if p.MeshFile == "sin_x2.msh"
        sin_x2_transform!(mesh,1.0)
    elseif p.MeshFile == "sin_wall.msh"
        sin_wall_transform!(mesh,1.0)
    end

    return precipMesh(sMesh,fMesh)
end