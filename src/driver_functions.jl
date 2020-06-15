function runSimulation(p::precipParam)
  doInitialPrinting(p)
  if p.runID == 1
    mainDriver(p::precipParam)
  elseif p.runID == 2
    pureDiffusionOneChemicalTestDriver(p)
  elseif p.runID == 3
    pureAdvectionOneChemicalTestDriver(p)
  elseif p.runID == 4
    AdvectionDiffusionOneChemicalTestDriver(p)
  elseif p.runID == 5
    MultiphaseBrinkmanTestDriver(p)
  elseif p.runID == 6
    DiffusionAndReactionTwoChemicalTestDriver(p)
  elseif p.runID == 7
    pureReactionTestDriver(p)
  else
    println("runID of ",p.runID," not recognized.")
    println("runID 1-7 required")
  end

  nothing
end

function mainDriver(p::precipParam)
  timer, mesh, soln, tsoln, prob = initializeData(p)

  # initial brinkman
  mpbparam = BrinkmanMPParam(1.0,mpbFriction(soln.θf.u),soln.θf.u)
  sol = solve(prob.fluid,mesh.fMesh,mpbparam)
  tsoln.fluid.u = sol.u
  tsoln.fluid.v = sol.v
  tsoln.fluid.p = sol.p

  intializeTimeSteppingTimer!(timer)
  doPrinting!(timer,p,mesh,tsoln)
  while timer.simulationTime < timer.finalTime
    #Ninterp = 80
    #interpVals = zeros(Ninterp,6)
    #surf = mySurface([0.5,0.37],[0.5,-0.37])

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

    # scrub scalars so they don't go below zero.
    for i=1:length(soln.ψA.u)
      if soln.ψA.u[i] < 0
        soln.ψA.u[i] = 0.0
      end
      if soln.ψB.u[i] < 0
        soln.ψB.u[i] = 0.0
      end
      if soln.ψC.u[i] < 0
        soln.ψC.u[i] = 0.0
      end
    end

    # brinkman
    # normalize thetaf
    mpbparam = BrinkmanMPParam(1.0,mpbFriction(θfn),θfn)
    sol = solve(prob.fluid,mesh.fMesh,mpbparam)
    soln.fluid.u = sol.u
    soln.fluid.v = sol.v
    soln.fluid.p = sol.p

    # inteprolate values
    #xvals = zeros(length(mesh.sMesh.xy))
    #yvals = zeros(length(mesh.sMesh.xy))
    #for i=1:length(mesh.sMesh.xy)
    #  xvals[i] = mesh.sMesh.xy[i].x
    #  yvals[i] = mesh.sMesh.xy[i].y
    #end
    #vmag = sqrt.(soln.fluid.u.^2 + soln.fluid.v.^2)
    #interpVals[:,1] = SurfaceInterp(mesh.sMesh,xvals,surf,Ninterp-1)
    #interpVals[:,2] = SurfaceInterp(mesh.sMesh,yvals,surf,Ninterp-1)
    #interpVals[:,3] = SurfaceInterp(mesh.sMesh,soln.ψA.u,surf,Ninterp-1)
    #interpVals[:,4] = SurfaceInterp(mesh.sMesh,soln.ψB.u,surf,Ninterp-1)
    #interpVals[:,5] = SurfaceInterp(mesh.sMesh,soln.θs.u,surf,Ninterp-1)
    #interpVals[:,6] = SurfaceInterp(mesh.fMesh,vmag,surf,Ninterp-1)

    tsoln = deepcopy(soln)
    if isPrintIndex(timer)
      # save interpolated values
      #savefile = string(p.SaveFolder,"/interpVals_",timer.printIndex,".jld2")
      #save(savefile, "time", timer.simulationTime, "interpVals", interpVals)

      # ordinary printing
      doPrinting!(timer,p,mesh,soln)
    end
    updateSimulationTime!(timer)
    updateStepIndex!(timer)
  end

  printSimulationFinishedMessage(timer)

  nothing
end

function normalize_thetaf(θf::Vector{T}) where T<:Real
  N=length(θf)
  θf_normalized = zeros(T,N)
  TOL = 0.1

  for i=1:N
    if θf[i] < TOL
      θf_normalized[i] = TOL
    else
      θf_normalized[i] = θf[i]
    end
  end

  return θf_normalized
end

function pureDiffusionOneChemicalTestDriver(p::precipParam)
    timer, mesh, soln, tsoln, prob = initializeData(p)

    ψAParam = PoissonParam(p.κA)

    for i=1:length(soln.θf.u)
      soln.θf.u[i] = 0.1
    end

    M = generateMassMatrixWScalar(mesh.sMesh,soln.θf.u)
    StiffA = GenerateSystem(mesh.sMesh,prob.ψA,ψAParam)
    ApplyBC!(StiffA,mesh.sMesh,prob.ψA,ψAParam,prob.ψA.OperatorType)

    intializeTimeSteppingTimer!(timer)
    doPrinting!(timer,p,mesh,soln)
    while timer.simulationTime < timer.finalTime
      if timer.stepIndex < 5
        soln.ψA.u = updateChemical(M,StiffA,tsoln.ψA.u,p;solver=:backwardeuler)
      else
        soln.ψA.u = updateChemical(M,StiffA,tsoln.ψA.u,p;solver=p.integrator)
      end

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

function pureAdvectionOneChemicalTestDriver(p::precipParam)
  println("pure advection test driver isn't programmed yet")
  nothing
end

function AdvectionDiffusionOneChemicalTestDriver(p::precipParam)
  println("advection-diffusion test driver isn't programmed yet")
  nothing
end

function MultiphaseBrinkmanTestDriver(p::precipParam)
  function computeNorms(N)
    start = time()
    # load mesh
    mesh = squareMeshFluid([-2,2,-1,1],N)

    α1(x,y) = 1.0
    α2(x,y) = 0.5+0.4*cos(3*x*y)
    α3(x,y) = 0.5+0.4*sin(7*x*y)

    xm = [i.x for i in mesh.xy]
    ym = [i.y for i in mesh.xy]
    α1arr = [α1(xm[i],ym[i]) for i=1:length(mesh.xy)]
    α2arr = [α2(xm[i],ym[i]) for i=1:length(mesh.xy)]
    α3arr = [α3(xm[i],ym[i]) for i=1:length(mesh.xy)]

    param  = BrinkmanMPParam(α1arr,α2arr,α3arr)

    OperatorType = :BrinkmanMP2D

    dNodes = Dirichlet(:left,:top,:bottom,:right)
    Nodes = [dNodes]

    u(x,y) = -x^4*y^2
    v(x,y) = 4x^3*y^3/3

    dUBCf = Dirichlet( (x,y) -> u(x,y) )
    dVBCf = Dirichlet( (x,y) -> v(x,y) )

    d2dudx(x,y) = -12*x^2*y^2; d2dudy(x,y) = -2*x^4
    d2dvdx(x,y) = 8*x*y^3;     d2dvdy(x,y) = 8*x^3*y
    lapu(x,y) = d2dudx(x,y) + d2dudy(x,y)
    lapv(x,y) = d2dvdx(x,y) + d2dvdy(x,y)

    dpdx(x,y) = 3*x^2*y^3; dpdy(x,y) = 3*x^3*y^2

    Fx(x,y) = -α1(x,y)*lapu(x,y) + α2(x,y)*u(x,y) + α3(x,y)*dpdx(x,y)
    Fy(x,y) = -α1(x,y)*lapv(x,y) + α2(x,y)*v(x,y) + α3(x,y)*dpdy(x,y)

    ffx   = Forcing(Fx)
    ffy   = Forcing(Fy)
    bcfun = [dUBCf,dVBCf,ffx,ffy]

    prob = Problem(mesh,Nodes,bcfun,OperatorType)
    sol = solve(prob,mesh,param)

    elapsed = time()-start

    # compute condition number of operator matrix
    #LinOp = GenerateSystem(mesh,prob,param)
    #ApplyBC!(LinOp,mesh,prob,param,OperatorType)
    #κ = cond(Array(LinOp.Op),2)
    κ = 0.5

    # compute mesh h
    h = hCalc(mesh)

    # generate exact solution
    Uexact(x,y) = dUBCf.f(x,y)
    Vexact(x,y) = dVBCf.f(x,y)
    Pexact(x,y) = x^3*y^3
    UexactArr = [Uexact.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    VexactArr = [Vexact.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    PexactArr = [Pexact.(mesh.xyp[i].x,mesh.xyp[i].y) for i=1:length(mesh.xyp)]

    # compute difference
    velErr = sqrt.((sol.u - UexactArr).^2 + (sol.v -VexactArr).^2)
    presErr = sol.p - PexactArr

    # compute Li norm
    L1v   = DomainNorm(mesh.xy,mesh.cm,velErr;normID="1")
    L2v   = DomainNorm(mesh.xy,mesh.cm,velErr;normID="2")
    Linfv = DomainNorm(mesh.xy,mesh.cm,velErr;normID="Inf")
    L1p   = DomainNorm(mesh.xyp,mesh.cmp,presErr;normID="1")
    L2p   = DomainNorm(mesh.xyp,mesh.cmp,presErr;normID="2")
    Linfp = DomainNorm(mesh.xyp,mesh.cmp,presErr;normID="Inf")

    # plot solution
    sD = ScalarData(sol.p,α1arr,α2arr,α3arr)
    sN = ScalarNames("pressure","alpha_1","alpha_2","alpha_3")
    vD = VectorData([sol.u,sol.v])
    vN = VectorNames("velocity")
    fn = Path("solution_var")
    vtksave(mesh,sD,sN,vD,vN,fn)

    # plot solution
    sD = ScalarData(PexactArr)
    sN = ScalarNames("pressure")
    vD = VectorData([UexactArr,VexactArr])
    vN = VectorNames("velocity")
    fn = Path("exact_var")
    vtksave(mesh,sD,sN,vD,vN,fn)

    return κ,h,L1v,L2v,Linfv,L1p,L2p,Linfp,elapsed
  end

  Narr     = [4,4,8,16,32,64]#,96]
  N        = length(Narr)
  harr     = zeros(N)
  L1arrV   = zeros(N)
  L2arrV   = zeros(N)
  LInfarrV = zeros(N)
  L1arrP   = zeros(N)
  L2arrP   = zeros(N)
  LInfarrP = zeros(N)
  timearr  = zeros(N)
  κarr     = zeros(N)

  for i=1:N
    n = Narr[i]
    κarr[i],harr[i],L1arrV[i],L2arrV[i],LInfarrV[i],
      L1arrP[i],L2arrP[i],LInfarrP[i],timearr[i] =
      computeNorms(n)
      println("completed N=$(n)")
  end

  κarr     = κarr[2:end]
  harr     = harr[2:end]
  L1arrV   = L1arrV[2:end]
  L2arrV   = L2arrV[2:end]
  LInfarrV = LInfarrV[2:end]
  L1arrP   = L1arrP[2:end]
  L2arrP   = L2arrP[2:end]
  LInfarrP = LInfarrP[2:end]
  timearr  = timearr[2:end]

  # export output to *.jld file
  run(`mkdir -p data`)

  save("data/TEMP_mpb2D_validation.jld", "harr",harr,
                                          "VL1arr", L1arrV,
                                          "VL2arr", L2arrV,
                                          "VLInfarr",LInfarrV,
                                          "PL1arr", L1arrP,
                                          "PL2arr", L2arrP,
                                          "PLInfarr",LInfarrP,
                                          "timearr",timearr,
                                          "κarr",κarr)


  nothing
end

function DiffusionAndReactionTwoChemicalTestDriver(p::precipParam)
  timer, mesh, soln, tsoln, prob = initializeData(p)

  ψAParam = PoissonParam(p.κA)
  ψBParam = PoissonParam(p.κB)
  ψCParam = PoissonParam(p.κC)

  M = generateMassMatrix(mesh.sMesh)

  StiffA = GenerateSystem(mesh.sMesh,prob.ψA,ψAParam)
  ApplyBC!(StiffA,mesh.sMesh,prob.ψA,ψAParam,prob.ψA.OperatorType)

  StiffB = GenerateSystem(mesh.sMesh,prob.ψB,ψBParam)
  ApplyBC!(StiffB,mesh.sMesh,prob.ψB,ψBParam,prob.ψB.OperatorType)

  StiffC = GenerateSystem(mesh.sMesh,prob.ψC,ψCParam)
  ApplyBC!(StiffC,mesh.sMesh,prob.ψC,ψCParam,prob.ψC.OperatorType)

  intializeTimeSteppingTimer!(timer)
  doPrinting!(timer,p,mesh,soln)
  while timer.simulationTime < timer.finalTime
    println("i1=",timer.stepIndex," ",computeTotalMass(p,mesh.sMesh,soln))
    if timer.stepIndex < 5
      soln.ψA.u = updateChemical(M,StiffA,tsoln.ψA.u,p;solver=:backwardeuler)
      soln.ψB.u = updateChemical(M,StiffB,tsoln.ψB.u,p;solver=:backwardeuler)
      soln.ψC.u = updateChemical(M,StiffC,tsoln.ψC.u,p;solver=:backwardeuler)
    else
      soln.ψA.u = updateChemical(M,StiffA,tsoln.ψA.u,p;solver=p.integrator)
      soln.ψB.u = updateChemical(M,StiffB,tsoln.ψB.u,p;solver=p.integrator)
      soln.ψC.u = updateChemical(M,StiffC,tsoln.ψC.u,p;solver=p.integrator)
    end
    println("i2=",timer.stepIndex," ",computeTotalMass(p,mesh.sMesh,soln))
    computeReaction!(p,soln)
    println("i3=",timer.stepIndex," ",computeTotalMass(p,mesh.sMesh,soln))
    println()

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

function pureReactionTestDriver(p::precipParam)
  timer, mesh, soln, tsoln, prob = initializeData(p)

  randomizeScalarSolutions!(soln)

  intializeTimeSteppingTimer!(timer)
  doPrinting!(timer,p,mesh,soln)
  while timer.simulationTime < timer.finalTime
    computeReaction!(p,soln)

    if isPrintIndex(timer)
      doPrinting!(timer,p,mesh,soln)
    end

    tsoln = deepcopy(soln)

    # step time
    updateSimulationTime!(timer)
    updateStepIndex!(timer)
  end

  printSimulationFinishedMessage(timer)

  nothing
end

function initializeData(p::precipParam)
  timer = initializeTimer(p)
  mesh = initializeMeshes(p)
  soln = initializeSolution(p)
  tsoln = initializeSolution(p)
  prob = defineProblems(p,mesh)

  return timer, mesh, soln, tsoln, prob
end

function mpbFriction(θf::Vector{T}) where T<:Real
  N = length(θf)
  friction = zeros(N)

  K=0.5
  n=2
  ξstar = 30.0
  θstar = 0.5
  h = 1000.0#computeH(ξstar,θstar,K,n)

  for i=1:N
    friction[i] = h*frictionKernel(θf[i],K,n)
  end

  return friction
end

function frictionKernel(θf,K,n)
  θs = (1.0-θf)
  return θs^n/(K^n + θs^n)
end

function computeH(ξstar,θstar,K,n)
  return ξstar/frictionKernel(θstar,K,n)
end
