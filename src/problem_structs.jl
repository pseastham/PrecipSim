import LinearAlgebra.cross

struct precipMesh
  sMesh::ScalarMesh
  fMesh::FluidMesh
end

struct precipSolution
  ψA::ScalarSolution
  ψB::ScalarSolution
  ψC::ScalarSolution
  θs::ScalarSolution
  θf::ScalarSolution
  fluid::FluidSolution      # velocity and pressure solution
end

struct precipProblem
  ψA::Problem
  ψB::Problem
  ψC::Problem
  fluid::Problem
end

function initializeSolution(p::precipParam)
  sMesh,fMesh = loadMesh(p)
  N = length(sMesh.xy); Np = length(fMesh.xyp)
  zs = zeros(N); zq = zeros(N); zp = zeros(Np)
  ψA = ScalarSolution(copy(zs))
  ψB = ScalarSolution(copy(zs))
  ψC = ScalarSolution(copy(zs))
  θs = ScalarSolution(copy(zs))
  θf = ScalarSolution(ones(N))
  fluidqp = FluidSolution(zq,copy(zq),zp)

  if p.runID == 6
    for i=1:N
      x = sMesh.xy[i].x; y = sMesh.xy[i].y
      ψA.u[i] = 5e-4*diffusionAndReactionIC(x,y,:ψA)
      ψB.u[i] = 5e-4*diffusionAndReactionIC(x,y,:ψB)
    end
  end

  return precipSolution(ψA,ψB,ψC,θs,θf,fluidqp)
end

function diffusionAndReactionIC(x,y,symbol)
  if symbol == :ψA
    #if inLeftArea(x,y)
    #  return 1.0
    #else
    #  return 0.0
    #end
    return abs(x)
  elseif symbol == :ψB
    #if inRightArea(x,y)
    #  return 1.0
    #else
    #  return 0.0
    #end
    return abs(y)
  end
  error("symbol not recognized")
end
function inLeftArea(x,y)
  if x < -0.5 && x > -1.5 && y < 0.5 && y > -0.5
    return true
  else
    return false
  end
end
function inRightArea(x,y)
  if x < 1.5 && x > 0.5 && y < 0.5 && y > -0.5
    return true
  else
    return false
  end
end

function randomizeScalarSolutions!(soln::precipSolution)
  N = length(soln.ψA.u)
  soln.ψA.u = 0.5*1e-2*rand(N)
  soln.ψB.u = 0.5*1e-2*rand(N)
  soln.ψC.u = zeros(N)
  soln.θs.u = zeros(N)
  soln.θf.u = 1 .- soln.θs.u
  nothing
end

function defineGMSHChemicalA(mesh::ScalarMesh)
  # define boundary conditions
  dNodes = Dirichlet(:inflow)
  nNodes = Neumann(:outflow,:wall)
  dBCf = Dirichlet((x,y) -> scalarBCf_ydom(x,y,:lower))
  nBCf = Neumann((x,y) -> 0.0)
  Nodes = [dNodes,nNodes]; bcfun = [dBCf,nBCf]

  # define problems
  chemAProb = Problem(mesh,Nodes,bcfun,:AdvDiff2D)
end

function defineGMSHChemicalB(mesh::ScalarMesh)
  # define boundary conditions
  dNodes = Dirichlet(:inflow)
  nNodes = Neumann(:outflow,:wall)
  dBCf = Dirichlet((x,y) -> scalarBCf_ydom(x,y,:upper))
  nBCf = Neumann((x,y) -> 0.0)
  Nodes = [dNodes,nNodes]; bcfun = [dBCf,nBCf]

  # define problems
  chemBProb = Problem(mesh,Nodes,bcfun,:AdvDiff2D)
end

function defineGMSHChemicalC(mesh::ScalarMesh)
  # define boundary conditions
  nNodes = Neumann(:inflow,:outflow,:wall)
  nBCf = Neumann((x,y) -> 0.0)
  Nodes = [nNodes]; bcfun = [nBCf]

  # define problems
  chemCProb = Problem(mesh,Nodes,bcfun,:AdvDiff2D)
end

function defineGMSHMPBrinkman(mesh::FluidMesh)
  dNodes = Dirichlet(:inflow,:wall)
  nNodes = Neumann(:outflow)
  dUBCf = Dirichlet(uBCf_ydom)
  dVBCf = Dirichlet(vBCf_ydom)
  nBCf = Neumann((x,y) -> 0.0)
  Nodes = [dNodes,nNodes]; bcfun = [dUBCf,dVBCf,nBCf]

  # define problems
  fluidProb = Problem(mesh,Nodes,bcfun,:BrinkmanMP2D)
end

function defineSquareChemical(mesh::ScalarMesh)
  # define boundary conditions
  nNodes = Neumann(:left,:right,:top,:bottom)
  nBCf = Neumann((x,y) -> 0.0)
  Nodes = [nNodes]; bcfun = [nBCf]

  # define problems
  chemAProb = Problem(mesh,Nodes,bcfun,:Poisson2D)
end

function defineSquareMPBrinkman(mesh::FluidMesh)
  dNodes = Dirichlet(:left,:right,:bottom,:top)
  dUBCf = Dirichlet((x,y) -> 0.0)
  dVBCf = Dirichlet((x,y) -> 0.0)
  Nodes = [dNodes]; bcfun = [dUBCf,dVBCf]

  # define problems
  fluidProb = Problem(mesh,Nodes,bcfun,:BrinkmanMP2D)
end

function defineProblems(p::precipParam,mesh::precipMesh)
  if p.runID == 1
    ψAprob = defineGMSHChemicalA(mesh.sMesh)
    ψBprob = defineGMSHChemicalB(mesh.sMesh)
    ψCprob = defineGMSHChemicalC(mesh.sMesh)
    fluidprob = defineGMSHMPBrinkman(mesh.fMesh)
  elseif p.runID == 2
    ψAprob = defineGMSHChemicalA(mesh.sMesh)
    ψBprob = defineGMSHChemicalB(mesh.sMesh)
    ψCprob = defineGMSHChemicalC(mesh.sMesh)
    fluidprob = defineGMSHMPBrinkman(mesh.fMesh)
  elseif p.runID == 3
    error("undefine problem. see problem_structs.jl")
  elseif p.runID == 4
    error("undefine problem. see problem_structs.jl")
  elseif p.runID == 5
    error("undefine problem. see problem_structs.jl")
  elseif p.runID == 6
    ψAprob = defineSquareChemical(mesh.sMesh)
    ψBprob = defineSquareChemical(mesh.sMesh)
    ψCprob = defineSquareChemical(mesh.sMesh)
    fluidprob = defineSquareMPBrinkman(mesh.fMesh)
  else
    error("you shouldn't get here")
  end

  return precipProblem(ψAprob,ψBprob,ψCprob,fluidprob)
end

function uBCf_ydom(x,y)
  # define relavant parameters for mesh
  θ = 0.7227342478134156
  α = pi/2 - θ
  sα = sin(α)
  sθ = sin(θ)
  cα = cos(α)
  cθ = cos(θ)
  IL = 1.0
  IW = 0.5
  OW = 1.5*IW
  mp = sqrt(IW^2-(0.5*OW)^2)

  # define (x,y) nodes defining inflow boundary
  x1  = -mp-IL*sα;  x2  = -IL*sα
  y1L = -IL*cα;  y2L = -IL*cα - OW/2
  y1U =  IL*cα;  y2U =  IL*cα + OW/2

  # determine if (x,y) sits on interface
  pTest = [x,y,0.0]
  p1L   = [x1,y1L,0.0]; p2L = [x2,y2L,0.0]
  p1U   = [x1,y1U,0.0]; p2U = [x2,y2U,0.0]

  v1L = pTest - p1L;  v2L = p2L - p1L
  v1U = pTest - p1U;  v2U = p2U - p1U

  cL = cross(v1L,v2L)[3]; cU = cross(v1U,v2U)[3]

  sectionName = :temp

  if abs(cL)<1e-8
    sectionName = :lower
  elseif abs(cU)<1e-8
    sectionName = :upper
  else
    sectionName = :wall
  end

  # if on inflow transformation = (x,y) -> (c*(-W+1)-SL, c*(-W+1))
  if sectionName == :lower
    xstart = x1; xend   = x2
    s      = linearInterp(x,xstart,xend,0,1)

    return 4*s*(1-s)*sα

  elseif sectionName == :upper  # top inflow
    xstart = x1; xend   = x2
    s      = linearInterp(x,xstart,xend,0,1)

    return 4*s*(1-s)*sα

  elseif sectionName == :wall
    return 0.0
  end
end

function vBCf_ydom(x,y)
  # define relavant parameters for mesh
  θ = 0.7227342478134156
  α = pi/2 - θ
  sα = sin(α)
  sθ = sin(θ)
  cα = cos(α)
  cθ = cos(θ)
  IL = 1.0
  IW = 0.5
  OW = 1.5*IW
  mp = sqrt(IW^2-(0.5*OW)^2)

  # define (x,y) nodes defining inflow boundary
  x1  = -mp-IL*sα;  x2  = -IL*sα
  y1L = -IL*cα;  y2L = -IL*cα - OW/2
  y1U =  IL*cα;  y2U =  IL*cα + OW/2

  # determine if (x,y) sits on interface
  pTest = [x,y,0.0]
  p1L   = [x1,y1L,0.0]; p2L = [x2,y2L,0.0]
  p1U   = [x1,y1U,0.0]; p2U = [x2,y2U,0.0]

  v1L = pTest - p1L;  v2L = p2L - p1L
  v1U = pTest - p1U;  v2U = p2U - p1U

  cL = cross(v1L,v2L)[3]; cU = cross(v1U,v2U)[3]

  sectionName = :temp

  if abs(cL)<1e-8
    sectionName = :lower
  elseif abs(cU)<1e-8
    sectionName = :upper
  else
    sectionName = :wall
  end

  # if on inflow transformation = (x,y) -> (c*(-W+1)-SL, c*(-W+1))
  if sectionName == :lower          # bottom inflow
    xstart = x1; xend = x2
    s      = linearInterp(x,xstart,xend,0,1)

    return 4*s*(1-s)*cα
  elseif sectionName == :upper      # top inflow
    xstart = x1; xend = x2
    s      = linearInterp(x,xstart,xend,0,1)

    return -4*s*(1-s)*cα
  elseif sectionName == :wall
    return 0.0
  end
end

function scalarBCf_ydom(x,y,ENTRANCE::Symbol)
  # define relavant parameters for mesh
  θ = 0.7227342478134156
  α = pi/2 - θ
  sα = sin(α)
  sθ = sin(θ)
  cα = cos(α)
  cθ = cos(θ)
  IL = 1.0
  IW = 0.5
  OW = 1.5*IW
  mp = sqrt(IW^2-(0.5*OW)^2)

  # define (x,y) nodes defining inflow boundary
  x1  = -mp-IL*sα;  x2  = -IL*sα
  y1L = -IL*cα;  y2L = -IL*cα - OW/2
  y1U =  IL*cα;  y2U =  IL*cα + OW/2

  # determine if (x,y) sits on interface
  pTest = [x,y,0.0]
  p1L   = [x1,y1L,0.0]; p2L = [x2,y2L,0.0]
  p1U   = [x1,y1U,0.0]; p2U = [x2,y2U,0.0]

  v1L = pTest - p1L;  v2L = p2L - p1L
  v1U = pTest - p1U;  v2U = p2U - p1U

  cL = cross(v1L,v2L)[3]; cU = cross(v1U,v2U)[3]

  sectionName = :temp

  if abs(cL)<1e-8
    sectionName = :lower
  elseif abs(cU)<1e-8
    sectionName = :upper
  else
    sectionName = :wall
  end

  # if on inflow transformation = (x,y) -> (c*(-W+1)-SL, c*(-W+1))
  if sectionName == ENTRANCE
    return 0.005
  else
    return 0.0
  end
end

"""
  linearInterp()

linear interpolates the point x0 which lies in (x1,x2) onto the domain (s1,s2)
"""
function linearInterp(x0,x1,x2,s1,s2)
  A = [x1 1; x2 1]; rhs=[s1,s2]
  mbtemp = A\rhs
  m = mbtemp[1]
  b = mbtemp[2]

  s0 = m*x0 + b

  return s0
end