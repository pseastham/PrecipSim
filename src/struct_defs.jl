struct precipMesh
    sMesh::eFEM.ScalarMesh
    fMesh::eFEM.FluidMesh
end

struct precipSolution
    ψA::eFEM.ScalarSolution
    ψB::eFEM.ScalarSolution
    ψC::eFEM.ScalarSolution
    θs::eFEM.ScalarSolution
    θf::eFEM.ScalarSolution
    fluid::eFEM.FluidSolution
end

struct precipProblem
    ψA::eFEM.Problem
    ψB::eFEM.Problem
    ψC::eFEM.Problem
    fluid::eFEM.Problem
end

@with_kw struct precipParam
    # --------------------
    # Physical Parameters
    # --------------------
    μ0::Float64  = 1.0  # viscosity
    μ1::Float64  = 1.0  # brinkman viscosity
    α0::Float64  = 0.0  # initial value of thetam/thetaf
    G::Float64   = 1e0  # pressure inflow
    ξϵ::Float64  = 1e0  # regularization of friction coefficient
  
    PeA::Float64 = 30    # Peclet number for chemical A
    PeB::Float64 = 30    # Peclet number for chemical B
    PeC::Float64 = 1     # Peclet number for chemical C
  
    ρf::Float64 = 1.0     # mass density of fluid
    ρs::Float64 = 4.1     # mass density of solid
  
    # computed parameters
    α::Float64 = 1
    h::Float64 = 1
  
    # ---------------------
    # Chemistry Parameters aA + bB -> cC
    # ---------------------
    a::Int = 1
    b::Int = 1
    c::Int = 1
    Ma::Float64 = 58.693            # molar mass of A
    Mb::Float64 = 17.008            # molar mass of A
    Mc::Float64 = 0.0               # molar mass of C
    
    # ----------------------------
    # Time-integration parameters
    # ----------------------------
    Δt::Float64        = 5e-2
    Tfinal::Float64    = 80.0
    integrator::Symbol = :cranknicholson   # :cranknicholson, :backwardeuler
  
    # ---------------
    # vtk Parameters
    # ---------------
    ψAname::String = "aqueous_A"
    ψBname::String = "aqueous_B"
    ψCname::String = "aqueous_C"
    θsname::String = "membrane_vf"
    θfname::String = "fluid_vf"
    printSkipIndex::Int = 10
    printBool::Bool     = true
    SaveFolder::String  = string(homedir(),"/data/6-7-2020-y-dom")
    SaveName::String    = "y-dom"
  
    # ----------------
    # Mesh Parameters
    # ----------------
    MeshFolder::String="/home/peastham/Dropbox/projects/project-precipitation---numerics/PrecipSim.jl/meshes/y-dom"
    MeshFile::String = "y-dom6.msh"
    MeshTransformFunc::Function = none
    MeshTransformParam::Float64 = 1.0
    
    # -------------
    # run ID (1-7)
    # -------------
    runID::Int = 1
  end
  

