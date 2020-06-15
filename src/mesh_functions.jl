# add example for if run case requires a square mesh (e.g. for testing)

function initializeMeshes(p::precipParam)
    sMesh,fMesh = loadMesh(p)

    mesh = precipMesh(sMesh,fMesh)
    if p.MeshFile == "sin_x2.msh"
        sin_x2_transform!(mesh,1.0)
    elseif p.MeshFile == "sin_wall.msh"
        sin_wall_transform!(mesh,1.0)
    end
    
    return precipMesh(sMesh,fMesh)
  end

function loadMesh(p::precipParam)
    fileDir = string(p.MeshFolder,"/",p.MeshFile)
    N = 50; square = [-2.0,2.0,-1.0,1.0]; order=2

    if p.runID == 1
        sMesh = Mesh(fileDir)
        fMesh = FluidMesh(fileDir)
    elseif p.runID == 2
        sMesh = Mesh(fileDir)
        fMesh = FluidMesh(fileDir)
    elseif p.runID == 3
        error("no mesh defined. see mesh_functions.jl")
    elseif p.runID == 4
        error("no mesh defined. see mesh_functions.jl")    
    elseif p.runID == 5
        error("no mesh defined. see mesh_functions.jl")   
    elseif p.runID == 6
        sMesh = squareMesh(square,N,order)
        fMesh = squareMeshFluid(square,N)
    elseif p.runID == 7
        sMesh = Mesh(fileDir)
        fMesh = FluidMesh(fileDir)
    else
        error("you shouldn't get this error")
    end

    return sMesh,fMesh
end

function sin_x2_transform!(mesh::precipMesh,ω::T) where T<:Real
    Ns = length(mesh.sMesh.xy)
    Nfu = length(mesh.fMesh.xy)
    Nfp = length(mesh.fMesh.xyp)
    amp = 0.1

    for ti=1:Ns
        xval = mesh.sMesh.xy[ti].x
        yval = mesh.sMesh.xy[ti].y
        
        # transform
        if xval >= 1.0
            mesh.sMesh.xy[ti].y += sin_x2_transform_func(xval,yval,ω,amp)
        end
    end

    for ti=1:Nfu
        xval = mesh.fMesh.xy[ti].x
        yval = mesh.fMesh.xy[ti].y
        
        # transform
        if xval >= 1.0
            mesh.fMesh.xy[ti].y += sin_x2_transform_func(xval,yval,ω,amp)
        end
    end

    for ti=1:Nfp
        xval = mesh.fMesh.xyp[ti].x
        yval = mesh.fMesh.xyp[ti].y
        
        # transform
        if xval >= 1.0
            mesh.fMesh.xyp[ti].y += sin_x2_transform_func(xval,yval,ω,amp)
        end
    end

    nothing
end

function sin_wall_transform!(mesh::precipMesh,ω::T) where T<:Real
    Ns = length(mesh.sMesh.xy)
    Nfu = length(mesh.fMesh.xy)
    Nfp = length(mesh.fMesh.xyp)
    amp = 0.1

    for ti=1:Ns
        xval = mesh.sMesh.xy[ti].x
        yval = mesh.sMesh.xy[ti].y
        
        # transform
        if xval >= 1.0
            mesh.sMesh.xy[ti].y += sin_wall_transform_func(xval,yval,ω,amp) 
        end
    end

    for ti=1:Nfu
        xval = mesh.fMesh.xy[ti].x
        yval = mesh.fMesh.xy[ti].y
        
        # transform
        if xval >= 1.0
            mesh.fMesh.xy[ti].y += sin_wall_transform_func(xval,yval,ω,amp) 
        end
    end

    for ti=1:Nfp
        xval = mesh.fMesh.xyp[ti].x
        yval = mesh.fMesh.xyp[ti].y
        
        # transform
        if xval >= 1.0
            mesh.fMesh.xyp[ti].y += sin_wall_transform_func(xval,yval,ω,amp) 
        end
    end

    nothing
end

sin_x2_transform_func(x,y,ω,amp) = amp*(cos(2*pi*ω*(x - 1.0)) - 1.0)

function sin_wall_transform_func(x,y,ω,amp) 
    valTop = amp*(cos(2*pi*ω*(x - 1.0)) - 1.0)
    valBottom = -3/8

    ytop = 3/8
    ybottom = -3/8
    
    return (y-ybottom)*valTop/(3/4)
end