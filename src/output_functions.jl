mutable struct simulationTimer
  initialTime::Float64
  wallTime::Float64
  timeStepStartTime::Float64
  simulationTime::Float64
  finalTime::Float64
  stepIndex::Int64
  printIndex::Int64
  printSkip::Int64
  Δt::Float64
end

function exportVTK!(timer::simulationTimer,p::precipParam,soln::precipSolution,mesh::precipMesh)
  # scalar data
  sData = ScalarData(soln.ψA.u,soln.ψB.u,soln.ψC.u,soln.θs.u,soln.θf.u)
  sNames = ScalarNames(p.ψAname,p.ψBname,p.ψCname,p.θsname,p.θfname)
  path = Path(p.SaveFolder,string(p.SaveName,"_scalars_",timer.printIndex))
  vtksave(mesh.sMesh,sData,sNames,path;time=timer.simulationTime)

  # fluids data
  sData = ScalarData(soln.fluid.p)
  sNames = ScalarNames("pressure")
  vData = VectorData([soln.fluid.u,soln.fluid.v])
  vNames = VectorNames("darcy_velocity")
  path = Path(p.SaveFolder,string(p.SaveName,"_fluid_",timer.printIndex))
  vtksave(mesh.fMesh,sData,sNames,vData,vNames,path;time=timer.simulationTime)
end

function doPrinting!(timer::simulationTimer,p::precipParam,mesh::precipMesh,soln::precipSolution)
  if p.printBool
    printProgressBar(timer)
    exportVTK!(timer,p,soln,mesh)
  end
  updatePrintIndex!(timer)
end

function initializeTimer(p::precipParam)
    return simulationTimer(time(),0.0,0.0,eps(),p.Tfinal,1,1,p.printSkipIndex,p.Δt)
end

function updateWallTime!(timer::simulationTimer)
    timer.wallTime = time() - timer.initialTime
end
function updateSimulationTime!(timer::simulationTimer)
    timer.simulationTime += timer.Δt
end
function updateStepIndex!(timer::simulationTimer)
    timer.stepIndex += 1
end
function updatePrintIndex!(timer::simulationTimer)
    timer.printIndex += 1
end
function intializeTimeSteppingTimer!(timer::simulationTimer)
  timer.timeStepStartTime::Float64
end

function isPrintIndex(timer::simulationTimer)
    if mod(timer.stepIndex,timer.printSkip)==0
        return true
    end
    return false
end

function printSimulationFinishedMessage(timer::simulationTimer)
  println()
  println("  Simulation done ... took ",sToHMS(timer.wallTime))
end

function printProgressBar(timer::simulationTimer)
  percentage = Int(ceil(timer.simulationTime/timer.finalTime*100))
  progressBarDisplay(percentage,timer)
end

function progressBarDisplay(percentage,timer::simulationTimer)
  TERMINALLENGTH = 80
  
  percent = string("\r  [",percentage,"%]")
  skip1  = 11-length(percent)
  exportFigure = string("exported fig. ",timer.printIndex)
  skip2  = 31-(length(percent)+skip1+length(exportFigure))
  updateWallTime!(timer)
  timeTaken = string("wall time: ",sToHMS(timer.wallTime),)
  skip3  = 48-(length(percent)+skip1+length(exportFigure)+skip2+length(timeTaken))
  tremain = (timer.finalTime - timer.simulationTime)*(timer.wallTime - timer.timeStepStartTime)/timer.simulationTime
  timeRemaining1 = string("remaining: ")
  timeRemaining2 = string(sToHMS(tremain),"\r")

  iterChar(TERMINALLENGTH," ")

  print(percent)

  space()
  iterChar(skip1,"-")
  space()
  print(exportFigure)

  space()
  iterChar(skip2,"-")
  space()
  print(timeTaken)

  space()
  iterChar(skip3,"-")
  space()
  print(timeRemaining1)
  if timer.printIndex > 4
    print(timeRemaining2)
  else
    print("...\r")
  end
end

function iterChar(iterations::Int,character::String)
  for i=1:iterations
    print(character)
  end
end

function space()
    print(" ")
end
  
function doInitialPrinting(p::precipParam)
  #if p.printBool
    println("  running ",returnRunIDName(p))
    println("  data being saved in: ",p.SaveFolder)
    
    
    run(`mkdir -p $(p.SaveFolder)`)
    exportParameters(p)
  #end
end

function returnRunIDName(p::precipParam)
  if p.runID == 1
    return "main program"
  elseif p.runID == 2
    return "pure diffusion one chemical test"
  elseif p.runID == 3
    return "pure advection one chemical test"
  elseif p.runID == 4
    return "diffusion + advection one chemical test"
  elseif p.runID == 5
    return "multiphase brinkman no chemical test"
  elseif p.runID == 6
    return "diffusion + reaction two chemical test"
  elseif p.runID == 7
    return "pure reaction test"
  else
    error("runID not recognized")
  end
end

"""
  sToHMS(s::Float64)

Converts Float64 in seconds to HH:MM:SS string format
"""
function sToHMS(s)
  minutes = 0.0
  hours   = 0.0
  timestr = ""

  # compute times
  hours   = Int(floor(s/3600))
  minutes = Int(floor(s/60 - 60*hours))
  seconds = Int(round(s - 3600*hours - 60*minutes,digits=0))

  hours < 10   ? hstr = "0" : hstr = ""
  minutes < 10 ? mstr = "0" : mstr = ""
  seconds < 10 ? sstr = "0" : sstr = ""

  if hours > 0
    timestr=string(timestr,hstr,hours,"h:")
  end
  if (minutes>0) || (hours > 0)
    timestr=string(timestr,mstr,minutes,"m:")
  end

  timestr = string(timestr,sstr,seconds,"s")

  return timestr
end

function exportParameters(p::precipParam)
  dict = type2dict(p)

  filename = string(p.SaveFolder,"/",p.SaveName,"_parameters.txt")
  open(filename, "w") do io
    for (key,value) in dict
      str=string(key,": ",value,"\n")
      write(io,str)
    end
  end

  nothing
end
