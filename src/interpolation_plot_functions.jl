# collection of functions to create movies/gifs/figures of quantities
# interpolated from 2D onto a line to make interpretation easier
using FileIO, Plots

function make_gif_of_files_in(file_folder::String,firstInt::Int,lastInt::Int)
    anim = @animate for j=firstInt:lastInt
        filename = string(file_folder,"/interpVals_$(j).jld2")
        time, interpVals = load(filename, "time", "interpVals")
        N,M = size(interpVals)

        nameArr = ["x","y","psi_A","psi_B","theta_s","|v|"]
        pL = plot()
        pR = plot()
        printTime = round(time,digits=2)
        for i=3:4
            plot!(pL,interpVals[:,2],interpVals[:,i],label=nameArr[i],title="t = $(printTime)",ylims=(0,0.005))
            plot!(pR,interpVals[:,2],interpVals[:,i+2],label=nameArr[i+2],title="t = $(printTime)",ylims=(0,1.4))
        end
        l = @layout [a{0.5w} b{0.5w}]
        plot(pL,pR,layout=l)
    end every 1

    gif(anim, "fpstest.gif", fps = 15)

    nothing
end