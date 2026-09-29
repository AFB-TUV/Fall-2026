using Plots

R = [0.156, 0.072, 0.135, 0.198]
β = [2.1, 0.7, 1.75, 2.8]

βrplot = plot(β, R, xlabel = "β", ylabel = "Return", title = "Expected Return vs. β", legend = false)
#definitely a straight line
png(βrplot, "hw5_plot")