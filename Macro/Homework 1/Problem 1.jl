using CSV, DataFrames, Plots, Parameters

UData = CSV.read("fredgraph.csv", DataFrame)
#remove the missing 10/1/25 data
dropmissing!(UData)

date_vec = UData[:, "observation_date"]
u_vec = UData[:, "UNEMPLOY"]./UData[:, "CLF16OV"]
us_vec = UData[:, "UEMPLT5"]./UData[:, "CLF16OV"]
rec_vec = UData[:, "USREC"]

rec_start = findall(diff(rec_vec) .> .1)
rec_end = findall(diff(rec_vec) .< -.1)
#recs = findall(rec_vec .== 1)
#rec_breaks = findall(diff(recs) .> 1)

function add_recs()
   for i in 1:length(rec_start)
        plot!([date_vec[rec_start[i]], date_vec[rec_end[i]]], [0,0], fillrange = (0,0.7), fillalpha = 0.3, fillcolor = :gray, label = "")
   end 
end
#rec_breaks = diff(recs) .> .1
#rec_start = recs[[true; rec_breaks]]
#rec_end = recs[[rec_breaks; true]]

f_vec = zeros(length(u_vec) - 1)
s_vec = zeros(length(f_vec))

for i in 1:length(f_vec)
    f_vec[i] = 1 - (u_vec[i+1] - us_vec[i+1])/u_vec[i]
    s_vec[i] = (u_vec[i+1] - (1-f_vec[i])*u_vec[i])/(1-u_vec[i])
end

uss_vec = s_vec ./ (s_vec .+ f_vec)

f_bar, s_bar = sum(f_vec)/length(f_vec), sum(s_vec)/length(s_vec)
u_bar, uss_bar = sum(u_vec)/length(u_vec), sum(uss_vec)/length(uss_vec)
mcf = 1 - s_bar - f_bar
hl = log(2)/log(mcf)

hazard_plot = plot(date_vec[1:length(date_vec) - 1],[f_vec,s_vec], ylim = (0,0.7), title = "Hazard Rates since 1948", label = ["Job Finding" "Separation"], legend = :bottomleft)
add_recs()
readhazard_plot = plot(date_vec[1:length(date_vec) - 1],[f_vec,s_vec.*10], ylim = (0,0.7), title = "Hazard Rates since 1948 (magnified)", label = ["Job Finding" "Separation (x10)"], legend = :bottomleft)
add_recs()

f_plot = plot(date_vec[1:length(date_vec) - 1],f_vec, ylim = (0,0.7), title = "Job-Finding Rates since 1948", label = "Job Finding", legend = :bottomleft)
add_recs()
s_plot = plot(date_vec[1:length(date_vec) - 1],s_vec, ylim = (0,0.05), title = "Separation Rates since 1948", label = "Separation", legend = :bottomleft)
add_recs()

u_plot = plot(date_vec[1:length(date_vec) - 1], [u_vec[1:length(u_vec)-1] uss_vec], ylim = (0.0, 0.10), title = "Unemployment since 1948", label = ["Unemployment" "Steady-State Unemployment"], legend = :bottomleft)
add_recs()

png(hazard_plot, "hazardplot")
png(readhazard_plot, "readablehazardplot")
png(u_plot, "unemploymentplot")
png(f_plot, "fplot")
png(s_plot, "splot")