using CSV, DataFrames, Plots

UData = CSV.read("fredgraph.csv", DataFrame)
#remove the missing 10/1/25 data
dropmissing!(UData)

date_vec = UData[:, "observation_date"]
u_vec = UData[:, "UNEMPLOY"]./UData[:, "CLF16OV"]
us_vec = UData[:, "UEMPLT5"]./UData[:, "CLF16OV"]
rec_vec = UData[:, "USREC"]
l_vec = UData[:, "CLF16OV"]

f_vec = zeros(length(u_vec) - 1)
s_vec = zeros(length(f_vec))
for i in 1:length(f_vec)
    f_vec[i] = 1 - (u_vec[i+1] - us_vec[i+1])/u_vec[i]
    s_vec[i] = (u_vec[i+1] - (1-f_vec[i])*u_vec[i])/(1-u_vec[i])
end

f_bar, s_bar = sum(f_vec)/length(f_vec), sum(s_vec)/length(s_vec)

hazard_plot = plot(date_vec[1:length(date_vec) - 1],[f_vec,s_vec])
readhazard_plot = plot(date_vec[1:length(date_vec) - 1],[f_vec,s_vec.*10])