using Plots

rev1 = max.(range(0.0, 200.0, 201) .- 100, 0)
profit1 = rev1 .- 10*(1 + 0.1)
plot_rp1 = plot([rev1 profit1], label=["Revenue of Call Option" "Profit of Call Option"])
png(plot_rp1, "rev_prof_plot_1")

rev2 = max.(100.0 .- (range(0.0, 200.0, 201)), 0)
profit2 = rev2 .- 7*(1 + 0.1)
plot_rp2 = plot([rev2 profit2], label=["Revenue of Put Option" "Profit of Put Option"])
png(plot_rp2, "rev_prof_plot_2")

brange = range(0.0, 100.0, 101)
S_0 = 50.0
price_option = 10.0
R = 1.1

payoff_stock = brange
profit_stock = payoff_stock .- S_0*R

sh_payoff_stock = -1*payoff_stock
sh_profit_stock = -1*profit_stock

payoff_call = max.(brange .- S_0, 0.0)
payoff_put = max.(S_0 .- brange, 0.0)
sh_payoff_call = -1*payoff_call
sh_payoff_put = -1*payoff_put

profit_call = payoff_call .- price_option*R
profit_put = payoff_put .- price_option*R 
sh_profit_call = -1*profit_call 
sh_profit_put = -1*profit_put 

price_bond = 50.0
payoff_bond = (price_bond*R).*ones(length(brange))
profit_bond = (price_bond*(R - 1.0)).*ones(length(brange))
sh_payoff_bond = -1*payoff_bond
sh_profit_bond = -1*profit_bond


stock_plot = plot([payoff_stock profit_stock sh_payoff_stock sh_profit_stock], title = "Stock Payoff and Profit", labels = ["Payoff of Stock (Long)" "Profit of Stock (Long)" "Payoff of Stock (Short)" "Profit of Stock (Short)"])
call_plot = plot([payoff_call profit_call sh_payoff_call sh_profit_call], title = "Call Payoff and Profit", labels = ["Payoff of Call (Long)" "Profit of Call (Long)" "Payoff of Call (Short)" "Profit of Call (Short)"])
put_plot = plot([payoff_put profit_put sh_payoff_put sh_profit_put], title = "Put Payoff and Profit", labels = ["Payoff of Put (Long)" "Profit of Put (Long)" "Payoff of Put (Short)" "Profit of Put (Short)"])
bond_plot = plot([payoff_bond profit_bond sh_payoff_bond sh_profit_bond], title = "Bond Payoff and Profit", labels = ["Payoff of Bond (Long)" "Profit of Bond (Long)" "Payoff of Bond (Short)" "Profit of Bond (Short)"])

png(stock_plot, "stock_plot")
png(call_plot, "call_plot")
png(put_plot, "put_plot")
png(bond_plot, "bond_plot")

#plot([payoff_stock profit_stock payoff_call profit_call payoff_put profit_put payoff_bond profit_bond], title = "Payoffs and Profits", labels = ["Payoff of Stock" "Profit of Stock" "Payoff of Call" "Profit of Call" "Payoff of Put" "Profit of Put" "Payoff of Bond" "Profit of Bond"])