clear

cd "C:\Users\Alejandro\Box Sync\S19\Quant Macro\Homework\Homework 3\extracts"
use rscfp2016.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 2016

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)

save 2016.dta, replace
clear

use rscfp2013.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 2013

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 2013.dta, replace
clear

use rscfp2010.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 2010

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 2010.dta, replace
clear

use rscfp2007.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 2007

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 2007.dta, replace
clear

use rscfp2004.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 2004

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 2004.dta, replace
clear

use rscfp2001.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 2001

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 2001.dta, replace
clear

use rscfp1998.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 1998

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 1998.dta, replace
clear

use rscfp1995.dta


keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 1995

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 1995.dta, replace
clear

use rscfp1992.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 1992

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 1992.dta, replace
clear

use rscfp1989.dta

keep wgt income wageinc houses asset NH_MORT ccbal

egen sumwgt = sum(wgt)

gen wgt_alt = wgt/sumwgt
gen alt_inc = wgt_alt*income
gen alt_winc = wgt_alt*wageinc
gen alt_house = wgt_alt*houses
gen alt_asset = wgt_alt*houses
gen alt_mort = wgt_alt*NH_MORT
gen alt_bal = wgt_alt*ccbal
gen year = 1989

egen mean_inc = mean(alt_inc)
egen med_inc = median(alt_inc)

egen mean_winc = mean(alt_winc)
egen med_winc = median(alt_winc)

egen mean_house = mean(alt_house)
egen med_house = median(alt_house)

egen mean_asset = mean(alt_asset)
egen med_asset = median(alt_asset)

egen mean_mort = mean(alt_mort)
egen med_mort = median(alt_mort)

egen mean_bal = mean(alt_bal)
egen med_bal = median(alt_bal)
save 1989.dta, replace

clear

use 2016.dta
append using 2013.dta
append using 2010.dta
append using 2007.dta
append using 2004.dta
append using 2001.dta
append using 1998.dta
append using 1995.dta
append using 1992.dta
append using 1989.dta
save scf.dta, replace

twoway connected mean_inc med_inc year
graph export inc.png, replace
twoway connected mean_winc med_winc year
graph export winc.png, replace
twoway connected mean_house med_house year
graph export house.png, replace
twoway connected mean_asset med_asset year
graph export asset.png, replace
twoway connected mean_mort med_mort year
graph export mort.png, replace
twoway connected mean_bal med_bal year
graph export bal.png, replace


