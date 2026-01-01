//ssc install mat2txt

/********************************************************************/
/* QQ Correlation Heatmap with Manual Block Bootstrap (Avik Sinha)  */
/********************************************************************/

clear
set more off
set linesize 255

/********************************************************************/
/* 1. Load data and check variables                                 */
/********************************************************************/
use data_core.dta, clear

ds
local vars "`r(varlist)'"
if wordcount("`vars'") != 2 {
    display as error "ERROR: Exactly two variables are required."
    exit 198
}

local x : word 1 of `vars'
local y : word 2 of `vars'
display "Variables: `x', `y'"

/********************************************************************/
/* 2. Create quantiles                                              */
/********************************************************************/
local Q = 19
xtile qx = `x', n(`=`Q'+1')
xtile qy = `y', n(`=`Q'+1')

/********************************************************************/
/* 3. Indicator variables (REVISED: highest quantile = 19 valid)    */
/********************************************************************/

forvalues i = 1/`Q' {
    gen Ix_`i' = (qx <= `i')
    gen Iy_`i' = (qy <= `i')
}

/********************************************************************/
/* 4. Compute QQ correlation matrix                                 */
/********************************************************************/
matrix QQ = J(`Q',`Q',.)
forvalues i = 1/`Q' {
    forvalues j = 1/`Q' {
        quietly corr Ix_`i' Iy_`j'
        matrix QQ[`i',`j'] = r(rho)
    }
}

/********************************************************************/
/* 5. Manual moving block bootstrap                                 */
/********************************************************************/
local B = 1000
local T = _N
local L = floor(_N^(1/3) + 0.9999)

tempname QQboot
matrix `QQboot' = J(`B', `Q'*`Q', .)

forvalues b = 1/`B' {
    preserve
    gen bootid = .
    local pos = 1
    while `pos' <= `T' {
        local start = runiformint(1, `T' - `L' + 1)
        forvalues j = 0/`=`L'-1' {
            if `pos' <= `T' {
                replace bootid = `start' + `j' in `pos'
                local ++pos
            }
        }
    }
    sort bootid
    keep if bootid <= `T'

    tempname M
    matrix `M' = J(`Q',`Q',.)
    forvalues i = 1/`Q' {
        forvalues j = 1/`Q' {
            quietly corr Ix_`i' Iy_`j'
            matrix `M'[`i',`j'] = r(rho)
        }
    }

    local k = 1
    forvalues i = 1/`Q' {
        forvalues j = 1/`Q' {
            matrix `QQboot'[`b',`k'] = `M'[`i',`j']
            local ++k
        }
    }
    restore
}

/********************************************************************/
/* 6. Compute bootstrap p-values                                    */
/********************************************************************/
matrix Pboot = J(`Q',`Q',.)

local k = 1
forvalues i = 1/`Q' {
    forvalues j = 1/`Q' {
        scalar rho0 = QQ[`i',`j']
        scalar cnt = 0
        forvalues b = 1/`B' {
            if abs(`QQboot'[`b',`k']) >= abs(rho0) {
                scalar cnt = cnt + 1
            }
        }
        matrix Pboot[`i',`j'] = cnt / `B'
        local ++k
    }
}

/********************************************************************/
/* 7. Export matrices for MATLAB                                    */
/********************************************************************/
mat2txt, matrix(QQ)    saving(QQ_matrix.csv) replace
mat2txt, matrix(Pboot) saving(P_matrix.csv)  replace

/********************************************************************/
/* 8. Call MATLAB                                                   */
/********************************************************************/

file open fh using "varnames.txt", write replace
file write fh "`x'" _n "`y'"
file close fh

shell cmd /c ""matlab.exe" -nosplash -nodesktop -r "cd('`c(pwd)''); drawQQheatmap; exit""

drop Ix_* Iy_* qx qy
