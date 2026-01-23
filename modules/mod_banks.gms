# mod_banks
*
* Short description
*activate with --mod_banks=1

#____________
# REFERENCES
* -
#=========================================================================
*   ///////////////////////       SETTING      ///////////////////////
#=========================================================================



##  CONF
#_________________________________________________________________________
* Definition of the global flags and settings specific to the module
$ifthen.ph %phase%=='conf'

$setglobal savings 'flexible'


## SETS
#_________________________________________________________________________
* In the phase SETS you should declare all your sets, or add to the existing
* sets the element that you need.
$elseif.ph %phase%=='sets'



## INCLUDE DATA
#_________________________________________________________________________
* In the phase INCLUDE_DATA you should declare and include all exogenous parameters.
* You can also modify the data loaded in data.gms
* Best practice : - create a .gdx containing those and to loading it
*                 - this is the only phase where we should have numbers...
$elseif.ph %phase%=='include_data'


Scalars
    bankkappa  'Bank leverage ratio'            / 4     /
    bankgamma  'Bank share of profits'          / 0.875 /
    invadjcost 'Investment adjustment cost'     / 0.5   /
    bankzeta_fixed 'Fixed bank zeta scaling factor' / 0.03  /
;

* Add NWRATIO to convergence check
vcheck('NWRATIO') = yes;

Parameters
        securities_tfp(t,n) 'Securities TFP'
        networth_tfp(t,n) 'Net worth TFP'
        securities_price_tfp(t,n) 'Securities price TFP'
        rk_tfp(t,n) 'Return to capital TFP'
        ri_tfp(t,n) 'Real interest rate TFP'
        cpc_tfp(t,n) 'Consumption per capita TFP'
        d_tfp(t,n) 'Deposits calibration (STOCK, not flow)'
        bank_zeta_param(t,n) 'Calibrated new bank entry equity'
        firm_s_tfp(t,n) 'Calibrated investment rate (I/Y)'
;



##  COMPUTE DATA
#_________________________________________________________________________
* In the phase COMPUTE_DATA you should declare and compute all the parameters
* that depend on the data loaded in the previous phase.
$elseif.ph %phase%=='compute_data'

* K/Y overrides: countries with high K/Y ratios cause solver instability
* CRO (4.59), SLO (5.26), BGR (3.95) need K/Y adjustment for convergence
k_tfp(t,'cro') = ykali(t,'cro') * 3.5;
k_tfp(t,'slo') = ykali(t,'slo') * 3.5;
k_tfp(t,'bgr') = ykali(t,'bgr') * 3.5;

cpc_tfp(t,n) = [ykali(t,n) - i_tfp(t,n)] / (pop(t,n)/1000);

loop((t,tp1)$pre(t,tp1),
   ri_tfp(t,n) = (1+prstp)*(cpc_tfp(tp1,n)/cpc_tfp(t,n))**(elasmu/tstep) - 1;
);

securities_price_tfp(t,n)$tnolast(t) = 1 / [ 1 - invadjcost*(i_tfp(t,n)/k_tfp(t,n) - dk(n)) ] ;

rk_tfp(t,n)$tnolast(t) = prodshare('capital',n)*ykali(t,n)/k_tfp(t,n) + (1-dk(n))*securities_price_tfp(t,n);

* Banking calibration (assumes leverage binds)
* SECURITIES(t) determines K(t+1) via eq_securities, so calibrate to k_tfp(t+1)
loop((t,tp1)$(pre(t,tp1) and tnolast(t)),
    securities_tfp(t,n) = k_tfp(tp1,n);
);
* Calibrate N and D from leverage constraint: S×Q = κ×N, so N = S×Q/κ
networth_tfp(t,n)$tnolast(t) = securities_tfp(t,n) * securities_price_tfp(t,n) / bankkappa;
d_tfp(t,n)$tnolast(t) = securities_tfp(t,n) * securities_price_tfp(t,n) * (bankkappa - 1) / bankkappa;

* Bank zeta: new equity proportional to net worth (about 3% of N)
bank_zeta_param(t,n)$tnolast(t) = bankzeta_fixed * networth_tfp(t,n);

* Calibrated investment rate for proper initialization
firm_s_tfp(t,n)$(ykali(t,n) gt 0) = i_tfp(t,n) / ykali(t,n);

* Diagnostic calibration parameters
Parameters
    calib_k_y_ratio(t,n) 'K/Y ratio from calibration'
    calib_i_k_ratio(t,n) 'I/K ratio from calibration'
    calib_net_k_growth(t,n) 'I/K - dk: positive means K grows in calibration'
    calib_n_y_ratio(t,n) 'N/Y ratio from calibration'
;
calib_k_y_ratio(t,n)$(ykali(t,n) gt 0) = k_tfp(t,n) / ykali(t,n);
calib_i_k_ratio(t,n)$(k_tfp(t,n) gt 0) = i_tfp(t,n) / k_tfp(t,n);
calib_net_k_growth(t,n)$(k_tfp(t,n) gt 0) = calib_i_k_ratio(t,n) - dk(n);
calib_n_y_ratio(t,n)$(ykali(t,n) gt 0 and tnolast(t)) = networth_tfp(t,n) / ykali(t,n);


##  DECLARE VARIABLES
#_________________________________________________________________________
* In the phase DECLARE VARS, you can DECLARE new variables for your module.
* Remember that by modifying sets, you already have some variables for free.
$elseif.ph %phase%=='declare_vars'

VARIABLES
    FIRM_PROFITS(t,n)   'Firm profits'
    WAGE(t,n)      'Wage'
    RK(t,n)             'Return to capital'
    PASSTHROUGH(t,n)    'Passthrough rate of climate policy to households'
    HANDOUTS(t,n)       'Government handouts'
    ETAX(t,n,ghg)       'Emission tax rate'

    EIND_FN(t,n,ghg)    'EIND without the YGROSS component'

    D(t,n)                 'Deposits of banks (D)'
    NETWORTH(t,n)          'Net worth of banks (N)'
    SECURITIES(t,n)        'Securities of banks (S)'
    SECURITIES_PRICE(t,n)  'Price of securities of banks (Q)'
    BANK_SHADOWVALUE(t,n)  'Shadow value of banks (phi)'
    BANK_LAMBDA(t,n)       'Lagrange multiplier on the bank leverage constraint (lambda)'
    BANK_DF(t,n)           'Discount factor of banks (Theta)'
    BANK_DIVIDENDS(t,n)    'Dividends of banks to households (Xi)'

    FIRM_I(t,n)         'Investment of firms'
    FIRM_S(t,n)         'Investment rate of firms'
;

* Post-solve diagnostic parameters
**PARAMETERS
**    leverage_slack(t,n) 'Post-solve: Leverage slack = kappa*N - S*Q'
**    leverage_binding(t,n) 'Post-solve: 1 if leverage binding, 0 if slack'
**    d_ratio(t,n) 'Post-solve: D/Y ratio'
**    bank_profit_margin(t,n) 'Post-solve: Bank profit margin = RK*κ - (1+RI)*(kappa-1)'
**    credit_spread(t,n) 'Post-solve: Credit spread RK - (1+RI)'
**    leverage_utilization(t,n) 'Post-solve: S×Q / (κN) - 1 if binding'
**;

FIRM_PROFITS.l(t,n) = 0 ;
RK.l(t,n) = rk_tfp(t,n);
SECURITIES_PRICE.l(t,n) = securities_price_tfp(t,n) ;
NETWORTH.l(t,n) = networth_tfp(t,n) ;
D.l(t,n) = d_tfp(t,n) ;
SECURITIES.l(t,n) = securities_tfp(t,n) ;
BANK_DIVIDENDS.l(t,n) = networth_tfp(t,n) / bankgamma * (1-bankgamma) ;
FIRM_S.l(t,n) = firm_s_tfp(t,n) ;
FIRM_I.l(t,n) = i_tfp(t,n) ;



##  COMPUTE VARIABLES
#_________________________________________________________________________
* In the phase COMPUTE_VARS, you fix starting points and bounds.
* DO NOT put VAR.l here! (use the declare_vars phase)
$elseif.ph %phase%=='compute_vars'

* Fix initial period values
YNET.fx(tfirst,n) = ykali(tfirst,n) ;
FIRM_PROFITS.fx(tfirst,n) = FIRM_PROFITS.l(tfirst,n) ;
RK.fx(tfirst,n) = RK.l(tfirst,n) ;
NETWORTH.fx(tfirst,n) = NETWORTH.l(tfirst,n);

* Variable bounds
RI.lo(t,n) = 1e-3 ;
RI.up(t,n)$(not t5last(t)) = 0.5 ;
SECURITIES_PRICE.lo(t,n) = 1e-3 ;
BANK_LAMBDA.lo(t,n) = 0 ;
BANK_SHADOWVALUE.fx(t,n)$tlast(t) = 1 ;
NETWORTH.lo(t,n)$(not tfirst(t)) = 1e-4 ;
D.lo(t,n) = 1e-4 ;
FIRM_S.lo(t,n) = 0.001 ;
FIRM_S.up(t,n) = 1 ;


#=========================================================================
*   ///////////////////////     OPTIMIZATION    ///////////////////////
#=========================================================================

##  EQUATION LIST
#_________________________________________________________________________
* List of equations
* One per line.
$elseif.ph %phase%=='eql'

eq_firmprofits
eq_firmFOC_L
eq_firmFOC_K
eq_firm_passthrough
eq_firmFOC_MIU
eq_EIND_fn
eq_govbudget

eq_bankbalancesheet # Eq. (4)
eq_bankFOC # Eq. (9) - Bank first-order condition
eq_bankdf # Bank discount factor
eq_bankshadowvalue # Eq. (12) - Shadow value of net worth
eq_banknetworth # Eq. (14) - Net worth evolution
eq_bankdividends
eq_securities
eq_leverage # Leverage constraint (always binding)
eq_household_budget # Eq. (2): Household budget constraint
eq_kfirmFOC # Eq. (26) - Investment adjustment cost pricing
eq_firmI
eq_goods_clearing # Goods market clearing: Y = C + I (CLIFF Section 2.6)


##  EQUATIONS
#_________________________________________________________________________
* In the phase EQS, you can include new equations to the model.
* The equations are always included.
* Best practice : - condition your equation to be able to do a run with tfix(t)
$elseif.ph %phase%=='eqs'


##### FIRM ---------------------------------------------------------
eq_firmprofits(tm1,t,n)$(reg(n) and pre(tm1,t) and tperiod(t) gt 1)..   
        FIRM_PROFITS(t,n) =E=  YGROSS(t,n)
                                - WAGE(t,n)*(pop(t,n)/1000)
                                - RK(t,n)*SECURITIES_PRICE(tm1,n)*K(t,n)
                                + (1-dk(n))*SECURITIES_PRICE(t,n)*K(t,n)      
                                - sum(ghg,  ABATECOST(t,n,ghg) + ETAX(t,n,ghg) * EIND(t,n,ghg)) ;

eq_firmFOC_L(t,n)$(reg(n))..  
        WAGE(t,n) =E= (prodshare('labour',n))*YGROSS(t,n)/(pop(t,n)/1000) * PASSTHROUGH(t,n) ;

eq_firmFOC_K(tm1,t,n)$(reg(n) and pre(tm1,t) and tperiod(t) gt 1)..
        RK(t,n) =E= { [prodshare('capital',n)*YGROSS(t,n)/K(t,n)] * PASSTHROUGH(t,n)  + (1-dk(n))*SECURITIES_PRICE(t,n) }/SECURITIES_PRICE(tm1,n) ;

eq_firm_passthrough(t,n)$(reg(n))..   
        PASSTHROUGH(t,n) =E= 1 - sum(ghg, ABATECOST_FN(t,n,ghg) + ETAX(t,n,ghg)*EIND_FN(t,n,ghg)) ;

eq_firmFOC_MIU(t,n,ghg)$(reg(n))..    
        ETAX(t,n,ghg) =E= MAC(t,n,ghg) * convy_ghg(ghg) ;

eq_EIND_fn(t,n,ghg)$(reg(n))..    
        EIND_FN(t,n,ghg) =E= sigma(t,n,ghg) * convq_ghg(ghg) * (1-MIU(t,n,ghg));



##### BANK ---------------------------------------------------------
* Eq. (4) - Balance sheet: S×Q = D + N
eq_bankbalancesheet(t,n)$(reg(n))..
        SECURITIES(t,n) =E= [D(t,n) + NETWORTH(t,n)]/SECURITIES_PRICE(t,n) ;

* Leverage constraint (always binding): S×Q = κ×N
eq_leverage(t,n)$(reg(n))..
        SECURITIES(t,n) * SECURITIES_PRICE(t,n) =E= bankkappa * NETWORTH(t,n) ;

* Securities determine next period capital
eq_securities(t,tp1,n)$(reg(n) and pre(t,tp1))..
        K(tp1,n) =E= SECURITIES(t,n) ;

* Eq. (14) - Net worth evolution
eq_banknetworth(t,tp1,n)$(reg(n) and pre(t,tp1))..
        NETWORTH(tp1,n) =E= bankgamma * [ RK(tp1,n)*SECURITIES_PRICE(t,n)*K(tp1,n) - (RI(t,n)+1)*D(t,n) ]
                           + bank_zeta_param(t,n) ;

* Eq. (9) - Bank first-order condition determines λ
eq_bankFOC(t,tp1,n)$(reg(n) and pre(t,tp1))..
        BANK_LAMBDA(t,n) =E= BANK_DF(tp1,n)*[RK(tp1,n) - (RI(t,n)+1)] * bankkappa ;

* Bank discount factor: Θ = (1/R) × [1 - γ + γ×φ]
eq_bankdf(t,tp1,n)$(reg(n) and pre(t,tp1))..
        BANK_DF(tp1,n) =E= {1/[RI(t,n)+1]} * [1 - bankgamma + bankgamma*BANK_SHADOWVALUE(tp1,n)] ;

* Eq. (12) - Shadow value of net worth
eq_bankshadowvalue(t,tp1,n)$(reg(n) and pre(t,tp1))..
        BANK_SHADOWVALUE(t,n) =E= [ BANK_LAMBDA(t,n) * SECURITIES(t,n)*SECURITIES_PRICE(t,n) ]
                                   / [ bankkappa * NETWORTH(t,n) ]
                                   + BANK_DF(tp1,n)*(RI(t,n)+1) ;

* Bank dividends to households
eq_bankdividends(t,tp1,n)$(reg(n) and pre(t,tp1))..
        BANK_DIVIDENDS(tp1,n) =E= (1-bankgamma) * [ RK(tp1,n)*SECURITIES_PRICE(t,n)*K(tp1,n) - (RI(t,n)+1)*D(t,n) ] ;

* CLIFF eq. (2) - Household budget constraint: C + D = wL + R×D_{-1} + Ξ + Π + T
eq_household_budget(t,tm1,n)$(reg(n) and pre(tm1,t) and not tfirst(t))..
        C(t,n) + D(t,n) =E= WAGE(t,n) * (pop(t,n)/1000)
                          + (1 + RI(tm1,n)) * D(tm1,n)
                          + BANK_DIVIDENDS(t,n)
                          + FIRM_PROFITS(t,n)
                          + HANDOUTS(t,n) ;

##### CAPITAL PRODUCERS --------------------------------------------
* Eq. (26) - Securities price determined by investment adjustment costs
eq_kfirmFOC(t,n)$(reg(n))..
        SECURITIES_PRICE(t,n) =E= 1 / { 1 - invadjcost*[FIRM_I(t,n)/K(t,n) - dk(n)] } ;

eq_firmI(t,n)$(reg(n))..
        FIRM_I(t,n) =E= FIRM_S(t,n) * YGROSS(t,n) ;

* CLIFF Section 2.6 - Goods market clearing: Y = C + I + Ω
eq_goods_clearing(t,n)$(reg(n))..
        C(t,n) + FIRM_I(t,n) =E= Y(t,n)
$if set mod_adaptation         - sum(g, I_ADA(g,t,n))
$if set mod_natural_capital    - sum(type, NAT_INV(type,t,n))
;


##### GOVERNMENT ---------------------------------------------------
eq_govbudget(t,n)$(reg(n))..   HANDOUTS(t,n) =E= sum[ghg, ETAX(t,n,ghg)*EIND(t,n,ghg)]
;


##  FIX VARIABLES
#_________________________________________________________________________
$elseif.ph %phase%=='fix_variables'








##  BEFORE SOLVE
#_________________________________________________________________________
* In the phase BEFORE_SOLVE, you can update parameters (fixed
* variables, ...) inside the nash loop and right before solving the
* model. This is typically done for externalities, spillovers, ...
* Best practice: record the variable that you update across iterations.
* Remember that you are inside the nash loop, so you cannot declare
* parameters, ...
$elseif.ph %phase%=='before_solve'






##  AFTER SOLVE
#_________________________________________________________________________
* In the phase AFTER_SOLVE, you compute what must be propagated across the
* regions after one bunch of parallel solving.
$elseif.ph %phase%=='after_solve'

* Post-solve diagnostics
**leverage_slack(t,n)$(reg(n)) = bankkappa * NETWORTH.l(t,n)
**                             - SECURITIES.l(t,n) * SECURITIES_PRICE.l(t,n);
**leverage_binding(t,n)$(reg(n) and leverage_slack(t,n) ne 0) =
**    1$(leverage_slack(t,n) lt 1e-6) + 0$(leverage_slack(t,n) ge 1e-6);
**d_ratio(t,n)$(reg(n) and Y.l(t,n) gt 0) = D.l(t,n) / Y.l(t,n);
**bank_profit_margin(t,n)$(reg(n)) =
**    RK.l(t,n) * bankkappa - (RI.l(t,n) + 1) * (bankkappa - 1);

* Track NWRATIO for convergence check
viter(iter,'NWRATIO',t,n)$(nsolve(n) and Y.l(t,n) gt 0) = NETWORTH.l(t,n) / Y.l(t,n);

* Additional diagnostics for analyzing N growth vs K growth
**credit_spread(t,n)$(reg(n)) = RK.l(t,n) - (RI.l(t,n) + 1);
**leverage_utilization(t,n)$(reg(n) and NETWORTH.l(t,n) gt 0) =
**    SECURITIES.l(t,n) * SECURITIES_PRICE.l(t,n) / (bankkappa * NETWORTH.l(t,n));




#===============================================================================
*     ///////////////////////     REPORTING     ///////////////////////
#===============================================================================
##  REPORT
#_________________________________________________________________________
* Post-solve evaluate report measures
$elseif.ph %phase%=='report'

* Social Cost of Carbon for banking module
* Uses eq_goods_clearing.m instead of eq_cc.m (both are goods market clearing)
scc(t,n,ghg)$(nsolve(n) and year(t) le 2200 and Y.l(t,n) gt 0) =
    -1e3*sum(nn$nsolve(nn), div0(eq_e.m(t,nn,ghg) , eq_goods_clearing.m(t,nn)) );


##  GDX ITEMS
#_________________________________________________________________________
* List the items to be kept in the final gdx
$elseif.ph %phase%=='gdx_items'

FIRM_PROFITS
WAGE
RK
PASSTHROUGH
HANDOUTS
ETAX

EIND_FN

D
NETWORTH
SECURITIES
SECURITIES_PRICE
BANK_SHADOWVALUE
BANK_LAMBDA
BANK_DF
BANK_DIVIDENDS

FIRM_I
FIRM_S

d_tfp
bank_zeta_param
firm_s_tfp
bankkappa
bankgamma
**credit_spread
**leverage_utilization
**
**calib_k_y_ratio
**calib_i_k_ratio
**calib_net_k_growth
**calib_n_y_ratio
*bank_profit_margin
dk

$endif.ph
