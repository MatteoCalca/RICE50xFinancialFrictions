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
$setglobal tmax_minus_10 '126'



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
    bankkappa  'Bank leverage ratio'            / 5     /
    bankgamma  'Bank share of profits'          / 0.9   /
    invadjcost 'Investment adjustment costl'    / 0     /
    ri_kali    'Real interest rate value for calibration' / 0.015 /
;



##  COMPUTE DATA
#_________________________________________________________________________
* In the phase COMPUTE_DATA you should declare and compute all the parameters
* that depend on the data loaded in the previous phase.
$elseif.ph %phase%=='compute_data'





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

    NETWORTH(t,n)          'Net worth of banks (N)'
    SECURITIES(t,n)        'Securities of banks (S)'
    SECURITIES_PRICE(t,n)  'Price of securities of banks (Q)'
    BANK_SHADOWVALUE(t,n)  'Shadow value of banks (phi)'
    BANK_LAMBDA(t,n)       'Lagrange multiplier on the bank leverage constraint (lambda)'
    BANK_DF(t,n)           'Discount factor of banks (Theta)'
    BANK_DIVIDENDS(t,n)    'Dividends of banks to households (Xi)'
;

FIRM_PROFITS.l(t,n) = 0 ;
RK.l(t,n) = prodshare('capital',n)*ykali(t,n)/k_tfp(t,n)  + (1-dk);
SECURITIES_PRICE.l(t,n) = 1 / [ 1 - invadjcost*(i_tfp(t,n)/k_tfp(t,n) - dk) ] ;
NETWORTH.l(t,n) = bankgamma * [RK.l(t,n)*SECURITIES_PRICE.l(t,n)*k_tfp(t,n) - (ri_kali+1)*i_tfp(t,n) ] ; 
BANK_DIVIDENDS.l(t,n) = (1-bankgamma) * [RK.l(t,n)*SECURITIES_PRICE.l(t,n)*k_tfp(t,n) - (ri_kali+1)*i_tfp(t,n) ] ;



##  COMPUTE VARIABLES
#_________________________________________________________________________
* In the phase COMPUTE_VARS, you fix starting points and bounds.
* DO NOT put VAR.l here! (use the declare_vars phase)
$elseif.ph %phase%=='compute_vars'

*ETAX.fx(t,n,ghg) = 0;

#RK.lo(t,n) = 0.08; # Needed because of the transversality condition

YNET.fx(tfirst,n) = ykali(tfirst,n) ;
FIRM_PROFITS.fx(tfirst,n) = FIRM_PROFITS.l(tfirst,n) ;
RK.fx(tfirst,n) = RK.l(tfirst,n) ;
SECURITIES_PRICE.fx(tfirst,n) = SECURITIES_PRICE.l(tfirst,n) ;
NETWORTH.fx(tfirst,n) = NETWORTH.l(tfirst,n) ;
BANK_DIVIDENDS.fx(tfirst,n) = BANK_DIVIDENDS.l(tfirst,n) ;


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
eq_bankFOC # Eq. (9)
eq_bankdf
eq_bankshadowvalue # Eq. (12)
eq_banknetworth # Eq. (14)
eq_bankdividends
eq_securities
eq_kfirmFOC # Eq. (26)


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
                                + (1-dk)*SECURITIES_PRICE(t,n)*K(t,n)      
                                - sum(ghg,  ABATECOST(t,n,ghg) + ETAX(t,n,ghg) * EIND(t,n,ghg)) ;

eq_firmFOC_L(t,n)$(reg(n))..  
        WAGE(t,n) =E= (prodshare('labour',n))*YGROSS(t,n)/(pop(t,n)/1000) * PASSTHROUGH(t,n) ;

eq_firmFOC_K(tm1,t,n)$(reg(n) and pre(tm1,t) and tperiod(t) gt 1)..
        RK(t,n) =E= { [prodshare('capital',n)*YGROSS(t,n)/K(t,n)] * PASSTHROUGH(t,n)  + (1-dk)*SECURITIES_PRICE(t,n) }/SECURITIES_PRICE(tm1,n) ;

eq_firm_passthrough(t,n)$(reg(n))..   
        PASSTHROUGH(t,n) =E= 1 - sum(ghg, ABATECOST_FN(t,n,ghg) + ETAX(t,n,ghg)*EIND_FN(t,n,ghg)) ;

eq_firmFOC_MIU(t,n,ghg)$(reg(n))..    
        ETAX(t,n,ghg) =E= MAC(t,n,ghg) * convy_ghg(ghg) ;

eq_EIND_fn(t,n,ghg)$(reg(n))..    
        EIND_FN(t,n,ghg) =E= sigma(t,n,ghg) * convq_ghg(ghg) * (1-MIU(t,n,ghg));



##### BANK ---------------------------------------------------------
eq_securities(t,tp1,n)$(reg(n) and pre(t,tp1)).. 
        K(tp1,n) =E= SECURITIES(t,n) ;
*Eq. (4)
eq_bankbalancesheet(t,n)$(reg(n))..  
        SECURITIES(t,n) =E= [I(t,n) + NETWORTH(t,n)]/SECURITIES_PRICE(t,n) ;
* Eq. (14)
eq_banknetworth(t,tp1,n)$(reg(n) and pre(t,tp1))..   
        NETWORTH(tp1,n) =E= bankgamma * [ RK(tp1,n)*SECURITIES_PRICE(t,n)*K(tp1,n) - (RI(t,n)+1)*I(t,n) ] ;
*Eq. (9)
eq_bankFOC(t,tp1,n)$(reg(n) and pre(t,tp1))..  
        BANK_LAMBDA(t,n) =E= BANK_DF(tp1,n)*[RK(tp1,n) - (RI(t,n)+1)] * bankkappa ;

eq_bankdf(t,tp1,n)$(reg(n) and pre(t,tp1))..  
        BANK_DF(t,n) =E= {1/[RI(t,n)+1]} * [1 - bankgamma + bankgamma*BANK_SHADOWVALUE(tp1,n)] ;
*Eq. (12)
eq_bankshadowvalue(t,tp1,n)$(reg(n) and pre(t,tp1))..   
        BANK_SHADOWVALUE(t,n) =E= [ BANK_LAMBDA(t,n) * SECURITIES(t,n)*SECURITIES_PRICE(t,n) ]
                                        /  [ bankkappa * NETWORTH(t,n) ]
                                        + BANK_DF(tp1,n)*(RI(t,n)+1) ;

eq_bankdividends(t,tp1,n)$(reg(n) and pre(t,tp1))..  
        BANK_DIVIDENDS(tp1,n) =E= (1-bankgamma) * [ RK(tp1,n)*SECURITIES_PRICE(t,n)*K(tp1,n) - (RI(t,n)+1)*I(t,n) ] ;



##### CAPITAL PRODUCERS --------------------------------------------
* Eq. (26)
eq_kfirmFOC(t,n)$(reg(n))..   SECURITIES_PRICE(t,n) =E= 1 / { 1 - invadjcost*[I(t,n)/K(t,n) - dk] } ; # invadjcost most likely around 1.5 or 2, start with 0



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







#===============================================================================
*     ///////////////////////     REPORTING     ///////////////////////
#===============================================================================
##  REPORT
#_________________________________________________________________________
* Post-solve evaluate report measures
$elseif.ph %phase%=='report'


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

NETWORTH
SECURITIES
SECURITIES_PRICE
BANK_SHADOWVALUE
BANK_LAMBDA
BANK_DF
BANK_DIVIDENDS



$endif.ph
