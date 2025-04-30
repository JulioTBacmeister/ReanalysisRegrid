#!/usr/bin/env python
# Import packages 
import sys
import os


print( "from inside GenRegrid.py")
print( sys.path ) 

import argparse as arg



import xarray as xr
import numpy as np
import pandas as pd

try:
    import ESMF as E
except ImportError:
    import esmpy as E

import importlib
import glob
import copy
import time


from Regridder import GlobalVarClass
from Regridder.GlobalVarClass import Gv

from Regridder import scripGen as SG
from Regridder import esmfRegrid as erg
from Utils import MyConstants as Con
from Utils import GridUtils as GrU
from Utils import MakePressures as MkP
# "ChatGPI version" --- 
from Regridder import VertRegridFlexLL as vrg

from Utils import humiditycalcs as hum



# Reload local packages that are under
# development
importlib.reload( erg )
importlib.reload( vrg )
importlib.reload( SG )
importlib.reload( MkP )
importlib.reload( hum )
importlib.reload( GrU )
importlib.reload( Con )
#importlib.reload( Gv )

#importlib.reload( Globals2 )

# Physical Constants
Rdry = Con.Rdry() # 


#Gv = Globals2.VariableContainer()

#-------------------------------------------------------------
#  Naming conventions
#-------------------------------------------------------------
# aaa_{CAM,ERA} 
# Indicates the immediate provenance of a variable, e.g.,
#      phis_CAM ==> phis from CAM on CAM's grid
#      phis_ERA ==> phis from ERA on ERA's grid
# lower case 'phis' indicates this is ndarray-like 
#
# aaa_{CAM,ERA}_x{ERA,CAM}
# Indicates variable has been remapped horizontall. So, e.g.
#      phis_ERA_xCAM ==> ERA phis remapped horizontally to the CAM grid 
#
# aaa_{CAM,ERA}_xz{ERA,CAM}
# Indicates variable has been remapped horizontally AND vertically. So, e.g.
#      te_ERA_xzCAM ==> ERA temperature remapped horizontally to the CAM horizontal grid 
#                       and then also vertically interpoated to the CAM vertical grid
#-------------------------------------------------------------
def xRegrid( ExitAfterTemperature=False , 
             HorzInterpLnPs=False , 
             Use_ps_ERA_xCAM_in_vert=True ):


    StartTime = time.asctime( time.localtime(time.time()) )
    tic_overall = time.perf_counter()
    print( f"starting xRegrid {Gv.MySrc} _x_ {Gv.MyDst} at {StartTime} ")
    #-----------------------------------------
    # Horz remap of ERA PHIS to CAM horz-grid 
    # - No time dimension
    #-----------------------------------------
    Gv.phis_ERA_xCAM = erg.HorzRG( aSrc = Gv.phis_ERA , 
                                regrd = Gv.regrd , 
                                srcField=Gv.srcf , 
                                dstField=Gv.dstf , 
                                srcGridkey=Gv.srcHkey ,
                                dstGridkey=Gv.dstHkey )
    
    toc_here = time.perf_counter()
    pTime = f"Finished phis Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)
    
    #-----------------------------------------
    # Calculate difference between phis's
    #-----------------------------------------
    Dphis = Gv.phis_ERA_xCAM - Gv.phis_CAM
    
        
    #-----------------------------------------
    # Horz remap of ERA PS to CAM horz grid
    # log-exp bit is to reproduce W&O
    #-----------------------------------------
    if (HorzInterpLnPs==True):
        xfld_ERA = np.log( Gv.ps_ERA )
    else:
        xfld_ERA = Gv.ps_ERA
        
    Gv.ps_ERA_xCAM    = erg.HorzRG( aSrc = xfld_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTHkey ,
                                 dstGridkey=Gv.dstHkey )

    if (HorzInterpLnPs==True):
        Gv.ps_ERA_xCAM = np.exp( Gv.ps_ERA_xCAM ) 

    toc_here = time.perf_counter()
    pTime = f"Finished ps Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    #-----------------------------------------------------------------
    # Make 3D all ERA pressure fields on ERA's full (horz+vert) grid
    #-----------------------------------------------------------------
    Gv.pmid_ERA, Gv.pint_ERA, Gv.delp_ERA \
        = MkP.Pressure (am=Gv.amid_ERA ,
                        bm=Gv.bmid_ERA ,
                        ai=Gv.aint_ERA ,
                        bi=Gv.bint_ERA ,
                        ps=Gv.ps_ERA ,
                        p_00=Gv.p_00_ERA ,
                        Gridkey = Gv.srcTZHkey )

    #-----------------------------------------------------
    # Horz remap of ERA temperature onto CAM horz grid
    #------------------------------------------------------
    Gv.te_ERA_xCAM    = erg.HorzRG( aSrc = Gv.te_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTZHkey ,
                                 dstGridkey= Gv.dstHkey )
    
    # Save off a copy of te before any funny business 
    # happens, e.g., Williamson&Olson
    #------------------------------------------------
    Gv.te_ERA_xCAM_00 = copy.deepcopy( Gv.te_ERA_xCAM  )
    toc_here = time.perf_counter()
    pTime = f"Finished te_ERA Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    #-----------------------------------------------------------
    # Find "T_bot" and "P_bot" defined in Williamson & Olson 
    # as occurring at "the first level above 150m ... " 
    # This uses ERA temperature and surface pressure remapped  
    # to CAM horz grid (adapted from Isla's function)
    #-----------------------------------------------------------
    
    tic_150 = time.perf_counter()

    Gv.te_150, Gv.pmid_150, Gv.L150  =\
                    MkP.Pressure_TandP150(
                    am=Gv.amid_ERA ,
                    bm=Gv.bmid_ERA ,
                    ai=Gv.aint_ERA ,
                    bi=Gv.bint_ERA ,
                    ps=Gv.ps_ERA_xCAM ,
                    te=Gv.te_ERA_xCAM , 
                    p_00=Gv.p_00_ERA , 
                    Gridkey = Gv.dstTZHkey )


    toc_150 = time.perf_counter()
    pTime = f"Finding Te150 and P150 took  {toc_150 - tic_150:0.4f} seconds"
    print(pTime)

    
    #-------------------------------------------------------------------
    #                    "CAM surface pressure"
    #-------------------------------------------------------------------
    # We don't actually have ps from CAM, so we make a guess based on
    # ps_ERA_xCAM and te_bot as described above.  In a sense this is a 
    # vertical remapping, so we could call this variable ps_ERA_xzCAM, 
    # but we'll just call it ps_CAM ...
    #
    # Note Gv.ps_CAM is on CAM horz grid
    #-------------------------------------------------------------------
    
    
    Gv.ps_CAM = vrg.PsAdjust( phis=Gv.phis_ERA_xCAM, 
                            phis_CAM=Gv.phis_CAM, 
                            ps= Gv.ps_ERA_xCAM , 
                            pm150=Gv.pmid_150 , 
                            te150=Gv.te_150 , 
                            Gridkey=Gv.dstTZHkey  )

   
    #-----------------------------------------------------------------------------------------------
    # Now we creat full 4(3)D pressure fields on the ERA and CAM vertical grids. These are used for 
    # vertical interpolation below. Not clear what surface pressure to use when building ERA vertical
    # grid, i.e., ps_CAM or ps_ERA_xCAM.
    #
    # The WO2015 document seems to suggest ps_CAM but W&O code definitely uses ps_ERA_xCAM
    #-----------------------------------------------------------------------------------------------
    tic_P3D = time.perf_counter()
    
    # Surface pressure: Choose wisely
    #---------------------------------
    if (Use_ps_ERA_xCAM_in_vert == True):
        ps_FLD = Gv.ps_ERA_xCAM
    else:
        ps_FLD = Gv.ps_CAM


    #------------------------------------------------------------------------
    # These 3D pressure outputs are on ERA's vert grid and CAM's horz grid.
    #------------------------------------------------------------------------
    Gv.pmid_CAM_zERA, Gv.pint_CAM_zERA, Gv.delp_CAM_zERA \
        = MkP.Pressure (am=Gv.amid_ERA ,
                        bm=Gv.bmid_ERA ,
                        ai=Gv.aint_ERA ,
                        bi=Gv.bint_ERA ,
                        ps=ps_FLD , # What to use here seems key: ps_CAM or ps_ERA_xCAM
                        p_00=Gv.p_00_ERA ,
                        Gridkey = Gv.dstTZHkey )

    #-----------------------------------------------------------------------
    # Now we need to calculate the 3d pressure field for CAM on to which variables 
    # will be interpolated
    #---------------------------------------------------------------------------------
    # These are 3D pressure fields for CAM calculated assuming hybrid coordntaes
    #     P = a * P_00 + b * Ps  
    # i.e., NOT MPAS
    #---------------------------------------------------------------------------------
    Gv.pmid_CAM,Gv.pint_CAM,Gv.delp_CAM \
        = MkP.Pressure (am=Gv.amid_CAM ,
                        bm=Gv.bmid_CAM ,
                        ai=Gv.aint_CAM ,
                        bi=Gv.bint_CAM ,
                        ps=Gv.ps_CAM ,
                        p_00=Gv.p_00_CAM , 
                        Gridkey = Gv.dstTZHkey )

    #-----------------------------------------------------------
    # Log-pressure is preferred for vertical interpolation
    # per Williamson&Olson
    #-----------------------------------------------------------
    p_00 = 100_000. # Here we just use the sensible value of p_00
    lnpint_CAM = -7_000. * np.log( Gv.pint_CAM / p_00 )
    lnpmid_CAM = -7_000. * np.log( Gv.pmid_CAM / p_00 )
    lnpmid_CAM_zERA = -7_000. * np.log( Gv.pmid_CAM_zERA /p_00 )

    toc_P3D = time.perf_counter()
    pTime = f"Creating 3D P-fields etc., took   {toc_P3D - tic_P3D:0.4f} seconds"
    print(pTime)
    
    
    if ( Gv.doWilliamsonOlson == True ):
        tic_WO = time.perf_counter()
        print( "WilliamsonOlson surface " )
        #----------------------------------------------------------
        # Calculate extrapolated surface temperature using 
        # Williamson & Olson standard lapse rate approach
        #----------------------------------------------------------
        Gv.ts_extrap = vrg.TsExtrap( ps = Gv.ps_CAM ,
                                  pm150 = Gv.pmid_150 ,
                                  te150 = Gv.te_150 )

                              
        #--------------------------------------------------------
        # If Williamson & Olson treatment of surface layer 
        # is selected then correct te_ERA_xzCAM between 
        # pmid_150 and ps_CAM
        #-------------------------------------------------------
        Gv.te_WO  =    vrg.TeWO( te = Gv.te_ERA_xCAM ,
                              pmid = Gv.pmid_CAM_zERA ,
                              te150 = Gv.te_150 ,
                              pm150 = Gv.pmid_150 ,
                              ts = Gv.ts_extrap,
                              ps = Gv.ps_CAM, 
                              L150 = Gv.L150 ,
                              Gridkey = Gv.dstTZHkey )

        Gv.te_ERA_xCAM = copy.deepcopy( Gv.te_WO )
        toc_WO = time.perf_counter()
        pTime = f"Williamson Olson surface took  {toc_WO - tic_WO:0.4f} seconds"
        print(pTime)                       
    

    print(" going into vertical regrid of T " )
    Gv.te_ERA_xzCAM = vrg.VertRG( a_x  = Gv.te_ERA_xCAM ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM ,
                               Gridkey =Gv.dstTZHkey ,
                               kind = 'quadratic' ) #linea
    

    #-------------------------------------------------------------
    # return statement to assist in debugging and analysis
    #-------------------------------------------------------------
    if ( ExitAfterTemperature == True ):
        return


    #-------------------------------------------------------------
    # Continue on to regridding Q, U, V, W ...
    #-------------------------------------------------------------

    #--------------------
    #  Regridding of Q
    #---------------------
    print(" going into horz+vertical regrid of Q " )
    Gv.q_ERA_xzCAM , Gv.q_ERA_xCAM =\
                        fullRegrid( a_ERA = Gv.q_ERA ,
                                    zSrc = lnpmid_CAM_zERA ,
                                    zDst = lnpmid_CAM )
        
    Gv.q_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.q_ERA_xzCAM ,
                                  a_zERA = Gv.q_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )

    qx = SaturateQ( q=Gv.q_ERA_xzCAM , 
                    te=Gv.te_ERA_xzCAM ,
                    p=Gv.pmid_CAM, 
                    Gridkey = Gv.dstTZHkey )

    # Looks like Q can get very negative in spikes ... 
    # Maybe best fixed in qsat ... 
    # qx = np.where( qx >= 0. , qx, 0. )
    
    Gv.q_ERA_xzCAM =  copy.deepcopy(qx)




    #--------------------
    #  Regridding of U
    #---------------------
    print(" going into horz+vertical regrid of U " )
    Gv.u_ERA_xzCAM, Gv.u_ERA_xCAM = fullRegrid ( a_ERA = Gv.u_ERA ,
                                           zSrc = lnpmid_CAM_zERA ,
                                           zDst = lnpmid_CAM )

    Gv.u_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.u_ERA_xzCAM ,
                                  a_zERA = Gv.u_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )
    
    #--------------------
    #  Regridding of V
    #---------------------
    print(" going into horz+vertical regrid of V " )
    Gv.v_ERA_xzCAM, Gv.v_ERA_xCAM = fullRegrid ( a_ERA = Gv.v_ERA ,
                                           zSrc = lnpmid_CAM_zERA ,
                                           zDst = lnpmid_CAM )

    Gv.v_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.v_ERA_xzCAM ,
                                  a_zERA = Gv.v_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )
    
    #--------------------
    #  Regridding of W
    #---------------------
    print(" going into horz+vertical regrid of W " )
    Gv.w_ERA_xzCAM, Gv.w_ERA_xCAM = fullRegrid ( a_ERA = Gv.w_ERA ,
                                          zSrc = lnpmid_CAM_zERA ,
                                          zDst = lnpmid_CAM )

    Gv.w_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.w_ERA_xzCAM ,
                                  a_zERA = Gv.w_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )
    

    
    toc = time.perf_counter()
    pTime = f"Overall time in this function  {toc - tic_overall:0.4f} seconds"
    print(pTime)
        
    rcode =1 
    return rcode

######################################################################
# xRegrid_mpas_phys will generate MPAS files for physics side nudging
######################################################################
def xRegrid_mpas_phys( ExitAfterTemperature=False , 
             HorzInterpLnPs=False , 
             Use_ps_ERA_xCAM_in_vert=True ):


    

    StartTime = time.asctime( time.localtime(time.time()) )
    tic_overall = time.perf_counter()
    print( f"starting xRegrid_mpas_phys {Gv.MySrc} _x_ {Gv.MyDst} at {StartTime} ")

    #--------------------------------------
    # Find size of source (ERA) data 
    #--------------------------------------
    nt,nz,ny,nx = np.shape( Gv.te_ERA )
    #--------------------------------------
    # Find size of destintation (CAM_mpas) grid 
    #--------------------------------------
    nzp_dst,ncol_dst = np.shape( Gv.zgrid_CAM )
    # tile in time to conform to ERA
    zgrid_CAM = np.tile( Gv.zgrid_CAM , (nt,1) ).reshape( nt, nzp_dst, ncol_dst )
    Gv.zgrid_CAM=zgrid_CAM
    print( f' .. size of tiled Gv.zgrid_CAM {np.shape(Gv.zgrid_CAM)}')

    #-----------------------------------------------------------------
    # zgrid_CAM is on interfaces, so we need to calculate mid-level values
    #-----------------------------------------------------------------
    Gv.zgrido_CAM = np.zeros( ( nt, nzp_dst-1, ncol_dst ) )
    for z in np.arange( nzp_dst-1 ):
        Gv.zgrido_CAM[:,z,:] = 0.5*( zgrid_CAM[:,z,:] + zgrid_CAM[:,z,:] )         
    
    #-----------------------------------------------------------------
    # Make 3D all ERA pressure fields on ERA's full (horz+vert) grid
    #-----------------------------------------------------------------
    Gv.pmid_ERA, Gv.pint_ERA, Gv.delp_ERA \
        = MkP.Pressure (am=Gv.amid_ERA ,
                        bm=Gv.bmid_ERA ,
                        ai=Gv.aint_ERA ,
                        bi=Gv.bint_ERA ,
                        ps=Gv.ps_ERA ,
                        p_00=Gv.p_00_ERA ,
                        Gridkey = Gv.srcTZHkey )

    
    #-----------------------------------------------------------------
    # Make ERA geopotential height fields on ERA's full (horz+vert) grid
    #-----------------------------------------------------------------
    topo_ERA = Gv.phis_ERA / Con.grav()
    
    Gv.z3e_ERA, Gv.z3o_ERA \
         = MkP.GeopHeight ( te=Gv.te_ERA , pint=Gv.pint_ERA ,  pmid=Gv.pmid_ERA , qv=Gv.q_ERA , topo=topo_ERA , 
                        Gridkey = Gv.srcTZHkey )

    
    #-----------------------------------------
    # Horz remap of ERA PHIS to CAM horz-grid 
    # - No time dimension
    #-----------------------------------------
    Gv.phis_ERA_xCAM = erg.HorzRG( aSrc = Gv.phis_ERA , 
                                regrd = Gv.regrd , 
                                srcField=Gv.srcf , 
                                dstField=Gv.dstf , 
                                srcGridkey=Gv.srcHkey ,
                                dstGridkey=Gv.dstHkey )
    
    toc_here = time.perf_counter()
    pTime = f"Finished phis Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)
    #-----------------------------------------
    # Calculate difference between phis's
    #-----------------------------------------
    Dphis = Gv.phis_ERA_xCAM - Gv.phis_CAM


    #-----------------------------------------
    # Set up tiled surface topo
    #-----------------------------------------
    htopo_ERA_xCAM =Gv.phis_ERA_xCAM / Con.grav()
    htopo_ERA_xCAM = np.tile( htopo_ERA_xCAM , (nt,1) ).reshape( nt, ncol_dst )
    htopo_CAM = Gv.phis_CAM / Con.grav()
    htopo_CAM = np.tile( htopo_CAM , (nt,1) ).reshape( nt, ncol_dst )



    
    #-----------------------------------------
    # Horz remap of ERA PS to CAM horz grid
    # log-exp bit is to reproduce W&O
    #-----------------------------------------
    if (HorzInterpLnPs==True):
        xfld_ERA = np.log( Gv.ps_ERA )
    else:
        xfld_ERA = Gv.ps_ERA
        
    Gv.ps_ERA_xCAM    = erg.HorzRG( aSrc = xfld_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTHkey ,
                                 dstGridkey=Gv.dstHkey )

    if (HorzInterpLnPs==True):
        Gv.ps_ERA_xCAM = np.exp( Gv.ps_ERA_xCAM ) 

    toc_here = time.perf_counter()
    pTime = f"Finished ps Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)





    
    #-----------------------------------------------------
    # Horz remap of ERA geopotential height (mid) onto 
    # CAM horz grid
    #------------------------------------------------------
    Gv.z3o_ERA_xCAM    = erg.HorzRG( aSrc = Gv.z3o_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTZHkey ,
                                 dstGridkey= Gv.dstHkey )

    toc_here = time.perf_counter()
    pTime = f"Finished z3o Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    #-----------------------------------------------------
    # Horz remap of ERA temperature onto CAM horz grid
    #------------------------------------------------------
    Gv.te_ERA_xCAM    = erg.HorzRG( aSrc = Gv.te_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTZHkey ,
                                 dstGridkey= Gv.dstHkey )

    toc_here = time.perf_counter()
    pTime = f"Finished TE Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)


    #----------------------------------------------------
    # At this point we have ERA Temp on CAM's horz grid



    
    print(" going into vertical regrid of T " )
    Gv.te_ERA_xzCAM = vrg.VertRG( a_x  = Gv.te_ERA_xCAM ,
                               zSrc = Gv.z3o_ERA_xCAM  ,
                               zDst = Gv.zgrido_CAM,
                               Gridkey =Gv.dstTZHkey ,
                               kind = 'quadratic' ) #linea
    

    toc_here = time.perf_counter()
    pTime = f"Finished TE Vert Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    Gv.te_ERA_xzCAM = vrg.BottomFillZgrid ( htopo_ERA_xCAM= htopo_ERA_xCAM, 
                                            htopo_CAM= htopo_CAM,
                                            zgrido_CAM=Gv.zgrido_CAM, 
                                            a_CAM=Gv.te_ERA_xzCAM , 
                                            Method='lapse_extrapolate', 
                                            lapse_rate=-6.5 * 1.e-3, fill_value=300., Gridkey=Gv.dstTZHkey )

    
    toc_here = time.perf_counter()
    pTime = f"Finished TE bottom fill {toc_here - tic_overall:0.4f} seconds"
    print(pTime)



    
    
    #-------------------------------------------------------------
    # Continue on to regridding Q, U, V, W ...
    #-------------------------------------------------------------

    #--------------------
    #  Regridding of Q
    #---------------------
    print(" going into horz+vertical regrid of Q " )
    Gv.q_ERA_xzCAM , Gv.q_ERA_xCAM =\
                        fullRegrid( a_ERA = Gv.q_ERA ,
                                    zSrc = Gv.z3o_ERA_xCAM  ,
                                    zDst = Gv.zgrido_CAM )
        

    Gv.q_ERA_xzCAM = vrg.BottomFillZgrid (  htopo_ERA_xCAM= htopo_ERA_xCAM, 
                                            htopo_CAM= htopo_CAM,
                                            zgrido_CAM=Gv.zgrido_CAM, 
                                            a_CAM=Gv.q_ERA_xzCAM , 
                                            Method='lapse_extrapolate', 
                                            lapse_rate=0.0, fill_value=300., Gridkey=Gv.dstTZHkey )



    print( f' Digressing to make CAM pressures ' , flush=True )
    #----------------------------------------------------------------
    # Digression to make ps,pmid,pint,delp on (Time,Vert,Horz)_CAM 
    # using provisional unsaturated Q's. SHould be OK ....
    #----------------------------------------------------------------
    # First adjust PS
    #----------------------------------------------------------------
    Gv.ps_CAM = vrg.PsAdjust_simple(phis_ERA_xCAM= Gv.phis_ERA_xCAM,
                                    phis_CAM=Gv.phis_CAM, 
                                    ps_ERA_xCAM=Gv.ps_ERA_xCAM, 
                                    te_ERA_xzCAM=Gv.te_ERA_xzCAM, 
                                    q_ERA_xzCAM=Gv.q_ERA_xzCAM,
                                    Gridkey=Gv.dstTZHkey  )

    #-----------------------------------------------------------------
    # Make (Time x Vert x Horz)_CAM pressure fields on CAM's grid
    #-----------------------------------------------------------------
    Gv.pmid_CAM, Gv.pint_CAM, Gv.delp_ERA \
        = MkP.Press_from_ZTQ( te=Gv.te_ERA_xzCAM, z3e=Gv.zgrid_CAM, qv=Gv.q_ERA_xzCAM, ps=Gv.ps_CAM,  Gridkey=Gv.dstTZHkey )

    toc_here = time.perf_counter()
    pTime = f"Finished making pressures for CAM {toc_here - tic_overall:0.4f} seconds"
    print(pTime)


    print( f' Saturating Q ' , flush=True )
    qx = SaturateQ( q=Gv.q_ERA_xzCAM , 
                    te=Gv.te_ERA_xzCAM ,
                    p=Gv.pmid_CAM, 
                    Gridkey = Gv.dstTZHkey )

    Gv.q_ERA_xzCAM =  copy.deepcopy(qx)

    toc_here = time.perf_counter()
    pTime = f"Finished Q Horz and Vert Rgrd and fill and saturate {toc_here - tic_overall:0.4f} seconds"
    print(pTime , flush=True )


    #--------------------
    #  Regridding of U
    #---------------------
    print(" going into horz+vertical regrid of U " )
    Gv.u_ERA_xzCAM , Gv.u_ERA_xCAM =\
                        fullRegrid( a_ERA = Gv.u_ERA ,
                                    zSrc = Gv.z3o_ERA_xCAM  ,
                                    zDst = Gv.zgrido_CAM )
        

    Gv.u_ERA_xzCAM = vrg.BottomFillZgrid (  htopo_ERA_xCAM= htopo_ERA_xCAM, 
                                            htopo_CAM= htopo_CAM,
                                            zgrido_CAM=Gv.zgrido_CAM, 
                                            a_CAM=Gv.u_ERA_xzCAM , 
                                            Method='extrapolate_to_zero', 
                                            lapse_rate=0.0, fill_value=300., Gridkey=Gv.dstTZHkey )
    toc_here = time.perf_counter()
    pTime = f"Finished U Horz and Vert Rgrd and fill {toc_here - tic_overall:0.4f} seconds"
    print(pTime , flush=True )


    #--------------------
    #  Regridding of V
    #---------------------
    print(" going into horz+vertical regrid of V " )
    Gv.v_ERA_xzCAM , Gv.v_ERA_xCAM =\
                        fullRegrid( a_ERA = Gv.v_ERA ,
                                    zSrc = Gv.z3o_ERA_xCAM  ,
                                    zDst = Gv.zgrido_CAM )
        

    Gv.v_ERA_xzCAM = vrg.BottomFillZgrid (  htopo_ERA_xCAM= htopo_ERA_xCAM, 
                                            htopo_CAM= htopo_CAM,
                                            zgrido_CAM=Gv.zgrido_CAM, 
                                            a_CAM=Gv.v_ERA_xzCAM , 
                                            Method='extrapolate_to_zero', 
                                            lapse_rate=0.0, fill_value=300., Gridkey=Gv.dstTZHkey )
    toc_here = time.perf_counter()
    pTime = f"Finished V Horz and Vert Rgrd and fill {toc_here - tic_overall:0.4f} seconds"
    print(pTime , flush=True )

    #--------------------
    #  Regridding of W
    #---------------------
    print(" going into horz+vertical regrid of W " )
    Gv.w_ERA_xzCAM , Gv.w_ERA_xCAM =\
                        fullRegrid( a_ERA = Gv.w_ERA ,
                                    zSrc = Gv.z3o_ERA_xCAM  ,
                                    zDst = Gv.zgrido_CAM )
        

    Gv.w_ERA_xzCAM = vrg.BottomFillZgrid (  htopo_ERA_xCAM= htopo_ERA_xCAM, 
                                            htopo_CAM= htopo_CAM,
                                            zgrido_CAM=Gv.zgrido_CAM, 
                                            a_CAM=Gv.w_ERA_xzCAM , 
                                            Method='extrapolate_to_zero', 
                                            lapse_rate=0.0, fill_value=300., Gridkey=Gv.dstTZHkey )
    toc_here = time.perf_counter()
    pTime = f"Finished W Horz and Vert Rgrd and fill {toc_here - tic_overall:0.4f} seconds"
    print(pTime , flush=True )



    rcode=1
    return rcode
    


#########################################################################
def fullRegrid( a_ERA,  zSrc ,  zDst , kind='linear', ReturnVars=2 ):
    
    print("Horz RG in fullRegrid " )
    a_ERA_xCAM    = erg.HorzRG( aSrc = a_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTZHkey,
                                 dstGridkey=Gv.dstHkey )

    print("Vert RG in fullRegrid " )
    a_ERA_xzCAM = vrg.VertRG( a_x  = a_ERA_xCAM ,
                              zSrc = zSrc ,
                              zDst = zDst ,
                              Gridkey=Gv.dstTZHkey,
                              kind = kind )

    if (ReturnVars==1):
        return a_ERA_xzCAM
    if (ReturnVars==2):
        return a_ERA_xzCAM,a_ERA_xCAM

def SaturateQ ( q , te , p, Gridkey ):

    qsat = hum.qsat( p=p, T=te )
    qx=np.minimum( q , qsat )

    """
    print( "UGLY temporary Kluge: Don't saturate/ZERO top 3 levels of model " )
    if (Gridkey == 'tzyx' ):
        qx[:,0:2,:,:]= 0*q[:,0:2,:,:]
    if (Gridkey == 'tzc' ):
        qx[:,0:2,:]= 0*q[:,0:2,:]
    if (Gridkey == 'zyx' ):
        qx[0:2,:,:]= 0*q[0:2,:,:]
    if (Gridkey == 'zc' ):
        qx[0:2,:]= 0*q[0:2,:]
    
    """
    
    return qx

