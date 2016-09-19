#!/bin/bash
#===============================================================================
# Processes the file structure, extracts the necessary variables,
# and generates climatologies necessary for GISS-driven GEOS-Chem simulations 
# from a standard GISS run directory.
#
# Script may be intiated from anywhere, as it uses the ModelE configuration 
# file (~/.modelErc) to find the run directory.
# 
# Arguments: 
#  (1) Required: GISS simulation $RUNID (i.e., rundeck filename = $RUNID.R)
#  (2) Optional: 1 year or 2-year range to loop over to process and generate 
#       climatologies; should exclude any GCM spin-up years [default 0000-4999]
#
# Note: If the user has NetCDF4 installed, they can achieve large compression
# benefits by setting nccopy_opts below.
#
# Examples: ./prepare_gcap2.sh GCAP_RCP45 2000 2100
#           ./prepare_gcap2.sh GCAP_PD    2000
#           ./prepare_gcap2.sh GCAP_LGM
#
# Pre-Requisite Tools for script to run
#  (1) Working GISS ModelE configuration and post-processing tools 
#      (e.g., scaleacc)
#  (2) NetCDF library (http://www.unidata.ucar.edu/software/netcdf/)
#  (3) NetCDF Operators (NCO; http://nco.sourceforge.net)
#  (4) Climate Data Operators (CDO; https://code.zmaw.de/projects/cdo)
#  (5) R with ncdf4 library (R; http://r-project.org)
#
# Intial Version 2015-03-02 by Lee T. Murray (ltmurray@post.harvard.edu)
#===============================================================================
# Check arguments
if [ $# -eq 0 ]; then
    echo "You must pass the RUNID of the simulation as the first argument"
    echo "Optionally, you may pass a year or range of years to loop over,"
    echo "   otherwise default slowly loops over 0000-4999"
    echo "e.g., ./process_output.sh GCAP 2000 2005"
    exit 1;
fi
RUNID=$1
yr1=0; yr2=4999
if [ $# -eq 3 ]; then
yr1=$2; yr2=$3
fi
if [ $# -eq 2 ]; then
yr1=$2; yr2=$2
fi

# Get ModelE configuration variables
source ~/.modelErc
RUNDIR=$SAVEDISK/$RUNID
if [ ! -d $RUNDIR ]; then
    echo "Simulation with RUNID=${RUNID} not found; aborting"
    exit 1;
else
    echo "GISS run directory found at $RUNDIR"
fi

# Make sure GISS input files are linked in the directory for GEOS-Chem to find
cd $RUNDIR
./${RUNID}ln   #just to keep a log of what the codes have done, lshen 
################################################################################
## Organize the SUBDD files and ACC files
################################################################################
echo "processing daily files begins (lshen)"
# Use NetCDF classic format, no compression, greater compatability 
nccopy_opts=''        
# Use NetCDF4 format for compression if installed; requires HDF5 library
#nccopy_opts='-k3 -d9' 

# Loop through time
for year in $(seq $yr1 $yr2); do
    YYYY=$(printf '%04i' $year);
    echo "Checking ${YYYY}"
    months=(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC)
    for month in {1..12}; do #..12}; do
	MM=$(printf '%02i' $month); MMM=${months[$(($month-1))]}
	for day in {1..31}; do
	    DD=$(printf '%02i' $day);
               echo 'loop at: ' $YYYY $MM $DD	    
	    # If the file exists
	    if [ -f ${YYYY}${MM}${DD}.subdd${RUNID}.nc ]; then
		echo 'Processing: ' $YYYY $MM $DD
		
		# Unpack subdd file
		$SAVEDISK/$RUNID/mk_diags/scaleacc ${YYYY}${MM}${DD}.subdd${RUNID}.nc all
		
		# Copy and compress files to destination directory
		if [ ! -d $YYYY/$MM ]; then mkdir -p $YYYY/$MM; fi
		
		nccopy $nccopy_opts ${YYYY}${MM}${DD}.aijh12i${RUNID}.nc \
		    $YYYY/$MM/${YYYY}${MM}${DD}.aijh12i${RUNID}.nc
		nccopy $nccopy_opts ${YYYY}${MM}${DD}.aijh6${RUNID}.nc \
		    $YYYY/$MM/${YYYY}${MM}${DD}.aijh6${RUNID}.nc
		nccopy $nccopy_opts ${YYYY}${MM}${DD}.aijlh12${RUNID}.nc \
		    $YYYY/$MM/${YYYY}${MM}${DD}.aijlh12${RUNID}.nc
		
		# Remove original and intermediaries
		rm ${YYYY}${MM}${DD}.aijh12i${RUNID}.nc
		rm ${YYYY}${MM}${DD}.aijh6${RUNID}.nc
		rm ${YYYY}${MM}${DD}.aijlh12${RUNID}.nc
		#rm ${YYYY}${MM}${DD}.subdd${RUNID}.nc
		
	    fi
	    
	done # Day Loop
	
	# Copy and compress acc file
	if [ -f ${MMM}${YYYY}.acc${RUNID}.nc ]; then 
	    if [ ! -d $YYYY/acc ]; then mkdir -p $YYYY/acc; fi
	    nccopy $nccopy_opts ${MMM}${YYYY}.acc${RUNID}.nc \
		$YYYY/acc/${MMM}${YYYY}.acc${RUNID}.nc
	   # rm ${MMM}${YYYY}.acc${RUNID}.nc
	fi
	
    done # Month Loop
done     # Year Loop

echo "processing daily files ends"
echo "processing monthly files begins"
################################################################################
## Process the ACC files to generate monthly mean climatologies
################################################################################

OVERWRITE=0 # 1 = overwrite; 0 = do not overwrite

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Required for GEOS-Chem data directories 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ACC Variable / Long Name
# grnd_alb - GROUND ALBEDO IN VISIBLE RANGE (%)
# ptrop    - TROPOPAUSE PRESSURE (WMO) (mb)
# sst      - SEA SURFACE TEMPERATURE (C)
# tgrnd    - GROUND TEMPERATURE (C)
# prec     - PRECIPITATION (mm/day)

# Get ModelE configuration values
source ~/.modelErc
#source ~/.bashrc  #lshen added it, Mar 29
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Have the acc files already been processed?
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if [[ ! -f $SAVEDISK/$RUNID/$RUNID.acc.ymonmean.nc || \
      ! -f $SAVEDISK/$RUNID/$RUNID.acc.yearmean.nc || \
         $OVERWRITE -eq 1 ]]; then

# Create and change to working directory for processing GISS acc files
wrkDir=$SAVEDISK/$RUNID/process
if [ -d $wrkDir ]; then 
    rm -rf $wrkDir; 
fi
mkdir -p $wrkDir; cd $wrkDir
    
# Initialize time index
t=0

for yr in $(seq $yr1 $yr2); do
    YYYY=$(printf '%04i' $(( 10#$yr )) ); 
    accDir=$SAVEDISK/$RUNID/$YYYY/acc
    if [ ! -d $accDir ]; then continue; fi # Skip if directory does not exist
    echo $accDir
    
    months=(JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC)
    for M in {1..12}; do
	MM=$(printf '%02i' $M); MMM=${months[$(($M-1))]}
           echo 'Loop acc at: $yr $M'	
	if [ -f $accDir/${months[$(($M-1))]}${yr}.acc${RUNID}.nc ]; then
	    echo 'Processing acc file for: ' $yr ${months[$(($M-1))]}
	    
	    if [ $t -eq 0 ]; then # If first month
		YYYY0=$YYYY; MM0=$MM
	    fi
	    
	    $SAVEDISK/$RUNID/mk_diags/scaleacc $accDir/${MMM}${YYYY}.acc${RUNID}.nc aij
	    ncks -O -v prsurf,grnd_alb,ptrop,sst,tgrnd,prec \
		${MMM}${yr}.aij${RUNID}.nc $RUNID.$YYYY.$MM.nc
            # Define missing value flag for certain tracers
	    ncatted -a _FillValue,grnd_alb,o,f,-1e30 $RUNID.$YYYY.$MM.nc
      	    ncatted -a _FillValue,sst,o,f,-1e30 $RUNID.$YYYY.$MM.nc
            # Create degenerate record dim named "time" and populate with t
	    ncecat -O -u time $RUNID.$YYYY.$MM.nc $RUNID.$YYYY.$MM.nc
	    ncap2 -O -s "time=array($t,1,\$time); time@units=\"months since ${YYYY0}-${MM0}-15 00:00\"" \
		$RUNID.$YYYY.$MM.nc $RUNID.$YYYY.$MM.nc
	    
            # aij file no longer needed
	    rm ${MMM}${YYYY}.aij${RUNID}.nc
	   
	    t=$(($t+1))  

	fi
    done
done

# Concatenate files in time
ncrcat $RUNID.*.nc $RUNID.nc

# Create climatologies
cdo ymonmean $RUNID.nc $RUNID.acc.ymonmean.nc
cdo yearmean $RUNID.acc.ymonmean.nc $RUNID.acc.yearmean.nc

# Move climatologies
mv $RUNID.acc.ymonmean.nc $SAVEDISK/$RUNID
mv $RUNID.acc.yearmean.nc $SAVEDISK/$RUNID

# Delete working directory and remaining files
rm -rf $wrkDir

fi
echo "processing monthly files ends"
################################################################################
## Calculate mean pressure edges for vertical regridding
################################################################################
cd $RUNDIR
echo "lshen test the C++ code. It begins here"
echo "it works on : $DECKS_REPOSITORY/$RUNID.R"
# Parse the rundeck file for model resolution
if [ -f $DECKS_REPOSITORY/$RUNID.R ]; then
    E_RES=`grep AIC.RES_ $DECKS_REPOSITORY/$RUNID.R | awk -F'.' '{print $2}' | awk -F'_' '{print $2}'`
else
    echo "No rundeck found for RUNID=${RUNID}"
    exit 1;
fi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Generate mean pressure file for regridding routines
if [[ ! -f pressure.$RUNID.$E_RES.nc  || $OVERWRITE == 1 ]]; then
  echo 'Preparing ModelE pressure edges for re-gridding scripts...'

RSCRIPT=`which Rscript`
if [ $RSCRIPT == "" ]; then 
    echo "Error: prepare_gcap2 cannot find R installation."
    exit 1;
fi
cat << EOF > make_pres.R
#!$RSCRIPT
.libPaths(c('~/R', .libPaths()))  #lshen added this
require(ncdf4)
E_RES <- "$E_RES"
RUNID <- "$RUNID"
SAVEDISK <- "$SAVEDISK"

if ( E_RES == "F40" ) {

id <- nc_open(sprintf("%s/%s/%s.acc.ymonmean.nc",SAVEDISK,RUNID,RUNID))
lat <- ncvar_get( id, "lat" ); lon <- ncvar_get( id, "lon" )

prsurf <- ncvar_get( id, "prsurf" )

IM=144; JM=90; LM=40; LS1=24
PSF=984.e0; PTOP = 150.e0; PMTOP = .1e0
PSFMPT=PSF-PTOP; PSTRAT=PTOP-PMTOP

CALC_VERT_AMP_F40 <- function(PS) {
# !@sum  CALC_VERT_AMPK calculates air mass and pressure vertical arrays
# !@auth Jean Lerner/Gavin Schmidt
# !@ver  1.0
      # USE CONSTANT, only : bygrav
      # USE MODEL_COM, only : lm,ls1,dsig,sig,sige,ptop,psfmpt,lm_req
     # *     ,req_fac,req_fac_m,req_fac_d,pmtop
      # IMPLICIT NONE

      # REAL*8, INTENT(IN) :: P0 !@var P0 surface pressure (-PTOP) (mb)
      # INTEGER, INTENT(IN) :: LMAX !@var LMAX max level for calculation
# !@var AM mass at each level (kg/m2)
# !@var PDSIG pressure interval at each level (mb)
# !@var PMID mid-point pressure (mb)
      # REAL*8, INTENT(OUT), DIMENSION(LMAX) :: AM,PDSIG,PMID,PL
# !@var PEDN edge pressure (top of box) (mb)
      # REAL*8, INTENT(OUT), DIMENSION(LMAX+1) :: PEDN
      # INTEGER :: L  !@var L  loop variables

# C**** Calculate air mass, layer pressures
# C**** Note that only layers LS1 and below vary as a function of surface
# C**** pressure.
# C**** Note Air mass is calculated in (kg/m^2)

P0 <- PS-PTOP

grav = 9.80665e0
BYGRAV = 1e0/grav

IM=144; JM=90; LM=40; LS1=24
PSF=984.e0; PTOP = 150.e0; PMTOP = .1e0
PSFMPT=PSF-PTOP; PSTRAT=PTOP-PMTOP

PLbot <- c( PSF,   964e0, 942e0, 917e0, 890e0, 860e0, 825e0, # Pbot L=1,..
          785e0, 740e0, 692e0, 642e0, 591e0, 539e0, 489e0,   #      L=...
          441e0, 396e0, 354e0, 316e0, 282e0, 251e0, 223e0,
          197e0, 173e0,
          PTOP,                                              #      L=LS1
          128e0, 108e0,  90e0,  73e0,  57e0,  43e0,  31e0,   #      L=...
          20e0,   10e0,  5.62e0,  3.16e0,  1.78e0,   1.e0,
          .562e0,  .316e0,  .178e0,  PMTOP )                 #      L=..,LM+1

      SIGE = (PLbot-PTOP)/PSFMPT
      SIG  = (SIGE[1:LM]+SIGE[2:(LM+1)])*0.5
      DSIG =  SIGE[1:LM]-SIGE[2:(LM+1)]
      byDSIG  =  1./DSIG

          AM <- array(dim=c(LM))
          PDSIG <- array(dim=c(LM))
          PMID <- array(dim=c(LM))
          PL <- array(dim=c(LM))
          PEDN <- array(dim=c(LM+1))

      for ( L in 1:(LS1-1) ) {
        PL[L]   = P0
        PDSIG[L]= P0*DSIG[L]
        PMID[L] = SIG[L]*P0+PTOP
        PEDN[L] = SIGE[L]*P0+PTOP
        AM  [L] = PDSIG[L]*1e2*BYGRAV
      }
      for ( L in LS1:LM ) {
        PL[L]   = PSFMPT
        PDSIG[L]= PSFMPT*DSIG[L]
        PMID[L] = SIG[L]*PSFMPT+PTOP
        PEDN[L] = SIGE[L]*PSFMPT+PTOP
        AM  [L] = PDSIG[L]*1e2*BYGRAV
      }
      PEDN[LM+1] = SIGE[LM+1]*PSFMPT+PTOP

return(list(PMID=PMID,PEDN=PEDN,AM=AM))
}

GET_PEDGE <- function( PS ) {
        IM=144; JM=90; LM=40
        PEDGE <- array( dim=c(IM,JM,LM+1) )
        for ( i in 1:IM ) { for ( j in 1:JM ) {
                PEDGE[i,j,] <- CALC_VERT_AMP_F40(PS[i,j])\$PEDN # kg m-2
        }}
return(PEDGE)
}

pedge <- array( dim=c(144,90,41,12) )
for ( t in 1:12 ) {
pedge[,,,t] <- GET_PEDGE( prsurf[,,t] ) # mb
}
pmid <- ( pedge[,,1:LM,] + pedge[,,2:(LM+1),] ) / 2.0

dimX <- ncdim_def( "lon", "degrees_east", lon )
dimY <- ncdim_def( "lat", "degrees_north", lat )
dimZ <- ncdim_def( "lev", "level", 1:LM )
dimZ2 <- ncdim_def( "leve", "level", 1:(LM+1) )
dimT <- ncdim_def( "time", "month", 1:12 )

var1 <- ncvar_def( "pmid",  "mb", list(dimX, dimY, dimZ, dimT ), -1e30 )
var2 <- ncvar_def( "pedge", "mb", list(dimX, dimY, dimZ2, dimT ), -1e30 )

ncout <- nc_create( sprintf("pressure.%s.%s.nc",RUNID,E_RES), list( var1, var2 ) )
ncvar_put( ncout, var1, pmid )
ncvar_put( ncout, var2, pedge )
nc_close( ncout )

} else { print( "make_pres.R needs to be edited to accommodate resolutions other than F40" ) }
EOF
chmod 755 make_pres.R
./make_pres.R
fi

exit 0;
