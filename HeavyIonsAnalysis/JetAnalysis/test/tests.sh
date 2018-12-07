#!/bin/bash

helpmessage() {
   echo -e "usage: tests [options]\n"
   echo -e "   -h, --help     show (this) help message"
   echo -e "   -i, --input    directory for inputs"
   echo -e "   -j, --jobs     number of jobs in parallel"
   echo -e "   -o, --output   keep output of tests"
   echo -e "   -s, --sample   keep samples"
   echo -e "   -t, --travis   customisation for travis build"
}

ARGS=()

while [ $# -gt 0 ]; do
   case "$1" in
      -h|--help)     helpmessage; exit 0 ;;
      -i|--input)    inputdir="$2"; shift 2 ;;
      --input=*)     inputdir="${1#*=}"; shift ;;
      -j|--jobs)     njobs="$2"; shift 2 ;;
      --jobs=*)      njobs="${1#*=}"; shift ;;
      -o|--output)   keepoutput=1; shift ;;
      -s|--sample)   keepsample=1; shift ;;
      -t|--travis)   travis=1; shift ;;
      -*)            echo -e "invalid option: $1\n"; exit 1 ;;
      *)             ARGS+=("$1"); shift ;;
   esac
done

set -- "${ARGS[@]}"

timestamp=$(date +"%Y-%m-%d_%H_%M_%S")

# setup testing area
area=$(pwd)/forest-tests-$timestamp
mkdir $area

# samples
samples=(
   https://github.com/bi-ran/samples/raw/master/step2_PbPb_1030pre6_HIRECO.root
   https://github.com/bi-ran/samples/raw/master/step2_PbPb_1030pre6_RECO.root
   https://github.com/bi-ran/samples/raw/master/step2_PbPb_1030_t0streamer_RECO.root
)

[ $inputdir ] && {
   sampledir=$(readlink -f $inputdir)
   echo -e "\n  using samples from:\n$inputdir"
} || {
   sampledir=$area/samples
   echo -e "\n  downloading samples to:\n$sampledir\n"

   mkdir -p $sampledir
   wget -nv -P $sampledir/ ${samples[@]}

   [ $keepsample ] && echo -e "\n  samples will be kept:\n\E[34m$sampledir\E[0m"
}

# event-plane calibration
calibrations=(
   $CMSSW_BASE/src/HeavyIonsAnalysis/EventAnalysis/data/HeavyIonRPRcd_PbPb2018_offline.db
)

for c in ${calibrations[@]}; do
   cp $c $area/
done

# setup foresting configs
configs=(
   runForestAOD_HI_MB_103X.py
   runForestAOD_HI_MIX_103X.py
   runForestAOD_pponAA_MB_103X.py
   runForestAOD_pponAA_MIX_103X.py
   runForestAOD_pponAA_JEC_103X.py
   runForestAOD_pponAA_DATA_103X.py
)

[ $travis ] && {
   exclude=(
      runForestAOD_HI_MB_103X.py
      runForestAOD_HI_MIX_103X.py
      runForestAOD_pponAA_JEC_103X.py
   )

   for e in ${exclude[@]}; do
      configs=( ${configs[@]/$e} )
   done
}

for c in ${configs[@]}; do
   if [[ $c == *"DATA"* ]]; then
      input=$sampledir/$(basename ${samples[2]})
   elif [[ $c == *"HI"* ]]; then
      input=$sampledir/$(basename ${samples[0]})
   else
      input=$sampledir/$(basename ${samples[1]})
   fi
   export input="fileNames = cms.untracked.vstring('file:"${input}"'),"

   awk '
      !f && /fileNames/ { f=1 }
      f && /)/ { f=0 ; p=1 ; next }
      p { system ("echo $input") ; p=0 }
      !f
   ' $c | sed "s|HiForestAOD.root|${c%.*}.root|g" > $area/test_$c
   unset input
done

# cmsRun configs
pushd $area > /dev/null
mkdir logs

echo -ne "\n  running configs"

children=()
exitc=()
for c in ${configs[@]}; do
   cmsRun test_$c &> logs/test_${c%.*}.log &
   children+=("$!")

   if [ $njobs ] && [ ${#children[@]} -eq $njobs ]; then
      wait ${children[0]}; exitc+=("$?"); printf .
      unset children[0]; children=( "${children[@]}" )
   fi
done

for c in "${!children[@]}"; do
   wait ${children[$c]}
   exitc+=("$?")
   printf .
done

echo

for c in "${!exitc[@]}"; do
   config=${configs[$c]}

   if [ ${exitc[$c]} -eq 0 ]; then
      echo -e "\E[32mPASS: $config\E[0m"
   else
      fail=1
      mkdir -p fail
      echo -e "\E[31mFAIL: $config\E[0m"
      mv test_$config logs/test_${config%.*}.log fail/
   fi
done

if [ $fail ]; then
   echo -e "\n  please check the following file(s):"
   for l in $(ls fail/); do
      echo -e "  $area/fail/$l"
   done
fi

# finish, clean up
popd > /dev/null

[ ! $keepsample ] && [ ! $inputdir ] && rm -r $sampledir
[ ! $keepsample ] && [ ! $keepoutput ] && [ ! $fail ] && rm -r $area

[ $keepoutput ] && echo -e "\n  output kept:\n\E[34m$area\E[0m"

echo -e "\n  done\n"

exit $fail
